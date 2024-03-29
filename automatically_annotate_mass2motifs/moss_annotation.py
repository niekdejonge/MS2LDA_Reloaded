import os
import tempfile
from matchms import Spectrum
from typing import List, Tuple, Optional
from tqdm import tqdm
import numpy as np
import pandas as pd
import subprocess
from automatically_annotate_mass2motifs.mass2motif import Mass2Motif
from automatically_annotate_mass2motifs.scores_matrix import ScoresMatrix
from automatically_annotate_mass2motifs.annotation import Annotation, AnnotationSettings
from automatically_annotate_mass2motifs.utils import store_as_json, return_non_existing_file_name


def create_moss_input_file(smiles_matching_mass2motif: List[str],
                           smiles_not_matching_mass2motif: List[str],
                           output_csv_file_location):
    """Creates an input file that can be recognized by MOSS"""
    assert not os.path.isfile(output_csv_file_location), "The output file already exists"
    category = np.append(np.zeros(len(smiles_matching_mass2motif), dtype=np.int8),
                         np.ones(len(smiles_not_matching_mass2motif), dtype=np.int8))
    ids = np.arange(0, len(category))
    smiles_df = pd.DataFrame(
        {"ids": ids, "category": category, "smiles": smiles_matching_mass2motif + smiles_not_matching_mass2motif})
    smiles_df.to_csv(output_csv_file_location, index=False, header=False)


def run_moss(smiles_file_name,
             output_file_name,
             minimal_support: int,
             maximal_support_complement: int):
    """Runs MoSS to detect substructures"""
    assert os.path.isfile(smiles_file_name), "The smiles file does not exist"
    assert not os.path.isfile(output_file_name), "The output file already exists"
    moss_executable = os.path.join(os.path.dirname(__file__), "MolecularSubstructureMiner/moss.jar")
    assert os.path.isfile(moss_executable)
    assert 0 <= minimal_support <= 100 and 0 <= maximal_support_complement <= 100, \
        "The support should be specified as percentage"
    assert isinstance(minimal_support, int) and isinstance(maximal_support_complement, int), \
        "Expected an percentage for the minimal and maximal support"
    subprocess.run(["java", "-cp",
                    moss_executable,
                    "moss.Miner",
                    f"-s{minimal_support}", f"-S{maximal_support_complement}", smiles_file_name, output_file_name],
                   check=True)


def run_moss_wrapper(smiles_matching_mass2motif: List[str],
                     smiles_not_matching_mass2motif: List[str],
                     minimal_relative_support: int,
                     maximal_relative_support_complement: int):
    """Runs Moss, by creating a temp dir the in between files are removed after running"""
    # Creates a temporary directory for all the intermediate files
    with tempfile.TemporaryDirectory() as temp_dir_name:
        # Creates output file names
        moss_input_file_name = os.path.join(temp_dir_name, "moss.smiles")
        moss_output_file_name = os.path.join(temp_dir_name, "moss.sub")
        create_moss_input_file(smiles_matching_mass2motif, smiles_not_matching_mass2motif, moss_input_file_name)

        run_moss(moss_input_file_name, moss_output_file_name, minimal_relative_support, maximal_relative_support_complement)

        moss_results = load_moss_results(moss_output_file_name)
    return moss_results


def get_moss_annotation(mass2motif: Mass2Motif,
                        scores_matrix: ScoresMatrix,
                        annotation_settings: AnnotationSettings) -> Optional[Annotation]:
    """A wrapper to add a moss annotation to a Mass2Motif"""
    # Selects the spectra that match with the mass2motif
    matching_spectra, not_matching_spectra = scores_matrix.select_matching_spectra(annotation_settings.minimal_similarity,
                                                                                   mass2motif)
    smiles_matching_mass2motif, smiles_not_matching_mass2motif = select_unique_matching_and_non_matching_smiles(
        matching_spectra, not_matching_spectra)

    if len(smiles_matching_mass2motif) < annotation_settings.minimal_number_of_matching_spectra:
        # Small number of smiles are not reliable and in addition can result in very long runtimes,
        # since two almost duplicate large molecules will have an exponentially growing number of possible substructures.
        print(f"The number of matching spectra is {len(matching_spectra)}, "
              f"which is less than {annotation_settings.minimal_number_of_matching_spectra}")
        return Annotation(None, annotation_settings, len(smiles_matching_mass2motif),
                          len(smiles_not_matching_mass2motif))
    moss_results = run_moss_wrapper(smiles_matching_mass2motif, smiles_not_matching_mass2motif,
                                    annotation_settings.moss_minimal_relative_support,
                                    annotation_settings.moss_maximal_relative_support_complement)

    return Annotation(moss_results, annotation_settings, len(smiles_matching_mass2motif),
                      len(smiles_not_matching_mass2motif))


def get_multiple_moss_annotations(mass2motifs: List[Mass2Motif],
                                  scores_matrix: ScoresMatrix,
                                  similarity_thresholds: Tuple[float, ...],
                                  minimal_relative_support: int,
                                  maximal_relative_support_complement: int,
                                  save_in_between_file=None):
    # pylint: disable=too-many-arguments
    if save_in_between_file is not None:
        save_in_between_file = return_non_existing_file_name(save_in_between_file)

    for mass2motif in tqdm(mass2motifs, desc="Annotating mass2motifs"):
        for similarity_threshold in similarity_thresholds:
            if mass2motif.check_if_annotation_exists(similarity_threshold,
                                                     minimal_relative_support,
                                                     maximal_relative_support_complement):
                print(f"Mass2Motif {mass2motif.motif_name} was already annotated with similarity threshold: {similarity_threshold}")
            else:
                annotation_settings = AnnotationSettings(similarity_threshold,
                                                         minimal_relative_support,
                                                         maximal_relative_support_complement,
                                                         minimal_number_of_matching_spectra=10)
                annotation = get_moss_annotation(mass2motif, scores_matrix, annotation_settings)
                mass2motif.add_moss_annotation(annotation)
                if save_in_between_file is not None:
                    store_as_json(mass2motifs, save_in_between_file)
    return mass2motifs


def select_unique_matching_and_non_matching_smiles(matching_spectra: List[Spectrum],
                                                   not_matching_spectra: List[Spectrum]):
    """Selects all inchikeys not in smiles, and returns 1 smile for each"""
    inchikey_smiles_dict = create_inchikeys_smiles_dict(matching_spectra+not_matching_spectra)

    matching_inchikeys = []
    for spectrum in matching_spectra:
        inchikey = spectrum.get("inchikey")[:14]
        matching_inchikeys.append(inchikey)

    matching_inchikeys = set(matching_inchikeys)
    all_inchikeys = set(inchikey_smiles_dict.keys())
    not_matching_inchikeys = [inchikey for inchikey in all_inchikeys if inchikey not in matching_inchikeys]

    matching_smiles = [inchikey_smiles_dict[inchikey] for inchikey in matching_inchikeys]
    not_matching_smiles = [inchikey_smiles_dict[inchikey] for inchikey in not_matching_inchikeys]
    return matching_smiles, not_matching_smiles


def create_inchikeys_smiles_dict(spectra):
    """Creates a dictionary with all inchikeys with one of the smiles that represents this inchikey"""
    # todo, change to the most common smile for the inchikey
    inchikey_dict = {}
    for spectrum in spectra:
        inchikey = spectrum.get("inchikey")[:14]
        if inchikey not in inchikey_dict:
            inchikey_dict[inchikey] = spectrum.get("smiles")
    return inchikey_dict


def load_moss_results(file_name: str) -> Optional[pd.DataFrame]:
    """Loads in a moss results file and returns a pd.DataFrame"""
    with open(file_name, "r", encoding="utf-8") as file:
        lines = file.readlines()
        if len(lines) == 0:
            # Moss returns an empty dataframe, when no annotations are found.
            return None
        assert lines[0] == "id,description,nodes,edges,s_abs,s_rel,c_abs,c_rel\n", \
            "Expected the header id,description,nodes,edges,s_abs,s_rel,c_abs,c_rel in a moss output file"

    with open(file_name, "r", encoding="utf-8") as file:
        moss_results = pd.read_csv(file)
    assert isinstance(moss_results, pd.DataFrame)
    assert list(moss_results.columns) == ["id", "description", "nodes", "edges", "s_abs", "s_rel", "c_abs", "c_rel"], \
        "Expected different moss_annotations, use load_moss_results to load the results from a file"
    moss_results.rename(columns={"description": "smiles"}, inplace=True)
    moss_results.set_index("smiles", inplace=True)
    moss_results["diff_s_rel_and_c_rel"] = moss_results["s_rel"] - moss_results["c_rel"]
    moss_results = moss_results[["s_abs", "c_abs", "s_rel", "c_rel", "diff_s_rel_and_c_rel"]]
    return moss_results
