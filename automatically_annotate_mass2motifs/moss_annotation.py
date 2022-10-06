import os
import tempfile
from typing import List
import numpy as np
import pandas as pd
from tqdm import tqdm
import subprocess
from automatically_annotate_mass2motifs.utils import return_non_existing_file_name
from automatically_annotate_mass2motifs.mass2motif import Mass2Motif
from automatically_annotate_mass2motifs.search_matching_spectra_for_mass2motif import SelectSpectraContainingMass2Motif
from automatically_annotate_mass2motifs.annotation import load_moss_results, Annotation


def create_moss_input_file(smiles_matching_mass2motif,
                           smiles_not_matching_mass2motif,
                           output_csv_file_location):
    assert not os.path.isfile(output_csv_file_location), "The output file already exists"
    category = np.append(np.zeros(len(smiles_matching_mass2motif), dtype=np.int8),
                         np.ones(len(smiles_not_matching_mass2motif), dtype=np.int8))
    ids = np.arange(0, len(category))
    smiles_df = pd.DataFrame(
        {"ids": ids, "category": category, "smiles": smiles_matching_mass2motif + smiles_not_matching_mass2motif})
    smiles_df.to_csv(output_csv_file_location, index=False, header=False)


def run_moss(smiles_file_name,
             output_file_name,
             minimal_support,
             maximal_support_complement):
    print(os.path.dirname(__file__))
    assert os.path.isfile(smiles_file_name), "The smiles file does not exist"
    assert not os.path.isfile(output_file_name), "The output file already exists"
    subprocess.run(["java", "-cp",
                    os.path.join(os.path.dirname(__file__), "MolecularSubstructureMiner/moss.jar"),
                    "moss.Miner",
                    f"-s{minimal_support}", f"-S{maximal_support_complement}", smiles_file_name, output_file_name])


def get_annotations(matching_spectra_selector: SelectSpectraContainingMass2Motif,
                    similarity_threshold,
                    minimal_relative_support,
                    maximal_relative_support_complement,
                    output_directory = None) -> List[Mass2Motif]:
    """Creates a file containing the smiles matching and not matching a mass2motif readable by MOSS

    :param matching_spectra_selector:
    :param similarity_threshold:
    :param minimal_relative_support:
    :param maximal_relative_support_complement:
    :param output_directory:

    :return: A list of mass2motifs with Annotation added.
    """
    # Creates a temporary directory for all the intermediate files
    temp_dir = None
    if output_directory is None:
        temp_dir = tempfile.TemporaryDirectory()
        output_directory = temp_dir.name
    else:
        assert not os.path.isfile(output_directory), "A folder is expected, but a file is given"
        if not os.path.isdir(output_directory):
            os.mkdir(output_directory)

    matrix_of_matching_spectra = matching_spectra_selector.select_matching_spectra(similarity_threshold)
    mass2motifs = matching_spectra_selector.mass2motifs
    for i, mass2motif in enumerate(tqdm(mass2motifs,
                                        desc="Annotating mass2motifs")):
        list_of_matching_spectra = matrix_of_matching_spectra[i]
        smiles_matching_mass2motif, smiles_not_matching_mass2motif = matching_spectra_selector.select_non_matching_smiles(
            list_of_matching_spectra)

        # Creates output file names
        base_file_name = os.path.join(output_directory, f"mass2motif_{mass2motif.motif_name}_min_{similarity_threshold}")
        moss_input_file_name = return_non_existing_file_name(base_file_name + ".smiles")
        moss_output_file_name = return_non_existing_file_name(base_file_name + ".sub")

        create_moss_input_file(smiles_matching_mass2motif, smiles_not_matching_mass2motif, moss_input_file_name)
        run_moss(moss_input_file_name,
                 output_file_name=moss_output_file_name,
                 minimal_support=minimal_relative_support,
                 maximal_support_complement=maximal_relative_support_complement)
        moss_results = load_moss_results(moss_output_file_name)
        mass2motif.manual_annotation = Annotation(moss_results,
                                                  similarity_threshold,
                                                  minimal_relative_support,
                                                  maximal_relative_support_complement)
    # Deletes the temporary directory
    if temp_dir is not None:
        temp_dir.cleanup()
    return mass2motifs


if __name__ == "__main__":
    moss_input_file = os.path.join("../data/moss/example2.smiles")
    moss_output_file = os.path.join("../data/moss/example2.sub")
    with open(moss_input_file, "w") as file:
        file.writelines(['0,0,CCC=O\n', '1,0,CCC\n', '2,1,CCCC\n', '3,1,CC\n'])
    run_moss(moss_input_file, moss_output_file, 0, 100)

    with open(moss_output_file, "r") as output:
        result = output.readlines()
        print(result)