import os
import pandas as pd
from matchms import Fragments, Spectrum
import numpy as np
from automatically_annotate_mass2motifs.search_matching_spectra_for_mass2motif import (overlap_in_fragments,
                                                                                       SelectSpectraContainingMass2Motif,                                                                                       similarity_mass2motif_and_spectrum)
from automatically_annotate_mass2motifs.mass2motif import Mass2Motif
from tests.test_automatically_annotate_mass2motifs.generate_test_data import binned_spectra_005
from automatically_annotate_mass2motifs.moss_annotation import create_moss_input_file, run_moss, get_annotations


def test_create_moss_input_file(tmp_path):
    file_name = os.path.join(tmp_path, 'moss_input_file.smiles')
    create_moss_input_file(smiles_matching_mass2motif=["CCC=O", "CCC"],
                           smiles_not_matching_mass2motif= ["CCCC", "CC"],
                           output_csv_file_location=file_name)
    assert os.path.isfile(file_name), "expected file to be created"
    with open(file_name, "r") as file:
        result = file.readlines()
        assert result == ['0,0,CCC=O\n', '1,0,CCC\n', '2,1,CCCC\n', '3,1,CC\n'], \
            "expected a different output file"

def test_run_moss(tmp_path):
    moss_input_file = os.path.join(tmp_path, 'moss_input_file.smiles')
    moss_output_file = os.path.join(tmp_path, 'moss_output_file.sub')
    with open(moss_input_file, "w") as file:
        file.writelines(['0,0,CCC=O\n', '1,0,CCC\n', '2,1,CCCC\n', '3,1,CC\n'])
    run_moss(moss_input_file, moss_output_file, 0, 100)
    assert os.path.isfile(moss_output_file), "expected file to be created"

    with open(moss_output_file, "r") as output:
        result = output.readlines()
        assert result == ['id,description,nodes,edges,s_abs,s_rel,c_abs,c_rel\n',
                          '1,O=C-C-C,4,3,1,50.0,0,0.0\n',
                          '2,C(-C)-C,3,2,2,100.0,1,50.0\n']


def test_select_matching_smiles():
    spectra = binned_spectra_005()
    mass2motifs = [Mass2Motif(fragments=[100.025], fragment_probabilities= [0.5],
                            losses=[128.075, 200.725], loss_probabilities=[0.1, 0.25],
                            bin_size=0.05),
                   Mass2Motif(fragments=[50.025], fragment_probabilities=[0.6],
                              losses=[80.025, 128.075], loss_probabilities=[0.5, 0.1],
                              bin_size=0.05)]
    spectrum_selector = SelectSpectraContainingMass2Motif(spectra, mass2motifs)
    expected_result = pd.DataFrame([[0.55, 0.05],
                                    [0.00, 0.60]])
    assert np.all(expected_result == spectrum_selector.scores_matrix), "Expected different scores matrix"
    spectra_per_mass2motif = spectrum_selector.select_matching_spectra(0.05)
    expected_spectra_per_mass2motif = [[spectra[0], spectra[1]], [spectra[1]]]
    assert np.all(spectra_per_mass2motif == expected_spectra_per_mass2motif), "Different spectra were expected"
    smiles_matching, smiles_not_matching = spectrum_selector.select_unique_matching_and_non_matching_smiles(
        spectra_per_mass2motif[0])
    assert set(smiles_matching) == {'CN=C=O', 'C1CCCCC1'}
    assert smiles_not_matching == []
    smiles_matching, smiles_not_matching = spectrum_selector.select_unique_matching_and_non_matching_smiles(
        spectra_per_mass2motif[1])
    assert smiles_matching == ["C1CCCCC1"]
    assert smiles_not_matching == ["CN=C=O"]


def test_get_annotations():
    spectra = binned_spectra_005()
    mass2motifs = [Mass2Motif(fragments=[100.025], fragment_probabilities= [0.5],
                              losses=[128.075, 200.725], loss_probabilities=[0.1, 0.25],
                              bin_size=0.05),
                   Mass2Motif(fragments=[50.025], fragment_probabilities=[0.6],
                              losses=[80.025, 128.075], loss_probabilities=[0.5, 0.1],
                              bin_size=0.05)]
    spectrum_selector = SelectSpectraContainingMass2Motif(spectra, mass2motifs)
    annotated_mass2motifs = get_annotations(spectrum_selector,
                                            0.05, 0.1, 0.8)
    assert annotated_mass2motifs == mass2motifs # This does currently not check if the added moss annotations actually match
    assert len(annotated_mass2motifs[0].moss_annotations) == 1, "Expected only one annotation to be added"
    annotation1 = annotated_mass2motifs[0].moss_annotations[0]
    assert annotation1.minimal_similarity == 0.05
    assert annotation1.moss_minimal_relative_support == 0.1
    assert annotation1.moss_maximal_relative_support_complement == 0.8

    assert annotation1.moss_annotations.to_dict() == \
           {'s_abs': {'C': 2, 'O=C=N-C': 1, 'C1-C-C-C-C-C-1': 1},
            'c_abs': {'C': 0, 'O=C=N-C': 0, 'C1-C-C-C-C-C-1': 0},
            's_rel': {'C': 100.0, 'O=C=N-C': 50.0, 'C1-C-C-C-C-C-1': 50.0},
            'c_rel': {'C': 0.0, 'O=C=N-C': 0.0, 'C1-C-C-C-C-C-1': 0.0},
            'diff_s_rel_and_c_rel': {'C': 100.0, 'O=C=N-C': 50.0, 'C1-C-C-C-C-C-1': 50.0}}
