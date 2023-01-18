import os
from typing import List
import pandas as pd
from matchms import Fragments, Spectrum
import numpy as np
from automatically_annotate_mass2motifs.search_matching_spectra_for_mass2motif import \
    (overlap_in_fragments, ScoresMatrix, similarity_mass2motif_and_spectrum, create_similarity_matrix)
from automatically_annotate_mass2motifs.mass2motif import Mass2Motif
from tests.test_automatically_annotate_mass2motifs.generate_test_data import binned_spectra_005
from automatically_annotate_mass2motifs.moss_annotation import (create_moss_input_file, run_moss,
                                                                run_moss_wrapper, get_moss_annotation)
from automatically_annotate_mass2motifs.annotation import Annotation


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

def test_run_moss_wrapper(tmp_path):
    result = run_moss_wrapper(smiles_matching_mass2motif=["CCC=O", "CCC"],
                              smiles_not_matching_mass2motif=["CCCC", "CC"],
                              minimal_relative_support=30,
                              maximal_relative_support_complement=100)
    expected_result = pd.DataFrame([[1, 0 ,50.0, 0.0, 50.0],
                                    [2, 1, 100.0, 50.0, 50.0]],
                                   index=["O=C-C-C", "C(-C)-C"],
                                   columns=["s_abs", "c_abs", "s_rel", "c_rel", "diff_s_rel_and_c_rel"])
    expected_result.index.name = "smiles"
    pd.testing.assert_frame_equal(expected_result, result)

def test_add_moss_annotations():
    spectra = binned_spectra_005()
    mass2motif = Mass2Motif(fragments=[100.025], fragment_probabilities= [0.5],
                            losses=[128.075, 200.725], loss_probabilities=[0.1, 0.25],
                            bin_size=0.05)
    scores_matrix = create_similarity_matrix([mass2motif], spectra)
    spectrum_selector = ScoresMatrix(spectra, scores_matrix)
    annotation = get_moss_annotation(mass2motif, spectrum_selector, 0.05, 10, 80)
    assert isinstance(annotation, Annotation)
    expected_result = pd.DataFrame([[2, 0 ,100.0, 0.0, 100.0],
                                    [1, 0, 50.0, 0.0, 50.0],
                                    [1, 0, 50.0, 0.0, 50.0]],
                                   index=["C", "O=C=N-C", "C1-C-C-C-C-C-1"],
                                   columns=["s_abs", "c_abs", "s_rel", "c_rel", "diff_s_rel_and_c_rel"])
    expected_result.index.name = "smiles"
    pd.testing.assert_frame_equal(expected_result, annotation.moss_annotations)

if __name__ == "__main__":
    pass
