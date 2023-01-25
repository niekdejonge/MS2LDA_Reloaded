import os
import pandas as pd
import json
from automatically_annotate_mass2motifs.annotation import Annotation, load_annotations_from_dict, AnnotationSettings
from automatically_annotate_mass2motifs.moss_annotation import load_moss_results
from automatically_annotate_mass2motifs.utils import store_as_json
from tests.test_automatically_annotate_mass2motifs.generate_test_data import generate_moss_results_file, generate_annotation


def test_create_annotation(tmp_path):
    file_name = os.path.join(tmp_path, "moss_results")
    generate_moss_results_file(file_name)
    moss_annotations = load_moss_results(file_name)
    annotation_settings = AnnotationSettings(minimal_similarity=0.3,
                                             moss_minimal_relative_support=40.0,
                                             moss_maximal_relative_support_complement=70.0,
                                             minimal_number_of_matching_spectra=10)
    annotation = Annotation(moss_annotations=moss_annotations,
                            annotation_settings=annotation_settings,
                            nr_of_spectra_matching_mass2motif=10,
                            nr_of_spectra_not_matching_mass2motif=100)
    result = annotation.moss_annotations
    assert isinstance(result, pd.DataFrame)
    assert list(result.columns) == ["s_abs", "c_abs", "s_rel", "c_rel", "diff_s_rel_and_c_rel"], \
        "different column names were expected"
    assert result.index.name == "smiles"


def test_create_annotation_with_none():
    annotation_settings = AnnotationSettings(minimal_similarity=0.3,
                                             moss_minimal_relative_support=40.0,
                                             moss_maximal_relative_support_complement=70.0,
                                             minimal_number_of_matching_spectra=10)
    annotation = Annotation(moss_annotations=None,
                            annotation_settings=annotation_settings,
                            nr_of_spectra_matching_mass2motif=10,
                            nr_of_spectra_not_matching_mass2motif=100)
    result = annotation.moss_annotations
    assert isinstance(result, pd.DataFrame)
    assert list(result.columns) == ["s_abs", "c_abs", "s_rel", "c_rel", "diff_s_rel_and_c_rel"], \
        "different column names were expected"
    assert result.index.name == "smiles"


def test_store_and_load_annotation_json(tmp_path):
    annotation = generate_annotation()
    file_name = os.path.join(tmp_path, "stored_annotation.json")
    store_as_json(annotation,
                  file_name)
    with open(file_name, 'r') as file:
        json_data = json.load(file)
        loaded_annotation = load_annotations_from_dict(json_data)
    assert len(loaded_annotation) == 1
    annotation.assert_equal(loaded_annotation[0])
