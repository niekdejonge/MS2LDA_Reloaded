import os
import pandas as pd

from automatically_annotate_mass2motifs.annotation import load_moss_results, Annotation
from tests.test_automatically_annotate_mass2motifs.generate_test_data import generate_moss_results_file


def test_load_moss_file(tmp_path):
    file_name = os.path.join(tmp_path, "moss_results")
    generate_moss_results_file(file_name)
    result = load_moss_results(file_name)
    assert isinstance(result, pd.DataFrame)
    assert list(result.columns) == ["id", "description", "nodes", "edges", "s_abs", "s_rel", "c_abs", "c_rel"], \
            "Expected different moss_annotations, use load_moss_reulsts to load the results from a file"


def test_load_empty_moss_file(tmp_path):
    file_name = os.path.join(tmp_path, "empty_moss_file")
    # Generate an empty file
    with open(file_name, "w") as moss_results_file:
        moss_results_file.write("")
    result = load_moss_results(file_name)
    assert result is None


def test_create_annotation(tmp_path):
    file_name = os.path.join(tmp_path, "moss_results")
    generate_moss_results_file(file_name)
    moss_annotations = load_moss_results(file_name)
    annotation = Annotation(moss_annotations=moss_annotations,
                            minimal_similarity=0.3,
                            moss_minimal_relative_support=40.0,
                            moss_maximal_relative_support_complement=70.0)
    result = annotation.moss_annotations
    assert isinstance(result, pd.DataFrame)
    assert list(result.columns) == ["s_abs", "c_abs", "s_rel", "c_rel", "diff_s_rel_and_c_rel"], \
        "different column names were expected"
    assert result.index.name == "smiles"


def test_annotation_to_dict(tmp_path):
    file_name = os.path.join(tmp_path, "moss_results")
    generate_moss_results_file(file_name)
    moss_annotations = load_moss_results(file_name)
    annotation = Annotation(moss_annotations=moss_annotations,
                            minimal_similarity=0.3,
                            moss_minimal_relative_support=40.0,
                            moss_maximal_relative_support_complement=70.0)
    assert annotation.to_dict() == \
           {'minimal_similarity': 0.3, 'moss_minimal_relative_support': 40.0,
            'moss_maximal_relative_support_complement': 70.0,
            'moss_annotations': {'s_abs': {'O(-C-C)-C(-C)-C-C': 6, 'O(-C-C)-C(-C)-C-C-C': 5},
                                 'c_abs': {'O(-C-C)-C(-C)-C-C': 5737, 'O(-C-C)-C(-C)-C-C-C': 4266},
                                 's_rel': {'O(-C-C)-C(-C)-C-C': 54.545456, 'O(-C-C)-C(-C)-C-C-C': 45.454544},
                                 'c_rel': {'O(-C-C)-C(-C)-C-C': 28.162584, 'O(-C-C)-C(-C)-C-C-C': 20.941534},
                                 'diff_s_rel_and_c_rel': {'O(-C-C)-C(-C)-C-C': 26.382872000000003,
                                                          'O(-C-C)-C(-C)-C-C-C': 24.513009999999998}}}
