import pandas as pd
from typing import List, Optional
import numpy as np


class Annotation:
    def __init__(self,
                 moss_annotations: Optional[pd.DataFrame],
                 minimal_similarity: float,
                 moss_minimal_relative_support: float,
                 moss_maximal_relative_support_complement: float,
                 nr_of_spectra_matching_mass2motif: int,
                 nr_of_spectra_not_matching_mass2motif: int,
                 minimal_number_of_matching_spectra: int = 0):

        # Settings
        self.minimal_similarity = minimal_similarity
        self.moss_minimal_relative_support = moss_minimal_relative_support
        self.moss_maximal_relative_support_complement = moss_maximal_relative_support_complement
        if moss_annotations is None:
            self.moss_annotations = pd.DataFrame(data=None,
                                                 columns=["s_abs", "c_abs", "s_rel", "c_rel", "diff_s_rel_and_c_rel"])
            self.moss_annotations.index.name = "smiles"
        else:
            self.moss_annotations = moss_annotations
        self.assert_correct_moss_annotations()
        # sort the annotations
        self.moss_annotations.sort_values("diff_s_rel_and_c_rel", inplace=True, ascending=False)
        self.nr_of_matching_spectra = nr_of_spectra_matching_mass2motif
        self.nr_of_not_matching_spectra = nr_of_spectra_not_matching_mass2motif
        self.minimal_number_of_matching_spectra = minimal_number_of_matching_spectra

    def assert_correct_moss_annotations(self):
        assert isinstance(self.moss_annotations, pd.DataFrame)
        assert list(self.moss_annotations.columns) == ["s_abs", "c_abs", "s_rel", "c_rel", "diff_s_rel_and_c_rel"], \
            "Expected different moss_annotations, use load_moss_results to load the results from a file"
        assert np.all(self.moss_annotations["s_rel"] >= self.moss_minimal_relative_support), \
            "moss minimal relative support given is wrong, since a lower s_rel was found"
        assert np.all(self.moss_annotations["c_rel"] <= self.moss_maximal_relative_support_complement),\
            "moss maximal relative support complement given is wrong, since a higher c_rel was found"

    def assert_equal(self, other):
        assert isinstance(other, Annotation)
        assert self.minimal_similarity == other.minimal_similarity
        assert self.moss_minimal_relative_support == other.moss_minimal_relative_support
        assert self.moss_maximal_relative_support_complement == other.moss_maximal_relative_support_complement
        pd.testing.assert_frame_equal(self.moss_annotations, other.moss_annotations)
        assert self.nr_of_matching_spectra == other.nr_of_matching_spectra
        assert self.nr_of_not_matching_spectra == other.nr_of_not_matching_spectra

    def __str__(self):
        pd.set_option('display.precision', 1)
        return str(self.moss_annotations.head())

    def to_dict(self):
        """Converts the object to a json storable format"""
        class_dict = self.__dict__.copy()
        annotations = self.moss_annotations.__deepcopy__()
        class_dict["moss_annotations"] = annotations.to_dict()
        return class_dict


def load_annotations_from_dict(dictionaries_from_json: List[dict])-> List[Annotation]:
    annotations = []
    for annotation_dict in dictionaries_from_json:
        dataframe_with_annotations = pd.DataFrame.from_dict(annotation_dict["moss_annotations"])
        dataframe_with_annotations.index.name="smiles"
        annotation = Annotation(
            moss_annotations=dataframe_with_annotations,
            minimal_similarity=annotation_dict["minimal_similarity"],
            moss_minimal_relative_support=annotation_dict["moss_minimal_relative_support"],
            moss_maximal_relative_support_complement=annotation_dict["moss_maximal_relative_support_complement"],
            nr_of_spectra_matching_mass2motif=annotation_dict["nr_of_matching_spectra"],
            nr_of_spectra_not_matching_mass2motif=annotation_dict["nr_of_not_matching_spectra"]
        )
        annotations.append(annotation)
    return annotations
