from typing import List, Optional
import numpy as np
from matchms import Fragments
import json

from automatically_annotate_mass2motifs.annotation import Annotation, load_annotations_from_dict


class Mass2Motif:
    """Stores Mass2Motif information

    """

    def __init__(self,
                 fragments: List[float],
                 fragment_probabilities: List[float],
                 losses: List[float],
                 loss_probabilities: List[float],
                 bin_size: float,
                 motif_name: Optional[str] = None,
                 motif_set_name: Optional[str] = None,
                 manual_annotation: Optional[str] = None,
                 moss_annotations: Optional[List[Annotation]] = None):
        assert all(fragments[i] <= fragments[i+1] for i in range(len(fragments) - 1)), "Expected sorted fragments"
        assert all(losses[i] <= losses[i+1] for i in range(len(losses) - 1)), "Expected sorted losses"

        self.fragments = Fragments(np.array(fragments), np.array(fragment_probabilities))
        self.losses = Fragments(np.array(losses), np.array(loss_probabilities))
        self.bin_size = bin_size
        self._assert_correct_bin_size()
        self.motif_name = motif_name
        self.motif_set_name = motif_set_name
        self.manual_annotation = manual_annotation
        if moss_annotations is None:
            self.moss_annotations: List[Annotation] = []
        else:
            self.moss_annotations = moss_annotations
        self._assert_correct_types()

    def _assert_correct_types(self):
        if self.manual_annotation is not None:
            assert isinstance(self.manual_annotation, str), \
                f"Expected a string for manual annotations, got {type(self.manual_annotation)}"
        if self.moss_annotations is not None:
            assert isinstance(self.moss_annotations, list), \
                f"Expected a list with Annotation object, got {type(self.moss_annotations)}"
            for annotation in self.moss_annotations:
                assert isinstance(annotation, Annotation), \
                    f"Expected a list with Annotation objects, got a list with {type(annotation)}"

    def _assert_correct_bin_size(self):
        assert isinstance(self.bin_size, float), "bin size is expected to be float"
        nr_of_decimals = len(str(self.bin_size).split(".", 1)[1])
        assert nr_of_decimals < 10, "An bin_size with this many decimals is unexpected"
        for mz in self.fragments.mz:
            assert round(mz%self.bin_size, 10) == self.bin_size*0.5, \
                f"Incorrect bin size of {self.bin_size} for mz of {mz}"
        for mz_loss in self.losses.mz:
            assert round(mz_loss%self.bin_size, 10) == self.bin_size*0.5, \
                f"Incorrect bin size of {self.bin_size} for mz of {mz_loss}"

    def __str__(self):
        return str({"fragments": self.fragments.mz,
                    "fragment_probabilities": self.fragments.intensities,
                    "losses": self.losses.mz,
                    "loss_probabilities": self.losses.intensities})

    def to_dict(self) -> dict:
        """Return a dictionary representation of a spectrum, to make it possible to store as json"""
        class_dict = self.__dict__.copy()
        del class_dict["fragments"]
        del class_dict["losses"]
        class_dict["fragment_mz"] = self.fragments.mz.tolist()
        class_dict["fragment_probabilities"] = self.fragments.intensities.tolist()
        class_dict["loss_mz"] = self.losses.mz.tolist()
        class_dict["loss_probabilities"] = self.losses.intensities.tolist()
        class_dict["moss_annotations"] = [moss_annotation.to_dict() for moss_annotation in self.moss_annotations]
        return class_dict

    def __eq__(self, other):
        return \
            self.fragments == other.fragments and \
            self.losses == other.losses and \
            self.bin_size == other.bin_size and \
            self.motif_name == other.motif_name and \
            self.motif_set_name == other.motif_set_name and \
            self.manual_annotation == other.manual_annotation and \
            self.moss_annotations == other.moss_annotations

    def assert_equal(self, other: "Mass2Motif"):
        assert isinstance(other, Mass2Motif)
        assert self.fragments == other.fragments
        assert self.losses == other.losses
        assert self.bin_size == other.bin_size
        assert self.motif_name == other.motif_name
        assert self.motif_set_name == other.motif_set_name
        assert self.manual_annotation == other.manual_annotation
        for i, annotation in enumerate(self.moss_annotations):
            annotation.assert_equal(other.moss_annotations[i])

    def check_if_annotation_exists(self,
                                   minimal_similarity,
                                   moss_minimal_relative_support,
                                   moss_maximal_relative_support_complement) -> bool:
        """Checks if there is already an annotation with these settings stored."""
        for annotation in self.moss_annotations:
            if (annotation.minimal_similarity == minimal_similarity and
                    annotation.moss_minimal_relative_support == moss_minimal_relative_support and
                    annotation.moss_maximal_relative_support_complement == moss_maximal_relative_support_complement):
                return True
        return False

    def add_moss_annotation(self, new_annotation: Annotation):
        assert isinstance(new_annotation, Annotation), "Expected type Annotation"
        for annotation in self.moss_annotations:
            assert annotation.minimal_similarity != new_annotation.minimal_similarity or \
                   annotation.moss_minimal_relative_support != new_annotation.moss_minimal_relative_support or \
                   annotation.moss_maximal_relative_support_complement != new_annotation.moss_minimal_relative_support, \
                "The annotation with these settings is already stored in this Mass2Motif"
        self.moss_annotations.append(new_annotation)


def load_mass2motifs_json(file_name) -> List[Mass2Motif]:
    with open(file_name, 'r') as file:
        mass2motifs = []
        for spectrum_dict in json.load(file):
            mass2motif = Mass2Motif(
                fragments=spectrum_dict["fragment_mz"],
                fragment_probabilities=spectrum_dict["fragment_probabilities"],
                losses=spectrum_dict["loss_mz"],
                loss_probabilities=spectrum_dict["loss_probabilities"],
                bin_size=spectrum_dict["bin_size"],
                motif_name=spectrum_dict["motif_name"],
                motif_set_name=spectrum_dict["motif_set_name"],
                manual_annotation=spectrum_dict["manual_annotation"],
                moss_annotations=load_annotations_from_dict(spectrum_dict["moss_annotations"])
            )
            mass2motifs.append(mass2motif)
    return mass2motifs
