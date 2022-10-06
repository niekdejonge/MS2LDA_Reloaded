import ast
from typing import List, Optional, Union
import numpy as np
from matchms import Fragments
import json
from automatically_annotate_mass2motifs.annotation import Annotation


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
        """Return a dictionary representation of a spectrum."""
        class_dict = self.__dict__
        class_dict["fragments"] = np.vstack((self.fragments.mz, self.fragments.intensities)).T.tolist()
        class_dict["losses"] = np.vstack((self.losses.mz, self.losses.intensities)).T.tolist()
        # todo add Annotation.
        return class_dict

    def __eq__(self, other):
        return \
            self.fragments == other.fragments and \
            self.losses == other.losses and \
            self.bin_size == other.bin_size and \
            self.motif_name == other.motif_name and \
            self.motif_set_name == other.motif_set_name and \
            self.manual_annotation == other.manual_annotation


def save_mass2motifs_json(mass2motifs: Union[List[Mass2Motif], Mass2Motif],
                          file_name):
    if not isinstance(mass2motifs, list):
        mass2motifs = [mass2motifs]
    json_str = []
    for mass2motif in mass2motifs:
        json_str.append(json.dumps(mass2motif.to_dict()))
    with open(file_name, "w", encoding="utf-8") as file:
        json.dump(json_str, file, indent=3)


def load_mass2motifs_json(file_name) -> List[Mass2Motif]:
    with open(file_name, 'rb') as file:
        mass2motifs = []
        for spectrum_dict in json.load(file):
            spectrum_dict = ast.literal_eval(spectrum_dict)
            mass2motifs.append(Mass2Motif(fragments=list(np.array(spectrum_dict["fragments"])[:, 0]),
                                          fragment_probabilities=list(np.array(spectrum_dict["fragments"])[:, 1]),
                                          losses=list(np.array(spectrum_dict["losses"])[:, 0]),
                                          loss_probabilities=list(np.array(spectrum_dict["losses"])[:, 1]),
                                          bin_size=spectrum_dict["bin_size"], motif_name=spectrum_dict["motif_name"],
                                          motif_set_name=spectrum_dict["motif_set_name"],
                                          manual_annotation=spectrum_dict["manual_annotation"]))
    return mass2motifs
