import ast
from typing import List, Optional, Union
import numpy as np
from matchms import Fragments
import json
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from matchms.plotting.spectrum_plots import plot_spectra_mirror, plot_spectrum

from automatically_annotate_mass2motifs.annotation import Annotation
from automatically_annotate_mass2motifs.utils import return_non_existing_file_name


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
                 manual_annotation: Optional[str] = None):
        assert all(fragments[i] <= fragments[i+1] for i in range(len(fragments) - 1)), "Expected sorted fragments"
        assert all(losses[i] <= losses[i+1] for i in range(len(losses) - 1)), "Expected sorted losses"

        self.fragments = Fragments(np.array(fragments), np.array(fragment_probabilities))
        self.losses = Fragments(np.array(losses), np.array(loss_probabilities))
        self.bin_size = bin_size
        self._assert_correct_bin_size()
        self.motif_name = motif_name
        self.motif_set_name = motif_set_name
        self.manual_annotation = manual_annotation
        self.moss_annotations: List[Annotation] = []
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
        class_dict["fragments"] = np.vstack((self.fragments.mz, self.fragments.intensities)).T.tolist()
        class_dict["losses"] = np.vstack((self.losses.mz, self.losses.intensities)).T.tolist()
        class_dict["moss_annotations"] = [moss_annotation.to_dict() for moss_annotation in self.moss_annotations]
        return class_dict

    def __eq__(self, other):
        return \
            self.fragments == other.fragments and \
            self.losses == other.losses and \
            self.bin_size == other.bin_size and \
            self.motif_name == other.motif_name and \
            self.motif_set_name == other.motif_set_name and \
            self.manual_annotation == other.manual_annotation

    def add_moss_annotation(self, new_annotation: Annotation):
        assert isinstance(new_annotation, Annotation), "Expected type Annotation"
        for annotation in self.moss_annotations:
            assert annotation.minimal_similarity != new_annotation.minimal_similarity or \
                   annotation.moss_minimal_relative_support != new_annotation.moss_minimal_relative_support or \
                   annotation.moss_maximal_relative_support_complement != new_annotation.moss_minimal_relative_support, \
                "The annotation with these settings is already stored in this Mass2Motif"
        self.moss_annotations.append(new_annotation)

    def visualize(self):
        fig = plt.figure()
        plot_mass2motif(self)
        return fig

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


def plot_mass2motif(mass2motif: Mass2Motif,
                    ax=None) -> plt.Axes:
    """Plots a mass2motif
    Code is largely taken from package "spectrum_utils".
    Parameters
    ----------
    mass2motif:
        The mass2motif that should be plotted
    ax:
        Axes instance on which to plot the Mass2motif. If None the current Axes
        instance is used.
    Returns
    -------
    plt.Axes
        The matplotlib Axes instance on which the spectra are plotted.
    """
    if ax is None:
        ax = plt.gca()

    # Top spectrum.
    plot_fragments(mass2motif.fragments, mirror_intensity=False, ax=ax, peak_color="darkblue")

    # Mirrored bottom spectrum.
    plot_fragments(mass2motif.losses, mirror_intensity=True, ax=ax, peak_color="red")
    ax.set_ylim(-1, 1)

    ax.axhline(0, color="#9E9E9E", zorder=10)
    if len(mass2motif.losses.mz) != 0 and len(mass2motif.fragments.mz) != 0:
        # Update axes so that both spectra fit.
        min_mz = max([0, np.floor(mass2motif.fragments.mz[0] / 100 - 1) * 100,
                      np.floor(mass2motif.losses.mz[0] / 100 - 1) * 100,])
        max_mz = max([np.ceil(mass2motif.fragments.mz[-1] / 100 + 1) * 100,
                      np.ceil(mass2motif.losses.mz[-1] / 100 + 1) * 100,])
    elif len(mass2motif.fragments.mz) != 0:
        min_mz = max(0, np.floor(mass2motif.fragments.mz[0] / 100 - 1) * 100)
        max_mz = np.ceil(mass2motif.fragments.mz[-1] / 100 + 1) * 100
    elif len(mass2motif.losses.mz) != 0:
        min_mz = max(0, np.floor(mass2motif.losses.mz[0] / 100 - 1) * 100)
        max_mz = np.ceil(mass2motif.losses.mz[-1] / 100 + 1) * 100
    else:
        assert False, "A Mass2motif needs to have at least one fragment or loss"
    ax.set_xlim(min_mz, max_mz)
    # makes sure the positive and negative axis are plotted correctly
    ax.yaxis.set_major_locator(mticker.AutoLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda x, pos: f"{abs(x):.2f}")
    )

    name1 = "Fragments"
    name2 = "Losses"

    x_text = 0.04 * (max_mz - min_mz)
    ax.text(x_text, 1, name1, ha="left", va="top", zorder=2, backgroundcolor="white")
    ax.text(x_text, -1, name2, ha="left", va="bottom", zorder=2, backgroundcolor="white")
    ax.set_title(f"Mass2Motif: {mass2motif.manual_annotation}")
    return ax


def plot_fragments(fragments: Fragments,
                   mirror_intensity=False,
                   ax=None,
                   grid=True,
                   peak_color="teal") -> plt.Axes:
    """
    Plots fragments
    """
    # pylint: disable=too-many-locals, too-many-arguments
    if ax is None:
        ax = plt.gca()

    def make_stems():
        """calculate where the stems of the spectrum peaks are going to be"""
        x = np.zeros([2, fragments.mz.size], dtype="float")
        y = np.zeros(x.shape)
        x[:, :] = np.tile(fragments.mz, (2, 1))
        y[1, :] = fragments.intensities
        return x, y

    x, y = make_stems()

    if mirror_intensity is True:
        y = -y
    ax.plot(x, y, color=peak_color, linewidth=1.0, marker="", zorder=5)

    y_max = 1.10
    ax.set_ylim(*(0, y_max) if not mirror_intensity else (-y_max, 0))

    ax.xaxis.set_minor_locator(mticker.AutoLocator())
    ax.yaxis.set_minor_locator(mticker.AutoLocator())
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    if grid in (True, "both", "major"):
        ax.grid(visible=True, which="major", color="#9E9E9E", linewidth=0.2)
    if grid in (True, "both", "minor"):
        ax.grid(visible=True, which="minor", color="#9E9E9E", linewidth=0.2)
    ax.set_axisbelow(True)

    ax.tick_params(axis="both", which="both", labelsize="small")
    y_ticks = ax.get_yticks()
    ax.set_yticks(y_ticks[y_ticks <= 1.0])

    ax.set_xlabel("m/z", style="italic")
    ax.set_ylabel("Probability")
    return ax
