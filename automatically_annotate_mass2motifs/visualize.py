import os
import numpy as np
from matchms import Fragments
from matchms.metadata_utils import is_valid_smiles
from matplotlib import pyplot as plt, ticker as mticker
from rdkit import Chem
from rdkit.Chem import Draw

from automatically_annotate_mass2motifs.annotation import Annotation
from automatically_annotate_mass2motifs.mass2motif import load_mass2motifs_json, Mass2Motif


def plot_fragments(fragments: Fragments,
                   mirror_intensity,
                   ax,
                   peak_color="teal") -> plt.Axes:
    """
    Plots fragments or losses, used in function plot_mass2motifs
    """
    # Calculates where the stems of the peaks are going to be
    x = np.zeros([2, fragments.mz.size], dtype="float")
    y = np.zeros(x.shape)
    x[:, :] = np.tile(fragments.mz, (2, 1))
    y[1, :] = fragments.intensities

    if mirror_intensity is True:
        y = -y
    ax.plot(x, y, color=peak_color, linewidth=1.0, marker="", zorder=5)

    # Add visualization of a grid
    ax.grid(visible=True, which="major", color="#9E9E9E", linewidth=0.2)
    ax.grid(visible=True, which="minor", color="#9E9E9E", linewidth=0.2)
    ax.set_axisbelow(True) # Makes sure the grid is the backgroung

    ax.set_xlabel("m/z", style="italic")
    ax.set_ylabel("Probability")
    return ax


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
    ax.set_ylim(-1.1, 1.1)

    # Adds a line between the top and bottom
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
    ax.xaxis.set_minor_locator(mticker.AutoLocator())
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, pos: f"{abs(x):.2f}"))

    # Adds fragments and losses to the plot
    x_text = 0.04 * (max_mz - min_mz)
    ax.text(x_text, 1, "Fragments", ha="left", va="top", zorder=2, backgroundcolor="white")
    ax.text(x_text, -1, "Losses", ha="left", va="bottom", zorder=2, backgroundcolor="white")

    ax.set_title(f"Mass2Motif: {mass2motif.manual_annotation}")
    return ax


def visualize_annotation(annotation: Annotation, nr_to_visualize):
    if nr_to_visualize > len(annotation.moss_annotations.index):
        smile_annotations = list(annotation.moss_annotations.index)
    else:
        smile_annotations = list(annotation.moss_annotations.index)[:nr_to_visualize]
    nr_of_colums = 3
    nr_of_rows = (len(smile_annotations)-1)//nr_of_colums+1
    fig = plt.figure(figsize=(5 * nr_of_colums, 5 * nr_of_rows))
    for i, smile in enumerate(smile_annotations):
        ax = fig.add_subplot(nr_of_rows, nr_of_colums, i + 1)
        if is_valid_smiles(smile):
            mol = Chem.MolFromSmiles(smile)
            im = Draw.MolToImage(mol)
            ax.imshow(im)
        ax.axis("off")
        ax.set_title(f"s_rel: {annotation.moss_annotations['s_rel'][smile]}\n"
                     f"c_rel: {annotation.moss_annotations['c_rel'][smile]}\n"
                     f"Smile: {smile}"
                     # f"s_abs: {self.moss_annotations['s_abs'][smile]}\n"
                     # f"c_abs: {self.moss_annotations['c_abs'][smile]}\n"
                     ,
                     )
    fig.suptitle(f"Minimal similarity: {annotation.minimal_similarity}\n"
                 f"Matching spectra: {annotation.nr_of_matching_spectra}\n"
                 f"Not matching spectra: {annotation.nr_of_not_matching_spectra}")
    fig.tight_layout()
    return fig

if __name__ == "__main__":
    base_directory = "../data/automatic_annotation/gnps_library_derived_mass2motifs/scores_matrix/"
    annotated_mass2motifs = load_mass2motifs_json(os.path.join(base_directory, "annotated_mass2motifs_gnps_library_derived_mass2motifs.json"))

    # fig = plt.figure()
    # plot_mass2motif(annotated_mass2motifs[0])
    # fig.show()
    for mass2motif in annotated_mass2motifs:
        fig =plt.figure()
        plot_mass2motif(mass2motif)
        fig.show()
        for annotation in mass2motif.moss_annotations:
            fig = visualize_annotation(annotation, 6)
            fig.show()