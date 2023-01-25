import os
from typing import Tuple

from automatically_annotate_mass2motifs.scores_matrix import ScoresMatrix
from automatically_annotate_mass2motifs.utils import load_pickled_file, store_as_json, store_pickled_file
from automatically_annotate_mass2motifs.moss_annotation import get_multiple_moss_annotations
from automatically_annotate_mass2motifs.mass2motif import load_mass2motifs_json
from automatically_annotate_mass2motifs.download_mass2motifs import download_motif_set_from_motifdb
from automatically_annotate_mass2motifs.clean_library_spectra import get_cleaned_and_binned_spectra
from automatically_annotate_mass2motifs.scores_matrix import create_similarity_matrix


def main(raw_library_spectra_file: str,
         motif_set_name: str,
         bin_size: float,
         ion_mode: str,
         directory_to_store_in_between_steps: str,
         similarity_thresholds: Tuple[float, ...],
         load_previous_annotation=False):
    """Pipeline for adding annotations to mass2motifs and storing in between files"""
    mass2motifs_file_name = os.path.join(directory_to_store_in_between_steps,
                                         f"mass2motifs_{motif_set_name}.json")
    annotated_mass2motifs_file_name = os.path.join(directory_to_store_in_between_steps,
                                                   f"annotated_mass2motifs_{motif_set_name}.json")
    if load_previous_annotation:
        assert os.path.isfile(annotated_mass2motifs_file_name), \
            "Annotated file does not yet exist, set load_previous_annotation to False"
        mass2motifs = load_mass2motifs_json(annotated_mass2motifs_file_name)
    elif os.path.isfile(mass2motifs_file_name):
        print("loading mass2motifs")
        mass2motifs = load_mass2motifs_json(mass2motifs_file_name)
    else:
        # Download mass2motifs from ms2lda.org and storing it
        mass2motifs = download_motif_set_from_motifdb(motif_set_name, bin_size)
        store_as_json(mass2motifs,
                      mass2motifs_file_name)
        print(f"mass2motifs are stored in: {mass2motifs_file_name}")

    spectra = get_cleaned_and_binned_spectra(raw_library_spectra_file,
                                             bin_size,
                                             ion_mode)

    scores_matrix_file_name = os.path.join(directory_to_store_in_between_steps,
                                           f"scores_matrix_{motif_set_name}.pickle")
    if os.path.isfile(mass2motifs_file_name):
        scores_matrix = load_pickled_file(scores_matrix_file_name)
        print(f"scores matrix is stored in: {scores_matrix_file_name}")
    else:
        print("Loading precalculated scores matrix and spectra")
        scores_matrix = create_similarity_matrix(mass2motifs,
                                                 spectra)
        store_pickled_file(scores_matrix, scores_matrix_file_name)

    spectra_selector = ScoresMatrix(spectra, scores_matrix)

    # Saves the annotations of the mass2motifs in the json containing the mass2motifs.
    mass2motifs = get_multiple_moss_annotations(mass2motifs, scores_matrix,
                                                similarity_thresholds,
                                                minimal_relative_support=60,
                                                maximal_relative_support_complement=80,
                                                save_in_between_file=annotated_mass2motifs_file_name)
    return spectra_selector, mass2motifs


if __name__ == "__main__":
    directory_to_store_in_between_steps = "../data/automatic_annotation/gnps_library_derived_mass2motifs/"

    main(raw_library_spectra_file=os.path.join(directory_to_store_in_between_steps,
                                               "../library_spectra/all_gnps.mgf"),
         motif_set_name="gnps library derived mass2motifs",
         bin_size=0.005,
         ion_mode="positive",
         directory_to_store_in_between_steps=directory_to_store_in_between_steps,
         similarity_thresholds=(0.1, 0.3, 0.6, 0.8, 0.9,),
         load_previous_annotation=True)
