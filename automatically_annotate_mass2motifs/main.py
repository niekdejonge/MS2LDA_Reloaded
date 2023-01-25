import os

from automatically_annotate_mass2motifs.scores_matrix import ScoresMatrix
from automatically_annotate_mass2motifs.utils import load_pickled_file
from automatically_annotate_mass2motifs.moss_annotation import get_multiple_moss_annotations
from automatically_annotate_mass2motifs.mass2motif import load_mass2motifs_json


def pipeline_from_scratch(raw_library_spectra_file,
                          motif_set_name,
                          directory_to_store_in_between_steps):
    pass
    # # Process library_spectra
    #
    # # Save the in between processed spectra.
    # # Save the binned spectra
    # binned_and_filtered_library_spectra =
    # # Load Mass2motifs
    # # Store the Mass2motifs as json
    # mass2motifs =
    #
    # # Create a similarity matrix
    # # Store the similarity matrix as pickled file
    # scores_matrix = create_similarity_matrix(mass2motifs,
    #                                          binned_and_filtered_library_spectra)
    # spectra_selector = ScoresMatrix(binned_and_filtered_library_spectra, scores_matrix)
    # store_pickled_file(spectra_selector.scores_matrix, scores_matrix_file)
    #
    # # Add annotations to mass2motifs
    # # Replace previously stored json with mass2motifs, with mass2motifs containing annotations.
    #
    # # Store visualizations of results


def pipeline_from_in_between_results(mass2motif_file,
                                     binned_spectra_file,
                                     scores_matrix_file,
                                     similarity_thresholds,
                                     results_file_name):
    print("loading mass2motifs")
    mass2motifs = load_mass2motifs_json(mass2motif_file)
    print("Loading precalculated scores matrix and spectra")
    scores_matrix = ScoresMatrix(load_pickled_file(binned_spectra_file),
                                 load_pickled_file(scores_matrix_file))
    mass2motifs = get_multiple_moss_annotations(mass2motifs, scores_matrix,
                                                similarity_thresholds,
                                                60,
                                                80,
                                                results_file_name)

if __name__ == "__main__":
    base_directory = "../data/automatic_annotation/gnps_library_derived_mass2motifs/scores_matrix/"

    pipeline_from_in_between_results(os.path.join(base_directory, "all_annotated_mass2motifs_gnps_library_derived_mass2motifs.json"),
                                     os.path.join(base_directory, "../../library_spectra/all_gnps_positive_ionmode_cleaned_filtered_binned_0.005.pickle"),
                                     os.path.join(base_directory, "mass2motif_scores_matrix_GNPS library derived Mass2Motifs(1).pickle"),
                                     [0.1, 0.3, 0.6, 0.8, 0.9],
                                     os.path.join(base_directory, "all_annotated_mass2motifs_gnps_library_derived_mass2motifs.json"))
    # annotated_mass2motifs = load_mass2motifs_json(os.path.join(base_directory, "annotated_mass2motifs_gnps_library_derived_mass2motifs.json"))
    # for mass2motif in annotated_mass2motifs:
    #     plt.Figure()
    #     fig = plot_mass2motif(mass2motif)
    #     fig.show()
    #     for annotation in mass2motif.moss_annotations:
    #         fig = annotation.visualize(6)
    #         fig.show()