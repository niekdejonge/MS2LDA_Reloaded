import os
from automatically_annotate_mass2motifs.download_mass2motifs import download_motif_set_from_motifdb
from automatically_annotate_mass2motifs.search_matching_spectra_for_mass2motif import SelectSpectraContainingMass2Motif
from automatically_annotate_mass2motifs.clean_library_spectra import FilterLibrarySpectra
from automatically_annotate_mass2motifs.utils import store_pickled_file, load_pickled_file
from automatically_annotate_mass2motifs.moss_annotation import get_annotations


def get_cleaned_and_binned_spectra(raw_library_spectra_file, bin_size):
    filtered_library_file = os.path.splitext(raw_library_spectra_file)[0] + "_positive_ionmode_cleaned_filtered.pickle"
    filtered_and_binned_library_file = os.path.splitext(raw_library_spectra_file)[
                                           0] + f"_positive_ionmode_cleaned_filtered_binned_{bin_size}.pickle"
    # Check if there are already files containing cleaned or filtered spectra
    if os.path.isfile(filtered_and_binned_library_file):
        print(f'Loading the file {filtered_and_binned_library_file} ')
        binned_and_filtered_library_spectra = load_pickled_file(filtered_and_binned_library_file)
    elif os.path.isfile(filtered_library_file):
        print(f'Loading the file {filtered_library_file} ')
        library_spectra = FilterLibrarySpectra(filtered_library_file, already_cleaned=True)
        library_spectra.bin_spectra(bin_size)
        library_spectra.save_pickled_spectra()
        binned_and_filtered_library_spectra = library_spectra.spectra
    else:
        # Load and filter library spectra
        library_spectra = FilterLibrarySpectra(raw_library_spectra_file, already_cleaned=False, ion_mode="positive")
        # Optional store the pickled spectra (for quicker processing)
        library_spectra.save_pickled_spectra()
        library_spectra.bin_spectra(bin_size)
        library_spectra.save_pickled_spectra()
        binned_and_filtered_library_spectra = library_spectra.spectra
    return binned_and_filtered_library_spectra


def select_smiles_for_mass2motifs(raw_library_spectra_file,
                                  motif_set_name,
                                  bin_size,
                                  scores_matrix_file):
    """Downloads mass2motif set and creates a moss input smile with matching spectra
    """
    # download mass2motifs
    mass2motifs = download_motif_set_from_motifdb(motif_set_name, bin_size)

    # Load in and filter and bin library spectra
    binned_and_filtered_library_spectra = get_cleaned_and_binned_spectra(raw_library_spectra_file, bin_size)

    # search mass2motif
    if os.path.isfile(scores_matrix_file):
        print("Loading precalculated scores matrix")
        scores_matrix = load_pickled_file(scores_matrix_file)
        spectra_selector = SelectSpectraContainingMass2Motif(binned_and_filtered_library_spectra,
                                                             mass2motifs,
                                                             assert_correct_spectra=False,
                                                             scores_matrix=scores_matrix)
    else:
        spectra_selector = SelectSpectraContainingMass2Motif(binned_and_filtered_library_spectra,
                                                             mass2motifs,
                                                             assert_correct_spectra=False)
        store_pickled_file(spectra_selector.scores_matrix, scores_matrix_file)

    mass2motifs = get_annotations(spectra_selector,
                                  similarity_thresholds=(0.8,0.6,0.4),
                                  minimal_relative_support=50,
                                  maximal_relative_support_complement=80,
                                  save_figures_dir="../data/test/figures/",
                                  output_directory=None)
    # Do a grid search
    # Save mass2motifs, including annotation
    return mass2motifs


if __name__ == "__main__":
    print("hello")
    # mass2motifs = download_motif_set_from_motifdb("GNPS library derived Mass2Motifs", 0.005)
    # for mass2motif in mass2motifs[:5]:
    #     mass2motif.visualize()

    mass2motifs = select_smiles_for_mass2motifs(raw_library_spectra_file="../data/all_gnps.mgf",
                                                motif_set_name="GNPS library derived Mass2Motifs",
                                                bin_size=0.005,
                                                scores_matrix_file="../data/moss/gnps_library_derived_mass2motifs/mass2motif_scores_matrix_GNPS library derived Mass2Motifs.pickle")
