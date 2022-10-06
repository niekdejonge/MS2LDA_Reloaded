import os
import matplotlib.pyplot as plt
import pandas as pd
from automatically_annotate_mass2motifs.download_mass2motifs import get_annotations
from automatically_annotate_mass2motifs.download_mass2motifs import download_motif_set_from_motifdb
from automatically_annotate_mass2motifs.search_mass2motif import SelectSpectraContainingMass2Motif
from automatically_annotate_mass2motifs.clean_library_spectra import FilterLibrarySpectra
from automatically_annotate_mass2motifs.utils import convert_file_to_matchms_spectrum_objects
from automatically_annotate_mass2motifs.utils import store_pickled_file, load_pickled_file
from automatically_annotate_mass2motifs.main import get_cleaned_and_binned_spectra
from automatically_annotate_mass2motifs.visualize_results import load_moss_results


def select_smiles_for_mass2motifs(raw_library_spectra_file,
                                  motif_set_name,
                                  bin_size,
                                  moss_files_folder
                                  ):
    """Downloads mass2motif set and creates a moss input smile with matching spectra
    """
    similarity_scores_matrix = load_pickled_file("/lustre/BIF/nobackup/jonge094/ms2lda_reloaded/data/moss/gnps_library_derived_mass2motifs/mass2motif_scores_matrix_four_motifs_GNPS library derived Mass2Motifs.pickle")[1:]
    print(similarity_scores_matrix.shape)
    # download mass2motifs
    mass2motifs = download_motif_set_from_motifdb(motif_set_name, bin_size)[1:4]
    print(len(mass2motifs))

    # Load in and filter and bin library spectra
    binned_and_filtered_library_spectra = get_cleaned_and_binned_spectra(raw_library_spectra_file, bin_size)

    # search mass2motif
    spectra_selector = SelectSpectraContainingMass2Motif(binned_and_filtered_library_spectra,
                                                         mass2motifs,
                                                         scores_matrix=similarity_scores_matrix,
                                                         assert_correct_spectra=False)
    print( "Loaded everything")
    for minimal_mass2motif_matching_score in [0.05, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8, 0.9]:
        spectra_selector.create_all_moss_files(output_folder=moss_files_folder,
                                               minimal_score=minimal_mass2motif_matching_score)


def frequencies_for_different_thresholds(mass2motif_dir,
                                         mass2motif_name):
    result_dict = {}
    for file_name in sorted(os.listdir(mass2motif_dir)):
        if os.path.splitext(file_name)[1] == ".sub":
            expected_start = "mass2motif_" + mass2motif_name+ "_min_"
            if file_name.startswith(expected_start):
                threshold = float(os.path.splitext(file_name)[0][len(expected_start):])
                result = load_moss_results(os.path.join(mass2motif_dir, expected_start + str(threshold) + ".sub"))
                result_dict[threshold] = result[["description", "s_abs", "s_rel", "c_rel", "diff_support"]]
    smiles_highest_threshold = result_dict[max(result_dict.keys())]["description"]
    five_highest_smiles = list(smiles_highest_threshold[:5])
    smiles_dict = {}
    for highest_smile in five_highest_smiles:
        s_rel_list = []
        for threshold, result in result_dict.items():
            try:
                s_rel_list.append(float(result.loc[result["description"]==highest_smile]["s_rel"]))
            except:
                s_rel_list.append(0.0)
        smiles_dict[highest_smile] = s_rel_list

    thresholds = list(result_dict.keys())
    plot_dependency_threshold(thresholds, smiles_dict)
    print(smiles_dict)


def plot_dependency_threshold(thresholds,
                              smiles_dict):
    for smiles in smiles_dict:
        plt.plot(thresholds, smiles_dict[smiles], label=smiles)
        plt.legend(loc="upper left")
        plt.xlabel("Threshold similarity")
        plt.ylabel("Support (%)")

    plt.show()

if __name__ == "__main__":
    frequencies_for_different_thresholds("../data/moss/gnps_library_derived_mass2motifs/different_thresholds/moss_results/",
                                       "gnps_motif_68.m2m")
    select_smiles_for_mass2motifs(raw_library_spectra_file="../data/all_gnps.mgf",
                                  motif_set_name="GNPS library derived Mass2Motifs",
                                  bin_size = 0.005,
                                  moss_files_folder="../data/moss/gnps_library_derived_mass2motifs/different_thresholds",)