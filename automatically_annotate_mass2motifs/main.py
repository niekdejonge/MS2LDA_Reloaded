import os
import pickle
from matchms.metadata_utils import is_valid_smiles

from automatically_annotate_mass2motifs.download_mass2motifs import download_motif_set_from_motifdb
from automatically_annotate_mass2motifs.search_mass2motif import SelectSpectraContainingMass2Motif
from automatically_annotate_mass2motifs.clean_library_spectra import FilterLibrarySpectra
from automatically_annotate_mass2motifs.utils import convert_file_to_matchms_spectrum_objects
from automatically_annotate_mass2motifs.utils import store_pickled_file, load_pickled_file

def get_cleaned_and_binned_spectra(raw_library_spectra_file, bin_size):
    filtered_library_file = os.path.splitext(raw_library_spectra_file)[0] + "_positive_ionmode_cleaned_filtered.pickle"
    filtered_and_binned_library_file = os.path.splitext(raw_library_spectra_file)[
                                           0] + f"_positive_ionmode_cleaned_filtered_binned_{bin_size}.pickle"
    # Check if there are already files containing cleaned or filtered spectra
    if os.path.isfile(filtered_and_binned_library_file):
        binned_and_filtered_library_spectra = load_pickled_file(filtered_and_binned_library_file)
        print(f'Using the file {filtered_and_binned_library_file} ')
    elif os.path.isfile(filtered_library_file):
        library_spectra = FilterLibrarySpectra(filtered_library_file, already_cleaned=True)
        library_spectra.bin_spectra(bin_size)
        library_spectra.save_pickled_spectra()
        binned_and_filtered_library_spectra = library_spectra.spectra
        print(f'Using the file {filtered_library_file} ')
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
                                  moss_files_folder,
                                  minimal_mass2motif_matching_score):
    """Downloads mass2motif set and creates a moss input smile with matching spectra
    """
    # download mass2motifs
    mass2motifs = download_motif_set_from_motifdb(motif_set_name, bin_size)
    print(len(mass2motifs))

    # Load in and filter and bin library spectra
    binned_and_filtered_library_spectra = get_cleaned_and_binned_spectra(raw_library_spectra_file, bin_size)

    # search mass2motif
    spectra_selector = SelectSpectraContainingMass2Motif(binned_and_filtered_library_spectra, mass2motifs)
    spectra_selector.create_all_moss_files(output_folder=moss_files_folder,
                                           minimal_score=minimal_mass2motif_matching_score)

if __name__ == "__main__":
    select_smiles_for_mass2motifs(raw_library_spectra_file="../data/all_gnps.mgf",
                                  motif_set_name="GNPS library derived Mass2Motifs",
                                  bin_size = 0.005,
                                  moss_files_folder="../data/moss/moss/gnps_library_derived_mass2motifs",
                                  minimal_mass2motif_matching_score=0.1)
