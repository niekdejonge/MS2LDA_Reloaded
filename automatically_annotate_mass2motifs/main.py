import pickle
from matchms.metadata_utils import is_valid_smiles

from automatically_annotate_mass2motifs.download_mass2motifs import get_single_motif_set, convert_json_to_mass2motif
from automatically_annotate_mass2motifs.search_mass2motif import SelectSpectraContainingMass2Motif
from automatically_annotate_mass2motifs.clean_library_spectra import FilterLibrarySpectra
from automatically_annotate_mass2motifs.utils import convert_file_to_matchms_spectrum_objects
from automatically_annotate_mass2motifs.utils import store_pickled_file, load_pickled_file
def select_smiles_for_mass2motifs(raw_library_spectra_file, motif_set_name, bin_size):
    # Clean and filter library spectra
    library_spectra = FilterLibrarySpectra(raw_library_spectra_file)
    library_spectra.save_pickled_spectra()
    # Bin library spectra to match with mass2motifs
    library_spectra.bin_spectra(bin_size)
    binned_library_file_name = library_spectra.save_pickled_spectra()
    # download mass2motifs
    json_mass_2motif_set = get_single_motif_set(motif_set_name)
    mass2motifs = convert_json_to_mass2motif(json_mass_2motif_set)
    # search mass2motif
    spectra_selector = SelectSpectraContainingMass2Motif(library_spectra.spectra, mass2motifs)
    selected_smiles = spectra_selector.select_smiles_mass2motif(0.1)
    print(spectra_selector.scores_matrix)
    return selected_smiles


if __name__ == "__main__":
    library_spectra_file_name = "../data/all_gnps_cleaned_filtered_binned_0.005.pickle"
    spectra = convert_file_to_matchms_spectrum_objects(library_spectra_file_name)
    # download mass2motifs
    json_mass_2motif_set = get_single_motif_set("Urine derived Mass2Motifs 2")
    mass2motifs = convert_json_to_mass2motif(json_mass_2motif_set)
    # search mass2motif
    spectra_selector = SelectSpectraContainingMass2Motif(spectra, mass2motifs)
    store_pickled_file(spectra_selector, "../data/select_spectra-containing_mass2motifs_1_motif")
    spectra_selector = load_pickled_file("../data/select_spectra-containing_mass2motifs_1_motif")
    selected_smiles = spectra_selector.select_smiles_mass2motif(0.1)
    print(len(selected_smiles[0]))
    selected_smiles = spectra_selector.select_smiles_mass2motif(0.2)
    print(len(selected_smiles[0]))
    selected_smiles = spectra_selector.select_smiles_mass2motif(0.3)
    print(len(selected_smiles[0]))
    selected_smiles = spectra_selector.select_smiles_mass2motif(0.4)
    print(len(selected_smiles[0]))
    selected_smiles = spectra_selector.select_smiles_mass2motif(0.5)
    print(len(selected_smiles[0]))
    selected_smiles = spectra_selector.select_smiles_mass2motif(0.6)
    print(len(selected_smiles[0]))