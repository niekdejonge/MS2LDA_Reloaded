import os
from tqdm import tqdm
from typing import Optional
from matchms import exporting
import matchms.filtering as msfilters
from matchms.logging_functions import set_matchms_logger_level
from matchms.metadata_utils import is_valid_smiles
from automatically_annotate_mass2motifs.bin_spectra import bin_spectra
from automatically_annotate_mass2motifs.utils import (store_pickled_file,
                                                      return_non_existing_file_name,
                                                      convert_file_to_matchms_spectrum_objects)


class FilterLibrarySpectra:
    """Filters spectra and can bin spectra for comparison to mass2motifs

    """
    def __init__(self,
                 spectra_file_name,
                 ion_mode: str = "positive",
                 already_cleaned = False):
        assert ion_mode in ("negative", "positive"), "ion_mode should be 'postive' or 'negative'"
        self.ion_mode = ion_mode
        self.spectra_file_name = spectra_file_name
        # store the filtering steps applied for naming stored files
        self.processing_log = ""
        set_matchms_logger_level("ERROR")
        self.spectra = convert_file_to_matchms_spectrum_objects(self.spectra_file_name, filter_metadata=True)

        if not already_cleaned:
            # clean and select spectra
            self.select_spectra_in_ion_mode()
            self.apply_metadata_filters()
            self.remove_spectra_missing_smiles_or_inchi()


    def apply_metadata_filters(self):
        cleaned_spectra = []
        for spectrum in tqdm(self.spectra, desc="Filting metadata of spectra"):
            # Default filters
            spectrum = msfilters.default_filters(spectrum)

            # Here, undefiend entries will be harmonized (instead of having a huge variation of None,"", "N/A" etc.)
            spectrum = msfilters.harmonize_undefined_inchikey(spectrum)
            spectrum = msfilters.harmonize_undefined_inchi(spectrum)
            spectrum = msfilters.harmonize_undefined_smiles(spectrum)

            # The repair_inchi_inchikey_smiles function will correct misplaced metadata (e.g. inchikeys entered as inchi etc.) and harmonize the entry strings.
            spectrum = msfilters.repair_inchi_inchikey_smiles(spectrum)

            # Where possible (and necessary, i.e. missing): Convert between smiles, inchi, inchikey to complete metadata.
            # This is done using functions from rdkit.
            spectrum = msfilters.repair_inchi_inchikey_smiles(spectrum)
            spectrum = msfilters.derive_inchi_from_smiles(spectrum)
            spectrum = msfilters.derive_smiles_from_inchi(spectrum)
            spectrum = msfilters.derive_inchikey_from_inchi(spectrum)
            cleaned_spectra.append(spectrum)
        self.spectra = cleaned_spectra
        self.processing_log += "_cleaned"


    def remove_spectra_missing_smiles_or_inchi(self):
        fully_annotated_spectra = []
        for spectrum in tqdm(self.spectra, desc="Removing spectra missing smiles or inchikeys"):
            inchikey = spectrum.get("inchikey")
            if inchikey is not None and len(inchikey) > 13:
                smiles = spectrum.get("smiles")
                inchi = spectrum.get("inchi")
                if smiles is not None and len(smiles) > 0 and is_valid_smiles(smiles):
                    if inchi is not None and len(inchi) > 0:
                        fully_annotated_spectra.append(spectrum)
        self.spectra = fully_annotated_spectra
        self.processing_log += "_filtered"

    def select_spectra_in_ion_mode(self):
        spectra_in_correct_ion_mode = []
        for spectrum in tqdm(self.spectra, desc=f"Selecting spectra in {self.ion_mode} ion mode"):
            if spectrum.get("ionmode") == self.ion_mode:
                spectra_in_correct_ion_mode.append(spectrum)
        self.spectra = spectra_in_correct_ion_mode
        self.processing_log += f"_{self.ion_mode}_ionmode"

    def bin_spectra(self, bin_width):
        # bin spectra
        self.spectra = bin_spectra(self.spectra, bin_width)
        # Normalize binned spectra
        self.spectra = [msfilters.normalize_intensities(spectrum) for spectrum in tqdm(self.spectra,
                                                                                       desc="Normalizing intensities of binned spectra")]
        self.processing_log += f"_binned_{bin_width}"

    def save_spectra_as_mgf(self, out_file_name: Optional[str] = None):
        if out_file_name is None:
            file_name_base = os.path.splitext(self.spectra_file_name)[0]
            out_file_name = file_name_base + self.processing_log + ".mgf"
        out_file_name = return_non_existing_file_name(out_file_name)
        exporting.save_as_mgf(self.spectra, out_file_name)
        return out_file_name

    def save_pickled_spectra(self, out_file_name: Optional[str] = None):
        if out_file_name is None:
            file_name_base = os.path.splitext(self.spectra_file_name)[0]
            out_file_name = file_name_base + self.processing_log + ".pickle"
        store_pickled_file(self.spectra, out_file_name)
        return out_file_name


if __name__ == "__main__":
    library_spectra = FilterLibrarySpectra("../data/all_gnps_cleaned_and_filtered.pickle", already_cleaned=True)
    library_spectra.save_pickled_spectra()
    library_spectra.bin_spectra(0.005)
    library_spectra.save_pickled_spectra()
