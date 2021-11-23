import os
from tqdm import tqdm
from typing import Union, List
from matchms import importing
from matchms.Spectrum import Spectrum


def convert_file_to_matchms_spectrum_objects(file_name
                                             ) -> Union[List[Spectrum], None]:
    """Loads spectra from your spectrum file into memory as matchms Spectrum object

    The following file extensions can be loaded in with this function:
    "mzML", "json", "mgf", "msp", "mzxml" and "usi".
    A pickled file is expected to directly contain a list of matchms spectrum objects.

    Args:
    -----
    file_name:
        Path to file containing spectra, with file extension "mzML", "json", "mgf", "msp",
        "mzxml", "usi" or "pickle"
    """
    assert os.path.exists(file_name), f"The specified file: {file_name} does not exists"

    file_extension = os.path.splitext(file_name)[1].lower()
    if file_extension == ".mzml":
        return list(tqdm(importing.load_from_mzml(file_name),
                    desc="Loading in spectra"))
    if file_extension == ".json":
        return list(tqdm(importing.load_from_json(file_name),
                    desc="Loading in spectra"))
    if file_extension == ".mgf":
        return list(tqdm(importing.load_from_mgf(file_name),
                    desc="Loading in spectra"))
    if file_extension == ".msp":
        return list(tqdm(importing.load_from_msp(file_name),
                    desc="Loading in spectra"))
    if file_extension == ".mzxml":
        return list(tqdm(importing.load_from_mzxml(file_name),
                    desc="Loading in spectra"))
    if file_extension == ".usi":
        return list(tqdm(importing.load_from_usi(file_name),
                    desc="Loading in spectra"))
    print(f"File extension of file: {file_name} is not recognized")
    return None


def add_unknown_charges_to_spectra(spectrum_list: List[Spectrum],
                                   charge_to_use: int = 1) -> List[Spectrum]:
    """Adds charges to spectra when no charge is known

    This is important for matchms to be able to calculate the parent_mass
    from the mz_precursor

    Args:
    ------
    spectrum_list:
        List of spectra
    charge_to_use:
        The charge set when no charge is known. Default = 1
    """
    for spectrum in spectrum_list:
        if spectrum.get("charge") is None:
            spectrum.set("charge", charge_to_use)
    return spectrum_list