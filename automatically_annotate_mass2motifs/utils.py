import json

import os
from typing import List
import pickle
from tqdm import tqdm
from matchms import importing
from matchms.Spectrum import Spectrum
from matchms.logging_functions import set_matchms_logger_level


def load_pickled_file(file_name):
    with open(file_name, "rb") as file:
        loaded_object = pickle.load(file)
    return loaded_object

def store_pickled_file(object_to_store,
                       file_name):
    return_non_existing_file_name(file_name)
    with open(file_name, "wb") as file:
        pickle.dump(object_to_store, file)

def return_non_existing_file_name(file_name):
    """Checks if a path already exists, otherwise creates a new filename with (1)"""
    if not os.path.exists(file_name):
        return file_name
    print(f"The file name already exists: {file_name}")
    file_name_base, file_extension = os.path.splitext(file_name)
    i = 1
    new_file_name = f"{file_name_base}({i}){file_extension}"
    while os.path.exists(new_file_name):
        i += 1
        new_file_name = f"{file_name_base}({i}){file_extension}"
    print(f"Instead the file will be stored in {new_file_name}")
    return new_file_name

def convert_file_to_matchms_spectrum_objects(file_name,
                                             filter_metadata: bool = False
                                             ) -> List[Spectrum]:
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
    # Setting the logger to error speeds up the loading process
    set_matchms_logger_level("ERROR")
    file_extension = os.path.splitext(file_name)[1].lower()
    if file_extension == ".mzml":
        return list(tqdm(importing.load_from_mzml(file_name, metadata_harmonization=filter_metadata),
                         desc="Loading in spectra"))
    if file_extension == ".json":
        return list(tqdm(importing.load_from_json(file_name, metadata_harmonization=filter_metadata),
                         desc="Loading in spectra"))
    if file_extension == ".mgf":
        return list(tqdm(importing.load_from_mgf(file_name, metadata_harmonization=filter_metadata),
                         desc="Loading in spectra"))
    if file_extension == ".msp":
        return list(tqdm(importing.load_from_msp(file_name, metadata_harmonization=filter_metadata),
                         desc="Loading in spectra"))
    if file_extension == ".mzxml":
        return list(tqdm(importing.load_from_mzxml(file_name, metadata_harmonization=filter_metadata),
                         desc="Loading in spectra"))
    if file_extension == ".usi":
        return list(tqdm(importing.load_from_usi(file_name, metadata_harmonization=filter_metadata),
                         desc="Loading in spectra"))
    if file_extension == ".pickle":
        spectra = load_pickled_file(file_name)
        assert isinstance(spectra, list), "Expected list of spectra"
        assert isinstance(spectra[0], Spectrum), "Expected list of spectra"
        return spectra
    assert False, f"File extension of file: {file_name} is not recognized"


def store_as_json(objects_to_store,
                  file_name):
    """Stores mass2motifs or annotations as json.

    The object needs the method .to_dict()"""
    if not isinstance(objects_to_store, list):
        objects_to_store = [objects_to_store]
    json_str = []
    for object_to_store in objects_to_store:
        json_str.append(object_to_store.to_dict())
    with open(file_name, "w", encoding="utf-8") as file:
        json.dump(json_str, file, indent=3)
