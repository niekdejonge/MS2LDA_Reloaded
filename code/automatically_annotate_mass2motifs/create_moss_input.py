from typing import List, Dict
import pandas as pd
import numpy as np
import os
import pickle
import json
import matchms.metadata_utils


def create_smiles_csv_for_moss(input_tsv_file_location: str,
                               output_csv_file_location: str,
                               all_inchikey_json_file
                               ) -> None:
    """From a tsv file with smiles and spectrum_id created a file readable by moss
    
    :param input_tsv_file_location:
    :param output_csv_file_location:
    """
    all_inchikey_dict = json.load(open(all_inchikey_json_file, "r"))
    smiles_df = pd.read_csv(input_tsv_file_location, delimiter="\t")
    smiles = list(smiles_df["Smiles"])
    smiles_matching_mass2motif, smiles_not_matching_mass2motif = select_smiles_not_in_library(smiles,
                                 all_inchikey_dict)
    category = np.append(np.zeros(len(smiles_matching_mass2motif), dtype=np.int8), np.ones(len(smiles_not_matching_mass2motif), dtype=np.int8))
    ids = np.arange(0, len(category))
    smiles_df = pd.DataFrame({"ids": ids, "category": category, "smiles": smiles_matching_mass2motif+smiles_not_matching_mass2motif})
    smiles_df.to_csv(output_csv_file_location, index=False, header=False)


def select_smiles_not_in_library(smiles: List[str],
                                 all_inchikey_dict: Dict[str, str]):
    """

    :param smiles:
    :param all_inchikey_dict:
    :return:
    """
    matching_inchikey_dict = {}
    for smile in smiles:
        inchi = matchms.metadata_utils.convert_smiles_to_inchi(smile)
        inchikey = matchms.metadata_utils.convert_inchi_to_inchikey(inchi)[:14]
        if inchikey not in matching_inchikey_dict:
            matching_inchikey_dict[inchikey] = smile

    smiles_not_matching_mass2motif = []
    for inchikey in all_inchikey_dict:
        if inchikey not in matching_inchikey_dict:
            smiles_not_matching_mass2motif.append(all_inchikey_dict[inchikey])

    smiles_matching_mass2motif = list(matching_inchikey_dict.values())
    return smiles_matching_mass2motif, smiles_not_matching_mass2motif


def select_inchikeys_and_smiles(list_of_spectra, output_file_name):
    inchikey_dict = {}
    for spectrum in list_of_spectra:
        inchikey = spectrum.get("inchikey")[:14]
        if inchikey not in inchikey_dict:
            inchikey_dict[inchikey] = spectrum.get("smiles")
    json.dump(inchikey_dict, open(output_file_name, "w"))


def load_pickled_file(filename: str):
    with open(filename, 'rb') as file:
        loaded_object = pickle.load(file)
    return loaded_object


if __name__ == "__main__":
    moss_dir = "/lustre/BIF/nobackup/jonge094/ms2lda_reloaded/data/moss"
    # spectra = load_pickled_file("C:/Users/jonge094/PycharmProjects/PhD_MS2Query/ms2query/data/libraries_and_models/gnps_15_12_2021/in_between_files/ALL_GNPS_15_12_2021_positive_annotated.pickle")
    # select_inchikeys_and_smiles(spectra, "C:/Users/jonge094/PycharmProjects/ms2lda_reloaded/data/moss/inchikeys_and_smiles.json")
    create_smiles_csv_for_moss(os.path.join(moss_dir, "test_data_gnps_motif_45.m2_massql_result.txt"),
                               os.path.join(moss_dir, "test_data_gnps_motif_45.m2_massql_result.csv"),
                               os.path.join(moss_dir, "inchikeys_and_smiles.json"))

