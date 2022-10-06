import os
from automatically_annotate_mass2motifs.moss_annotation import create_moss_input_file


def test_create_moss_input_file(tmp_path):
    file_name = os.path.join(tmp_path, 'moss_input_file.smiles')
    create_moss_input_file(smiles_matching_mass2motif=["CCC=O", "CCC"],
                           smiles_not_matching_mass2motif= ["CCCC", "CC"],
                           output_csv_file_location=file_name)
