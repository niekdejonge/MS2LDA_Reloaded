import os
from automatically_annotate_mass2motifs.moss_annotation import create_moss_input_file, run_moss


def test_create_moss_input_file(tmp_path):
    file_name = os.path.join(tmp_path, 'moss_input_file.smiles')
    create_moss_input_file(smiles_matching_mass2motif=["CCC=O", "CCC"],
                           smiles_not_matching_mass2motif= ["CCCC", "CC"],
                           output_csv_file_location=file_name)
    assert os.path.isfile(file_name), "expected file to be created"
    with open(file_name, "r") as file:
        result = file.readlines()
        assert result == ['0,0,CCC=O\n', '1,0,CCC\n', '2,1,CCCC\n', '3,1,CC\n'], \
            "expected a different output file"

def test_run_moss(tmp_path):
    moss_input_file = os.path.join(tmp_path, 'moss_input_file.smiles')
    moss_output_file = os.path.join(tmp_path, 'moss_output_file.sub')
    with open(moss_input_file, "w") as file:
        file.writelines(['0,0,CCC=O\n', '1,0,CCC\n', '2,1,CCCC\n', '3,1,CC\n'])
    run_moss(moss_input_file, moss_output_file, 0, 100)
    assert os.path.isfile(moss_output_file), "expected file to be created"

    with open(moss_output_file, "r") as output:
        result = output.readlines()
        assert result == ['id,description,nodes,edges,s_abs,s_rel,c_abs,c_rel\n',
                          '1,O=C-C-C,4,3,1,50.0,0,0.0\n',
                          '2,C(-C)-C,3,2,2,100.0,1,50.0\n']
