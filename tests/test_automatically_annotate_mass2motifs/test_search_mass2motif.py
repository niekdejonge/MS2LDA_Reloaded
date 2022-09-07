import os.path
from typing import List

import pandas as pd
from matchms import Fragments, Spectrum
import numpy as np
from automatically_annotate_mass2motifs.search_mass2motif import (overlap_in_fragments,
                                                                  SelectSpectraContainingMass2Motif,
                                                                  similarity_mass2motif_and_spectrum)
from automatically_annotate_mass2motifs.mass2motif import Mass2Motif
from tests.test_automatically_annotate_mass2motifs.generate_test_data import binned_spectra_005

def test_overlap_in_fragments():
    spectrum_fragment = Fragments(np.array([12.3, 14.5, 23.4, 26.7]), np.array([1.0, 0.2, 0.5, 0.3]))
    mass2motif_fragment = Fragments(np.array([12.3, 18.5, 26.7]), np.array([1.0, 0.2, 0.3]))

    result = overlap_in_fragments(mass2motif_fragment, spectrum_fragment)
    assert np.all(result.mz == np.array([12.3, 18.5, 26.7]))
    assert np.all(result.intensities == np.array([1.0, 0, 0.3]))

def test_similarity_mass2motif_and_spectrum():
    binned_spectrum = Spectrum(mz=np.array([100.0, 200.0], dtype="float"),
                         intensities=np.array([0.8, 0.6], dtype="float"),
                         metadata={})
    binned_spectrum.losses = Fragments(mz=np.array([200.0, 300.0], dtype="float"),
                                 intensities=np.array([0.4, 0.6], dtype="float"))
    mass2motif = Mass2Motif(['fragment_100.0', 'loss_200.0'],
                            [0.5, 0.25])
    score = similarity_mass2motif_and_spectrum(binned_spectrum, mass2motif)
    assert score == 0.5


def test_select_spectra_matching_mass2motif(tmp_path):
    spectra = binned_spectra_005()

    mass2motifs = [Mass2Motif(['fragment_100.025', 'loss_200.725', 'loss_128.075'],
                              [0.5, 0.25, 0.1], 0.05),
                   Mass2Motif(['fragment_50.025', 'loss_80.025', 'loss_128.075'],
                              [0.6, 0.5, 0.1], 0.05)
                   ]
    spectra_containing_mass2motif = SelectSpectraContainingMass2Motif(spectra, mass2motifs)
    expected_result = pd.DataFrame([[0.55, 0.05],
                                    [0.00, 0.60]])
    assert np.all(expected_result == spectra_containing_mass2motif.scores_matrix), "Expected different scores matrix"
    spectra_per_mass2motif = spectra_containing_mass2motif.select_spectra_matching_mass2motif(0.05)
    expected_spectra_per_mass2motif = [[spectra[0], spectra[1]], [spectra[1]]]
    assert np.all(spectra_per_mass2motif == expected_spectra_per_mass2motif), "Different spectra were expected"
    smiles_matching, smiles_not_matching = spectra_containing_mass2motif.select_non_matching_smiles(spectra_per_mass2motif[0])
    assert set(smiles_matching) == {'CN=C=O', 'C1CCCCC1'}
    assert smiles_not_matching == []
    smiles_matching, smiles_not_matching = spectra_containing_mass2motif.select_non_matching_smiles(spectra_per_mass2motif[1])
    assert smiles_matching == ["C1CCCCC1"]
    assert smiles_not_matching == ["CN=C=O"]
    spectra_containing_mass2motif.create_all_moss_files(tmp_path, 0.05)
    assert os.path.exists(os.path.join(tmp_path, "mass2motif_0.smiles")), "Moss smiles file was not created"
    assert os.path.exists(os.path.join(tmp_path, "mass2motif_1.smiles")), "Moss smiles file was not created"
    with open(os.path.join(tmp_path, "mass2motif_0.smiles"), "rb") as file_0:
        assert file_0.readlines() in ([b'1,0,CN=C=O\r\n', b'0,0,C1CCCCC1\r\n'], [b'0,0,C1CCCCC1\r\n', b'1,0,CN=C=O\r\n']), "Expected different output in moss file"
    with open(os.path.join(tmp_path, "mass2motif_1.smiles"), "rb") as file_0:
        assert file_0.readlines() == [b'0,0,C1CCCCC1\r\n', b'1,1,CN=C=O\r\n'], "Expected different output in moss file"


if __name__ == "__main__":
    pass