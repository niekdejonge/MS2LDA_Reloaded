from typing import List

import pandas as pd
from matchms import Fragments, Spectrum
import numpy as np
from automatically_annotate_mass2motifs.search_mass2motif import (overlap_in_fragments,
                                                                  SelectSpectraContainingMass2Motif,
                                                                  similarity_mass2motif_and_spectrum)
from automatically_annotate_mass2motifs.mass2motif import Mass2Motif
from tests.test_automatically_annotate_mass2motifs.generate_test_data import spectra_with_losses

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


def test_SelectSpectraContainingMass2Motif():
    spectra = spectra_with_losses()

    mass2motifs = [Mass2Motif(['fragment_100.025', 'loss_200.725', 'loss_128.075'],
                              [0.5, 0.25, 0.1], 0.05),
                   Mass2Motif(['fragment_50.025', 'loss_80.025', 'loss_128.075'],
                              [0.6, 0.5, 0.1], 0.05)
                   ]
    spectra_containing_mass2motif = SelectSpectraContainingMass2Motif(spectra, mass2motifs)
    expected_result = pd.DataFrame([[0.55, 0.05],
                                    [0.00, 0.60]], columns = ["CN=C=O", "C1CCCCC1"])
    assert np.all(expected_result == spectra_containing_mass2motif.scores_matrix)
    smiles_per_mass2motif = spectra_containing_mass2motif.select_smiles_mass2motif(0.05)
    expected_smiles = pd.Series([['CN=C=O', 'C1CCCCC1'], ['C1CCCCC1']])
    assert np.all(smiles_per_mass2motif == expected_smiles), "Different smiles were expected"

if __name__ == "__main__":
    pass