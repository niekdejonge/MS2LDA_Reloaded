import pandas as pd
from matchms import Fragments, Spectrum
import numpy as np
from automatically_annotate_mass2motifs.scores_matrix import (overlap_in_fragments,
                                                              ScoresMatrix,
                                                              similarity_mass2motif_and_spectrum,
                                                              create_similarity_matrix)
from automatically_annotate_mass2motifs.mass2motif import Mass2Motif
from tests.test_automatically_annotate_mass2motifs.generate_test_data import binned_spectra_005


def test_overlap_in_fragments():
    spectrum_fragment = Fragments(np.array([12.3, 14.5, 23.4, 26.7]), np.array([1.0, 0.2, 0.5, 0.3]))
    mass2motif_fragment = Fragments(np.array([12.3, 18.5, 26.7]), np.array([1.0, 0.2, 0.3]))

    result = overlap_in_fragments(mass2motif_fragment, spectrum_fragment)
    assert np.all(result.mz == np.array([12.3, 18.5, 26.7]))
    assert np.all(result.intensities == np.array([1.0, 0, 0.3]))


def test_similarity_mass2motif_and_spectrum():
    binned_spectrum = Spectrum(mz=np.array([100.0025, 200.0025], dtype="float"),
                         intensities=np.array([0.8, 0.6], dtype="float"),
                         metadata={})
    binned_spectrum.losses = Fragments(mz=np.array([200.0025, 300.0025], dtype="float"),
                                 intensities=np.array([0.4, 0.6], dtype="float"))
    mass2motif = Mass2Motif(fragments=[100.0025],
                            fragment_probabilities= [0.5],
                            losses=[200.0025],
                            loss_probabilities=[0.25],
                            bin_size=0.005)
    score = similarity_mass2motif_and_spectrum(binned_spectrum, mass2motif)
    assert score == 0.5


def test_generate_similarity_matrix():
    spectra = binned_spectra_005()
    mass2motifs = [Mass2Motif(fragments=[100.025], fragment_probabilities= [0.5],
                            losses=[128.075, 200.725], loss_probabilities=[0.1, 0.25],
                              motif_name="motif_1", motif_set_name="motif_set_1", bin_size=0.05),
                   Mass2Motif(fragments=[50.025], fragment_probabilities=[0.6],
                              losses=[80.025, 128.075], loss_probabilities=[0.5, 0.1],
                              motif_name="motif_2", motif_set_name="motif_set_1",
                              bin_size=0.05),
                   Mass2Motif(fragments=[50.025], fragment_probabilities=[0.6],
                              losses=[80.025, 128.075], loss_probabilities=[0.5, 0.1],
                              motif_name="motif_3", motif_set_name="motif_set_1",
                              bin_size=0.05)
                   ]
    similarity_matrix = create_similarity_matrix(mass2motifs, spectra)
    expected_result = pd.DataFrame([[0.55, 0.05],[0.00, 0.60],[0.00, 0.60]],
                                   index=["motif_set_1 motif_1", "motif_set_1 motif_2", "motif_set_1 motif_3"],
                                   columns=["ABCDEFGHIJKLAQ", "BBBBBBBBBBBBBB"])
    pd.testing.assert_frame_equal(expected_result, similarity_matrix)


def test_select_matching_spectra():
    spectra = binned_spectra_005()
    mass2motifs = [Mass2Motif(fragments=[100.025], fragment_probabilities= [0.5],
                            losses=[128.075, 200.725], loss_probabilities=[0.1, 0.25],
                              motif_name="motif_1", motif_set_name="motif_set_1", bin_size=0.05),
                   Mass2Motif(fragments=[50.025], fragment_probabilities=[0.6],
                              losses=[80.025, 128.075], loss_probabilities=[0.5, 0.1],
                              motif_name="motif_2", motif_set_name="motif_set_1",
                              bin_size=0.05),
                   Mass2Motif(fragments=[50.025], fragment_probabilities=[0.6],
                              losses=[80.025, 128.075], loss_probabilities=[0.5, 0.1],
                              motif_name="motif_3", motif_set_name="motif_set_1",
                              bin_size=0.05)
                   ]
    similarity_matrix = pd.DataFrame([[0.55, 0.05],[0.00, 0.60],[0.00, 0.60]],
                                     index=["motif_set_1 motif_1", "motif_set_1 motif_2", "motif_set_1 motif_3"],
                                     columns=["ABCDEFGHIJKLAQ", "BBBBBBBBBBBBBB"])
    scores_matrix = ScoresMatrix(spectra, similarity_matrix)
    result = scores_matrix.select_matching_spectra(0.6, mass2motifs[0])
    assert result[0] == [] and result[1] == spectra
    result = scores_matrix.select_matching_spectra(0.6, mass2motifs[1])
    assert result[0] == [spectra[1],] and result[1] == [spectra[0],]
    result = scores_matrix.select_matching_spectra(0, mass2motifs[0])
    assert result[0] == spectra and result[1] == []

if __name__ == "__main__":
    pass