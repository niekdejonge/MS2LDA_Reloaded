from typing import List
import pandas as pd
from tqdm import tqdm
from matchms import Spectrum, Fragments
from matchms.metadata_utils import is_valid_smiles
import numpy as np
from automatically_annotate_mass2motifs.mass2motif import Mass2Motif


class ScoresMatrix:
    """Stores a similarity matrix between mass2motifs and spectra"""
    def __init__(self,
                 binned_spectra: List[Spectrum],
                 scores_matrix: pd.DataFrame):
        assert scores_matrix.shape[1] == len(binned_spectra), "The wrong number of spectra is given"
        assert [spectrum.get('inchikey') for spectrum in binned_spectra] == list(scores_matrix.columns), \
            "The inchikeys of the spectra given do not match the (order of the) scores matrix"
        self.spectra = binned_spectra
        self.scores_matrix = scores_matrix

    def select_matching_spectra(self,
                                minimal_score: float,
                                mass2motif: Mass2Motif):
        row_corresponding_to_mass2motif = self.scores_matrix.loc[f"{mass2motif.motif_set_name} {mass2motif.motif_name}", :]
        matching_spectra_booleans = row_corresponding_to_mass2motif >= minimal_score
        matching_spectrum_indexes = [i for i, x in enumerate(matching_spectra_booleans) if x]
        not_matching_spectrum_indexes = [i for i, x in enumerate(matching_spectra_booleans) if not x]
        matching_spectra = [self.spectra[spectrum_index] for spectrum_index in matching_spectrum_indexes]
        not_matching_spectra = [self.spectra[spectrum_index] for spectrum_index in not_matching_spectrum_indexes]
        return matching_spectra, not_matching_spectra


def create_similarity_matrix(mass2motifs: List[Mass2Motif],
                             binned_spectra: List[Spectrum]) -> pd.DataFrame:
    """Creates a similarity matrix between mass2motifs and spectra

    For the index the motif_set_name and motif_name are combined and for the column names the inchikeys are used"""
    # Check on the first spectrum if all filtering steps were applied
    assert_correct_spectrum(binned_spectra[0])
    bin_size = mass2motifs[0].bin_size
    assert np.all([mass2motif.bin_size == bin_size for mass2motif in mass2motifs]), \
        f"Expected mass2motifs with the same binning method found {[mass2motif.bin_size for mass2motif in mass2motifs]}"
    # todo also check if it is the same binning as the spectra
    # Create the indexes of the mass2motifs
    mass2motif_indexes = [f"{mass2motif.motif_set_name} {mass2motif.motif_name}" for mass2motif in mass2motifs]
    assert len(mass2motif_indexes) == len(set(mass2motif_indexes)), "One of the Motifs has the same motif set name and motif name"
    spectra_indexes = [spectrum.get("inchikey") for spectrum in binned_spectra]
    scores_matrix = []

    for mass2motif in tqdm(mass2motifs,
                           "Calculating scores for mass2motif"):
        scores_mass2motif = []
        for spectrum in binned_spectra:
            scores_mass2motif.append(similarity_mass2motif_and_spectrum(spectrum, mass2motif))
        scores_matrix.append(scores_mass2motif)
    return pd.DataFrame(scores_matrix, index=mass2motif_indexes, columns=spectra_indexes)


def assert_correct_spectrum(spectrum: Spectrum):
    """Checks if the spectra are normalized have losses and if the spectra have valid smiles"""
    assert isinstance(spectrum, Spectrum), "Expected a list of matchms spectrum objects"
    assert is_valid_smiles(spectrum.get("smiles")), f"The smile {spectrum.get('smiles')} is not a valid smile"
    assert len(spectrum.losses.mz) > 0, "Losses have not yet been added"
    assert max(spectrum.peaks.intensities) == 1, "Spectra peaks were not normalized"
    assert max(spectrum.losses.intensities) == 1, "Spectra losses were not normalized"


def overlap_in_fragments(mass2motif_fragment: Fragments,
                         spectrum_fragment: Fragments):
    """Returns the intensities in spectrum_fragments for the mz values in mass2motif_fragments"""
    overlapping_intensities_spectrum = np.empty(mass2motif_fragment.mz.shape)
    for i_mass2motif, mz in enumerate(mass2motif_fragment.mz):
        if mz in spectrum_fragment.mz:
            i_spectrum = np.where(spectrum_fragment.mz == mz)
            spectrum_intensity = spectrum_fragment.intensities[i_spectrum]
            overlapping_intensities_spectrum[i_mass2motif] = spectrum_intensity
        else:
            overlapping_intensities_spectrum[i_mass2motif] = 0.0
    return Fragments(mass2motif_fragment.mz, overlapping_intensities_spectrum)


def similarity_mass2motif_and_spectrum(binned_spectrum: Spectrum,
                                       mass2motif: Mass2Motif):
    """Calculates the similarity between a spectrum and a mass2motif

       A dot product is calculated between the intensities for the peaks
       in the mass2motifs times the probabilities of the mass2motif"""
    spectrum_fragments = overlap_in_fragments(mass2motif.fragments, binned_spectrum.peaks)
    spectrum_losses = overlap_in_fragments(mass2motif.losses, binned_spectrum.losses)
    spectrum_intensities = np.concatenate((spectrum_fragments.intensities, spectrum_losses.intensities))
    mass2motif_probabilities = np.concatenate((mass2motif.fragments.intensities, mass2motif.losses.intensities))
    return np.dot(mass2motif_probabilities, spectrum_intensities)
