from typing import List, Optional
import pandas as pd
from tqdm import tqdm
from matchms import Spectrum, Fragments
from matchms.metadata_utils import is_valid_smiles
import numpy as np
from automatically_annotate_mass2motifs.mass2motif import Mass2Motif


class SelectSpectraContainingMass2Motif:
    """Creates a similarity matrix between spectra and Mass2motifs"""
    def __init__(self,
                 binned_spectra: List[Spectrum],
                 mass2motifs: List[Mass2Motif],
                 assert_correct_spectra: bool = True,
                 scores_matrix: Optional[pd.DataFrame] = None):
        self.mass2motifs = mass2motifs
        self.bin_size = \
            mass2motifs[0].bin_size
        assert np.all([mass2motif.bin_size == self.bin_size for mass2motif in mass2motifs]), \
            f"Expected mass2motifs with the same binning method found {[mass2motif.bin_size for mass2motif in mass2motifs]}"

        self.spectra = binned_spectra
        if assert_correct_spectra:
            self.assert_correct_spectra()
        if scores_matrix is not None:
            assert scores_matrix.shape == (len(self.mass2motifs), len(self.spectra)), \
                "A scores matrix with the wrong dimensions was given"
            self.scores_matrix = scores_matrix
        else:
            self.similarity_matrix_mass2motifs_and_spectra()
        self.inchikey_smiles_dict = self.create_inchikeys_smiles_dict()

    def assert_correct_spectra(self):
        """Checks if the spectra have the expected format"""
        assert isinstance(self.spectra, list), "Expected a list of spectra"
        for spectrum in self.spectra:
            assert isinstance(spectrum, Spectrum), "Expected a list of matchms spectrum objects"
            assert is_valid_smiles(spectrum.get("smiles")), f"The smile {spectrum.get('smiles')} is not a valid smile"
        # Checks for first spectrum (to not cost too much time)
        checked_spectrum = self.spectra[0]
        assert len(checked_spectrum.losses.mz) > 0, "Losses have not yet been added"
        assert max(checked_spectrum.peaks.intensities) == 1, "Spectra peaks were not normalized"
        assert max(checked_spectrum.losses.intensities) == 1, "Spectra losses were not normalized"

    def create_inchikeys_smiles_dict(self):
        """Creates a dictionary with all inchikeys with one of the smiles that represents this inchikey"""
        inchikey_dict = {}
        for spectrum in self.spectra:
            inchikey = spectrum.get("inchikey")[:14]
            if inchikey not in inchikey_dict:
                inchikey_dict[inchikey] = spectrum.get("smiles")
        return inchikey_dict

    def similarity_matrix_mass2motifs_and_spectra(self):
        scores_matrix = []
        for mass2motif in tqdm(self.mass2motifs, "Calculating scores for mass2motif"):
            scores_mass2motif = []
            for spectrum in self.spectra:
                scores_mass2motif.append(similarity_mass2motif_and_spectrum(spectrum, mass2motif))
            scores_matrix.append(scores_mass2motif)
        self.scores_matrix = pd.DataFrame(scores_matrix)

    def select_matching_spectra(self, minimal_score):
        spectrum_indexes_per_mass2motif = self.scores_matrix.apply(lambda row: row[row >= minimal_score].index.tolist(),
                                                                   axis=1)
        spectra_per_mass2motif = []
        for spectrum_indexes in spectrum_indexes_per_mass2motif:
            spectra_per_mass2motif.append([self.spectra[i] for i in spectrum_indexes])
        return spectra_per_mass2motif

    def select_unique_matching_and_non_matching_smiles(self,
                                   matching_spectra: List[Spectrum]):
        """Selects all inchikeys not in smiles, and returns 1 smile for each

        This output is needed for creating moss files"""
        matching_inchikeys = []
        for spectrum in matching_spectra:
            inchikey = spectrum.get("inchikey")[:14]
            matching_inchikeys.append(inchikey)

        matching_inchikeys = set(matching_inchikeys)
        all_inchikeys = set(self.inchikey_smiles_dict.keys())
        not_matching_inchikeys = [inchikey for inchikey in all_inchikeys if inchikey not in matching_inchikeys]

        matching_smiles = [self.inchikey_smiles_dict[inchikey] for inchikey in matching_inchikeys]
        not_matching_smiles = [self.inchikey_smiles_dict[inchikey] for inchikey in not_matching_inchikeys]
        return matching_smiles, not_matching_smiles

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
