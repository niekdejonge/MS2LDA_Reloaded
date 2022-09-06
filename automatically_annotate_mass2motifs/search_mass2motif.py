from typing import List
import pandas as pd
from tqdm import tqdm
from matchms import Spectrum, Fragments
from matchms.metadata_utils import is_valid_smiles
from matchms.filtering import normalize_intensities
import numpy as np
from automatically_annotate_mass2motifs.mass2motif import Mass2Motif
from automatically_annotate_mass2motifs.bin_spectra import bin_spectra


class SelectSpectraContainingMass2Motif:
    def __init__(self,
                 binned_spectra: List[Spectrum],
                 mass2motifs: List[Mass2Motif]):
        self.mass2motifs = mass2motifs
        self.bin_size = mass2motifs[0].bin_size
        assert np.all([mass2motif.bin_size == self.bin_size for mass2motif in mass2motifs]), \
            "Expected mass2motifs with the same binning method"

        self.spectra = binned_spectra
        # self.assert_correct_spectra()
        self.smiles_spectra = [spectrum.get("smiles") for spectrum in self.spectra]
        self.scores_matrix = self._similarity_matrix_mass2motifs_and_spectra()

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


    def _similarity_matrix_mass2motifs_and_spectra(self):
        scores_matrix = []
        for mass2motif in tqdm(self.mass2motifs, "Calculating scores for mass2motif"):
            scores_mass2motif = []
            for spectrum in self.spectra:
                scores_mass2motif.append(similarity_mass2motif_and_spectrum(spectrum, mass2motif))
            scores_matrix.append(scores_mass2motif)
        return pd.DataFrame(scores_matrix, columns=self.smiles_spectra)

    def select_smiles_mass2motif(self, minimal_score) -> pd.Series:
        """Returns the smiles of the spectra per Mass2Motif that have a high enough similarity score"""
        result = self.scores_matrix.apply(lambda row: row[row >= minimal_score].index.tolist(), axis=1)
        return result


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


if __name__ == "__main__":
    from clean_library_spectra import convert_file_to_matchms_spectrum_objects
    binned_library_spectra = convert_file_to_matchms_spectrum_objects("../data/all_gnps_cleaned_filtered_binned_0.005.pickle")
    # SelectSpectraContainingMass2Motif(binned_library_spectra, mass2motifs)

