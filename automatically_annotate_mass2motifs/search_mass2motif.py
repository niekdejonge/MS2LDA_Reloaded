import os
from typing import List
import pandas as pd
from tqdm import tqdm
from matchms import Spectrum, Fragments
from matchms.metadata_utils import is_valid_smiles
import numpy as np
from automatically_annotate_mass2motifs.mass2motif import Mass2Motif
from automatically_annotate_mass2motifs.utils import return_non_existing_file_name


class SelectSpectraContainingMass2Motif:
    def __init__(self,
                 binned_spectra: List[Spectrum],
                 mass2motifs: List[Mass2Motif],
                 assert_correct_spectra: bool = True):
        self.mass2motifs = mass2motifs
        self.bin_size = mass2motifs[0].bin_size
        assert np.all([mass2motif.bin_size == self.bin_size for mass2motif in mass2motifs]), \
            "Expected mass2motifs with the same binning method"

        self.spectra = binned_spectra
        if assert_correct_spectra:
            self.assert_correct_spectra()
        self.scores_matrix = self._similarity_matrix_mass2motifs_and_spectra()
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


    def _similarity_matrix_mass2motifs_and_spectra(self):
        scores_matrix = []
        for mass2motif in tqdm(self.mass2motifs, "Calculating scores for mass2motif"):
            scores_mass2motif = []
            for spectrum in self.spectra:
                scores_mass2motif.append(similarity_mass2motif_and_spectrum(spectrum, mass2motif))
            scores_matrix.append(scores_mass2motif)
        return pd.DataFrame(scores_matrix)

    def select_spectra_matching_mass2motif(self, minimal_score) -> List[List[Spectrum]]:
        """Returns the smiles of the spectra per Mass2Motif that have a high enough similarity score"""
        spectrum_indexes_per_mass2motif = self.scores_matrix.apply(lambda row: row[row >= minimal_score].index.tolist(), axis=1)
        spectra_per_mass2motif = []
        for spectrum_indexes in spectrum_indexes_per_mass2motif:
            spectra_per_mass2motif.append([self.spectra[i] for i in spectrum_indexes])
        return spectra_per_mass2motif

    def create_inchikeys_smiles_dict(self):
        """Creates a dictionary with all inchikeys with one of the smiles that represents this inchikey"""
        inchikey_dict = {}
        for spectrum in self.spectra:
            inchikey = spectrum.get("inchikey")[:14]
            if inchikey not in inchikey_dict:
                inchikey_dict[inchikey] = spectrum.get("smiles")
        return inchikey_dict

    def select_non_matching_smiles(self, spectra: List[Spectrum]):
        """Selects all inchikeys not in smiles, and returns 1 smile for each

        This is done to prevent overfitting to inchikeys containing many smiles"""
        matching_inchikeys = []
        for spectrum in spectra:
            inchikey = spectrum.get("inchikey")[:14]
            matching_inchikeys.append(inchikey)

        matching_inchikeys = set(matching_inchikeys)
        all_inchikeys = set(self.inchikey_smiles_dict.keys())
        not_matching_inchikeys = [inchikey for inchikey in all_inchikeys if inchikey not in matching_inchikeys]

        matching_smiles = [self.inchikey_smiles_dict[inchikey] for inchikey in matching_inchikeys]
        not_matching_smiles = [self.inchikey_smiles_dict[inchikey] for inchikey in not_matching_inchikeys]
        return matching_smiles, not_matching_smiles

    def create_moss_file(self, list_of_matching_spectra: List[Spectrum], output_csv_file_location):
        """Creates a file containing the smiles matching and not matching a mass2motif readable by MOSS

        :param list_of_smiles:
        :param output_csv_file_location:
        :return:
        """
        output_csv_file_location = return_non_existing_file_name(output_csv_file_location)
        smiles_matching_mass2motif, smiles_not_matching_mass2motif = self.select_non_matching_smiles(list_of_matching_spectra)
        category = np.append(np.zeros(len(smiles_matching_mass2motif), dtype=np.int8),
                             np.ones(len(smiles_not_matching_mass2motif), dtype=np.int8))
        ids = np.arange(0, len(category))
        smiles_df = pd.DataFrame(
            {"ids": ids, "category": category, "smiles": smiles_matching_mass2motif + smiles_not_matching_mass2motif})
        smiles_df.to_csv(output_csv_file_location, index=False, header=False)

    def create_all_moss_files(self, output_folder: str, minimal_score: float):
        """Creates moss files for all mass2motifs"""
        spectra_per_mass2motifs = self.select_spectra_matching_mass2motif(minimal_score)
        for i, spectra_per_mass2motif in enumerate(spectra_per_mass2motifs):
            output_file_name = os.path.join(output_folder, f"mass2motif_{i}.smiles")
            self.create_moss_file(spectra_per_mass2motif, output_file_name)


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
    pass
