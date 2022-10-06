import numpy as np
from typing import List
from tqdm import tqdm
from matchms import Spectrum, Fragments
from matchms.filtering import add_losses


class Binner:
    """Class for creating bins and storing settings"""
    def __init__(self, bin_width):
        self.bin_width = bin_width
        self.number_of_decimals = len(str(bin_width).split(".")[1])

    def return_mass_bin(self, mass):
        mass_bin = mass - mass % self.bin_width + 0.5 * self.bin_width
        return round(mass_bin, self.number_of_decimals + 1)

    def return_binned_masses(self, list_of_masses, list_of_intensities):
        binned_masses = {}
        for i, mass in enumerate(list_of_masses):
            mass_bin = self.return_mass_bin(mass)
            intensity = list_of_intensities[i]
            if mass_bin in binned_masses:
                binned_masses[mass_bin] += intensity
            else:
                binned_masses[mass_bin] = intensity
        return binned_masses


def bin_spectrum(spectrum: Spectrum,
                 binner: Binner) -> Spectrum:
    """Bins one spectrum"""
    if spectrum.losses is None:
        spectrum = add_losses(spectrum)
    binned_masses = binner.return_binned_masses(spectrum.mz, spectrum.intensities)
    binned_losses = binner.return_binned_masses(spectrum.losses.mz, spectrum.losses.intensities)
    fragments = Fragments(np.array(list(binned_masses.keys())), np.array(list(binned_masses.values())))
    losses = Fragments(np.array(list(binned_losses.keys())), np.array(list(binned_losses.values())))
    spectrum.peaks = fragments
    spectrum.losses = losses
    return spectrum


def bin_spectra(spectra: List[Spectrum],
                bin_width) -> List[Spectrum]:
    """Bins all spectra"""
    binner = Binner(bin_width)
    binned_spectra = []
    for spectrum in tqdm(spectra,
                         desc="Binning spectra"):
        binned_spectra.append(bin_spectrum(spectrum, binner))
    return binned_spectra
