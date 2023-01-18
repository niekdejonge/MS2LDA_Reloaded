import numpy as np
from typing import List
from tqdm import tqdm
from matchms import Spectrum, Fragments
from matchms.filtering import add_losses


def bin_spectra(spectra: List[Spectrum],
                bin_width) -> List[Spectrum]:
    """Bins all spectra"""
    binned_spectra = []
    for spectrum in tqdm(spectra,
                         desc="Binning spectra"):
        binned_spectra.append(bin_spectrum(spectrum, bin_width))
    return binned_spectra

def bin_spectrum(spectrum: Spectrum,
                 bin_width) -> Spectrum:
    """Bins one spectrum"""
    if spectrum.losses is None:
        spectrum = add_losses(spectrum)

    binned_masses = return_binned_masses(spectrum.mz, spectrum.intensities, bin_width)
    binned_losses = return_binned_masses(spectrum.losses.mz, spectrum.losses.intensities, bin_width)
    fragments = Fragments(np.array(list(binned_masses.keys())), np.array(list(binned_masses.values())))
    losses = Fragments(np.array(list(binned_losses.keys())), np.array(list(binned_losses.values())))
    spectrum.peaks = fragments
    spectrum.losses = losses
    return spectrum


def return_binned_masses(list_of_masses,
                         list_of_intensities,
                         bin_width):
    binned_masses = {}
    for i, mass in enumerate(list_of_masses):
        mass_bin = return_mass_bin(mass, bin_width)
        intensity = list_of_intensities[i]
        if mass_bin in binned_masses:
            binned_masses[mass_bin] += intensity
        else:
            binned_masses[mass_bin] = intensity
    return binned_masses

def return_mass_bin(mass, bin_width):
    number_of_decimals = len(str(bin_width).split(".")[1])
    mass_bin = mass - mass % bin_width + 0.5 * bin_width
    return round(mass_bin, number_of_decimals + 1)

def check_correct_bin_width(spectrum: Spectrum,
                            bin_width):
    for fragment in spectrum.peaks.mz:
        assert fragment == return_mass_bin(fragment, bin_width), \
            f"The fragment {fragment} does not match the bin width {bin_width}"
