from typing import List
import numpy as np
from matchms import Spectrum, Fragments

def spectra_with_losses() -> List[Spectrum]:
    spectrum1 = Spectrum(mz=np.array([100.0123918902183, 200.7213821], dtype="float"),
                         intensities=np.array([1.0, 0.5], dtype="float"),
                         metadata={"smiles": "CN=C=O"})
    spectrum1.losses = Fragments(mz=np.array([200.7213821, 300.12312321], dtype="float"),
                         intensities=np.array([0.2, 1.0], dtype="float"))

    spectrum2 = Spectrum(mz=np.array([50.0123918902183, 200.7213821], dtype="float"),
                         intensities=np.array([1.0, 0.5], dtype="float"),
                         metadata={"smiles": "C1CCCCC1"})
    spectrum2.losses = Fragments(mz=np.array([200.7213821, 400.12312321], dtype="float"),
                         intensities=np.array([0.2, 1.0], dtype="float"))
    return [spectrum1, spectrum2]

def binned_spectra_005() -> List[Spectrum]:
    """List of spectra binned at 0.05"""
    spectrum1 = Spectrum(mz=np.array([100.025, 200.725], dtype="float"),
                         intensities=np.array([1.0, 0.5], dtype="float"),
                         metadata={"smiles": "CN=C=O"})
    spectrum1.losses = Fragments(mz=np.array([200.725, 300.125], dtype="float"),
                         intensities=np.array([0.2, 1.0], dtype="float"))

    spectrum2 = Spectrum(mz=np.array([50.025, 200.725], dtype="float"),
                         intensities=np.array([1.0, 0.5], dtype="float"),
                         metadata={"smiles": "C1CCCCC1"})
    spectrum2.losses = Fragments(mz=np.array([200.725, 400.125], dtype="float"),
                         intensities=np.array([0.2, 1.0], dtype="float"))
    return [spectrum1, spectrum2]