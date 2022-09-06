from matchms import Spectrum, Fragments
import numpy as np

from automatically_annotate_mass2motifs.bin_spectra import bin_spectra, bin_spectrum, Binner
from tests.test_automatically_annotate_mass2motifs.generate_test_data import spectra_with_losses, binned_spectra_005

def test_bin_spectra():
    spectra = spectra_with_losses()
    result = bin_spectra(spectra, 0.05)
    expected_result = binned_spectra_005()
    for i, spectrum in enumerate(result):
        assert spectrum.__eq__(expected_result[i]), "Spectra are not binned correctly"

def test_bin_spectrum():
    spectrum = spectra_with_losses()[0]
    result = bin_spectrum(spectrum, Binner(0.1))
    expected_result = Spectrum(mz=np.array([100.05, 200.75], dtype="float"),
                         intensities=np.array([1.0, 0.5], dtype="float"),
                         metadata={"smiles": "CN=C=O"})
    expected_result.losses = Fragments(mz=np.array([200.75, 300.15], dtype="float"),
                         intensities=np.array([0.2, 1.0], dtype="float"))
    assert result.__eq__(expected_result)
