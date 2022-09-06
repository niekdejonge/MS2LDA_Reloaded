from typing import List
import numpy as np
from matchms import Spectrum


def create_test_spectra() -> List[Spectrum]:
    """ Returns a list with two test spectra
    """
    spectrum1 = Spectrum(mz=np.array([100.0, 200.0, 300.0, 300.01
                                      ], dtype="float"),
                         intensities=np.array([0.11106008, 0.12347332,
                                               0.16352988, 1.],
                                              dtype="float"),
                         metadata={'precursor_mz': 900.0,
                                   'charge': 1})
    spectrum2 = Spectrum(mz=np.array([100.0, 100.02, 500.0, 700.00], dtype="float"),
                         intensities=np.array([0.5, 0.3, 0.000001, 1.],
                                              dtype="float"),
                         metadata={'precursor_mz': 907.0,
                                   'charge': 1})
    return [spectrum1, spectrum2]

def test_create():
    assert True