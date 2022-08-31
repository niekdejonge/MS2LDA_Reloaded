import numpy as np
from automatically_annotate_mass2motifs.mass2motif import Mass2Motif


def test_creating_mass2motif():
    words = ['fragment_375.2225',
             'loss_80.0275',
             'loss_128.0625']
    probabilities = [0.0055,
                     0.0012,
                     0.0019]
    mass2motif = Mass2Motif(words, probabilities)
    assert mass2motif.fragments.mz == [375.2225]
    assert np.all(mass2motif.losses.mz == [80.0275, 128.0625])
    assert mass2motif.fragments.intensities == [0.0055]
    assert np.all(mass2motif.losses.intensities == [0.0012, 0.0019])

    assert mass2motif.bin_size == 0.005, "bin_size is not correctly set"


if __name__ == "__main__":
    test_creating_mass2motif()
