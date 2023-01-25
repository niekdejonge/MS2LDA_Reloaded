import os
import numpy as np
from automatically_annotate_mass2motifs.mass2motif import Mass2Motif, load_mass2motifs_json
from automatically_annotate_mass2motifs.utils import store_as_json
from automatically_annotate_mass2motifs.download_mass2motifs import convert_words_to_peaks
from automatically_annotate_mass2motifs.annotation import Annotation, AnnotationSettings
from automatically_annotate_mass2motifs.moss_annotation import load_moss_results
from tests.test_automatically_annotate_mass2motifs.generate_test_data import (generate_moss_results_file,
                                                                              generate_annotation)


def test_convert_words_to_peaks():
    words = ['fragment_375.2225',
             'loss_80.0275',
             'loss_128.0625']
    probabilities = [0.0055,
                     0.0012,
                     0.0019]
    fragments, fragment_probabilities, losses, loss_probabilities = convert_words_to_peaks(words, probabilities)
    assert fragments == [375.2225]
    assert np.all(losses == [80.0275, 128.0625])
    assert fragment_probabilities == [0.0055]
    assert np.all(loss_probabilities == [0.0012, 0.0019])


def test_convert_words_to_peaks_unordered():
    words = ['fragment_375.2225',
             'loss_128.0625',
             'loss_80.0275']
    probabilities = [0.0055,
                     0.0019,
                     0.0012,]
    fragments, fragment_probabilities, losses, loss_probabilities = convert_words_to_peaks(words, probabilities)
    assert fragments == [375.2225]
    assert np.all(losses == [80.0275, 128.0625])
    assert fragment_probabilities == [0.0055]
    assert np.all(loss_probabilities == [0.0012, 0.0019])


def test_creating_mass2motif_with_correct_mass_bin():
    mass2motif = Mass2Motif(fragments=[375.2225], fragment_probabilities=[0.0055], losses=[80.0275, 128.0625],
                            loss_probabilities=[0.0012, 0.0019], bin_size=0.005)
    assert mass2motif.fragments.mz == [375.2225]
    assert np.all(mass2motif.losses.mz == [80.0275, 128.0625])
    assert mass2motif.fragments.intensities == [0.0055]
    assert np.all(mass2motif.losses.intensities == [0.0012, 0.0019])
    assert mass2motif.bin_size == 0.005, "bin_size is not correctly set"


def test_creating_mass2motif_with_incorrect_mass_bin():
    error_raised = False
    try:
        Mass2Motif(fragments=[375.2225], fragment_probabilities=[0.0055], losses=[80.0275, 128.0625],
                   loss_probabilities=[0.0012, 0.0019], bin_size=0.5)
    except AssertionError:
        error_raised = True
    assert error_raised, "An assert should warn for a wrongly specified bin_size"

def test_save_and_load_mass2motif_as_json(tmp_path):
    annotation = generate_annotation()

    mass2motifs = [Mass2Motif(fragments=[100.025], fragment_probabilities= [0.5],
                            losses=[128.075, 200.725], loss_probabilities=[0.1, 0.25],
                            bin_size=0.05, moss_annotations=[annotation]),
                   Mass2Motif(fragments=[50.025], fragment_probabilities=[0.6],
                              losses=[], loss_probabilities=[],
                              bin_size=0.05)]

    file_name = os.path.join(tmp_path, "mass2motifs.json")
    store_as_json(mass2motifs, file_name)
    loaded_mass2motifs = load_mass2motifs_json(file_name)
    for i, mass2motif in enumerate(mass2motifs):
        mass2motif.assert_equal(loaded_mass2motifs[i])

if __name__ == "__main__":
    pass