import os
from tqdm import tqdm
from typing import List, Dict, Tuple
from collections import defaultdict

from spec2vec import SpectrumDocument
from matchms import Spectrum
from matchms.filtering import add_losses, default_filters, normalize_intensities, \
    reduce_to_number_of_peaks, require_minimum_number_of_peaks, select_by_mz, \
    select_by_relative_intensity
from utils import convert_file_to_matchms_spectrum_objects
from gensim import corpora
import logging

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)


def load_and_filter_spectra(file_name) -> List[SpectrumDocument]:
    """Loads in all spectra from a file, filters and returns SpectrumDocuments
    :param file_name: The file name of
    :return: A list with filtered spectrum documents
    """
    spectra = convert_file_to_matchms_spectrum_objects(file_name)
    filtered_spectra = tqdm([spectrum_processing(s) for s in spectra],
                            desc="Filtering spectra")
    spectrum_documents = list(tqdm([SpectrumDocument(filtered_spectrum, n_decimals=1)
                                    for filtered_spectrum in filtered_spectra],
                              desc="Creating spectrum documents"))
    return spectrum_documents


def spectrum_processing(spectrum: Spectrum):
    """Preprocessing of spectra"""
    spectrum = default_filters(spectrum)
    spectrum = normalize_intensities(spectrum)
    spectrum = select_by_relative_intensity(spectrum, intensity_from=0.0005, intensity_to=1.0)
    spectrum = reduce_to_number_of_peaks(spectrum, n_required=10, n_max=500)
    spectrum = select_by_mz(spectrum, mz_from=0, mz_to=1000)
    spectrum = add_losses(spectrum, loss_mz_from=10.0, loss_mz_to=1000.0)
    spectrum = require_minimum_number_of_peaks(spectrum, n_required=10)
    return spectrum


def spectrumdoc2bow(spectrum_document: SpectrumDocument,
                    token2id: Dict[str, int]) -> List[Tuple[int, int]]:
    """Converts a spectrum_document into the bag-of-words format (BoW)

    :param spectrum_document:
        A spectrum document for which a indexed bag of word should be created.
    :param token2id:

    :return:
    """
    assert isinstance(spectrum_document, SpectrumDocument), "spectrumdoc2bow expects a SpectrumDocument as input"
    words = spectrum_document.words
    weights = spectrum_document.weights
    assert len(words) == len(weights), "The words and weights are expected to be of equal length"
    # The defaultdict, automatically creates a new key when the key is not present and a value is added.
    bag_of_words = defaultdict(int)
    for i, word in enumerate(words):
        bag_of_words[token2id[word]] += int(round(weights[i] * 10000))
    # return tokenids, in ascending id order
    bag_of_words = sorted(bag_of_words.items())
    return bag_of_words


def create_dict_and_corpus(spectrum_docs: List[SpectrumDocument],
                           dir_for_corpus_and_dict):
    """Creates a dictionary and a corpus for a list of spectrum_documents

    :param spectrum_docs:
        List of spectrum_docs
    :param dir_for_corpus_and_dict:
        The directory in which the output files are stored.
    :return: Stores the dictionary and the corpus in seperate files in the dir_for_corpus_and_dict
    """
    # to do. This seems
    dictionary = corpora.Dictionary(spectrum_docs)
    dictionary.save(os.path.join(dir_for_corpus_and_dict, "dictionary"))
    my_corpus = list(tqdm([spectrumdoc2bow(spectrum_doc, dictionary.token2id)
                           for spectrum_doc in spectrum_docs],
                     desc="Creating corpus"))
    corpora.MmCorpus.serialize(os.path.join(dir_for_corpus_and_dict, "corpus.mm"), my_corpus)


if __name__ == "__main__":
    filtered_spectra = load_and_filter_spectra("../data/Brocadia-Excl1-POS-1.mzML")
    create_dict_and_corpus(filtered_spectra,
                           "../data/speed_up_ms2lda/test_corpus")
