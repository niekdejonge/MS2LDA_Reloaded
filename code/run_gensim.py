import os
import pickle
from typing import List
from spec2vec import SpectrumDocument
from matchms import Spectrum
from matchms.filtering import add_losses
from matchms.filtering import add_parent_mass
from matchms.filtering import default_filters
from matchms.filtering import normalize_intensities
from matchms.filtering import reduce_to_number_of_peaks
from matchms.filtering import require_minimum_number_of_peaks
from matchms.filtering import select_by_mz
from utils import convert_files_to_matchms_spectrum_objects, add_unknown_charges_to_spectra
from gensim.models.callbacks import PerplexityMetric
from gensim.corpora import Dictionary
from gensim.models import LdaMulticore
from gensim.models import CoherenceModel
from gensim import corpora
from tqdm import tqdm
import logging

from pprint import pprint


logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)

def spectrum_processing(s):
    """This is how one would typically design a desired pre- and post-
    processing pipeline."""
    s = default_filters(s)
    s = add_parent_mass(s)
    s = normalize_intensities(s)
    s = reduce_to_number_of_peaks(s, n_required=10, ratio_desired=0.5, n_max=500)
    s = select_by_mz(s, mz_from=0, mz_to=1000)
    s = add_losses(s, loss_mz_from=10.0, loss_mz_to=200.0)
    s = require_minimum_number_of_peaks(s, n_required=10)
    return s


def train_lda_model(dictionary,
                    corpus,
                    num_topics=10,
                    chunksize=20000,
                    passes=1,
                    iterations=400,
                    eval_every=None,
                    ):
    # Make a index to word dictionary.
    temp = dictionary[0]  # This is only to "load" the dictionary.
    id2word = dictionary.id2token

    model = LdaMulticore(
        corpus=corpus,
        id2word=id2word,
        chunksize=chunksize,
        alpha='symmetric',
        eta='auto',
        iterations=iterations,
        num_topics=num_topics,
        passes=passes,
        eval_every=eval_every)
    return model


def create_corpus(filtered_spectra: List[Spectrum],
                  dir_for_corpus_and_dict):
    docs = tqdm([SpectrumDocument(filtered_spectrum, n_decimals=1) for filtered_spectrum in filtered_spectra],
                desc="converting spectra into spectrum documents")
    print("Creating a dictionary")
    dict = Dictionary(docs)
    dict.save_as_text(os.path.join(dir_for_corpus_and_dict, "dictionary"))
    my_corpus = tqdm([dict.doc2bow(doc) for doc in docs],
                     desc="Creating corpus")
    corpora.MmCorpus.serialize(os.path.join(dir_for_corpus_and_dict, "corpus.mm"), my_corpus)


def train_model(corpus_dir):
    dict = Dictionary.load_from_text(os.path.join(corpus_dir, "dictionary"))
    corpus = corpora.MmCorpus(os.path.join(corpus_dir, "corpus.mm"))
    print(corpus)
    model = train_lda_model(dict, tqdm(corpus), num_topics=1000)
    model.save(os.path.join(corpus_dir, "lda_model"))

def analyze_model_results(corpus_dir):
    # dict = Dictionary.load_from_text(os.path.join(corpus_dir, "dictionary"))
    corpus = corpora.MmCorpus(os.path.join(corpus_dir, "corpus.mm"))
    model = LdaMulticore.load(os.path.join(corpus_dir, "lda_model"))
    print(model)
    results = model.print_topics(num_words=30)
    pprint(results)
    # top_topics = model.top_topics(corpus)
    # pprint(top_topics)

if __name__ == "__main__":
    # spectra = convert_files_to_matchms_spectrum_objects("C:/Users/jonge094/PycharmProjects/PhD_MS2Query/ms2query/data/test_dir/test_spectra/Brocadia-Excl1-POS-1.mzML")
    # filtered_spectra = pickle.load(open("../data/ALL_GNPS_filtered.pickle", "rb"))
    # filtered_spectra = [filtered_spectrum for filtered_spectrum in filtered_spectra if filtered_spectrum is not None]
    # create_corpus(filtered_spectra,
    #               "../data/corpus")
    # train_model("../data/corpus")
    analyze_model_results("../data/corpus")
    #
    # avg_topic_coherence = sum([t[1] for t in top_topics]) / 10
    # print('Average topic coherence: %.4f.' % avg_topic_coherence)
    #
    # print(len(top_topics))
    # pprint(top_topics)
    # print(CoherenceModel(model=model, corpus=my_corpus, coherence='u_mass').get_coherence())
