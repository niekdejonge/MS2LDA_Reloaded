import os
import logging
from gensim.models import LdaMulticore
from gensim import corpora


def train_lda_model(dictionary,
                    corpus,
                    num_topics=20,
                    chunksize=20000,
                    passes=5,
                    iterations=400,
                    ):
    """
    :param dictionary:
    :param corpus:
    :param num_topics:
    :param chunksize:
    :param passes:
    :param iterations:
    :param eval_every:
    :return:
    """
    # pylint: disable=too-many-arguments
    # Make a index to word dictionary.
    logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
    # token2id = dictionary.token2id  # This is only to "load" the dictionary.
    # id2word = {v: k for k, v in token2id.items()}
    model = LdaMulticore(
        corpus=corpus,
        id2word=dictionary,
        chunksize=chunksize,
        alpha=0.001,
        eta=0.001,
        iterations=iterations,
        num_topics=num_topics,
        passes=passes,
        eval_every=1)
    return model


def train_and_save_model(corpus_dir, overwrite_existing_file = False):

    save_location = os.path.join(corpus_dir, "lda_model")
    if os.path.isfile(save_location):
        assert overwrite_existing_file, "Save file name already exists"
    dictionary = corpora.Dictionary.load(os.path.join(corpus_dir, "dictionary"))
    corpus = corpora.MmCorpus(os.path.join(corpus_dir, "corpus.mm"))
    model = train_lda_model(dictionary, corpus)
    model.save(save_location)


if __name__ == "__main__":
    corpus_dir = "../data/speed_up_ms2lda/test_corpus"
    train_and_save_model(corpus_dir, True)
