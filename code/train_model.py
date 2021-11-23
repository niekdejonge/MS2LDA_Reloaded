import os
from tqdm import tqdm

from gensim.models import LdaMulticore, CoherenceModel
from gensim import corpora

from pprint import pprint


def train_lda_model(dictionary,
                    corpus,
                    num_topics=10,
                    chunksize=20000,
                    passes=1,
                    iterations=400,
                    eval_every=None,
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

def train_model(corpus_dir):
    dict = corpora.Dictionary.load_from_text(os.path.join(corpus_dir, "dictionary"))
    corpus = corpora.MmCorpus(os.path.join(corpus_dir, "corpus.mm"))
    print(corpus)
    model = train_lda_model(dict, tqdm(corpus), num_topics=1000)
    model.save(os.path.join(corpus_dir, "lda_model"))


def analyze_model_results(corpus_dir):
    # dict = corpora.Dictionary.load_from_text(os.path.join(corpus_dir, "dictionary"))
    corpus = corpora.MmCorpus(os.path.join(corpus_dir, "corpus.mm"))
    model = LdaMulticore.load(os.path.join(corpus_dir, "lda_model"))
    print(model)
    results = model.print_topics(num_words=30)
    pprint(results)
    # top_topics = model.top_topics(corpus)
    # pprint(top_topics)


if __name__ == "__main__":
    train_model("../data/corpus")
    analyze_model_results("../data/corpus")

    avg_topic_coherence = sum([t[1] for t in top_topics]) / 10
    print('Average topic coherence: %.4f.' % avg_topic_coherence)

    print(len(top_topics))
    pprint(top_topics)
    print(CoherenceModel(model=model, corpus=my_corpus, coherence='u_mass').get_coherence())
