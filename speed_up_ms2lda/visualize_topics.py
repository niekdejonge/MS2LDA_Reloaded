import os
from gensim import corpora
from gensim.models import LdaMulticore, CoherenceModel
from pprint import pprint


def analyze_model_results(corpus_dir):
    # dict = corpora.Dictionary.load(os.path.join(corpus_dir, "dictionary"))
    # print(dict)
    corpus = corpora.MmCorpus(os.path.join(corpus_dir, "corpus.mm"))
    model = LdaMulticore.load(os.path.join(corpus_dir, "lda_model"))
    results = model.print_topics(num_words=30)
    pprint(results)
    # top_topics = model.top_topics(corpus)
    # avg_topic_coherence = sum([t[1] for t in top_topics]) / len(top_topics)
    print(model)
    # print('Average topic coherence: %.4f.' % avg_topic_coherence)
    # print(len(top_topics))
    # pprint(top_topics)
    print(CoherenceModel(model=model, corpus=corpus, coherence='u_mass').get_coherence())

if __name__ == "__main__":
    corpus_dir = "../data/speed_up_ms2lda/corpus"
    analyze_model_results(corpus_dir)

    # logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
    # hdp = gensim.models.HdpModel(corpus, dict)
    # topic_info = hdp.print_topics(num_topics=1000, num_words=100)
    # print(sys.argv[0])
    # print(os.path.join(sys.argv[0], "../data/corpus_brocadia_sample"))
    # analyze_model_results("../data/corpus_brocadia_sample")
