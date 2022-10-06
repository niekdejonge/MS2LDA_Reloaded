import re
from typing import List, Tuple
import requests
from math import gcd
from functools import reduce
from automatically_annotate_mass2motifs.mass2motif import Mass2Motif


def download_motif_set_from_motifdb(motifset_name, bin_size) -> List[Mass2Motif]:
    """Downloads 1 motif_set from motifdb

    To find possible motifdb names check out: https://ms2lda.org/motifdb/list_motifsets/
    """

    server_url = "http://ms2lda.org/motifdb"
    # Get the list of motif sets
    output = requests.get(server_url + '/list_motifsets', timeout=60)
    motif_set_list = output.json()

    motif_set_id = motif_set_list[motifset_name]
    motif_set = requests.get(server_url + f"/get_motifset/{motif_set_id}", timeout=60).json()
    motif_set_metadata= requests.get(server_url + f"/get_motifset_metadata/{motif_set_id}", timeout=60).json()

    mass2motif_list = []
    for motif_name in motif_set:
        words = list(motif_set[motif_name].keys())
        probabilities = list(motif_set[motif_name].values())
        fragments, fragment_probabilities, losses, loss_probabilities = convert_words_to_peaks(words, probabilities)
        mass2motif = Mass2Motif(fragments, fragment_probabilities, losses, loss_probabilities, bin_size=bin_size,
                                motif_name=motif_name, motif_set_name=motifset_name,
                                manual_annotation=motif_set_metadata[motif_name]["annotation"])
        mass2motif_list.append(mass2motif)
    return mass2motif_list


def get_annotations(motifset_name: str, motif_name: str) -> str:
    server_url = "http://ms2lda.org/motifdb"
    # Get the list of motif sets
    output = requests.get(server_url + '/list_motifsets', timeout=60)
    motif_set_list = output.json()

    motif_set_id = motif_set_list[motifset_name]
    motif_set_metadata= requests.get(server_url + f"/get_motifset_metadata/{motif_set_id}", timeout=60).json()
    annotation = motif_set_metadata[motif_name]["annotation"]
    return annotation


def convert_words_to_peaks(words: List[str], probabilities: List[float]) -> \
        Tuple[List[float], List[float], List[float], List[float]]:
    """Converts words and probabilites into fragments, losses and probabilities
    """
    assert_correct_input(words, probabilities)
    fragments_dict = {}
    losses_dict = {}
    for i, word in enumerate(words):
        result = re.search("\\A(fragment|loss)_([0-9]+\\.[0-9]+)\\Z", word)
        is_loss = result.group(1) == "loss"
        fragment_mz = float(result.group(2))
        if not is_loss:
            assert fragment_mz not in fragments_dict, "No duplicated fragments are expected in 1 mass2motif"
            fragments_dict[fragment_mz] = probabilities[i]
        else:
            assert fragment_mz not in losses_dict, "No duplicated losses are expected in 1 mass2motif"
            losses_dict[fragment_mz] = probabilities[i]
    fragments = sorted(fragments_dict)
    losses = sorted(losses_dict)
    fragment_probabilities = [fragments_dict[fragment] for fragment in fragments]
    loss_probabilities = [losses_dict[loss] for loss in losses]
    return fragments, fragment_probabilities, losses, loss_probabilities


def assert_correct_input(words: List[str], probabilities: List[float]):
    """Checks if words and probabilities are in the expected format
    """
    assert isinstance(words, list)
    assert isinstance(probabilities, list)
    assert len(words) > 0, "More than 0 words are expected in a mass2motif"
    assert len(words) == len(probabilities), "An equal number of words and probabilities was expected"
    for word in words:
        assert isinstance(word, str), "A list of strings was expected as words"
        word_and_fragment = word.split("_")
        assert len(word_and_fragment) == 2, "Words should contain 1 underscore"
        assert word_and_fragment[0] == "fragment" or \
               word_and_fragment[0] == "loss", "Words should start with fragment or loss"
        fragment = word_and_fragment[1]
        whole_number_and_decimal = fragment.split(".")
        assert len(word_and_fragment) == 2, "The mz should be seperated by a '.'"
        whole_number = whole_number_and_decimal[0]
        decimal = whole_number_and_decimal[1]
        assert whole_number.isdigit(), "An unexpected word format was given"
        assert decimal.isdigit(), "An unexpected word format was given"
    for probability in probabilities:
        assert isinstance(probability, float)


def get_largest_possible_bin_size(words):
    """Determines the largest possible_bin_size from the fragment sizes"""
    list_of_decimals = []
    assert len(words) > 0, "More than 0 words are expected in a mass2motif"
    for word in words:
        decimals_str = word.split(".", 1)[1]
        decimals = int(decimals_str)
        list_of_decimals.append(decimals)

    nr_of_decimals = max(list_of_decimals)
    greatest_common_denominator = reduce(gcd, list_of_decimals)
    largest_possible_bin_size = greatest_common_denominator/(10**nr_of_decimals)*2
    return largest_possible_bin_size


def get_all_motifsets():
    pass

