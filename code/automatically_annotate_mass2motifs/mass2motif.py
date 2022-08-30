from typing import List
import numpy as np
from math import gcd
from functools import reduce
from matchms import Fragments

class Mass2Motif:
    """Stores Mass2Motif information

    """

    def __init__(self, words: List[str],
                 probabilities: List[float],
                 motif_name):
        assert words
        self.words = words
        self.probabilities = probabilities
        self.motif_name = motif_name
        fragments = Fragments(np.array(list(binned_masses.keys())), np.array(list(binned_masses.values())))

        self.bin_size = self.get_bin_size()
        self.assert_correct_mass2motif()

    def __str__(self):
        return str({"motif_name": self.motif_name, "words": self.words, "probabilities": self.probabilities})

    def assert_correct_mass2motif(self):
        """Checks if words and probabilities are in the expected format"""
        assert isinstance(self.words, list)
        assert isinstance(self.probabilities, list)
        assert len(self.words) > 0, "More than 0 words are expected in a mass2motif"
        assert len(self.words) == len(self.probabilities), "An equal number of words and probabilities was expected"
        for word in self.words:
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
        for probability in self.probabilities:
            assert isinstance(probability, float)

    def get_bin_size(self):
        """Determines the bin_size from the fragment sizes"""
        list_of_decimals = []
        for word in self.words:
            decimals_str = word.split(".", 1)[1]
            decimals = int(decimals_str)
            list_of_decimals.append(decimals)

        nr_of_decimals = len(decimals_str)
        greatest_common_denominator = reduce(gcd, list_of_decimals)
        bin_size = greatest_common_denominator/(10**nr_of_decimals)*2
        return bin_size
