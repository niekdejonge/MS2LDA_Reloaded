from typing import List
import numpy as np
from math import gcd
from functools import reduce
from matchms import Fragments
import re


class Mass2Motif:
    """Stores Mass2Motif information

    """

    def __init__(self, words: List[str],
                 probabilities: List[float],
                 motif_name: str = None):
        self.assert_correct_input(words, probabilities)
        self.motif_name = motif_name

        peaks, peak_probabilities, losses, loss_probabilities = self.convert_words_to_peaks(words, probabilities)
        self.fragments = Fragments(np.array(peaks), np.array(peak_probabilities))
        self.losses = Fragments(np.array(losses), np.array(loss_probabilities))

        self.bin_size = self.get_bin_size(words)

    @staticmethod
    def convert_words_to_peaks(words, probabilities):
        peaks = []
        peak_probabilities = []
        losses = []
        loss_probabilities = []
        for i, word in enumerate(words):
            result = re.search("\A(fragment|loss)_([0-9]+\.[0-9]+)\Z", word)
            is_loss = result.group(1) == "loss"
            fragment_mz = float(result.group(2))
            if not is_loss:
                peaks.append(fragment_mz)
                peak_probabilities.append(probabilities[i])
            else:
                losses.append(fragment_mz)
                loss_probabilities.append((probabilities[i]))
        return peaks, peak_probabilities, losses, loss_probabilities

    @staticmethod
    def assert_correct_input(words, probabilities):
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

    @staticmethod
    def get_bin_size(words):
        """Determines the bin_size from the fragment sizes"""
        list_of_decimals = []
        for word in words:
            decimals_str = word.split(".", 1)[1]
            decimals = int(decimals_str)
            list_of_decimals.append(decimals)

        nr_of_decimals = len(decimals_str)
        greatest_common_denominator = reduce(gcd, list_of_decimals)
        bin_size = greatest_common_denominator/(10**nr_of_decimals)*2
        return bin_size

    def __str__(self):
        return str({"motif_name": self.motif_name,
                    "fragments": self.fragments.mz,
                    "fragment_probabilities": self.fragments.intensities,
                    "losses": self.losses.mz,
                    "loss_probabilities": self.losses.intensities})



if __name__ == "__main__":
    mass2motif = Mass2Motif(['fragment_375.2225', 'loss_80.0275', 'loss_128.0625', 'loss_141.1025', 'loss_147.0525'],
                            [0.005584277434409811, 0.0012797822397661643, 0.0019930865131162247, 0.0026998426721283244, 0.0032813585667783463],
                            'Urine_MM_motif_319.m2m')
    print(mass2motif)