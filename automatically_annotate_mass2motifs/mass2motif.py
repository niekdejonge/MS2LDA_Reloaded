from typing import List, Optional
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
                 bin_size: Optional[float] = None,
                 motif_name: str = None,
                 motif_set_name:str = None,
                 annotation: str = None):
        self.assert_correct_input(words, probabilities)

        peaks, peak_probabilities, losses, loss_probabilities = self.convert_words_to_peaks(words, probabilities)
        self.fragments = Fragments(np.array(peaks), np.array(peak_probabilities))
        self.losses = Fragments(np.array(losses), np.array(loss_probabilities))
        if bin_size is None:
            self.bin_size = self.get_bin_size(words)
        else:
            self.bin_size = bin_size
            self.assert_correct_bin_size()
        self.motif_name = motif_name
        self.motif_set_name = motif_set_name
        self.annotation = annotation

    @staticmethod
    def convert_words_to_peaks(words, probabilities):
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
        assert len(words) > 0, "More than 0 words are expected in a mass2motif"
        for word in words:
            decimals_str = word.split(".", 1)[1]
            decimals = int(decimals_str)
            list_of_decimals.append(decimals)

        nr_of_decimals = len(decimals_str)
        greatest_common_denominator = reduce(gcd, list_of_decimals)
        bin_size = greatest_common_denominator/(10**nr_of_decimals)*2
        return bin_size

    def assert_correct_bin_size(self):
        assert isinstance(self.bin_size, float), "bin size is expected to be float"
        nr_of_decimals = len(str(self.bin_size).split(".", 1)[1])
        assert nr_of_decimals < 10, "An bin_size with this many decimals is unexpected"
        for mz in self.fragments.mz:
            assert round(mz%self.bin_size, 10) == self.bin_size*0.5, \
                f"Incorrect bin size of {self.bin_size} for mz of {mz}"
        for mz_loss in self.losses.mz:
            assert round(mz_loss%self.bin_size, 10) == self.bin_size*0.5, \
                f"Incorrect bin size of {self.bin_size} for mz of {mz_loss}"

    def __str__(self):
        return str({"fragments": self.fragments.mz,
                    "fragment_probabilities": self.fragments.intensities,
                    "losses": self.losses.mz,
                    "loss_probabilities": self.losses.intensities})


if __name__ == "__main__":
    pass
    # mass2motif = Mass2Motif(['fragment_375.2225', 'loss_80.0275', 'loss_128.0625', 'loss_141.1025', 'loss_147.0525'],
    #                         [0.005584277434409811, 0.0012797822397661643, 0.0019930865131162247, 0.0026998426721283244, 0.0032813585667783463],
    #                         'Urine_MM_motif_319.m2m')
    # print(mass2motif)