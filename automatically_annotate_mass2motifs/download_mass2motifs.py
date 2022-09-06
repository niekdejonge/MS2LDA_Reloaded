from typing import Dict
import requests
from automatically_annotate_mass2motifs.mass2motif import Mass2Motif


def get_single_motif_set(motifset_name):
    """Downloads 1 motif_set from motifdb

    To find possible motifdb names check out: https://ms2lda.org/motifdb/list_motifsets/
    """

    server_url = "http://ms2lda.org/motifdb"
    # Get the list of motif sets
    output = requests.get(server_url + '/list_motifsets', timeout=60)
    motif_set_list = output.json()

    motif_set_id = motif_set_list[motifset_name]
    motif_set = requests.get(server_url + f"/get_motifset/{motif_set_id}", timeout=60).json()
    return motif_set


def convert_json_to_mass2motif(mass2motif_set: Dict):
    mass2motif_list = []
    for motif_name in mass2motif_set:
        words = list(mass2motif_set[motif_name].keys())
        probabilities = list(mass2motif_set[motif_name].values())
        mass2motif = Mass2Motif(words, probabilities)
        mass2motif_list.append(mass2motif)
    return mass2motif_list


if __name__ == "__main__":
    json_mass_2motif_set = get_single_motif_set("Urine derived Mass2Motifs 2")
    mass2motifs = convert_json_to_mass2motif(json_mass_2motif_set)
    for mass2motif in mass2motifs:
        print(sum(mass2motif.fragments.intensities) + sum(mass2motif.losses.intensities))
