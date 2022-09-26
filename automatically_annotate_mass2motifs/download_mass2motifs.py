from typing import Dict, List
import requests
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
        mass2motif = Mass2Motif(words,
                                probabilities,
                                bin_size=bin_size,
                                motif_name=motif_name,
                                motif_set_name=motif_name,
                                annotation = motif_set_metadata[motif_name]["annotation"])
        mass2motif_list.append(mass2motif)
    return mass2motif_list

def get_all_motifsets():
    pass



if __name__ == "__main__":
    mass2motifs = download_motif_set_from_motifdb("Urine derived Mass2Motifs 2")
    for mass2motif in mass2motifs:
        print(sum(mass2motif.fragments.intensities) + sum(mass2motif.losses.intensities))
