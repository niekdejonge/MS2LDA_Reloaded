import os
import pandas as pd
from automatically_annotate_mass2motifs.download_mass2motifs import get_annotations

def load_moss_results(file_name):
    with open(file_name, "r") as file:
        moss_results = pd.read_csv(file, index_col="id")

    moss_results["diff_support"] = moss_results["s_rel"] - moss_results["c_rel"]
    moss_results.sort_values("diff_support", inplace=True, ascending=False)
    return moss_results


def print_substructures_and_annotation(mass2motif_dir,
                                       mass2motif_set_name):

    for file_name in sorted(os.listdir(mass2motif_dir)):
        if os.path.splitext(file_name)[1] == ".sub":
            mass2motif_name = file_name.split("_min_")[0][11:]
            annotation = get_annotations(mass2motif_set_name,
                                         mass2motif_name)
            result = load_moss_results(os.path.join(mass2motif_dir, file_name))
            pd.set_option('display.precision', 1)
            print(result[["description", "s_abs", "s_rel", "c_rel"]].head())
            print(annotation)

if __name__ == "__main__":
    print_substructures_and_annotation("../data/moss/gnps_library_derived_mass2motifs/different_thresholds/moss_results/",
                                       "GNPS library derived Mass2Motifs")
