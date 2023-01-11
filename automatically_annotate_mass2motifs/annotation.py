from typing import List, Optional
import re
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt


class Annotation:
    def __init__(self,
                 moss_annotations: pd.DataFrame,
                 minimal_similarity: float,
                 moss_minimal_relative_support: float,
                 moss_maximal_relative_support_complement: float):

        # Settings
        self.minimal_similarity = minimal_similarity
        self.moss_minimal_relative_support = moss_minimal_relative_support
        self.moss_maximal_relative_support_complement = moss_maximal_relative_support_complement

        self.moss_annotations = moss_annotations
        self.assert_correct_moss_annotations()
        self.reformat_moss_annotations()

    def reformat_moss_annotations(self):
        # rename columns to new names
        self.moss_annotations.rename(columns={"description": "smiles"}, inplace=True)
        self.moss_annotations.set_index("smiles", inplace=True)
        self.moss_annotations["diff_s_rel_and_c_rel"] = self.moss_annotations["s_rel"] - self.moss_annotations["c_rel"]
        self.moss_annotations.sort_values("diff_s_rel_and_c_rel", inplace=True, ascending=False)
        self.moss_annotations = self.moss_annotations[["s_abs", "c_abs", "s_rel", "c_rel", "diff_s_rel_and_c_rel"]]

    def assert_correct_moss_annotations(self):
        assert isinstance(self.moss_annotations, pd.DataFrame)
        assert list(self.moss_annotations.columns) == ["id", "description", "nodes", "edges", "s_abs", "s_rel", "c_abs", "c_rel"], \
            "Expected different moss_annotations, use load_moss_reulsts to load the results from a file"
        assert np.all(self.moss_annotations["s_rel"] >= self.moss_minimal_relative_support), \
            "moss minimal relative support given is wrong, since a lower s_rel was found"
        assert np.all(self.moss_annotations["c_rel"] <= self.moss_maximal_relative_support_complement),\
            "moss maximal relative support complement given is wrong, since a higher c_rel was found"

    def __str__(self):
        pd.set_option('display.precision', 1)
        return str(self.moss_annotations.head())

    def to_dict(self):
        """Converts the object to a json storable format"""
        class_dict = self.__dict__.copy()
        annotations = self.moss_annotations.__deepcopy__()
        class_dict["moss_annotations"] = annotations.to_dict()
        return class_dict

    def visualize(self):
        smile_annotations = list(self.moss_annotations.index)
        fig = plt.figure()
        for i, smile in enumerate(smile_annotations):
            mol = Chem.MolFromSmiles(smile)
            ax = plt.subplot2grid((1,len(smile_annotations)), (0,i))
            im = Draw.MolToImage(mol)
            ax.imshow(im)
            ax.axis("off")
            ax.set_title(f"s_rel: {self.moss_annotations['s_rel'][smile]}\n"
                         f"c_rel: {self.moss_annotations['s_rel'][smile]}\n"
                         f"s_abs: {self.moss_annotations['s_abs'][smile]}\n"
                         f"c_abs: {self.moss_annotations['c_abs'][smile]}\n",
                         y=-0.5)
        fig.suptitle(f"Minimal similarity: {self.minimal_similarity}\n"
                     f"Min rel supp: {self.moss_minimal_relative_support}\n"
                     f"Max rel supp complement: {self.moss_maximal_relative_support_complement}")
        fig.tight_layout()
        return fig


def load_moss_results(file_name: str) -> Optional[pd.DataFrame]:
    """Loads in a moss results file and returns a pd.DataFrame"""
    with open(file_name, "r") as file:
        lines = file.readlines()
        if len(lines) == 0:
            # Moss returns an empty dataframe, when no annotations are found.
            return None
        assert lines[0] == "id,description,nodes,edges,s_abs,s_rel,c_abs,c_rel\n", \
            "Expected the header id,description,nodes,edges,s_abs,s_rel,c_abs,c_rel in a moss output file"

    with open(file_name, "r") as file:
        moss_results = pd.read_csv(file)
    return moss_results


def parse_moss_file_name(file_name):
    """Parses the file name of a moss file"""

    match = re.match(r"\Amass2motif_(.*)_min_(.*).sub\Z", file_name)
    if match is None:
        return None
    mass2motif_name = match.group(0)
    mass2motif_similarity_threshold = match.group(1)
    return mass2motif_name, mass2motif_similarity_threshold
