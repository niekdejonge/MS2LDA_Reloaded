import os
from typing import List
import numpy as np
import pandas as pd
from matchms import Spectrum, Fragments
from automatically_annotate_mass2motifs.annotation import Annotation
from automatically_annotate_mass2motifs.moss_annotation import load_moss_results

def spectra_with_losses() -> List[Spectrum]:
    spectrum1 = Spectrum(mz=np.array([100.0123918902183, 200.7213821], dtype="float"),
                         intensities=np.array([1.0, 0.5], dtype="float"),
                         metadata={"smiles": "CN=C=O", "inchikey": "ABCDEFGHIJKLAQ"})
    spectrum1.losses = Fragments(mz=np.array([200.7213821, 300.12312321], dtype="float"),
                         intensities=np.array([0.2, 1.0], dtype="float"))

    spectrum2 = Spectrum(mz=np.array([50.0123918902183, 200.7213821], dtype="float"),
                         intensities=np.array([1.0, 0.5], dtype="float"),
                         metadata={"smiles": "C1CCCCC1", "inchikey": "BBBBBBBBBBBBBB"})
    spectrum2.losses = Fragments(mz=np.array([200.7213821, 400.12312321], dtype="float"),
                         intensities=np.array([0.2, 1.0], dtype="float"))
    return [spectrum1, spectrum2]


def binned_spectra_005() -> List[Spectrum]:
    """List of spectra binned at 0.05"""
    spectrum1 = Spectrum(mz=np.array([100.025, 200.725], dtype="float"),
                         intensities=np.array([1.0, 0.5], dtype="float"),
                         metadata={"smiles": "CN=C=O", "inchikey": "ABCDEFGHIJKLAQ"})
    spectrum1.losses = Fragments(mz=np.array([200.725, 300.125], dtype="float"),
                         intensities=np.array([0.2, 1.0], dtype="float"))

    spectrum2 = Spectrum(mz=np.array([50.025, 200.725], dtype="float"),
                         intensities=np.array([1.0, 0.5], dtype="float"),
                         metadata={"smiles": "C1CCCCC1", "inchikey": "BBBBBBBBBBBBBB"})
    spectrum2.losses = Fragments(mz=np.array([200.725, 400.125], dtype="float"),
                         intensities=np.array([0.2, 1.0], dtype="float"))
    return [spectrum1, spectrum2]


def generate_moss_results_file(file_name):
    with open(file_name, "w") as moss_results_file:
        moss_results_file.write("id,description,nodes,edges,s_abs,s_rel,c_abs,c_rel\n"
                                "1,O(-C-C)-C(-C)-C-C-C,8,7,5,45.454544,4266,20.941534\n"
                                "2,O(-C-C)-C(-C)-C-C,7,6,6,54.545456,5737,28.162584\n")
    return file_name


def generate_moss_annotation():
    moss_annotation = pd.DataFrame([[1, 0, 50.0, 0.0, 50.0],
                                    [2, 1, 100.0, 50.0, 50.0]],
                                   index=["O=C-C-C", "C(-C)-C"],
                                   columns=["s_abs", "c_abs", "s_rel", "c_rel", "diff_s_rel_and_c_rel"])
    moss_annotation.index.name = "smiles"
    return moss_annotation

def generate_annotation()->Annotation:
    moss_annotations = generate_moss_annotation()
    annotation = Annotation(moss_annotations=moss_annotations,
                            minimal_similarity=0.3,
                            moss_minimal_relative_support=40.0,
                            moss_maximal_relative_support_complement=70.0,
                            nr_of_spectra_matching_mass2motif=10,
                            nr_of_spectra_not_matching_mass2motif=100)
    return annotation
