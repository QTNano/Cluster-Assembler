import itertools
import os
import glob
from ase import Atoms
# from dscribe.descriptors import CoulombMatrix
from sklearn.preprocessing import StandardScaler
import tqdm
import numpy as np
import random

import core.tools as tools
import core.representatives as rep
import core.connectivity as net

def get_templates(params: dict) -> list:
    natoms = params["NUMELEM1"]+params["NUMELEM2"]
    combinations = itertools.combinations(range(natoms), params["NUMELEM1"])
    return list(combinations)

def gen_permutations(params: dict, files: str) -> None:
    """
    Generate all permutations of the alloys
    """
    atom1 = params["ELEM1"]
    atom2 = params["ELEM2"]
    natoms = params["NUMELEM1"]+params["NUMELEM2"]
    outfolder = params["MOD2"]["TMP_FOLDER"]+"/unfiltered/"
    maxgen_frame = params["MOD2"]["MAX_GEN_PER_FRAME"]
    pre_selection = params["MOD2"]["PRE_SELECTION"]
        # "RUN_MOD_ZERO": false,

    combinations = list(itertools.combinations(range(natoms), params["NUMELEM1"]))
    ids_combinations = list(range(len(combinations)))

    for file in tqdm.tqdm(files, desc='Processing Frames', position=0):
        pfiles = []
        X = []
        list_coords = []
        list_atoms = []
        coulomb = []
        base_name = os.path.basename(file)
        base_name = os.path.splitext(base_name)[0]
        natoms, atomtypes, coords, energy = tools.xyzRead(file)
        random.shuffle(ids_combinations)

        for index, id_combination in enumerate(ids_combinations):
            if index < maxgen_frame:
                for id_atom in range(len(atomtypes)):
                    if id_atom in combinations[id_combination]:
                        atomtypes[id_atom] = atom1
                    else:
                        atomtypes[id_atom] = atom2
                pfile = base_name+"_P"+str(index)
                if pre_selection:
                    pfiles.append(pfile[:])
                    list_coords.append(coords[:])
                    list_atoms.append(atomtypes)
                    coulomb.append(tools.eigenCoulomb(natoms, list_atoms[index], coords))
                else:
                    tools.generateXYZ(atomtypes, coords, 0.0, pfile, outfolder)

        if pre_selection:
            coulomb = np.array(coulomb)
            coulomb = StandardScaler().fit_transform(coulomb)
            sel_samples = rep.get_representatives(params["MOD2"], coulomb, pfiles)

            for idx in sel_samples:
                tools.generateXYZ(list_atoms[idx], list_coords[idx], 0.0, pfiles[idx], outfolder)

