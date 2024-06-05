import ase.neighborlist as ase_n
from ase.io import read
from scipy import sparse
import glob
import shutil
import tqdm
import numpy as np
import time

# from multiprocessing import Pool, Manager
from joblib import Parallel, delayed


def integrity_test(folder: str, threshold: float):
    files=glob.glob(folder+"/unfiltered/*.xyz")
    folder_out = folder+"/filtered"

    for file in tqdm.tqdm(files):
        mol=read(file)
        cutOff = ase_n.natural_cutoffs(mol)
        # print("DEBUG: Check cutoff at connectivity")
        cutOff = [dist*threshold for dist in cutOff]
        neighborList = ase_n.NeighborList(cutOff, self_interaction=False, bothways=True)
        neighborList.update(mol)
        matrix = neighborList.get_connectivity_matrix(sparse=False)
        n_components, component_list = sparse.csgraph.connected_components(matrix)

        if n_components==1:
            shutil.copy(file,folder_out+"/"+file.split("/")[-1])

def integrity_complexes(args):
    file, threshold, folder_out = args
    mol=read(file)
    cutOff = ase_n.natural_cutoffs(mol)
    rmax = np.max(cutOff)
    i, j, d = ase_n.neighbor_list('ijd', mol, rmax*2.0)
    file_ok = True
    for index_i, index_j, dist in zip(i, j, d):
        # if index_i in [3,4,5]:
        #     print(f"{index_i}:{index_j} - {cutOff[index_i]}+{cutOff[index_j]} -> {dist:.2f} Å")
        if dist*threshold < (cutOff[index_i]+cutOff[index_j]):
            file_ok = False
            break
    if file_ok:
        shutil.copy(file,folder_out+"/"+file.split("/")[-1])

def integrity_test_complexes_hpc(folder: str, threshold: float):
    start_time = time.time()
    files=glob.glob(folder+"/unfiltered/*.xyz")
    folder_out = folder+"/filtered"
    num_tasks = len(files)

    args = [(files[i], threshold, folder_out) for i in range(num_tasks)]

    Parallel(n_jobs=-1)(delayed(integrity_complexes)(n) for n in tqdm.tqdm(args, total=num_tasks))
    end_time = time.time()
    print(f"\t\t\tExecution time: {end_time - start_time} seconds")


def integrity_test_complexes(folder: str, threshold: float):
    start_time = time.time()
    files=glob.glob(folder+"/unfiltered/*.xyz")
    folder_out = folder+"/filtered"

    for file in tqdm.tqdm(files):
            # print("="*20)
            # print("FILE: ", file)
            mol=read(file)
            cutOff = ase_n.natural_cutoffs(mol)
            rmax = np.max(cutOff)
            i, j, d = ase_n.neighbor_list('ijd', mol, rmax*2.0)
            file_ok = True
            for index_i, index_j, dist in zip(i, j, d):
                # if index_i in [3,4,5]:
                #     print(f"{index_i}:{index_j} - {cutOff[index_i]}+{cutOff[index_j]} -> {dist:.2f} Å")
                if dist*threshold < (cutOff[index_i]+cutOff[index_j]):
                    file_ok = False
                    break
            if file_ok:
                shutil.copy(file,folder_out+"/"+file.split("/")[-1])
    end_time = time.time()
    print(f"\t\t\tExecution time: {end_time - start_time} seconds")
