from sklearn.preprocessing import StandardScaler
import numpy as np
import numpy.typing as npt

import tqdm
import glob
import math
import ast
import os

def xyzRead(fname: str):
    fin = open(fname, "r")
    line1 = fin.readline().split()
    natoms = int(line1[0])
    comments = fin.readline()[:-1]
    if comments: 
        try:
            energy = int(comments.split()[-1])
        except:
            energy = None
    else:
        energy = None
    coords = np.zeros([natoms, 3], dtype="float64")
    atomtypes = []
    for x in coords:
        line = fin.readline().split()
        atomtypes.append(line[0])
        x[:] = list(map(float, line[1:4]))
    return natoms, atomtypes, coords, energy

def replaceAtomSymbols(params: dict, fname: str, permut: tuple):
    with open(fname, 'r') as file:
        line = file.readline()
        modified = [line]
        line = file.readline()
        modified.append(line)
        lines = file.readlines()
        for at, lin in enumerate(lines):
            if at in permut:
                modified.append(lin.replace('At', params["ELEM1"]))
            else:
                modified.append(lin.replace('At', params["ELEM2"]))

    with open(fname, 'w') as file:
        file.writelines(modified)
    file.close()


def getCharge(element):
    f = open("mol.txt")
    atomicnum = [line.split()[1] for line in f if line.split()[0] == element]
    f.close()
    return int(atomicnum[0])

def coulombMatrix(natoms, atomtypes, coords):
    i=0 ; j=0    
    colM = np.zeros((natoms,natoms))
    chargearray = np.zeros((natoms,1))
    charge = [getCharge(symbol)  for symbol in atomtypes]
    for i in range(0,natoms):
        colM[i,i]=0.5*charge[i]**2.4   # Diagonal term described by Potential energy of isolated atom
        for j in range(i+1,natoms):
            dist= np.linalg.norm(coords[i,:] - coords[j,:])   
            colM[j,i] = charge[i]*charge[j]/dist   #Pair-wise repulsion 
            colM[i,j] = colM[j,i]
    return colM
 
def eigenCoulomb(natoms, atomtypes, coords):
    sCoulomb = coulombMatrix(natoms, atomtypes, coords)
    sCoulomb = sCoulomb.astype(int)
    eigValues = -np.sort(-np.linalg.eigvals(sCoulomb))
    if np.any(np.iscomplex(eigValues)) == False:
        return eigValues#[0:natoms]
    else:
        # print('\nWARNING: complex (coulomb matrix) engenvalues for ' + fname + ' employing np.linalg.eigvals function.\n') 
        eigValues = -np.sort(-np.linalg.eigvalsh(sCoulomb))
        if np.any(np.iscomplex(eigValues)):
            # print('\nWARNING: complex (coulomb matrix) engenvalues for ' + fname + ' employing np.linalg.eigvalsh function.\n WARNING: only the real part will be returned.\n')
            return eigValues.real#[0:natoms]
        else :
            return eigValues#[0:atoms]
        
def generateXYZ(list_atoms, list_coords, energy, pfile, outfolder):
    fname = outfolder+pfile+'.xyz'
    with open(fname, 'w') as file:
        file.write(f'{len(list_atoms)}\n')
        file.write(f"Energy = {energy}\n")
        for idx in range(len(list_atoms)):
            file.write(f" {list_atoms[idx]} "+" ".join(map(str, list_coords[idx])) + '\n')

def getCoulombEig(files):
    X = []
    list_atoms = []
    list_coords = []
    coulomb = []
    energies = []
    for file in files:
        natoms, atomtypes, coords, energy = xyzRead(file)
        coulomb.append(eigenCoulomb(natoms, atomtypes, coords))
        energies.append(energy)            
    coulomb = np.array(coulomb)
    coulomb = StandardScaler().fit_transform(coulomb)
    return coulomb, energies
