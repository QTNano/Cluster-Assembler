#!/usr/bin/python

#   Copyright 2023 Raphael Bühler, Max Schütz, Karla F. Andriani,
#   João Paulo A. de Mendonça, Vivianne K. Ocampo-Restrepo, Marcos G. Quiles, 
#   Christian Gemel, Juarez L. F. Da Silva, Roland A. Fischer
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

import os
import sys
import random
import json
import glob
import numpy as np
import tqdm

# Local packages
import core.connectivity as net
import core.representatives as rep
import core.permute as perm
import core.tools as tools
import core.frames as gen
import core.complexes as complexes



def start_message() -> bool:
	"""
	Open info about the program parameters
	"""

	print("""\n\n\n\t~?JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ?!     
	P&#P55555555555555555555555555555555555B&#:    
	P&P                                    ?&B:    
	P&P            ^YPPPPPPPP5!            J&B:    
	P&P           7#&5??????J#&J           J&B:    
	P&P          J&#7        ~B&5:         J&B:    
	P&P   ^!~~~~5&B~          :P&G!~~~~^   J&B:    
	P&P   JBBBBB##J            !#&BBBBBP.  J&B:    
	P&P   ......7##J          7##Y......   J&B:    
	P&P          ~B&5:      .J##?          J&B:    
	P&P           :P&BGGGGGGG&#!           J&B:    
	P&P            .!7777777?B#5.          J&B:    
	P&P                      ^G&G^         J&B:    
	P&P                       .?!.         ?&B:    
	P&B?????????????^ .^^^^. .?????????????P&#:    
	7PPPPPPPPPPPPPPG! .~~~~: :PPPPPPPPPPPPPPPY.   
	""")

	print("\n	QTNano Cluster Assembler \n\n")


	if len(sys.argv) != 2:
		print("\tUsage: \n")
		print(f"\t\t$ python {sys.argv[0]} parameters.json\n")
		print(f"\t\t\tparameters.json - a json file that contains all the simulation parameters")
		return False	
	return True

def load_parameters() -> dict:
	"""
	Read the simulation parameters from a json file
	"""
	try:
		with open(sys.argv[1]) as f:
			j = json.load(f)
	except FileNotFoundError:
		print("[ERROR] input json file not found")
		exit()
	except json.decoder.JSONDecodeError:
		print("[ERROR] input json format is invalid")
		exit()
	except Exception:
		print("[ERROR] an unexpected error happened trying to read the input")
		exit()

	return j

def M0_Selection(params: dict, adhoc: bool) -> None:
	print("\n\t\tModule 0: Representative Sample Selection:")

	if adhoc:
		inputfolder = params["INPUT_FOLDER"]
		outfolder = params["OUTPUT_FOLDER"]
	else:
		inputfolder = params["TMP_FOLDER"]+"/filtered/"
		outfolder = params["TMP_FOLDER"]+"/selected/"

	os.system(f"rm -rf {outfolder}")
	os.system(f"mkdir {outfolder}")

	print(f"\n\t\t\tLoading files from: {inputfolder}")

	files=sorted(glob.glob(f"{inputfolder}/*.xyz"))
	if not files:
		print(f"\n\t\t\tFail - Folder {inputfolder} is empty\n\n")
		exit()

	print(f"\n\t\t\tSample Selection - Files Loaded {len(files)}")
	if len(files) > 1:
		coulomb, energies = tools.getCoulombEig(files)
		sel_samples = rep.get_representatives(params, coulomb, energies, files)
	else:
		sel_samples = [0]

	for id in sel_samples:
		infile = files[id]
		outfile = f"{outfolder}/{os.path.basename(infile)}"
		os.system(f"cp {infile} {outfile}")
	inFiles = len(glob.glob(f"{outfolder}/*.xyz"))
	print(f"\n\t\t\tTotal of Selected Samples: {inFiles}")

	print(f"\n\t\t\tSelectec Files are available at {outfolder}")

	print("\n\t\tEnd of module 0 - Representative Selection")

def M1_frame_family(params: dict) -> None:
	print("\n\nModule 1: Frame Family:")

	tmpfolder = params["MOD1"]["TMP_FOLDER"]
	outfolder = params["MOD1"]["OUTPUT_FOLDER"]
	threshold = params["MOD1"]["INTEGRITY"]
	os.system("rm -rf log.txt")
	os.system(f"rm -rf {tmpfolder}")
	os.system(f"rm -rf {outfolder}")
	os.system(f"mkdir {tmpfolder}")
	os.system(f"mkdir {tmpfolder}/unfiltered")
	os.system(f"mkdir {tmpfolder}/filtered")
	os.system(f"mkdir {tmpfolder}/selected")
	os.system(f"mkdir {outfolder}")

	inFiles = 0
	if params["MOD1"]["INPUT_FOLDER"]:
		inputfolder = params["MOD1"]["INPUT_FOLDER"]
		files=sorted(glob.glob(f"{inputfolder}/*.xyz"))
		if files:
			print(f"\t\tLoading Unary geometries from: {inputfolder}")
			os.system(f"cp {inputfolder}/*.xyz {tmpfolder}/unfiltered/")

		subfolder_pattern = os.path.join(inputfolder, '*/')
		subfolders = glob.glob(f"{inputfolder}/*/")
		if subfolders:
			for subfolder in subfolders:
				files=sorted(glob.glob(f"{subfolder}/*.xyz"))
				if files:
					subn = os.path.basename(os.path.normpath(subfolder))
					for file in files:
						subf = os.path.basename(file)
						fout = f"{tmpfolder}/unfiltered/{subn}_{subf}"
						os.system(f"cp {file} {fout}")
		
		inFiles = len(glob.glob(f"{tmpfolder}/unfiltered/*.xyz"))
		print(f"\t\t\tTotal of Unary geometries from: {inputfolder}: {inFiles}")

	numsamples = params["MOD1"]["NUMGEN"]
	factor = params["MOD1"]["RADIUS_FACTOR"]
	gamma = params["MOD1"]["GAMMA"]

	natoms = 0
	elements = []
	for key, value in params.items():
		if key.startswith("NUMELEM"):
			natoms += value
		if key.startswith("ELEM"):
			elements.append(value)

	for atom in elements:
		print(f"\t\tGenerating new structures of {atom}{natoms}...")
		gen.genSamples(numsamples, "CUBE", atom, natoms, factor, gamma, f"{tmpfolder}/unfiltered/")
		gen.genSamples(numsamples, "SPHERE", atom, natoms, factor, gamma, f"{tmpfolder}/unfiltered/")

	inFiles = len(glob.glob(f"{tmpfolder}/unfiltered/*.xyz")) - inFiles
	print(f"\t\t\tTotal of Generated Frames: {inFiles}")

	print("\n\t\tTesting the integrity of the nanoclusters")
	net.integrity_test(tmpfolder, threshold) # /unfiltered to /filtered

	inFiles = len(glob.glob(f"{tmpfolder}/filtered/*.xyz"))
	print(f"\t\t\tTotal of Frames After the Integrity Test: {inFiles}")

	files=glob.glob(f"{tmpfolder}/filtered/*.xyz")
	if not files:
		print("\n\nERROR: Frames did not pass the integrity test\n\n")
		exit()

	runZero = params["MOD1"]["RUN_MOD_ZERO"]

	if runZero:
		print("\n\t\tSelecting representative frames (k-means)")
		M0_Selection(params["MOD1"], False)
		filfold = tmpfolder+"/selected"

	else:
		filfold = tmpfolder+"/filtered"

	# os.system(f"cp {filfold}/*.xyz {outfolder}/")

	print(f"\n\t\tCopying final files to the output folder")
	files=glob.glob(f"{filfold}/*.xyz")
	for file in tqdm.tqdm(files):
		os.system(f"cp {file} {outfolder}/")

	print(f"\n\t\tFrames are available at {outfolder}")

	print("\nEnd of module 01 - Frame Family")

def M2_core_family(params: dict) -> None:
	print("\n\nModule 2: Core Family:")

	inputfolder = params["MOD2"]["INPUT_FOLDER"]
	tmpfolder = params["MOD2"]["TMP_FOLDER"]
	outfolder = params["MOD2"]["OUTPUT_FOLDER"]
	threshold = params["MOD2"]["INTEGRITY"]
	maxgen_frame = params["MOD2"]["MAX_GEN_PER_FRAME"]

	os.system(f"rm -rf {tmpfolder}")
	os.system(f"rm -rf {outfolder}")
	os.system(f"mkdir {tmpfolder}")
	os.system(f"mkdir {tmpfolder}/unfiltered")
	os.system(f"mkdir {tmpfolder}/filtered")
	os.system(f"mkdir {tmpfolder}/selected")
	os.system(f"mkdir {outfolder}")

	print(f"\n\t\tLoading frames from: {inputfolder}")

	files=sorted(glob.glob(f"{inputfolder}/*.xyz"))
	if not files:
		print(f"\n\t\tFail - Folder {inputfolder} is empty\n\n")
		exit()

	combinations = perm.get_templates(params)
	print(f"\t\t\tTotal of files: {len(files)}")
	print(f"\t\t\tMax Permutations per file: {len(combinations)}")

	print(f"\n\t\tGenerating Cores:")
	if len(combinations) < maxgen_frame:
		print(f"\t\t\tWARNING: Generations ({maxgen_frame}) exceed the maximum permutations")
		print(f"\t\t\t         A total of {len(combinations)} cores will be generated per frame")
	else:
		print(f"\t\t\tA total of {maxgen_frame} cores will be generated per frame")

	perm.gen_permutations(params, files)

	inFiles = len(glob.glob(f"{tmpfolder}/unfiltered/*.xyz"))
	print(f"\t\t\tTotal of Cores: {inFiles}")

	print("\n\t\tTesting the integrity of the nanoclusters:")
	net.integrity_test(tmpfolder, threshold) # /unfiltered to /filtered

	inFiles = len(glob.glob(f"{tmpfolder}/filtered/*.xyz"))
	print(f"\t\t\tTotal of Cores After the Integrity Test: {inFiles}")

	files=glob.glob(f"{tmpfolder}/filtered/*.xyz")
	if not files:
		print("\n\nERROR: Cores did not pass the integrity test\n\n")
		exit()

	runZero = params["MOD2"]["RUN_MOD_ZERO"]

	if runZero:
		print("\n\t\tSelecting final representative cores (k-means)")
		M0_Selection(params["MOD2"], False)
		filfold = tmpfolder+"/selected"

	else:
		filfold = tmpfolder+"/filtered"

	# os.system(f"cp {filfold}/*.xyz {outfolder}/")
	print(f"\n\t\tCopying final files to the output folder")
	files=glob.glob(f"{filfold}/*.xyz")
	for file in tqdm.tqdm(files):
		os.system(f"cp {file} {outfolder}/")


	print(f"\n\t\tCores are available at {outfolder}")

	print("\nEnd of module 02 - Core Family")

def M3_add_ligants(params: dict) -> None:
	print("\n\nModule 3: Complexes Generation:")

	coresfolder = params["MOD3"]["CORES_FOLDER"]
	ligandsfolder = params["MOD3"]["LIGANDS_FOLDER"]
	tmpfolder = params["MOD3"]["TMP_FOLDER"]
	outfolder = params["MOD3"]["OUTPUT_FOLDER"]
	threshold = params["MOD3"]["INTEGRITY"]
	os.system(f"rm -rf {tmpfolder}")
	os.system(f"rm -rf {outfolder}")
	os.system(f"mkdir {tmpfolder}")
	os.system(f"mkdir {tmpfolder}/unfiltered")
	os.system(f"mkdir {tmpfolder}/filtered")
	os.system(f"mkdir {tmpfolder}/selected")
	os.system(f"mkdir {outfolder}")

	lig_distribution = params["MOD3"]["LIGANDS_DISTRIBUTION"]
	lig_orientation = params["MOD3"]["LIGANDS_ORIENTATION"]
	numsim = params["MOD3"]["N_SAMPLES"]
	ncores = params["MOD3"]["N_CORES"]
	nsites = np.sum(lig_distribution)
	ligdist = params["MOD3"]["LIGANDS_DISTANCE"]
	sites_ids = list(range(nsites))
	deformation = params["MOD3"]["DEFORMATION"]

	natoms = 0
	for key, value in params.items():
		if key.startswith("NUMELEM"):
			natoms += value

	print(f"\n\t\tLoading Cores from: {coresfolder}")
	cores=sorted(glob.glob(f"{coresfolder}/*.xyz"))
	if len(cores) < ncores:
		ncores = len(cores)
	print(f"\t\t\tTotal of Cores in the Folder: {len(cores)} - Up to {ncores} cores will be loaded")

	if not cores:
		print(f"\n\t\tFail - Folder {coresfolder} is empty\n\n")
		exit()

	print(f"\n\t\tLoading Ligands from: {ligandsfolder}")
	ligands=sorted(glob.glob(f"{ligandsfolder}/*.xyz"))
	print(f"\t\t\tTotal of Ligands: {len(ligands)} - {len(lig_distribution)} ligands will be loaded")

	if not ligands:
		print(f"\n\t\tFail - Folder {ligandsfolder} is empty\n\n")
		exit()

	if ncores < len(cores):
		idx_cores = random.sample(range(len(cores) + 1), ncores)
		cores = [core for (id, core) in enumerate(cores) if id in idx_cores]
	else:
		idx_cores = list(range(len(cores)))
	
	print(f"\n\t\tPositioning Ligands:")

	if natoms < np.sum(lig_distribution):
		print(f"\n\t\t\tFail - Number of ligands [{np.sum(lig_distribution)}] is superior to the number of atoms in the core [{natoms}]\n\n")
		exit()

	for core in tqdm.tqdm(cores):
		base_name = os.path.basename(core)
		base_name = os.path.splitext(base_name)[0]
		core_natoms, core_atomtypes, core_coords, _ = tools.xyzRead(core)
		core_coords = complexes.center_mol(core_coords)
		sites = complexes.fibonacci_sphere(nsites)
		sites = complexes.optimize_sites(sites) # Force-field approach to optimize the sites

		for sim in range(numsim):
			sites = complexes.rotate_atoms(sites) #random rotation of the sites
			decision = random.random()
			if decision < 0.5:
				core_sites, radius = complexes.adjust_sites(core_coords, sites, ligdist, deformation)
			else:
				core_sites, radius = complexes.adjust_sites_oriented(core_coords, sites, ligdist, deformation)
				
			idx_site = 0
			random.shuffle(sites_ids)
			complex_atoms = core_atomtypes[:]
			complex_coords = core_coords[:]
			for idx, numlig in enumerate(lig_distribution):
				lig_natoms, lig_atomtypes, lig_coords, _ = tools.xyzRead(ligands[idx])
				lig_coords = complexes.rotate_atoms(lig_coords) #random rotation of the ligands
				lig_coords = complexes.center_mol(lig_coords)
				for lig in range(numlig):
					dist = np.linalg.norm(core_sites[sites_ids[idx_site]])
					direction = np.array(core_sites[sites_ids[idx_site]])/dist
					to_append = complexes.positining_ligand(lig_coords, direction, dist, core_coords, lig_orientation[idx])
					complex_atoms.extend(lig_atomtypes)
					complex_coords = np.concatenate((complex_coords, to_append), axis=0)
					idx_site += 1

			pfile = "/"+base_name+"_C"+str(sim)
			tools.generateXYZ(complex_atoms, complex_coords, 0.0, pfile, tmpfolder+"/unfiltered")

	inFiles = len(glob.glob(f"{tmpfolder}/unfiltered/*.xyz"))
	print(f"\t\t\tTotal of Generated Complexes: {inFiles}")
	print("\n\t\tTesting the integrity of the complexes")
	net.integrity_test_complexes_hpc(tmpfolder, threshold)
	inFiles = len(glob.glob(f"{tmpfolder}/filtered/*.xyz"))
	print(f"\t\t\tTotal of Complexes After The Integrity Test:: {inFiles}")

	files=glob.glob(f"{tmpfolder}/filtered/*.xyz")
	if not files:
		print("\n\nERROR: Complexes did not pass the integrity test\n\n")
		exit()

	runZero = params["MOD3"]["RUN_MOD_ZERO"]
	if runZero:
		print("\n\t\tSelecting representative complexes (k-means)")
		M0_Selection(params["MOD3"], False)
		filfold = tmpfolder+"/selected"
	else:
		filfold = tmpfolder+"/filtered"

	print(f"\n\t\tCopying final files to the output folder")
	files=glob.glob(f"{filfold}/*.xyz")
	for file in tqdm.tqdm(files):
		os.system(f"cp {file} {outfolder}/")
	# os.system(f"cp {filfold}/*.xyz {outfolder}/")	

	print(f"\n\t\tComplexes are available at {outfolder}")
	print("\nEnd of module 03 - Complexes Generation")

def main() -> None:
	"""
	Main routines
	"""
	if start_message() == False:
		exit("---DONE---")
	
	params = load_parameters()

	if 0 in params["MODULES"]:
		M0_Selection(params["MOD0"], True)

	if 1 in params["MODULES"]:
		M1_frame_family(params)

	if 2 in params["MODULES"]:
		M2_core_family(params)

	if 3 in params["MODULES"]:
		M3_add_ligants(params)

	print("\n\nE N D !")

if __name__ == "__main__":
    main()

