![alt text](https://qtnano.iqsc.usp.br/wp-content/themes/qtnano/images/logo_foot2.png)
# `cluster_assembler` - QTNano's protected clusters (cluster-ligands) assembling procedure

[![DOI](https://zenodo.org/badge/520031576.svg)](https://zenodo.org/badge/latestdoi/520031576)

This packege is the workflow implementation used in [Living Libraries of Mixed Metal Clusters](https://doi.org/).
This package can be used to generate atomistic structures for metallic complexes that model the strutures presents in living libraries.
This is done via a commbination of DFT ab-initio calculations and data science techniques.


## Recomendation before use
Before running the code, the user should check the following:

1. The directory where you plan to run the workflow should be set up with the following files/folder:
 - `main.py`, the main python script to execute the workflow.
 - `mol.txt`, simple list of atomic numbers and simbles, used in the code. Can be adapted by the user if specific labels are used instead of traditional atomic species (e.g., `O1, O2, O3, ...` is used to represent diferet parametrizations for oxigen). This file MOST be in the same folder as the `main.py`
 - `core/`, folder containing geometries (XYZ files) for metalic cores, in case they are to be inported.  **DONT KNOW IF WE STILL ARE USING THE FILES THAT ARE THERE**
 - `ligands/`, folder containing the geometrie (XYZ files) for the ligands
 - `parameters.json`, input file for the workflow.
 - Opitionaly, files as `job_pbs` and `requirements.txt` should be provided to run `cluster_assembler` on HPC facilities within specified python modules. Sample submission job files for 
 personal computer (`job.sh`) and HPC clusters (`job.pbs`) are provided ere in the `jobs/` folder.
2. Check if you have a installation of [`Python 3`](https://python.org) working properly. In case you dont, we strongly incorage the use of [`Anaconda`](https://anaconda.org/).


## Editing `parameters.json`

The file `parameters.json` has to be set-up for each specific use case. This is done by changing the parameters according to the tasks that `cluster_assember` will perform. Here, a reference guide for the meaning of each input field in this file:

**... BELLOW NEEDS TO BE EXPLAINED**

**Main settings** 
* `MODULES` ... (sample value: `[1,2,3]`)
* `ELEM1` ... (sample value: `"Cu"`)
* `ELEM2` ... (sample value: `"Zn"`)
* `NUMELEM1` ... (sample value: `3`)
* `NUMELEM2` ... (sample value: `4`)

**MOD0**
* `COMMENTS` ... (sample value: `" "`)
* `KMEANS` ... (sample value: `[70]`)
* `RUN_KMEANS` ... (sample value: `10`)
* `MAXSAMPLES` ... (sample value: `1`)
* `INPUT_FOLDER` ... (sample value: `"./Mod02"`)
* `OUTPUT_FOLDER` ... (sample value: `"./Selected_kmeans_user"`) 

**MOD1**
* `COMMENTS` ... (sample value: `" "`)
* `RUN_MOD_ZERO` ... (sample value: `false`)
* `INTEGRITY` ... (sample value: `1.5`)
* `NUMGEN` ... (sample value: `120`)
* `KMEANS` ... (sample value: `[3, 9, 1]`)
* `RUN_KMEANS` ... (sample value: `10`)
* `MAXSAMPLES` ... (sample value: `1`)
* `RADIUS_FACTOR` ... (sample value: `0.5`)
* `GAMMA` ... (sample value: `0.2`)
* `INPUT_FOLDER` ... (sample value: `""`)
* `TMP_FOLDER` ... (sample value: `"./varmod1"`)
* `OUTPUT_FOLDER` ... (sample value: `"./Mod01"`)


**MOD2**
* `COMMENTS` ... (sample value: `""`)
* `MAX_GEN_PER_FRAME` ... (sample value: `35`)
* `PRE_SELECTION` ... (sample value: `false`)
* `RUN_MOD_ZERO` ... (sample value: `false`)
* `INTEGRITY` ... (sample value: `1.5`)
* `KMEANS` ... (sample value: `[3, 9, 1]`)
* `RUN_KMEANS` ... (sample value: `10`)
* `MAXSAMPLES` ... (sample value: `1`)
* `INPUT_FOLDER` ... (sample value: `"./Mod01"`)
* `TMP_FOLDER` ... (sample value: `"./varmod2"`)
* `OUTPUT_FOLDER` ... (sample value: `"./Mod02"`) 

**MOD3**
* `COMMENTS` ... (sample value: `"Ligands orientation ... (sample value: [-1] for random, [0 or 1] for diatomic molecules, [1,2,3] for selecting the atoms 1, 2 and 3 for orienting the ligant"`) 
* `RUN_MOD_ZERO` ... (sample value: `true`)
* `DEFORMATION` ... (sample value: `true`)
* `INTEGRITY` ... (sample value: `1.25`)
* `KMEANS` ... (sample value: `[6,20,3]`)
* `RUN_KMEANS` ... (sample value: `10`)
* `MAXSAMPLES` ... (sample value: `1`)
* `LIGANDS_DISTRIBUTION` ... (sample value: `[0,5]`)
* `LIGANDS_ORIENTATION` ... (sample value: `[[-1], [0,1,2]]`)
* `N_SAMPLES` ... (sample value: `100`)
* `N_CORES` ... (sample value: `100`)
* `LIGANDS_DISTANCE` ... (sample value: `2.25`)
* `CORES_FOLDER` ... (sample value: `"./Cu3Zn4_structures"`)
* `LIGANDS_FOLDER` ... (sample value: `"./ligands"`)
* `TMP_FOLDER` ... (sample value: `"./varmod3"`)
* `OUTPUT_FOLDER` ... (sample value: `"./Mod03"`)


We recommend that inexperienced users keep advanced parameters set to their `default` values. However, all parameters can be adjusted as needed. 

## Conda

## Running
After setting up the `parameters.json` file, the workflow can be started by running:

```
python main.py /adress/to/parameters.json
```

All folders and outputs are created in the same working directory where `python main.py` was called. Therefore, when the code or parameter file requires an address (e.g., in the `INPUT_FOLDER` field), the user only needs to type the file's name. However, if you want to use another folder,  providing the complete path for that folder will suffice.

## List of Modules
### Module 0 - Mod0
This module is a *k*-means clustering tool to select the representatives unary cores, binary cores, and complexes. The original version of the *k*-means clustering of nanoclusters and molecules we based this module on can be accessed [here](https://github.com/quiles/Adsorption\_Clus/blob/main/representative.py). **MORE DETAILS?**

*NOTE:* The use of this module has `pandas` as prerequisite. However, if you are using Anaconda, this package is tipicaly included in your base enviropment.


### Module 1 - Mod1.
In this module, unary metal cores are generated, processed, analyzed, and clustered. We included a series of user options that allow combining external sourcing and internally generating the unary cores as needed by the case of study. The cores are checked according to a integrity filter and clustered by *k*-means if the `RUN_MOD_ZERO` flag is activated.

### Module 2 - Mod2.
Here, permutations (i.e., all possible or user-specified) are performed in the atomic species of each of the selected unary cores. The binary cores are subsequently checked for integrity and clustered by *k*-means if the `RUN_MOD_ZERO` flag is activated.

### Module 3 - Mod3.

Module responsable to distribute the ligands around the selected binary cores. The ligands XYZ are imported from `ligands/` folder, particularly when the ligand is oriented distributed. An overlap filter iteratively excludes structures with overlapping atoms. The default threshold is set to t=**??????**, and can be manually tuned it in the `/filters/overlapping.py` file if n ecessary. **ITS NOT IN THE PARAMETERS FILES????**. The filtered complexes are again clustered by *k*-means if the `RUN_MOD_ZERO` flag is activated. Finally, the selected structures are the representatives within the living library family of structures, and can be externally optimized using methods such as Density Functional Theory.

*NOTE:* Energy calculations must be done externally and can be automated using appropriate quantum chemistry software such as FHI-aims, CREST, VASP, or ADF.




---

#**DIDN'T CHANGED AFTER THAT** 



## Sample inputs and outputs

CuZn2 + Cp*3


We provided a ZIP file (`Cu3Zn4Cp5.zip`, in two parts) with a folder that contains both the code, the inputs, and the outputs for the test case of Cu3Zn4Cp5 complex. The folders are already named as the default values at `main.py`. You can directly execute the `main.py` on that folder and follow the default options and the given INPUT_Cu3Zn4Cp5 to reproduce the data. 
In Step 1, the pre-generated structures can be found in the folder `all_xyz`. 

Additionally, we provided another ZIP file (`sample_step_9.zip`) for testing the advanced option for Step 9. The folder contains a sample of `data_2`and the `/final_representatives` for the test case of Cu3Zn4Cp5. The result for Step 10 for 15 representatives Cu3Zn4Cp5 structures is also available in the `\to_dft_tight` folder.  

The code was originally tested in Ubuntu 20.04 system but should be compatible with any Linux distribution. The time-consuming part of the workflow is the quantum chemistry calculation, i.e., Steps 7 and 10, which are related to the cluster nature, HPC facilities, and the used quantum chemistry program package.


*HAVE FUN and, if our codes were useful for you, remember to cite our work!*

