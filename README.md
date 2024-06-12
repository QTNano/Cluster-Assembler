![alt text](https://qtnano.iqsc.usp.br/wp-content/themes/qtnano/images/logo_foot2.png)
# `Cluster-Assembler`

[![DOI](https://zenodo.org/badge/520031576.svg)](https://zenodo.org/badge/latestdoi/520031576)

This package can be used to generate atomistic structures for metallic complexes that model the strutures presents in living libraries.
This is done via a commbination of DFT ab-initio calculations and data science techniques.


## Table of contents

[Recomendations](#recomendations)

[Installation](#installation)

[Execution](#execution)

[Modules](#modules)

[Example](#example)

[How to cite](#how-to-cite)

## Recomendations
Before running the code, the user should check the following:

1. The directory where you plan to run the workflow should be set up with the following files/folder:
 - `main.py`, the main python script to execute the workflow.
 - `mol.txt`, simple list of atomic numbers and simbles, used in the code. Can be adapted by the user if specific labels are used instead of traditional atomic species (e.g., `O1, O2, O3, ...` is used to represent diferet parametrizations for oxigen). This file MOST be in the same folder as the `main.py`
 - `core/`, folder containing the core code of the toolbox
 - `parameters.json`, input file for the workflow.

2. We strongly encourage the use of [`Anaconda`](https://anaconda.org/).

## Installation

This tool utilizes Anaconda to set up a virtual environment where all required dependencies are managed, although the system can also be configured manually.

To begin, clone or download this repository and install [Anaconda](https://www.anaconda.com/) for your specific operating system. Once installed, navigate to the root directory of this repository using your terminal and execute the following command:

```shell
$ conda env create -f environment.yml
```
Executing this command will generate a new Anaconda virtual environment equipped with all the dependencies needed to run the Cluster-Assembler toolbox.

---

## Execution

### Starting conda

To execute the toolbox, start by activating the `py396` Anaconda environment:
```shell
$ conda activate py396
```
Next, launch the program with Python:
```shell
$ python3 main.py parameters.json
```

The `parameters.json` contains all parameters for the execution. This file is explained next.

---

### Toolbox parameters

The file `parameters.json` has to be set-up for each specific use case. This is done by changing the parameters according to the tasks that `Cluster-Assember` will perform. Here, a reference guide for the meaning of each input field in this file:

**Main settings** 

| Parameter | Type | Description | Default Value |
|--|--|--|--|
| `MODULES` | list of integers | Select which module (or modules) will be run | `[1, 2, 3]`|
| `ELEM1` | string | The first species of the clusters | `"Cu"` |
| `ELEM2` | string | The second species of the clusters | `"Zn"` |
| `NUMELEM1` | integer | Number of atoms of the first species | `3` |
| `NUMELEM2` | integer | Number of atoms of the second species | `4` |

**Common Parameters - Modules 0, 1, 2, and/or 3**

| Parameter | Type | Description | Default Value |
|--|--|--|--|
| `KMEANS` | list of integers | Define the number of seeds (k) of the k-means algorithm. A single value indicate the exact number of seeds; otherwise, the user can define a range and the system will automatically set the best number of seeds following a Silhouette analysis. In this last case, the first and the second number defines the range (e.g., from 2 to 50), and the last number defines the evaluation step  | `[10]` or `[2, 50, 1]`|
| `RUN_KMEANS` | integer | Defines the number of k-means runs | `10`|
| `MAXSAMPLES` | integer | Number of samples selected per cluster | `1`|
| `INPUT_FOLDER` | string | Input data folder | `./input_data`|
| `OUTPUT_FOLDER` | string | Output data folder | `./output_data`|
| `TMP_FOLDER` | string | Temporary data folder - storage of intermediate steps | `./var_data`|
| `RUN_MOD_ZERO` | boolean | Define if the module 0 will be run | `false`|

**MOD1 - Module 1 - Frame Family**
| Parameter | Type | Description | Default Value |
|--|--|--|--|
| `INTEGRITY` | float | Define the tolerance regarding the covalent radius allowed to form a structure, e.g., when set to 1, the covalent radius of the species is the limit to form a link.
 | `1.5`|
| `NUMGEN` | integer | Number of generated frames | `100`|
| `RADIUS_FACTOR` | float | Parameter that controls the size of the box/sphere | `0.5`|
| `GAMMA` | float | Tolerance regarding the species covalent radius | `0.2`|


**MOD2 - Module 2 - Core Family**
| Parameter | Type | Description | Default Value |
|--|--|--|--|
| `MAX_GEN_PER_FRAME` | integer | Number of cores to be generated for each input frame | `100`|
| `INTEGRITY` | float | Define the tolerance regarding the covalent radius allowed to form a structure | `1.5`|


**MOD3 - Module 3 - Complexes Generation**
| Parameter | Type | Description | Default Value |
|--|--|--|--|
| `DEFORMATION` | boolean | If true, the sites will be adjusted to the surface of the core | `true`|
| `INTEGRITY` | float | Define the tolerance regarding the covalent radius allowed to form the complexes | `1.25`|
| `LIGANDS_DISTRIBUTION` | list of integers | Define the ligands and their quantities to generate the complexes. The first number indicates the number of times that ligand 1 will be added, the second number the second ligand, and so on. The ligands are selected by alphabetical order of files in the folder informed in `LIGANDS_FOLDER`  | `[3,2]`|
| `LIGANDS_ORIENTATION` | list of lists | Ligands orientation: sample value: [-1] for random, [0 or 1] for diatomic molecules, [1,2,3] for selecting the atoms 1, 2 and 3 for orienting the ligant | `[[-1], [0,1,2]]`|
| `N_SAMPLES` | integer | Number of complexes that will be generated for each input core | `100`|
| `N_CORES` | integer | Number of cores that will be used in the simulation. If the number of cores are superior to the number of files presented in `CORES_FOLDER`, the number of existant files will be used | `3`|
| `LIGANDS_DISTANCE` | float | Distance between the ligand and the core | `2.25`|
| `CORES_FOLDER` | string | Folder containing the cores (XYZ) files | `./cores`|
| `LIGANDS_FOLDER` | string | Folder containing the ligands (XYZ) files | `./ligands`|


We recommend that inexperienced users keep advanced parameters set to their `default` values. However, all parameters can be adjusted as needed. 

## Modules
### Module 0 - MOD0
This module is a *k*-means clustering tool to select the representatives unary cores, binary cores, and complexes. This module uses the eigenvalues of the Coulomb matrix of the provided structures to select a set of representative ones.


### Module 1 - MOD1.
In this module, unary cores are generated, processed, analyzed, and clustered. We included a series of user options that allow combining external sourcing and internally generating the unary cores as needed by the case of study. The cores are checked according to a integrity filter (covalent radius) and clustered by *k*-means if the `RUN_MOD_ZERO` flag is activated.

### Module 2 - MOD2.
Here, permutations (i.e., all possible or user-specified) are performed in the atomic species of each of the selected unary cores. The binary cores are subsequently checked for integrity and clustered by *k*-means if the `RUN_MOD_ZERO` flag is activated.

### Module 3 - MOD3.

Module responsable to distribute the ligands around the selected binary cores. The ligands XYZ are imported from `ligands/` folder. Several ligands might be used in the process. Their orientation regarding the core is controled by the user, which can select a random  or specific orientation (provided by a set of base atoms). An overlap filter iteratively excludes structures with overlapping atoms (integrity). The default threshold is set to t=`1.25`, and can be manually tuned. The filtered complexes are again clustered by *k*-means if the `RUN_MOD_ZERO` flag is activated. Finally, the selected structures are the representatives within the living library family of structures, and can be externally optimized using methods such as Density Functional Theory.

*NOTE:* Energy calculations must be done externally and can be automated using appropriate quantum chemistry software such as FHI-aims, CREST, VASP, or ADF.


## Example

---

#**DIDN'T CHANGED AFTER THAT** 



## Sample inputs and outputs

CuZn2 + Cp*3


We provided a ZIP file (`Cu3Zn4Cp5.zip`, in two parts) with a folder that contains both the code, the inputs, and the outputs for the test case of Cu3Zn4Cp5 complex. The folders are already named as the default values at `main.py`. You can directly execute the `main.py` on that folder and follow the default options and the given INPUT_Cu3Zn4Cp5 to reproduce the data. 
In Step 1, the pre-generated structures can be found in the folder `all_xyz`. 

Additionally, we provided another ZIP file (`sample_step_9.zip`) for testing the advanced option for Step 9. The folder contains a sample of `data_2`and the `/final_representatives` for the test case of Cu3Zn4Cp5. The result for Step 10 for 15 representatives Cu3Zn4Cp5 structures is also available in the `\to_dft_tight` folder.  

The code was originally tested in Ubuntu 20.04 system but should be compatible with any Linux distribution. The time-consuming part of the workflow is the quantum chemistry calculation, i.e., Steps 7 and 10, which are related to the cluster nature, HPC facilities, and the used quantum chemistry program package.


*HAVE FUN and, if our codes were useful for you, remember to cite our work!*

## How to cite