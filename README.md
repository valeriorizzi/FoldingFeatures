# README

This repository contains two Python scripts to prepare PLUMED[^1] files for OneOPES[^2] simulations as outlined in the following paper:

> PAPER LINK

Specifically `bioinspired_features.py` generates features based on two unbiased simulations in PLUMED syntax, while `oneopes_files_generation.py` prepares the input files for OneOPES simulations.

## Requirements
The script requires Python3 with
- mdtraj
- PLUMED

PLUMED needs to be installed and sourced when the scripts are run. The version of PLUMED used for the development and testing is the [2.9.0](https://github.com/plumed/plumed2/releases/tag/v2.9.0), which we suggest for running both the simulations and the scripts to avoid compatibility issues. A minimum working environment can be installed in conda with the following

```
conda env create -f environment_biofeat.yml
```
If requested, `bioinspired_features.py` also produces a PyMOL session with the relevant hydrogen bonds and side chain contacts features for visualization. A minimum working environment with PyMOL can be installed in conda with the following
```
conda env create -f environment_biogfeat_pymol.yml
```

## Basic usage of bioinspired_features.py
### Required flags
The script for feature extraction `bioinspired_features.py` requires a few mandatory flags that can be checked with
```
python bioinspired_features.py --h
```
Six files are requested, namely
- xtc trajectory of the protein folded (`-F`)
- xtc trajectory of the protein unfolded (`-U`)
- pdb file of the full solvated box (`-r`)
- pdb file of the folded protein (`-rp`)
- pdb file of the alpha carbons of the protein (`-rca`)
- mcfile containing the masses and the charges of the system (`-mc`). This can be generated with the PLUMED function [DUMPMASSCHARGE](https://www.plumed.org/doc-v2.9/user-doc/html/_d_u_m_p_m_a_s_s_c_h_a_r_g_e.html)

Note that the reference file for the protein (`-rp`) will be the atoms considered for the features. It can be only part of a bigger system, but must keep the same atom index as the topology of the trajectory. 

The script produces a series of PLUMED files to analyze the trajectories. As such, the folded and unfolded trajectories should be in different directories so that these files do not overwrite each other. An archetypical directory organization may look like the following
```
parent_directory
│   bioinspired_features_generation.py
│   folded_protein.pdb
│   folded_protein_CA.pdb
|   full_solvated_box.pdb
|   mcfile
|
└───folded_trajectory
│   │   folded.xtc
|
└───unfolded_trajectory
    │   unfolded.xtc
```
With respect to this organization, the basic command to run the script would look like the following
```
python bioinspired_features_generation.py -F ./folded_trajectory/folded.xtc -U ./unfolded_trajectory/unfolded.xtc -r full_solvated_box.pdb -rp folded_protein.pdb -rca folded_protein_CA.pdb -mc mcfile
```
An example of files and directory organization is provided in [example_TRP](https://github.com/heritiem/bioinspired_features/tree/main/example_TRP) for tryptophan-cage, which is the same as the one used in the publication. At the end of the analysis, the script sumarizes the main run details and the features in the `bioinspired_features.log` file. 

### Optional flags
The `bioinspired_features_generation.py` script has a number of optional flags, namely
- For more control over the filtering process, both the cutoff (default: 0.6 nm) and the lda threshold values (default: 0.3) can be changed with the flags `-c` and `-l`, respectively
- The flag `-py` produces a PyMOL session which contains all the detected features (h-bonds and contacts)
- The flag `-e` produces more explicit PLUMED files. These are easier to read (for humans) but result in a much slower simulation if used for running biased MD
- The flag `-y` run the script without asking any confirmation to the user
- The flag `-s` sets the stride for the analysis of the trajectories (defaul: 10). We recommend to include at least 1000 frames to obtain good statistics.

## Basic usage of oneopes_files_generation.py
### Required flags
The script for OneOPES files generation `oneopes_files_generation.py` requires a few mandatory flags that can be checked with
```
python oneopes_files_generation.py --h
```
This script is supposed to be used after running the features generation script `bioinspired_features.py`, as it depends on a few output files of this script, namely `plumed_final_noduplicate.dat`, generated where the feature generation script is run, and `COLVAR_diff` and `COLVAR_solvation` for both the folded and unfolded trajectory, generated inside the corresponding directories of the trajectories. It is suggested to use these scripts in the same directory.

The mandatory flags are the following
- the directory containing the file plumed_final_noduplicate.dat (`-s`)
- the directory containing COLVAR_diff and COLVAR_solvation for the folded trajectory (`-f`)
- the directory containing COLVAR_diff and COLVAR_solvation for the unfolded trajectory (`-u`)
- the seven minimum (`-minTs`) and maximum (`-maxTs`) temperatures for the replica exchange for replicas 1 to 7
- the PACE for the two main CVs, in units of simulation steps (`-p`)
- the BARRIER for the two main CVs in kJ/mol (`-b`)
- pdb file of the folded protein (`-rp`)
- pdb file of the alpha carbons of the protein (`-rca`)

Further information about the OPES-specific PLUMED flags can be found on the corresponding [PLUMED manual page](https://www.plumed.org/doc-v2.9/user-doc/html/_o_p_e_s__m_e_t_a_d__e_x_p_l_o_r_e.html) and in-depth interpretations are reported in the corresponding publications.[^2][^3][^4]

The output files are organised in the following tree
```
oneopes_files
│   folded_protein.pdb
│   folded_protein_CA.pdb
|
└───rep0
│   │   plumed.dat
└───rep1
│   │   plumed.dat
└───rep2
│   │   plumed.dat
└───rep3
│   │   plumed.dat
└───rep4
│   │   plumed.dat
└───rep5
│   │   plumed.dat
└───rep6
│   │   plumed.dat
└───rep7
    │   plumed.dat
```
The `plumed.dat` files for the eight replicas reference the `.pdb` files in the parent directory. Please note that the `plumed.dat` files are different and should not be mixed between directories. To run a simulation, once a `.tpr` file is compiled with GROMACS, a copy of it should be put in each `rep?` directory. The simulation can then be run for example with the following command
```
mpirun -n 8 gmx_mpi mdrun -deffnm prod -s prod.tpr -multidir ./rep? -plumed plumed.dat -replex 10000 -hrex
```
There is no default value for the `-replex` flag. We suggest to use twice the PACE set for the main CVs.

### Optional flags
The `oneopes_files_generation.py` script has a number of optional flags, namely
- The name of the output directory to contain the OneOPES files (`--outdir`, default: oneopes_files)
- The PACE for auxiliary CVs, in units of simulation steps (`--pace_minor`, default: 100000)
- The BARRIER for auxiliary CVs in kJ/mol (`--barrier_minor`, default: 3)
- The STRIDE for COLVAR writing in output (`--stride`, default: 500)
- The name of the output file containing the collective variables (`--colvar`, default: COLVAR)
- The PACE for OPES multithermal bias update, in units of simulation steps (`--opesx_pace`, default: 1000)
- The UPDATE_FROM for starting the OPES multithermal bias application, in units of simulation steps (`--opesx_update`, default: 0)
- The OBSERVATION_STEPS for starting the OPES multithermal observation for bias construction, in PACE units (`--opesx_obs`, default: 100)

Further information about the OPES multithermal PLUMED flags can be found on the corresponding [PLUMED manual page](https://www.plumed.org/doc-v2.9/user-doc/html/_o_p_e_s__e_x_p_a_n_d_e_d.html) and in-depth interpretations are reported in the corresponding publications.[^2][^3][^4][^5]

For reweighting the simulations, we suggest the [original script](https://github.com/invemichele/OPES-explore/blob/main/postprocessing/FES_from_Reweighting.py) from the OPES papers or the [modified version](https://github.com/obzehn/multithermal_fes) that supports also the multithermal reweighting.

## References
[^1]: Tribello, G. A., et al. "PLUMED 2: New feathers for an old bird." Computer physics communications 185.2 (2014): 604-613. [DOI:10.1016/j.cpc.2013.09.018](https://doi.org/10.1016/j.cpc.2013.09.018)
[^2]: Rizzi, V., et al. "OneOPES, a combined enhanced sampling method to rule them all." Journal of Chemical Theory and Computation 19.17 (2023): 5731-5742. [DOI:10.1021/acs.jctc.3c00254](https://doi.org/10.1021/acs.jctc.3c00254)
[^3]: Invernizzi, M., and M. Parrinello. "Rethinking metadynamics: from bias potentials to probability distributions." The journal of physical chemistry letters 11.7 (2020): 2731-2736. [DOI:10.1021/acs.jpclett.0c00497](https://doi.org/10.1021/acs.jpclett.0c00497)
[^4]: Invernizzi, M., and M. Parrinello. "Exploration vs convergence speed in adaptive-bias enhanced sampling." Journal of Chemical Theory and Computation 18.6 (2022): 3988-3996. [DOI:10.1021/acs.jctc.2c00152](https://doi.org/10.1021/acs.jctc.2c00152)
[^5]: Invernizzi, M., P. M. Piaggi, and M. Parrinello. "Unified approach to enhanced sampling." Physical Review X 10.4 (2020): 041034. [DOI:10.1103/PhysRevX.10.041034](https://doi.org/10.1103/PhysRevX.10.041034)
