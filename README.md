# Scripts and data for "Finding reaction pathways efficiently using reaction databases and active learning".

# Table of Contents
1. [Iterative training](#iterative-training)
2. [Model for finding SN2 reaction NEBs](#model-for-finding-sn2-reaction-nebs)
3. [Active learning on ethyl chloride and fluoride ion SN2](#active-learning-on-ethyl-chloride-and-fluoride-ion-sn2)
4. [Active learning on microsolvated ethyl chloride and fluoride ion SN2](#active-learning-on-microsolvated-ethyl-chloride-and-fluoride-ion-sn2)
5. [Aromatic substitution with AL and ReaxFF](#aromatic-substitution-with-al-and-reaxff)
6. [Figures](#figures)


## Iterative training

Directory `01_iterative_training` contains the initial training set (reactants only) `train.xyz` along with isolated atom energies. The test set `test.xyz` contains transition states that are transferred to the training set in small batches during iterative training.

The script `orca_eval.py` is used to evaluate the MP2 energies and forces and contains the ORCA calculation details.

Iterative training is run using the `train.py` script with supplied arguments. `train.py` with example arguments can be run using `train.sh`

## Model for finding SN2 reaction NEBs

The files `r.xyz`, `p.xyz`, and `ts.xyz` contain the respective reactant, product, and transition state structures of 1,778 reactions. The MP2 energies and forces are also present in the xyz files.

Use the `make_sets.py` script to randomly select reactions for the train, test and valid splits. In `make_sets.sh`, we give commands for the sets used here

For the rattled training sets, we need to make R, P, TS configurations with added random atomic displacements using the `rattle.py` script and recalculate their MP2 energies and forces using `orca_eval.py`

## Active learning on ethyl chloride and fluoride ion SN2

The reactant and product configurations are stored in `r.xyz` and `p.xyz`, respectively.

Use the `make_sets.py` script to make rattled configurations of the reactants and products, calculate their MP2 energies and forces, and make train, valid, test sets. `train.sh` has the example inputs for `train.py` to generate the sets.

To run the active learning, use `train.sh` which contains example input arguments for the `train.py` script

## Active learning on microsolvated ethyl chloride and fluoride ion SN2

The reactant and product configurations are stored in `r.xyz` and `p.xyz`, respectively.

Use the `make_sets.py` script to make rattled configurations of the reactants and products, calculate their MP2 energies and forces, and make train, valid, test sets. `train.sh` has the example inputs for `train.py` to generate the sets.

To run the active learning, use `train.sh` which contains example input arguments for the `train.py` script

## Aromatic substitution with AL and ReaxFF

In addition to active learning scripts, we provide the ReaxFF force field parameters used "ffield.reax.hcnos", the ReaxFF retraining code, and input files for running the ReaxFF NEB

## Figures

Scripts to make figures used in the manuscript