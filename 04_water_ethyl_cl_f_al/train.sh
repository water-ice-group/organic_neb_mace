#!/bin/bash

python train.py \
    --reactants_fname='r.xyz' \
    --products_fname='p.xyz' \
    --mace_dir='mace' \
    --orca_path='orca/orca' \
    --n_comm=2 \
    --num_new_configs=1 \
    --neb_steps=250 \
    --n_images=40 \
    --k=0.04 \
    --neb_fmax=0.05 \
    --qm_fmax=10 \
    --max_n_iterations=40 \
    --charge=-1 \
    --mult=1 \
