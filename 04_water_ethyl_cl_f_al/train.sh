#!/bin/bash

python train.py \
    --reactants_fname='r.xyz' \
    --products_fname='p.xyz' \
    --mace_dir='/my/mace/directory' \
    --orca_path='/path/to/my/orca/executable' \
    --num_models=2 \
    --num_new_configs=1 \
    --neb_steps=250 \
    --num_images=40 \
    --k=0.04 \
    --neb_fmax=0.05 \
    --qm_fmax=10 \
    --charge=-1 \
    --mult=1 \
