#!/bin/bash


python3 train.py \
    --reactants_fname='r.xyz' \
    --products_fname='p.xyz' \
    --mace_dir='/my/mace/directory' \
    --orca_path='/path/to/my/orca/executable' \
    --num_models=2 \
    --num_new_configs=5 \
    --neb_steps=250 \
    --num_images=21 \
    --k=0.2 \
    --neb_fmax=0.05 \
    --qm_fmax=10 \
    --max_n_iterations=50 \
    --charge=0 \
    --mult=1 \
    --extra_neb_steps \
    --use_idpp \
    --use_climb \
    --neb_method="improvedtangent" \
    --num_channels=64

