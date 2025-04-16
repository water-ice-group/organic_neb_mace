#!/bin/bash


python3 ../../../scripts/train.py \
    --reactants_fname='r.xyz' \
    --products_fname='p.xyz' \
    --train_num_rattle=30 \
    --valid_num_rattle=30 \
    --test_num_rattle=30 \
    --rattle_stdev=0.02
    --charge=0 \
    --mult=1 \
    --orca_path='/path/to/my/orca/executable' \
    --mace_dir='/my/mace/directory' \
    --num_models=2 \
    --num_new_configs=5 \
    --neb_steps=250 \
    --num_images=21 \
    --k=0.2 \
    --neb_fmax=0.05 \
    --qm_fmax=10 \
    --max_n_iterations=25
    --extra_neb_steps \
    --use_idpp \
    --use_climb \
    --neb_method='improvedtangent'
    --mace_hidden_irreps='64x0e + 64x1o'
    --E0s="{1: -13.600517822134075, 6: -1027.5647719688623, 7: -1482.8258826434405, 8: -2039.1187770138454, 16: -10821.988689368247}"