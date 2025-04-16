#!/bin/bash


python ../../scripts/train.py \
    --reactants_fname='r.xyz' \
    --products_fname='p.xyz' \
    --train_num_rattle=30 \
    --valid_num_rattle=30 \
    --test_num_rattle=30 \
    --rattle_stdev=0.02
    --charge=-1 \
    --mult=1 \
    --orca_path='/path/to/my/orca/executable' \
    --mace_dir='/my/mace/directory' \
    --num_models=2 \
    --num_new_configs=1 \
    --neb_steps=250 \
    --num_images=40 \
    --k=0.04 \
    --neb_fmax=0.05 \
    --qm_fmax=10 \
    --max_n_iterations=20
    --neb_method='aseneb'
    --mace_hidden_irreps='128x0e + 128x1o'
    --E0s="{1: -13.600517822134075, 6: -1027.5647719688623, 8: -2039.1187770138454, 9: -2711.1240911830328, 17: -12510.40325874303}"