#!/bin/bash

python /path_to_mace/mace/mace/cli/run_train.py \
    --name="MACE" \
    --train_file="../train.xyz" \
    --valid_file='../valid.xyz' \
    --test_file="../test.xyz" \
    --model="MACE" \
    --config_type_weights='{"Default":1.0}' \
    --hidden_irreps='128x0e + 128x1o' \
    --r_max=5.0 \
    --batch_size=32 \
    --max_num_epochs=1000 \
    --swa \
    --start_swa=500\
    --ema \
    --ema_decay=0.99 \
    --amsgrad \
    --restart_latest \
    --device=cuda \
    --save_cpu \
    --seed=1 \
    --forces_weight=10.0 \
    --energy_weight=1.0 \
    --swa_forces_weight=100.0 \
    --swa_energy_weight=1000.0 \
    --lr=0.001 \
    --forces_key="mp2_forces" \
    --energy_key="mp2_energy" \
