#!/bin/bash

python3 train.py \
  --train_fname='train.xyz' \
  --test_fname='test.xyz' \
  --train_eval_fname='train.xyz' \
  --test_eval_fname='test.xyz' \
  --mace_dir='~/my_mace_dir/mace' \