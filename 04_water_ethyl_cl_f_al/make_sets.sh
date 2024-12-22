#!/bin/bash

python3 make_sets.py \
  --reactants_input='r.xyz' \
  --products_input='p.xyz' \
  --E0s_file='e0s.xyz' \
  --ntrain=30 \
  --nvalid=20 \
  --ntest=20 \
  --orca_path='my_orca_path' \
  --output='.' \
  --charge=-1 \
  --mult=1 \
  --stdev=0.02