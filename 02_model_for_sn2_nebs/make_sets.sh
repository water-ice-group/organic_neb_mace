#!/bin/bash

python3 make_sets.py \
  --ntrain=1422 \
  --nvalid=178 \
  --ntest=178 \
  --reactants_input='r.xyz' \
  --products_input='p.xyz' \
  --ts_input='ts.xyz' \
  --E0s_file='e0s.xyz' \
  --output='no-rattle'

python3 make_sets.py \
  --ntrain=1422 \
  --nvalid=178 \
  --ntest=178 \
  --reactants_input='rattled-stdev-0.02/r-mp2.xyz' \
  --products_input='rattled-stdev-0.02/p-mp2.xyz' \
  --ts_input='rattled-stdev-0.02/ts-mp2.xyz' \
  --E0s_file='e0s.xyz' \
  --output='rattled-stdev-0.02/model'

python3 make_sets.py \
  --ntrain=1422 \
  --nvalid=178 \
  --ntest=178 \
  --reactants_input='rattled-stdev-0.05/r-mp2.xyz' \
  --products_input='rattled-stdev-0.05/p-mp2.xyz' \
  --ts_input='rattled-stdev-0.05/ts-mp2.xyz' \
  --E0s_file='e0s.xyz' \
  --output='rattled-stdev-0.05/model'