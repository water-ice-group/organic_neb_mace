import numpy as np
from ase.io import read, write
import argparse
import os

def main(args):
    ntrain = args.ntrain
    nvalid = args.nvalid
    ntest = args.ntest

    indices = np.arange(ntrain + nvalid + ntest)

    np.random.seed(123)

    """ Randomly select which reactions go into train, valid, test splits """

    train_id = np.random.choice(indices, size=ntrain, replace=False)
    indices = np.setdiff1d(indices, train_id)

    valid_id = np.random.choice(indices, size=nvalid, replace=False)
    indices = np.setdiff1d(indices, valid_id)

    test_id = np.random.choice(indices, size=ntest, replace=False)
    indices = np.setdiff1d(indices, test_id)

    ''' Generate training set, which may contain the rattled geometries '''

    e0s = read(args.E0s_file, ':')
    reactants = read(args.reactants_input, ':')
    products = read(args.products_input, ':')
    transition_states = read(args.ts_input, ':')

    train = e0s \
        + [reactants[i] for i in train_id] \
        + [transition_states[i] for i in train_id] \
        + [products[i] for i in train_id]

    """ Generate valid, and test sets, which contain the equilibrium geometries (no added random displacements) """

    reactants = read('r.xyz', ':')
    products = read('p.xyz', ':')
    transition_states = read('ts.xyz', ':')

    valid = [reactants[i] for i in valid_id] \
        + [transition_states[i] for i in valid_id] \
        + [products[i] for i in valid_id]
        
    test = [reactants[i] for i in test_id] \
        + [transition_states[i] for i in test_id] \
        + [products[i] for i in test_id]

    os.makedirs(args.output, exist_ok=True)

    write(f'{args.output}/train.xyz', train)
    write(f'{args.output}/valid.xyz', valid)
    write(f'{args.output}/test.xyz', test)

    """ Write test set reactants, products, and TS separately for NEB calculations """

    write('r-test.xyz', [reactants[i] for i in test_id])
    write('ts-test.xyz', [transition_states[i] for i in test_id])
    write('p-test.xyz', [products[i] for i in test_id])


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ntrain", default=1422, type=int)
    parser.add_argument("--nvalid", default=178, type=int)
    parser.add_argument("--ntest", default=178, type=int)
    parser.add_argument("--reactants_input", required=True)
    parser.add_argument("--products_input", required=True)
    parser.add_argument("--ts_input", required=True)
    parser.add_argument("--E0s_file", default='e0s.xyz')
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    main(args)