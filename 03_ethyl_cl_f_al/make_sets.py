from ase.io import read, write
import os
import argparse
from orca_eval import single_point

def main(args):
    reactants = read(args.reactants_input)
    products = read(args.products_input)
    e0s = read(args.E0s_file, ':')

    r_rattled = []
    p_rattled = []

    for i in range(args.ntrain + args.nvalid + args.ntest):
        atoms = reactants.copy()
        atoms.rattle(stdev=args.stdev, seed=2*i)
        atoms = single_point(orca_path=args.orca_path, config=atoms, charge=args.charge, mult=args.mult)
        r_rattled.append(atoms)

        atoms = products.copy()
        atoms.rattle(stdev=args.stdev, seed=2*i+1)
        atoms = single_point(orca_path=args.orca_path, config=atoms, charge=args.charge, mult=args.mult)
        p_rattled.append(atoms)

    train = e0s + r_rattled[:args.ntrain] + p_rattled[:args.ntrain]
    valid = r_rattled[args.ntrain:args.ntrain + args.nvalid] + p_rattled[args.ntrain:args.ntrain + args.nvalid]
    test = r_rattled[args.ntrain + args.nvalid:] + p_rattled[args.ntrain + args.nvalid:]
    
    os.makedirs(args.output, exist_ok=True)
    write(f'{args.output}/train.xyz', train)
    write(f'{args.output}/valid.xyz', valid)
    write(f'{args.output}/test.xyz', test)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reactants_input", required=True)
    parser.add_argument("--products_input", required=True)
    parser.add_argument("--E0s_file", default='e0s.xyz')
    parser.add_argument("--ntrain", default=30, type=int)
    parser.add_argument("--nvalid", default=20, type=int)
    parser.add_argument("--ntest", default=20, type=int)
    parser.add_argument("--orca_path", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--charge", default=-1, type=int)
    parser.add_argument("--mult", default=1, type=int)
    parser.add_argument("--stdev", default=0.02, type=float)
    args = parser.parse_args()
    main(args)