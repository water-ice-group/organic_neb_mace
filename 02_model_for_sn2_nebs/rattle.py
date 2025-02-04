from ase.io import read, write
import os
import argparse

def main(args):
    
    reactants = read(args.reactants_input, ':')
    products = read(args.products_input, ':')
    transition_states = read(args.ts_input, ':')

    """ If stdev is 0, simply copy the R, P, TS to the outpur dir """
    
    if args.stdev == 0:
        r_rattled = reactants
        p_rattled = products
        ts_rattled = transition_states
    else:
        r_rattled = []
        p_rattled = []
        ts_rattled = []

        seed = 0

        for at in reactants:
            at2 = at.copy()
            at2.rattle(stdev=args.stdev, seed=seed)
            r_rattled.append(at2)
            seed += 1

        for at in products:
            at2 = at.copy()
            at2.rattle(stdev=args.stdev, seed=seed)
            p_rattled.append(at2)
            seed += 1

        for at in transition_states:
            at2 = at.copy()
            at2.rattle(stdev=args.stdev, seed=seed)
            ts_rattled.append(at2)
            seed  += 1

    os.makedirs(args.output, exist_ok=True)
    write(f'{args.output}/r.xyz', r_rattled)
    write(f'{args.output}/p.xyz', p_rattled)
    write(f'{args.output}/ts.xyz', ts_rattled)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reactants_input", required=True)
    parser.add_argument("--products_input", required=True)
    parser.add_argument("--ts_input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--stdev", default=0.02, type=float)
    args = parser.parse_args()
    main(args)