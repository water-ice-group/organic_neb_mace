from ase.calculators.orca import ORCA, OrcaProfile
import argparse
from ase.io import read, write
import os
from ase.optimize import LBFGS

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="path to input configs", required=True)
    parser.add_argument("--output", help="path to save output configs to", required=True)
    parser.add_argument("--charge", help="charge", default=0, type=int)
    parser.add_argument("--mult", help="spin multiplicity", default=1, type=int)
    parser.add_argument("--task", help="task", required=True)
    parser.add_argument('--theory', help='level of theory', default='RI-MP2 6-311G(d) autoaux nofrozencore tightscf engrad')
    parser.add_argument("--index", help="which config to take", default=None)
    return parser.parse_args()

def get_calc(theory, label, charge=0, mult=1):
    calc = ORCA(
        directory=label,
        profile=OrcaProfile('/rds/user/dk584/hpc-work/programs/orca_5_0_3_linux_x86-64_shared_openmpi411/orca'),
        charge=charge, mult=mult,task='gradient',
        orcasimpleinput=theory,
        orcablocks='%pal nprocs 1 end'
    )
    return calc

def single_point(config, charge, mult, theory, prefix="mp2_"):
    calc = get_calc(charge, mult, theory)
    calc.calculate(atoms=config, properties={"energy", "forces"}, system_changes=None)
    os.system("rm orca_property.txt orca.*")
    config.info.update(
        {
            f"{prefix}energy": calc.results["energy"]
        }
    )
    config.arrays.update(
        {
            f"{prefix}forces": calc.results["forces"]
        }
    )
    return config

def optimize(config, charge, mult, out, theory, fmax=0.01, nsteps=100):
    calc = get_calc(theory, charge, mult)
    config.calc = calc
    opt = LBFGS(config)
    opt.run(fmax=fmax, steps=nsteps)
    write(out, opt.atoms)

def main():
    args = parse_args()
    index = args.index
    configs = read(args.input, ":")
    if args.task == "single_point":
        for config in configs:
            config = single_point(config, charge=args.charge, mult=args.mult, theory=args.theory)
            write(args.output, config, append=True)
        os.system("rm orca.* orca_property.txt")
    elif args.task == "optimize":
        optimize(
            configs[0], 
            charge=args.charge, 
            mult=args.mult,
            out=args.output,
            theory=args.theory
        )

if __name__ == "__main__":
    main()
