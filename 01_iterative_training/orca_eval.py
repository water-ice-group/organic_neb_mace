from ase.calculators.orca import ORCA, OrcaProfile
import argparse
from ase.io import read, write
import os
from ase import Atoms

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="path to input configs", required=True)
    parser.add_argument("--ids", help="which configuration(s) to take", default="0")
    parser.add_argument("--output", help="path to save output configs to", required=True)
    parser.add_argument("--charge", help="charge", default=0, type=int)
    parser.add_argument("--mult", help="spin multiplicity", default=1, type=int)
    parser.add_argument("--orca_path", help="path to orca installation", required=True)
    return parser.parse_args()

def get_calc(orca_path, charge=0, mult=1):
    calc = ORCA(
        profile=OrcaProfile(orca_path),
        label="orca",
        charge=charge, mult=mult,task='gradient',
        orcasimpleinput='RI-MP2 6-311G(d) autoaux nofrozencore TightSCF engrad',
        orcablocks='%pal nprocs 1 end \n %scf ConvForced=1 end'
    )
    return calc

def single_point(config, orca_path, charge, mult, prefix="mp2_"):
    calc = get_calc(orca_path, charge, mult)
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


def main():
    args = parse_args()
    configs = read(args.input, ":")
    orca_path = args.orca_path

    if isinstance(configs, list):
        for config in configs:
            config = single_point(config, orca_path, args.charge, args.mult)
        write(args.output, configs)
        os.system("rm orca.* orca_property.txt")

    elif isinstance(configs, Atoms):
        configs = single_point(configs, orca_path, args.charge, args.mult)
        write(args.output, configs)
        os.system("rm orca.* orca_property.txt")

if __name__ == "__main__":
    main()
