from ase.io import read, write
import numpy as np
from ase.optimize import FIRE, LBFGS
from ase.constraints import FixedLine, FixAtoms, FixInternals
from ase.calculators.orca import ORCA, OrcaProfile
import argparse

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="path to input configs", required=True)
    parser.add_argument("--id", help="path to input configs", required=True, type=int)

    return parser.parse_args()


args = parse_args()

atoms = read(args.input, args.id)

calc = ORCA(
    profile=OrcaProfile('/rds/user/dk584/hpc-work/programs/orca_5_0_3_linux_x86-64_shared_openmpi411/orca'),
    charge=-1, mult=1,task='gradient',
    orcasimpleinput='RI-MP2 6-311G(d) autoaux TightSCF nofrozencore engrad',
    orcablocks='%pal nprocs 1 end \n %scf ConvForced=1 end'
)

calc.calculate(atoms=atoms, properties=['energy', 'forces'], system_changes=[])

atoms.info.update({'mp2_energy': calc.results['energy']})
atoms.arrays.update({'mp2_forces': calc.results['forces']})

output_fname = f'mp2.xyz'
 

write(output_fname, atoms)
