import numpy as np
from ase.mep.neb import NEB, NEBOptimizer
from ase.optimize import FIRE, BFGS
from orca_dft import get_calc
from ase.io import read, write
import argparse
import ase

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--a", help="path to XYZ configurations", required=True)
    parser.add_argument("--b", help="path to XYZ configurations", required=False)
    parser.add_argument("--output", help="output path", default="neb.xyz")
    parser.add_argument(
        "--nimages", help="number of images", type=int, default=50
    )
    parser.add_argument(
        "--fmax", help="max force to stop neb", type=float, default=0.05
    )
    parser.add_argument(
        "--k", help="spring constant", type=float, default=1.0
    )
    parser.add_argument(
        "--nsteps", help="number of steps", type=int, default=100
    )
    parser.add_argument("--idpp", help="use idpp or not", type=bool, default=False)
    parser.add_argument('--theory', help='level of theory', default='RI-MP2 6-311G(d) autoaux nofrozencore tightscf engrad')


    return parser.parse_args()

def committee_error(dyn, selector=str, reg=0.2):
    if isinstance(dyn.atoms, NEB):
        images = dyn.atoms.images[1:-1]
        errors = []
    else:
        pass

    selector_types = ["energy", "forces"]

    if selector == "energy":
        errors = [np.sqrt(image.calc.results["energy_var"]) for image in images]
    elif selector == "forces":
        for i, image in enumerate(images):
            f_err = np.sqrt(np.sum(np.var(image.calc.results["forces_comm"], axis=0), axis=1))
            f_mean = np.linalg.norm(image.calc.results["forces"], axis=1)
            errors[i] = np.max(f_err / (f_mean + reg))
    else:
        raise ValueError(f"Incorrect selector value. Must be one of {selector_types}.")
    dyn.atoms.errors = np.array(errors)
    print(np.max(dyn.atoms.errors))
def write_neb(dyn, out_file):
    write(out_file, dyn.atoms.images)
""" User inputs """
args = parse_args()

nimages = args.nimages
nsteps = args.nsteps
fmax = args.fmax
k = args.k

a = read(args.a, ":")
if len(a) == 1:
    a = a[0]
if isinstance(a, list):
    images = a
    neb = NEB(images, k=k, climb=False, parallel=True, remove_rotation_and_translation=True)
elif isinstance(a, ase.Atoms):
    b = read(args.b)
    images = [a]
    images += [a.copy() for i in range(nimages-2)]
    images += [b]
    neb = NEB(images, k=k, climb=False, parallel=True, remove_rotation_and_translation=True)
    # Interpolate the potisions of the three middle images:
    if args.idpp:
        neb.interpolate("idpp")
    else:
        neb.interpolate()



# Set calculators:

for i, image in enumerate(images[1:-1]):
    image.calc = get_calc(args.theory, label=f'orca{i}', charge=-1, mult=1)
# Optimize:
#dyn = FIRE(neb)
dyn = BFGS(neb)
#dyn.attach(committee_error, interval=1, dyn=dyn, selector="energy")
dyn.attach(write_neb, interval=1, dyn=dyn, out_file=args.output)
dyn.run(fmax=0.5, steps=nsteps)

dyn.atoms.climb=True
dyn.run(fmax=fmax, steps=nsteps)
