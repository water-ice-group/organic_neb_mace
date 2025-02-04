import numpy as np
from neb import NEB, NEBOptimizer
from ase.optimize import FIRE
from mace.calculators.mace import MACECalculator
from ase.io import read, write
import argparse
import ase
from ase.parallel import world

def write_neb(dyn, fname):
    write(fname, dyn.atoms.images)


rxn_ids = np.arange(178)
nimages = 21
nsteps = 1000
fmax = 0.05
k = 0.02
model_path = '../../model/MACE_swa.model'
reactants_fname = 'r-test.xyz'
products_fname = 'p-test.xyz'

for rxn in rxn_ids:
    r = read(reactants_fname, rxn)
    p = read(products_fname, rxn)
    images = [r]
    images += [r.copy() for i in range(nimages-2)]
    images += [p]
    for image in images[1:-1]:
        image.calc = MACECalculator(model_paths=model_path, default_dtype='float64', device='cuda')
    neb = NEB(images, k=k, remove_rotation_and_translation=True, use_gpu=True, calc=images[1].calc)
    neb.interpolate()
    dyn = FIRE(neb)
    dyn.attach(write_neb, interval=1, dyn=dyn, fname=f'nebs/neb-{rxn}.xyz')
    dyn.run(steps=nsteps, fmax=fmax)
