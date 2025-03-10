from ase.io import read, write
import numpy as np
from ase.mep import NEB
import glob


def get_fmax(images, k=0.02):
    neb = NEB(images, k=k)
    neb.get_forces()
    return neb.get_residual()

images = read('./no-rattle/nebs/nebs/neb-0.xyz', ':')
print(get_fmax(images))

fmax_0 = np.array([get_fmax(read(fname, ':')) for fname in glob.glob('./no-rattle/nebs/nebs/neb-*.xyz')])
fmax_1 = np.array([get_fmax(read(fname, ':')) for fname in glob.glob('./rattled-stdev-0.02/nebs/nebs/neb-*.xyz')])
fmax_2 = np.array([get_fmax(read(fname, ':')) for fname in glob.glob('./rattled-stdev-0.05/nebs/nebs/neb-*.xyz')])

print(np.sum(fmax_0>0.05))
print(np.sum(fmax_1>0.05))
print(np.sum(fmax_1>0.05))