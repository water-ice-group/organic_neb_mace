
from ase.io import read, write
import numpy as np


# Load Final AL path

images = read('../../mace-neb-al-250steps/iter18/neb.xyz', ':')
energies = [atoms.get_potential_energy() for atoms in images]

ts_idx = np.argmax(energies)

r = images[ts_idx-2]
p = images[ts_idx+2]

write('r.xyz', r)
write('p.xyz', p)