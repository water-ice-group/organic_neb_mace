from ase.io import read, write
from ase import Atoms, Atom
from ase import units

images = []

z_dict = {1:"C", 2:"H", 3:"N", 4:"O", 5:"S"}

for i in range(32):
    fname = f'final-{i+1}'
    image = Atoms(None)
    with open(fname) as lines:
        for line in lines:
            items = line.strip().split(" ")
            at = Atom(
                symbol=z_dict[int(items[1])], 
                position=(float(items[2]), float(items[3]), float(items[4]))
                )
            image += at
    images.append(image)

with open('log.lammps') as f:
    for line in f:
        pass
    last_line = line
    items = line.strip().split()
    for i in range(32):
        images[i].info["energy"] = float(items[10+2*i]) * units.kcal / units.mol

write('neb.xyz', images)

