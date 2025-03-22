from ase.io import read, write
from ase.io.lammpsdata import write_lammps_data
from ase.mep import NEB
import argparse


def make_final(config):

    lines = [str(len(config)) + '\n']

    for i, pos in enumerate(config.get_positions()):
        line = [str(i+1)] + pos.astype(str).tolist()
        lines.append(" ".join(line) + '\n')

    return lines


parser = argparse.ArgumentParser()
parser.add_argument('--nimages', type=int, default=32)
parser.add_argument('--r', type=str)
parser.add_argument('--p', type=str)
parser.add_argument('--idpp', action='store_true')


args = parser.parse_args()

if args.p is None:
    print('p not supplied, building configs using r.')
    images = read(args.r, ':')
    if len(images) < 2:
        print('p was not supplied and r does not contain enough images to build bot end points')
    
    for atoms in images:
        atoms.center()
        atoms.translate([5, 5, 5])
        atoms.set_cell([10, 10, 10])
    
    write_lammps_data('r.txt', images[0], units='real', atom_style='full', masses=True)

    for i, atoms in enumerate(images[1:]):
        lines = make_final(atoms)
        
        with open(f"coords.initial.{i+1}", 'w') as f:
            f.writelines(lines)

else:
    r = read(args.r)
    r.center()
    r.translate([5, 5, 5])
    r.set_cell([10, 10, 10])
    write_lammps_data(args.r.replace('xyz', '') + 'txt', r, units='real', atom_style='full', masses=True)

    p = read(args.p)
    p.center()
    p.translate([5, 5, 5])
    p.set_cell([10, 10, 10])

    lines = make_final(p)

    with open('p.txt', 'w') as f:
        f.writelines(lines)

    if args.idpp:
        images = [r] + [r.copy() for i in range(args.nimages-2)] + [p]
        neb = NEB(images, k=args.nimages*0.01, remove_rotation_and_translation=True)
        neb.interpolate('idpp')

        for i, image in enumerate(neb.images[1:]):
            lines = make_final(image)
            
            with open(f"coords.initial.{i+1}", 'w') as f:
                f.writelines(lines)




