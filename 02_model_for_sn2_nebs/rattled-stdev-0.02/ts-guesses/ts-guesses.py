from ase.io import read, write
import numpy as np

for i in range(178):
    images = read(f'nebs/nebs/neb-{i}.xyz', '1:-1')
    energies = [at.get_potential_energy() for at in images]
    ts_guess = images[np.argmax(energies)]
    write('ts-guesses.xyz', ts_guess, append=True)
