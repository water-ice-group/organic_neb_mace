from ase.io import read, write
import numpy as np
from matplotlib import pyplot as plt
import ase
from mace.calculators.mace import MACECalculator
from matplotlib.cm import get_cmap
from matplotlib import rc
import matplotlib as mpl


def d1_d2(config=ase.Atoms, n=int, l=int, c=int):
    pos = config.get_positions()
    d1 = np.linalg.norm(pos[c,:] - pos[n,:])
    d2 = np.linalg.norm(pos[c,:] - pos[l,:])
    return d1, d2


plt.rcParams.update({
    "text.usetex": False,
    "font.family": "serif"
})

s = 12
m = 18
l = 14

plt.rc('font', size=s)          # controls default text sizes
plt.rc('xtick', labelsize=m)    # fontsize of the tick labels
plt.rc('ytick', labelsize=m)    # fontsize of the tick labels
plt.rc('axes', labelsize=m)    # fontsize of the x and y labels

display_mp2 = np.array([0,2,4,5,6,7,8,9,10,11,12,13,14,15,17,19])

# Create a figure and a gridspec with two subplots
fig, axs = plt.subplots(1, 2, figsize=(10,6), gridspec_kw={'width_ratios': [3, 2]})

cmap = get_cmap('viridis')
cm = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.9]
cm.reverse()



lines = []

for i in [0,1,2,3,4,5,6,7,8,9]:
    neb = read(f'../03_ethyl_cl_f_al/iter{i}/neb.xyz', ':')

    d1d2 = np.array([list(d1_d2(atoms,  n=3, l=4, c=0)) for atoms in neb])
    d1 = d1d2[:,0]
    d2 = d1d2[:,1]
    ls = '-'
    alpha=1
    if i < 9:
        ls = '--'
        alpha=0.7
    
    line = axs[0].plot(d2,d1, label=f'AL iter {i}', zorder=0, c=cmap(cm[i]), ls=ls, alpha=alpha, marker='.')
    lines += line



for i in [0,9]:
    neb = read(f'../03_ethyl_cl_f_al/iter{i}/neb.xyz', ':')
    mace_energies = []
    calc = MACECalculator(model_paths=f'../03_ethyl_cl_f_al/iter{i}/MACE_*_swa.model', device='cpu', default_dtype='float64')
    for atoms in neb:
        atoms.calc = calc
        mace_energies.append(atoms.get_potential_energy())
    mace_energies = np.array(mace_energies)
    mace_energies -= np.min(mace_energies)

    d1d2 = np.array([list(d1_d2(atoms,  n=3, l=4, c=0)) for atoms in neb])
    d1 = d1d2[:,0]
    d2 = d1d2[:,1]
    ls = '-'
    alpha=1
    if i < 9:
        ls = '--'
        alpha=0.7
    
    #axs[0].scatter(d2,d1, zorder=0, c=cmap(cm[i]))
    if i == 9 or i == 0:

        axs[1].plot((d2-d1), mace_energies, c=cmap(cm[i]), ls=ls, marker='.', zorder=1)
    #axs[1].scatter((d2-d1), mace_energies, c=cmap(cm[i]), marker='.')

    if i == 0:
        mp2_energies = np.array([atoms.info['mp2_energy'] for atoms in read('../03_ethyl_cl_f_al/iter0/neb-mp2.xyz', ':')])
        mp2_energies -= np.min(mp2_energies)
        axs[1].scatter((d2-d1)[display_mp2], mp2_energies[display_mp2], c=cmap(cm[i]), marker='D',s=16, zorder=2, label=f'MP2 AL iter {i}')

    if i == 9:
        mp2_energies = np.array([atoms.info['mp2_energy'] for atoms in read('../03_ethyl_cl_f_al/iter9/neb-mp2.xyz', ':')])
        mp2_energies -= np.min(mp2_energies)
        axs[1].scatter((d2-d1)[display_mp2], mp2_energies[display_mp2], c=cmap(cm[i]), marker='D',s=16, zorder=2, label=f'MP2 AL iter {i}')



axs[0].set_ylabel('$d_1\ /\ \mathrm{\AA}$')
axs[0].set_xlabel('$d_2\ /\ \mathrm{\AA}$')
axs[0].axis("equal")

axs[1].set_ylabel('Energy / eV')
axs[1].set_xlabel('$d_2 - d_1 /\ \mathrm{\AA}$')

fname = "../03_ethyl_cl_f_al/train.xyz"
n_atom = 5
n_train = 60

configs = read(fname, ":")
initial = configs[n_atom : n_atom  + n_train]
neb_configs = configs[n_atom + n_train :]

d1i = []
d2i = []

for c in initial:
    d1, d2 = d1_d2(c, n=3, l=4, c=0)
    d1i.append(d1)
    d2i.append(d2)

d1n = []
d2n = []

for c in neb_configs:
    d1, d2 = d1_d2(c, n=3, l=4, c=0)
    d1n.append(d1)
    d2n.append(d2)

s1 = axs[0].scatter(d2i, d1i, label="Initial data", marker=".")
s2 = axs[0].scatter(d2n, d1n, label="NEB active learning", marker="D", c='red', s=12)

legend1 = axs[0].legend(handles=lines, loc=1, title='NEB pathways')
axs[0].add_artist(legend1)

legend2 = axs[0].legend(handles=[s1,s2], loc=2, title='Training data')
#legend2 = axs[0].legend(handles=[s1], loc=2, title='Training data')
axs[0].add_artist(legend2)

#axs[0].legend()

axs[0].axis('equal')
axs[0].set_ylabel('$d_1\ /\ \mathrm{\AA}$')
axs[0].set_xlabel('$d_2\ /\ \mathrm{\AA}$')
axs[1].legend()

plt.savefig('ethyl_cl_f_al.pdf')
plt.tight_layout()
plt.show()