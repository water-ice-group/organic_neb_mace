from ase.io import read, write
import numpy as np
from matplotlib import pyplot as plt
from mace.calculators.mace import MACECalculator
from ase.geometry import distance
from ase.build import minimize_rotation_and_translation
from ase.visualize import view
from ase.mep import NEB
from matplotlib.ticker import ScalarFormatter
from matplotlib.spines import Spine
from matplotlib.gridspec import GridSpec


true_ts_0 = read('../02_model_for_sn2_nebs/ts-test-mp2.xyz', ':')
guess_ts_0 = read('../02_model_for_sn2_nebs/no-rattle/ts-guesses/ts-guesses-mp2.xyz', ':')
true_ts_mp2_0 = np.array([at.info['mp2_energy'] for at in true_ts_0])
guess_ts_mp2_0 = np.array([at.info['mp2_energy'] for at in guess_ts_0])
true_ts_mp2_guess_ts_mp2_err_0 = guess_ts_mp2_0 - true_ts_mp2_0

true_ts_1 = read('../02_model_for_sn2_nebs/ts-test-mp2.xyz', ':')
guess_ts_1 = read('../02_model_for_sn2_nebs/rattled-stdev-0.02/ts-guesses/ts-guesses-mp2.xyz', ':')
true_ts_mp2_1 = np.array([at.info['mp2_energy'] for at in true_ts_1])
guess_ts_mp2_1 = np.array([at.info['mp2_energy'] for at in guess_ts_1])
true_ts_mp2_guess_ts_mp2_err_1 = guess_ts_mp2_1 - true_ts_mp2_1

true_ts_2 = read('../02_model_for_sn2_nebs/ts-test-mp2.xyz', ':')
guess_ts_2 = read('../02_model_for_sn2_nebs/rattled-stdev-0.05/ts-guesses/ts-guesses-mp2.xyz', ':')
true_ts_mp2_2 = np.array([at.info['mp2_energy'] for at in true_ts_2])
guess_ts_mp2_2 = np.array([at.info['mp2_energy'] for at in guess_ts_2])
true_ts_mp2_guess_ts_mp2_err_2 = guess_ts_mp2_2 - true_ts_mp2_2

plt.rcParams['axes.titlesize'] = 14
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['legend.title_fontsize'] = 16


fig,(ax2,ax1)=plt.subplots(2, 1, figsize=(10,6))

#ax.set_facecolor('#eeeeee')
bins = 40
bins = np.histogram(np.hstack((true_ts_mp2_guess_ts_mp2_err_0, true_ts_mp2_guess_ts_mp2_err_1)), bins=bins)[1]

bin_width = bins[1] - bins[0]

ax1.hist(true_ts_mp2_guess_ts_mp2_err_0, bins, width=bin_width/3, color='blue', zorder=2, edgecolor='#666666', label='MACE-sn2')
ax1.hist(true_ts_mp2_guess_ts_mp2_err_1, bins+bin_width/3, width=bin_width/3, color='red', zorder=2, edgecolor='#666666', label='MACE-sn2-rattled-A')

ax2.hist(true_ts_mp2_guess_ts_mp2_err_0, bins, width=bin_width/3, color='blue', zorder=2, edgecolor='#666666', label='MACE-sn2')
ax2.hist(true_ts_mp2_guess_ts_mp2_err_1, bins+bin_width/3, width=bin_width/3, color='red', zorder=2, edgecolor='#666666', label='MACE-sn2-rattled-A')

ax1.set_ylim(0, 11)
ax2.set_ylim(105, 198)

#ax1.set_xlim(-1.5, 5)
#ax2.set_xlim(-1.5, 5)


ax1.spines['top'].set(alpha=0.3, ls='--')
#ax2.xaxis.set_ticks([])
ax2.tick_params(labeltop=False)
ax2.spines['bottom'].set(alpha=0.3, ls='--')
ax2.xaxis.tick_top()

#ax.hist(true_ts_mp2_guess_ts_mp2_err_2, bins+bin_width*2/3, alpha=0.5, width=bin_width/3, color='green', zorder=2, edgecolor='green', label='MACE-sn2-rattled-stdev-0.02')
ax2.legend(loc=2)
err_below_x_0 = [np.sum((np.absolute(true_ts_mp2_guess_ts_mp2_err_0)<i).astype(int))/178 for i in np.linspace(0,4.47,100000)]
err_below_x_1 = [np.sum((np.absolute(true_ts_mp2_guess_ts_mp2_err_1) <i).astype(int))/178 for i in np.linspace(0,4.47,100000)]

#err_below_x_2 = [np.sum((np.absolute(true_ts_mp2_guess_ts_mp2_err_2) <i).astype(int)) for i in np.linspace(0,3.5,50)]

#ax.grid(True, which='both', color='white', zorder=0)
ax1.set_xlabel('Energy error on TS / eV')
ax2.set_ylabel('Frequency')
ax2.yaxis.set_label_coords(-0.08, 0)
ax2.set_xticks([])

left, bottom, width, height = [0.55, 0.65, 0.4, 0.25]
ax3 = fig.add_axes([left, bottom, width, height])
ax3.plot(np.linspace(0,4.47,100000), err_below_x_0, c='blue', zorder=2)
ax3.plot(np.linspace(0,4.47,100000), err_below_x_1, c='red', zorder=2)
ax3.set_xscale('log')
ax3.set_xticks([0.0001, 0.001, 0.01, 0.1, 1])
ax3.set_yticks([0, 1])
ax3.set_ylim(0,1)

ax3.set_xlim(0.0001, 4.470)
#formatter = ScalarFormatter()
#formatter.set_powerlimits((-1, 1))
#formatter.set_scientific(True)
#ax3.xaxis.set_major_formatter(formatter)
ax3.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
#ax1.plot(np.linspace(0,3.5,50), err_belox_x_2, c='green', zorder=2)
ax3.set_xlabel('Absolute error / eV')
ax3.set_ylabel('Fraction')
ax3.set_title('Fraction of TS with error below given value')
ax3.grid(zorder=0, which='both', axis='y')
ax3.grid(zorder=0, axis='x')

d=0.015

kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax1.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax2.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

left, bottom, width, height = [0.55, 0.25, 0.4, 0.25]
ax4 = fig.add_axes([left, bottom, width, height])
bins2 = np.linspace(-0.05, 0.05, 22)

bin_width = 0.1 / 22

ax4.hist(true_ts_mp2_guess_ts_mp2_err_0, bins2, rwidth=1/3, color='blue', zorder=2, edgecolor='#666666', label='MACE-sn2', align='mid')
ax4.hist(true_ts_mp2_guess_ts_mp2_err_1, bins2 + bin_width/3, rwidth=1/3, color='red', zorder=2, edgecolor='#666666', label='MACE-sn2-r', align='mid')
ax4.set_ylabel('Frequency')

ax3.text(0.05, 0.95, 'a', transform=ax3.transAxes, fontsize=24, verticalalignment='top', horizontalalignment='left')
ax4.text(0.05, 0.95, 'b', transform=ax4.transAxes, fontsize=24, verticalalignment='top', horizontalalignment='left')

#ax[0].set_title('MACE trained with SN2 R+TS+P')
#ax[1].set_title('MACE trained with rattled SN2 R+TS+P')
plt.tight_layout()

plt.savefig('ts-guesses.pdf')


# Print out some statistics

from ase.mep import NEB
import glob


def get_fmax(images, k=0.02):
    neb = NEB(images, k=k)
    neb.get_forces()
    return neb.get_residual()


def rmse(arr):
    return np.sqrt(np.mean(arr**2))


datapath = '../02_model_for_sn2_nebs/'

fmax_0 = np.array([get_fmax(read(f'{datapath}no-rattle/nebs/nebs/neb-{i}.xyz', ':')) for i in range(178)])
fmax_1 = np.array([get_fmax(read(f'{datapath}rattled-stdev-0.02/nebs/nebs/neb-{i}.xyz', ':')) for i in range(178)])
fmax_2 = np.array([get_fmax(read(f'{datapath}rattled-stdev-0.05/nebs/nebs/neb-{i}.xyz', ':')) for i in range(178)])
print(fmax_0)

print(np.sum(fmax_0>0.05))
print(np.sum(fmax_1>0.05))
print(np.sum(fmax_2>0.05))

true_ts_0 = read(datapath+'ts-test-mp2.xyz', ':')
guess_ts_0 = read(datapath+'no-rattle/ts-guesses/ts-guesses-mp2.xyz', ':')

true_ts_mp2_0 = np.array([at.info['mp2_energy'] for at in true_ts_0])
guess_ts_mp2_0 = np.array([at.info['mp2_energy'] for at in guess_ts_0])
true_ts_mp2_guess_ts_mp2_err_0 = guess_ts_mp2_0 - true_ts_mp2_0

true_ts_1 = read(datapath+'ts-test-mp2.xyz', ':')
guess_ts_1 = read(datapath+'rattled-stdev-0.02/ts-guesses/ts-guesses-mp2.xyz', ':')

true_ts_mp2_1 = np.array([at.info['mp2_energy'] for at in true_ts_1])
guess_ts_mp2_1 = np.array([at.info['mp2_energy'] for at in guess_ts_1])
true_ts_mp2_guess_ts_mp2_err_1 = guess_ts_mp2_1 - true_ts_mp2_1

true_ts_2 = read(datapath+'ts-test-mp2.xyz', ':')
guess_ts_2 = read(datapath+'rattled-stdev-0.05/ts-guesses/ts-guesses-mp2.xyz', ':')

true_ts_mp2_2 = np.array([at.info['mp2_energy'] for at in true_ts_2])
guess_ts_mp2_2 = np.array([at.info['mp2_energy'] for at in guess_ts_2])
true_ts_mp2_guess_ts_mp2_err_2 = guess_ts_mp2_2 - true_ts_mp2_2

print(rmse(true_ts_mp2_guess_ts_mp2_err_0[fmax_0<0.05]))
print(rmse(true_ts_mp2_guess_ts_mp2_err_1[fmax_1<0.05]))
print(rmse(true_ts_mp2_guess_ts_mp2_err_2[fmax_2<0.05]))