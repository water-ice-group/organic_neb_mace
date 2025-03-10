import numpy as np
from matplotlib import pyplot as plt
from mace.calculators.mace import MACECalculator
from ase.io import read
from matplotlib import colors
import matplotlib.ticker as ticker
import matplotlib.pyplot as mpl

number_iters = 11
name = "../01_iterative_training/iterative_training_models"
data = []

mpl.rcParams['font.size'] = 14

fig = plt.figure(constrained_layout=False, figsize=(13,5))
gs = fig.add_gridspec(2, 3)
linestyles = ["-", "--", ":", "-."]

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0])
ax3 = fig.add_subplot(gs[0:, 1:])

rmse_barrier = []

for i in range(number_iters):
    ts_dict = {}
    r_list = []
    ts_list = []
    test_configs = read("{}/iter{}/test_eval.xyz".format(name, i), ":")
    train_configs = read("{}/iter{}/train_eval.xyz".format(name, i), "7:")

    for config in test_configs:
        ts_dict.update({config.info["hash"]:config})
    for config in train_configs:
        if config.info['config'] != 'ts':
            if config.info['transition_state_hash'] in ts_dict:
                r_list.append(config)
                ts_list.append(ts_dict[config.info['transition_state_hash']])
    #ts_list = [ts_dict[config.info['transition_state_hash']] for config in r_list]
    mace_barriers = [ts.info['MACE_energy'] - r.info['MACE_energy'] for r, ts in zip(r_list, ts_list)]
    mp2_barriers = [ts.info['mp2_energy'] - r.info['mp2_energy'] for r, ts in zip(r_list, ts_list)]
    rmse_barrier.append(np.sqrt(np.mean((np.array(mace_barriers) - np.array(mp2_barriers))**2)))
    if i == 0:
        hist1 = ax1.hist2d(mace_barriers, mp2_barriers, bins=30, cmap='viridis', norm=colors.LogNorm())
    if i == 10:
        hist2 = ax2.hist2d(mace_barriers, mp2_barriers, bins=30, cmap='viridis', norm=colors.LogNorm())



# Set the same axis limits for all subplots
xlim = (-2, 4)
ylim = (-2, 4)

formatter = ticker.LogFormatter()

ax1.set_xlim(xlim)
ax1.set_ylim(ylim)
ax1.set_aspect('equal')
# Add a red dashed line for y=x equation
ax1.plot(ylim, ylim, 'r--')

# Add axis ticks only at integer values
ax1.set_xticks(np.arange(np.ceil(xlim[0]), np.floor(xlim[1])+1, 1))
ax1.set_yticks(np.arange(np.ceil(ylim[0]), np.floor(ylim[1])+1, 1))
ax2.set_xlim(xlim)
ax2.set_ylim(ylim)
ax2.set_aspect('equal')
# Add a red dashed line for y=x equation
ax2.plot(ylim, ylim, 'r--')

# Add axis ticks only at integer values
ax2.set_xticks(np.arange(np.ceil(xlim[0]), np.floor(xlim[1])+1, 1))
ax2.set_yticks(np.arange(np.ceil(ylim[0]), np.floor(ylim[1])+1, 1))

ax2.set_xlabel('Reference barrier / eV')
ax2.set_ylabel('MACE barrier / eV')

cbar = fig.colorbar(hist2[3], ax=[ax1, ax2], format=formatter)
cbar.set_label('Counts')

ax1.text(-6.8, 0.8, "0% TS", fontsize=14)
ax2.text(-6.8, 0.8, "10% TS", fontsize=14)

ax3.plot(np.arange(0, 11), np.array(rmse_barrier)*1000, linestyle='--', zorder=1, linewidth=2)
print(rmse_barrier)
ax3.set_xlim(0, 10)
ax3.set_ylim(0,np.max(rmse_barrier)*1000+10)
ax3.set_aspect(0.03)
#ax.set_xticks(xval)
ax3.set_xlabel("% of TS in training set")
ax3.set_ylabel("RMSE in barrier height / meV")

ax3.plot([0, 10], [43, 43], linestyle="--", c="black", zorder=3)
text = ax3.text(x=1.6, y=40, s="43 meV", bbox=dict(facecolor='none', edgecolor='black', boxstyle='round'), zorder=4)
text.set_backgroundcolor("white")
#ax3.set_xticks(np.arange(6)*2)
#plt.tight_layout()
plt.savefig("iters.pdf", bbox_inches='tight')