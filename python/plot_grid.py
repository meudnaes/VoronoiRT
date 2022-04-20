import h5py

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from plot_line import wavelength, lambda0
from plot_searchlight import font_size

plt.rcParams['text.usetex'] = True

font_size()

center = np.argmin(np.abs(wavelength - lambda0))

f = h5py.File("../data/voronoi_ul7n12_3e6.h5", "r")
p_cont = np.array(f["positions"])
f.close()

f = h5py.File("../data/ionised_hydrogen_3000000.h5", "r")
p_ion = np.array(f["positions"])
f.close()

f = h5py.File("../data/destruction_3e6.h5", "r")
p_des = np.array(f["positions"])
f.close()

tau_unity = np.load("./sourcedata/tau_unity.npz")
tau_unity = tau_unity.reshape(tau_unity.shape[0], -1)


"""
fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)

for i in range(len(wavelength)):
    sns.histplot(data=tau_unity[i,:]/1e6, element="step", fill=False, stat="frequency",
                 ax=ax, lw=0.05, color="k")
    
ax.set_xlabel(r"$\textrm{Height [Mm]}$")
ax.set_ylabel(r"$\textrm{Line Formation Height}$")

plt.savefig("../img/compare_line/tau_unity.pdf")
"""

z_cont = p_cont[:,0]
z_ion = p_ion[:,0]
z_des = p_des[:,0]

fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)

sns.histplot(data=z_cont/1e6, element="step", fill=False, stat="frequency",
             ax=ax, label=r"$\alpha^c$")
sns.histplot(data=z_ion/1e6, element="step", fill=False, stat="frequency",
             ax=ax, label=r"$\textrm{N}_{HII}$")
sns.histplot(data=z_des/1e6, element="step", fill=False, stat="frequency",
             ax=ax, label=r"$\varepsilon$")

y_vals = ax.get_yticks()
ax.set_yticklabels([r"$%.1f$" %(x*1e-6) for x in y_vals])

ax.legend()

#ax.hist(z/1e6, bins=500, histtype="step")
ax.set_xlabel(r"$\textrm{Height [Mm]}$")
ax.set_ylabel(r"$\textrm{Column Density of Sites [1e6/Mm]}$")

ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis

ax2.set_ylabel(r"$\textrm{Height of Line formation}~h(\tau=1)$")  # we already handled the x-label with ax1
sns.histplot(data=tau_unity[center,:]/1e6, element="step", fill=False, stat="frequency",
                 ax=ax2, color="k")

plt.savefig("../img/sites_histogram.pdf")