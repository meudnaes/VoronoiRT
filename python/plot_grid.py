import h5py

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from plot_line import wavelength, lambda0
from plot_searchlight import font_size

plt.rcParams['text.usetex'] = True

font_size()

center = np.argmin(np.abs(wavelength - lambda0))
blue_wing = center - 10

f = h5py.File("../data/voronoi_ul7n12_3e6.h5", "r")
p_cont = np.array(f["positions"])
f.close()

f = h5py.File("../data/ionised_hydrogen_3000000.h5", "r")
p_ion = np.array(f["positions"])
f.close()

f = h5py.File("../data/destruction_3e6.h5", "r")
p_des = np.array(f["positions"])
f.close()

z_regular = np.load("./sourcedata/half_res_ul7n12_z.npy")

tau_unity = np.load("./sourcedata/tau_unity.npz")
tau_unity = tau_unity.reshape(tau_unity.shape[0], -1)

tau_mean = np.mean(tau_unity[0:len(wavelength),:], axis=0)

"""
fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)

sns.displot(data=tau_unity[0:len(wavelength),:]/1e6, y=z_regular, x=wavelength, element="step", fill=False, stat="frequency",
                 ax=ax, lw=0.05, color="k")

ax.set_xlabel(r"$\textrm{Height [Mm]}$")
ax.set_ylabel(r"$\textrm{Line Formation Height}$")

plt.savefig("../img/compare_line/tau_unity.pdf")
"""

z_cont = p_cont[:,0]
z_ion = p_ion[:,0]
z_des = p_des[:,0]

fig, ax = plt.subplots(figsize=(6,5), constrained_layout=True)

sns.histplot(data=z_cont/1e6, element="step", fill=False, stat="frequency",
             ax=ax, label=r"$\alpha^c$", color="blue", ls="solid")
sns.histplot(data=z_ion/1e6, element="step", fill=False, stat="frequency",
             ax=ax, label=r"$\textrm{N}_\textrm{HII}$", color="dodgerblue", ls="solid")
sns.histplot(data=z_des/1e6, element="step", fill=False, stat="frequency",
             ax=ax, label=r"$\varepsilon$", color="cyan", ls="solid")

y_vals = ax.get_yticks()
ax.set_yticklabels([r"$%.1f$" %(x*1e-6) for x in y_vals])

ax.legend(loc="upper left")
ax.yaxis.label.set_color("blue")
ax.tick_params(axis="y", colors="blue")

#ax.hist(z/1e6, bins=500, histtype="step")
ax.set_xlabel(r"$\textrm{Height [Mm]}$")
ax.set_ylabel(r"$\textrm{Density of Sites per height [1e6/Mm]}$")

ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
# ax2.set_zorder(-1)

sns.histplot(data=tau_unity[center, :]/1e6, element="step", fill=True, stat="density",
                 ax=ax2, color="maroon", label=r"$\textrm{line core}$", alpha=0.2)
sns.histplot(data=tau_unity[blue_wing, :]/1e6, element="step", fill=True, stat="density",
                 ax=ax2, color="orangered", label=r"$\textrm{line wing}$", alpha=0.2)

ax2.set_ylabel(r"$\textrm{Density of Line formation per height}$")
ax2.yaxis.label.set_color("#980002")
ax2.tick_params(axis="y", colors="#980002")
ax2.legend(loc="upper right")

plt.savefig("../img/sites_histogram.pdf")
