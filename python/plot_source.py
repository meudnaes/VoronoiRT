import numpy as np
import matplotlib as mpl
# import matplotlib.cm as cm
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm
from plot_searchlight import get_intensity, font_size, iunits
from plot_line import wavelength, lambda0

plt.rcParams['text.usetex'] = True

centre = np.argmin(np.abs(wavelength - lambda0))

PATH = "./sourcedata/"

font_size()

source_half_res = np.load("./sourcedata/half_res_ul7n12_source.npy")[:,:,:,:]
source_irregular = np.load("./sourcedata/voronoi_ul7n12_3e6_2_source.npy")[:,:,:,:]

print(source_half_res.shape)
print(source_irregular.shape)

S_diff = np.abs(1 - source_irregular/source_half_res)
S_diff = np.max(S_diff, axis=0)

regular_z = np.load("./sourcedata/half_res_ul7n12_z.npy")
irregular_z = np.load("./sourcedata/voronoi_ul7n12_3e6_2_z.npy")

print(np.sum(regular_z - irregular_z))

# source_half_res = source_half_res.reshape(source_half_res.shape[0], -1)
# source_irregular = source_irregular.reshape(source_irregular.shape[0], -1)

"""
fig, ax = plt.subplots(1, 2, figsize=(10, 6), constrained_layout=True, sharey=True)

ax[0].plot(regular_z/1e6,
           source_half_res[:, :],
           color='k',
           lw=0.005)
ax[0].set_xlabel(r"$\textrm{Height [Mm]}$")
ax[0].set_ylabel(r"$S_\lambda~\left[\textrm{kW}\,\textrm{m}^{-2}\,\textrm{nm}^{-1}\right]$")
ax[0].set_title(r"$\textrm{Regular Grid}$")

ax[1].plot(irregular_z/1e6,
           source_irregular[:, ::4],
           color='k',
           lw=0.005)
ax[1].set_xlabel(r"$\textrm{Height [Mm]}$")
ax[1].set_title(r"$\textrm{Irregular Grid}$")

fig.suptitle(r"$\textrm{Source Function}$")
plt.savefig("../img/compare_line/source_function.png", dpi=300)
"""

"""
S_diff = S_diff.reshape(S_diff.shape[0], S_diff.shape[1], -1)
S_diff = np.mean(S_diff, axis=2)

print(S_diff.shape)

fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
im=ax.imshow(S_diff.T,
             rasterized=True,
             origin="lower",
             norm=LogNorm(vmin=1e-3, vmax=1e5),
             aspect="auto",
             cmap="cividis",
             extent=[regular_z[0]/1e6, regular_z[-1]/1e6, wavelength[0], wavelength[-1]])
fig.colorbar(im, label=r"$\textrm{Average relative error}$")

ax.set_xlabel(r"$\textrm{height [Mm]}$")
ax.set_ylabel(r"$\textrm{wavelength [nm]}$")

plt.savefig("../img/compare_line/source_diff.pdf")
"""

S_diff = S_diff.reshape(S_diff.shape[0], -1)
S_median = np.median(S_diff, axis=1)
print(S_median.shape)

fig, ax = plt.subplots(figsize=(6, 5.5), constrained_layout=True, sharey=True)

ax.plot(regular_z/1e6,
        S_diff,
        color='k',
        lw=0.0035,
        alpha=1.0,
        rasterized=True)

ax.plot(regular_z/1e6,
        S_median,
        color="cyan",
        label="median")

ax.legend(loc="upper right")

ax.set_xlabel(r"$\textrm{Height [Mm]}$")
ax.set_ylabel(r"$\textrm{Max rel. diff.,}~\max_\lambda\left(1 - S_\textrm{irregular}/S_\textrm{regular}\right)$")
ax.set_yscale("log")
ax.set_ylim(1e-4, 1e4)

# fig.suptitle(r"$\textrm{Difference in Source Function}$")
plt.savefig("../img/compare_line/source_diff.pdf", dpi=400)
