import numpy as np
import matplotlib as mpl
# import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from plot_searchlight import get_intensity

#plt.rc('text.latex', preamble=r'\usepackage{cmbright}')
#plt.rc('text', usetex=False)

plt.rcParams['text.usetex'] = True

PATH = "../data/LTE/"

# Empirical values, gotten from simulation
CMAX = 60.52755753644262
CMIN = 20.169713390042126

intensity_100_000 = get_intensity("I_irregular_100000.npy", PATH)
intensity_250_000 = get_intensity("I_irregular_250000.npy", PATH)
intensity_500_000 = get_intensity("I_irregular_500000.npy", PATH)
intensity_1_000_000 = get_intensity("I_irregular_1000000.npy", PATH)
intensity_2_500_000 = get_intensity("I_irregular_2500000.npy", PATH)
intensity_5_000_000 = get_intensity("I_irregular_5000000.npy", PATH)
intensity_10_000_000 = get_intensity("I_irregular_10000000.npy", PATH)

intensity_quarter = get_intensity("I_regular_quarter.npy", PATH)
intensity_half = get_intensity("I_regular_half.npy", PATH)
intensity_full = get_intensity("I_regular_full.npy", PATH)

intensity_250_000_temp = get_intensity("I_irregular_250000_temp.npy", PATH)
intensity_500_000_temp = get_intensity("I_irregular_500000_temp.npy", PATH)
intensity_1_000_000_temp = get_intensity("I_irregular_1000000_temp.npy", PATH)

intensity_250_000_uniform = get_intensity("I_irregular_250000_uniform.npy", PATH)
intensity_500_000_uniform = get_intensity("I_irregular_500000_uniform.npy", PATH)
intensity_1_000_000_uniform = get_intensity("I_irregular_1000000_uniform.npy", PATH)

"""
fig, ax = plt.subplots(1, 3, figsize=(9,3), constrained_layout=True)

ax[0].imshow(intensity_100_000,
               cmap="magma",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[0].axis(False)
ax[0].set_title("$100\,000$ Sites")

# Line:
x = np.load(PATH+"x_irregular_100000.npy")
pix2Mm = (x.max() - x.min())*1e-6/len(x)
ax[0].hlines(y=30, xmin=70, xmax=70 + 2/pix2Mm, lw=2, color='w',
             path_effects=[pe.Stroke(linewidth=3, foreground="black"),pe.Normal()])

# Text:
ax[0].text(70, 50, r"\textbf{2 Mm}", color='w', fontsize=12,
           path_effects=[pe.Stroke(linewidth=1, foreground="black"),pe.Normal()])

ax[1].imshow(intensity_1_000_000,
               cmap="magma",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[1].axis(False)
ax[1].set_title("$1\,000\,000$ Sites")

im = ax[2].imshow(intensity_10_000_000,
               cmap="magma",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[2].axis(False)
ax[2].set_title("$10\,000\,000$ Sites")


plt.colorbar(im, fraction=0.05, pad=0.06)
plt.savefig("../img/compare_continuum/LTE500_irregular/LTEmaps.png", dpi=300)
plt.close()
"""

"""
fig, ax = plt.subplots(1, 3, figsize=(9,3), constrained_layout=True)

ax[0].imshow(intensity_quarter,
               cmap="magma",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[0].axis(False)
ax[0].set_title("Quarter Resolution")

# Line:
x = np.load(PATH+"x_regular_quarter.npy")
pix2Mm = (x.max() - x.min())*1e-6/len(x)
ax[0].hlines(y=8, xmin=10, xmax=10 + 2/pix2Mm, lw=2, color='w',
             path_effects=[pe.Stroke(linewidth=3, foreground="black"),pe.Normal()])

# Text:
ax[0].text(10, 10, r"\textbf{2 Mm}", color='w', fontsize=12,
           path_effects=[pe.Stroke(linewidth=1, foreground="black"),pe.Normal()])


ax[1].imshow(intensity_half,
               cmap="magma",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[1].axis(False)
ax[1].set_title("Half Resolution")

im = ax[2].imshow(intensity_full,
                  cmap="magma",
                  origin="lower",
                  vmax=CMAX,
                  vmin=CMIN)
ax[2].axis(False)
ax[2].set_title("Full Resolution")

plt.colorbar(im, fraction=0.05, pad=0.06)
plt.savefig("../img/compare_continuum/LTE500_irregular/LTEmaps_regular_grid.png", dpi=300)
plt.close()
"""


fig, ax = plt.subplots(1, 3, constrained_layout=True, figsize=(9,3))

ax[0].imshow(intensity_1_000_000,
               cmap="magma",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[0].axis(False)
ax[0].set_title("Sites from extinction")

# Line:
x = np.load(PATH+"x_irregular_100000.npy")
pix2Mm = (x.max() - x.min())*1e-6/len(x)
ax[0].hlines(y=30, xmin=70, xmax=70 + 2/pix2Mm, lw=2, color='w',
             path_effects=[pe.Stroke(linewidth=3, foreground="black"),pe.Normal()])

# Text:
ax[0].text(70, 50, r"\textbf{2 Mm}", color='w', fontsize=12,
           path_effects=[pe.Stroke(linewidth=1, foreground="black"),pe.Normal()])

ax[1].imshow(intensity_1_000_000_temp,
               cmap="magma",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[1].axis(False)
ax[1].set_title(r"Sites from $dT/dz$")

im = ax[2].imshow(intensity_1_000_000_uniform,
               cmap="magma",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[2].axis(False)
ax[2].set_title(r"Sites from $U$")

plt.colorbar(im, fraction=0.05, pad=0.06)
plt.savefig("../img/compare_continuum/LTE500_irregular/LTEmaps_compare.png", dpi=300)
plt.close()
