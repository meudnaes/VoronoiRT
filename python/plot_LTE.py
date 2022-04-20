import numpy as np
import matplotlib as mpl
# import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from plot_searchlight import get_intensity, iunits, font_size

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
intensity_15_000_000 = get_intensity("I_irregular_15000000.npy", PATH)

intensity_quarter = get_intensity("I_regular_quarter.npy", PATH)
intensity_third = get_intensity("I_regular_third.npy", PATH)
intensity_half = get_intensity("I_regular_half.npy", PATH)
intensity_full = get_intensity("I_regular_full.npy", PATH)

intensity_250_000_temp = get_intensity("I_irregular_250000_temp.npy", PATH)
intensity_500_000_temp = get_intensity("I_irregular_500000_temp.npy", PATH)
intensity_1_000_000_temp = get_intensity("I_irregular_1000000_temp.npy", PATH)

intensity_250_000_uniform = get_intensity("I_irregular_250000_uniform.npy", PATH)
intensity_500_000_uniform = get_intensity("I_irregular_500000_uniform.npy", PATH)
intensity_1_000_000_uniform = get_intensity("I_irregular_1000000_uniform.npy", PATH)

intensity_100_000_ionised_hydrogen = get_intensity("I_irregular_100000_ionised_hydrogen.npy", PATH)
intensity_250_000_ionised_hydrogen = get_intensity("I_irregular_250000_ionised_hydrogen.npy", PATH)
intensity_500_000_ionised_hydrogen = get_intensity("I_irregular_500000_ionised_hydrogen.npy", PATH)
intensity_1_000_000_ionised_hydrogen = get_intensity("I_irregular_1000000_ionised_hydrogen.npy", PATH)
intensity_2_500_000_ionised_hydrogen = get_intensity("I_irregular_2500000_ionised_hydrogen.npy", PATH)
intensity_5_000_000_ionised_hydrogen = get_intensity("I_irregular_5000000_ionised_hydrogen.npy", PATH)
intensity_10_000_000_ionised_hydrogen = get_intensity("I_irregular_10000000_ionised_hydrogen.npy", PATH)
intensity_15_000_000_ionised_hydrogen = get_intensity("I_irregular_15000000_ionised_hydrogen.npy", PATH)


font_size()

fig, ax = plt.subplots(1, 3, figsize=(9,3.1), constrained_layout=True)
ax[0].imshow(intensity_250_000_ionised_hydrogen,
               cmap="magma",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[0].axis(False)
ax[0].set_title(r"$2.5\cdot 10^5~\rm{Sites}$")

# Line:
x = np.load(PATH+"x_irregular_100000.npy")
pix2Mm = (x.max() - x.min())*1e-6/len(x)
ax[2].hlines(y=30, xmin=70, xmax=70 + 1/pix2Mm, lw=2, color='w',
             path_effects=[pe.Stroke(linewidth=3, foreground="black"),pe.Normal()])

# Text:
ax[2].text(70, 50, r"\textbf{1 Mm}", color='w', fontsize=12,
           path_effects=[pe.Stroke(linewidth=1, foreground="black"),pe.Normal()])

ax[1].imshow(intensity_2_500_000_ionised_hydrogen,
               cmap="magma",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[1].axis(False)
ax[1].set_title(r"$2.5\cdot 10^6~\rm{Sites}$")

im = ax[2].imshow(intensity_10_000_000_ionised_hydrogen,
               cmap="magma",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[2].axis(False)
ax[2].set_title(r"$10^7~\rm{Sites}$")


plt.colorbar(im, fraction=0.046, pad=0.04, label=iunits)
# plt.savefig("../img/compare_continuum/LTE500_irregular/LTEmaps.pdf")
# plt.show()
plt.close()


fig, ax = plt.subplots(1, 3, figsize=(11,4), constrained_layout=True)

ax[0].imshow(intensity_250_000,
               cmap="magma",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[0].axis(False)
ax[0].set_title(r"$2.5\cdot 10^5~\rm{Sites}$")

ax[1].imshow(intensity_2_500_000,
               cmap="magma",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[1].axis(False)
ax[1].set_title(r"$2.5\cdot 10^6~\rm{Sites}$")

im = ax[2].imshow(intensity_10_000_000,
               cmap="magma",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[2].axis(False)
ax[2].set_title(r"$10^7~\rm{Sites}$")

x = np.load(PATH+"x_irregular_100000.npy")
pix2Mm = (x.max() - x.min())*1e-6/len(x)

# Line:
ax[0].hlines(y=32, xmin=40, xmax=40 + 1/pix2Mm, lw=3, color='w',
             path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])
# Text:
ax[0].text(40, 40, r"\textbf{1 Mm}", color='w', fontsize=14,
           path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])



# Line:
ax[1].hlines(y=32, xmin=40, xmax=40 + 1/pix2Mm, lw=3, color='w',
             path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])

# Text:
ax[1].text(40, 40, r"\textbf{1 Mm}", color='w', fontsize=14,
           path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])

# Line:
ax[2].hlines(y=32, xmin=40, xmax=40 + 1/pix2Mm, lw=3, color='w',
             path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])
# Text:
ax[2].text(40, 40, r"\textbf{1 Mm}", color='w', fontsize=14,
           path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])

plt.colorbar(im, fraction=0.046, pad=0.04, label=iunits)

fig.suptitle(r"$\textbf{Disk-centre Intensity 500\,nm, Irregular Grid}$")
plt.savefig("../img/compare_continuum/LTE500_irregular/LTEmaps.pdf")
plt.close()


fig, ax = plt.subplots(1, 3, figsize=(11,4), constrained_layout=True)

ax[0].imshow(intensity_third,
               cmap="magma",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[0].axis(False)
ax[0].set_title(r"$\rm{One-Third~Resolution}$")

ax[1].imshow(intensity_half,
               cmap="magma",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[1].axis(False)
ax[1].set_title(r"$\rm{Half~Resolution}$")

im = ax[2].imshow(intensity_full,
                  cmap="magma",
                  origin="lower",
                  vmax=CMAX,
                  vmin=CMIN)
ax[2].axis(False)
ax[2].set_title(r"$\rm{Full~Resolution}$")

# Line:
x = np.load(PATH+"x_regular_third.npy")
pix2Mm = (x.max() - x.min())*1e-6/len(x)
ax[0].hlines(y=16/3, xmin=20/3, xmax=20/3 + 1/pix2Mm, lw=3, color='w',
             path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])
# Text:
ax[0].text(20/3, 20/3, r"\textbf{1 Mm}", color='w', fontsize=14,
           path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])



# Line:
x = np.load(PATH+"x_regular_half.npy")
pix2Mm = (x.max() - x.min())*1e-6/len(x)
ax[1].hlines(y=8, xmin=10, xmax=10 + 1/pix2Mm, lw=3, color='w',
             path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])

# Text:
ax[1].text(10, 10, r"\textbf{1 Mm}", color='w', fontsize=14,
           path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])

# Line:
x = np.load(PATH+"x_regular_full.npy")
pix2Mm = (x.max() - x.min())*1e-6/len(x)
ax[2].hlines(y=16, xmin=20, xmax=20 + 1/pix2Mm, lw=3, color='w',
             path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])
# Text:
ax[2].text(20, 20, r"\textbf{1 Mm}", color='w', fontsize=14,
           path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])

fig.suptitle(r"$\textbf{Disk-centre Intensity 500\,nm, Regular Grid}$")
plt.colorbar(im, fraction=0.046, pad=0.04, label=iunits)
plt.savefig("../img/compare_continuum/LTE500_irregular/LTEmaps_regular_grid.pdf")
plt.close()


fig, ax = plt.subplots(1, 3, constrained_layout=True, figsize=(11,4))

im = ax[2].imshow(intensity_1_000_000,
               cmap="magma",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[2].axis(False)
ax[2].set_title(r"$\rm{Sites~from~extinction}~\alpha_{500}$")

ax[1].imshow(intensity_1_000_000_temp,
               cmap="magma",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[1].axis(False)
ax[1].set_title(r"$\rm{Sites~from}~\rm{d} T/\rm{d} z$")

ax[0].imshow(intensity_1_000_000_uniform,
               cmap="magma",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[0].axis(False)
ax[0].set_title(r"$\rm{Sites~from}~U$")

x = np.load(PATH+"x_irregular_100000.npy")
pix2Mm = (x.max() - x.min())*1e-6/len(x)

# Line:
ax[0].hlines(y=32, xmin=40, xmax=40 + 1/pix2Mm, lw=3, color='w',
             path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])
# Text:
ax[0].text(40, 40, r"\textbf{1 Mm}", color='w', fontsize=14,
           path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])



# Line:
ax[1].hlines(y=32, xmin=40, xmax=40 + 1/pix2Mm, lw=3, color='w',
             path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])

# Text:
ax[1].text(40, 40, r"\textbf{1 Mm}", color='w', fontsize=14,
           path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])

# Line:
ax[2].hlines(y=32, xmin=40, xmax=40 + 1/pix2Mm, lw=3, color='w',
             path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])
# Text:
ax[2].text(40, 40, r"\textbf{1 Mm}", color='w', fontsize=14,
           path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])

plt.colorbar(im, fraction=0.046, pad=0.04, label=iunits)

fig.suptitle(r"$\textbf{Disk-centre Intensity 500\,nm, Irregular Grid}$")
plt.savefig("../img/compare_continuum/LTE500_irregular/LTEmaps_compare.pdf")
plt.close()

"""
fig, ax = plt.subplots(1, 2, constrained_layout=True)

ax[0].imshow(intensity_10_000_000,
               cmap="magma",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[0].axis(False)
ax[0].set_title("$10\,000\,000$ Sites")

# Line:
x = np.load(PATH+"x_irregular_100000.npy")
pix2Mm = (x.max() - x.min())*1e-6/len(x)
ax[0].hlines(y=30, xmin=70, xmax=70 + 1/pix2Mm, lw=2, color='w',
             path_effects=[pe.Stroke(linewidth=3, foreground="black"),pe.Normal()])

# Text:
ax[0].text(70, 50, r"\textbf{1 Mm}", color='w', fontsize=12,
           path_effects=[pe.Stroke(linewidth=1, foreground="black"),pe.Normal()])

im = ax[1].imshow(intensity_15_000_000,
               cmap="magma",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[1].axis(False)
ax[1].set_title("$15\,000\,000$ Sites")

plt.colorbar(im, fraction=0.05, pad=0.06)
plt.show()
"""
