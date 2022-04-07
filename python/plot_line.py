import numpy as np
import matplotlib as mpl
# import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from plot_searchlight import get_intensity

#plt.rc('text.latex', preamble=r'\usepackage{cmbright}')
#plt.rc('text', usetex=False)


lambda0 = 121.56841096386111 # nm
wavelength = np.array([120.85647513019845, 121.04863120292787, 121.18861407155109,
                       121.29060823835265, 121.36494181786958, 121.41913498583673,
                       121.45866338283498, 121.48751396885412, 121.50858975582949,
                       121.52400450419977, 121.53529729933265, 121.54358879034517,
                       121.54969495175679, 121.55420991638432, 121.55756628818231,
                       121.56007905763141, 121.56197757770408, 121.5634288464184,
                       121.56455445948667, 121.56544295399095, 121.56615879614279,
                       121.56674892551068, 121.56724752004678, 121.56767946562972,
                       121.56806288233183, 121.56841096386111, 121.56875904539042,
                       121.56914246209253, 121.56957440767549, 121.57007300221157,
                       121.57066313157947, 121.57137897373131, 121.5722674682356,
                       121.57339308130386, 121.57484435001817, 121.57674287009084,
                       121.57925563953995, 121.58261201133794, 121.58712697596546,
                       121.5932331373771,  121.6015246283896,  121.6128174235225,
                       121.62823217189276, 121.64930795886815, 121.67815854488727,
                       121.71768694188553, 121.77188010985269, 121.84621368936962,
                       121.94820785617118, 122.0881907247944,  122.2803467975238]) # nm


plt.rcParams['text.usetex'] = True

PATH = "./linedata/"

# Empirical values, gotten from simulation
CMAX = 150
CMIN = 1

intensity_half = get_intensity("regular_half_disk_centre.npy", PATH)
intensity_third = get_intensity("regular_third_disk_centre.npy", PATH)
intensity_quarter = get_intensity("regular_quarter_disk_centre.npy", PATH)

intensity_5e5 = get_intensity("voronoi_5e5_disk_centre.npy", PATH)
intensity_5e5_1dot5 = get_intensity("voronoi_ul7n12_5e5_disk_centre_1dot5.npy", PATH)
intensity_1e6 = get_intensity("voronoi_ul7n12_1e6_disk_centre_1.npy", PATH)
intensity_1e6_1dot5 = get_intensity("voronoi_ul7n12_1e6_disk_centre_1dot5.npy", PATH)
intensity_ext_5e5 = get_intensity("total_ext_5e5_disk_centre_1.npy", PATH)
intensity_des_5e5 = get_intensity("destruction_5e5_disk_centre_1.npy", PATH)
intensity_density_5e5 = get_intensity("density_5e5_disk_centre_1.npy", PATH)

convergence_quarter = np.load(PATH+"regular_ul7n12_quarter.npy")
convergence_half = np.load(PATH+"regular_ul7n12_half.npy")
convergence_third = np.load(PATH+"regular_ul7n12_third.npy")
convergence_5e5 = np.load(PATH+"voronoi_ul7n12_5e5.npy")


center = np.argmin(np.abs(wavelength - lambda0))
left_wing = center - 10
right_wing = center + 10

"""
fig, ax = plt.subplots(1, 3, figsize=(9,3), constrained_layout=True)

ax[0].imshow(intensity_half[left_wing, :, :],
               cmap="gist_gray_r",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[0].axis(False)
wl = "{idx}".format(idx=wavelength[left_wing])
ax[0].set_title(r"$\textrm{Blue Wing}$"+wl)

# Line:
x = np.load("../data/LTE/x_regular_half.npy")
pix2Mm = (x.max() - x.min())*1e-6/len(x)
ax[0].hlines(y=8, xmin=10, xmax=10 + 1/pix2Mm, lw=2, color='w',
             path_effects=[pe.Stroke(linewidth=3, foreground="black"),pe.Normal()])

# Text:
ax[0].text(10, 10, r"\textbf{1 Mm}", color='w', fontsize=12,
           path_effects=[pe.Stroke(linewidth=1, foreground="black"),pe.Normal()])

ax[1].imshow(intensity_half[center, :, :],
               cmap="gist_gray_r",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[1].axis(False)
wl = "{idx}".format(idx=wavelength[center])
ax[1].set_title(r"$\textrm{Line Centre}$"+wl)

im = ax[2].imshow(intensity_half[right_wing, :, :],
               cmap="gist_gray_r",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[2].axis(False)
wl = "{idx}".format(idx=wavelength[right_wing])
ax[2].set_title(r"$\textrm{Red Wing}$"+wl)

fig.colorbar(im, fraction=0.05, pad=0.06)

plt.savefig("../img/compare_line/disk_centre_half.pdf")
"""

"""
fig, ax = plt.subplots(1, 3, figsize=(9,3), constrained_layout=True)

ax[0].imshow(intensity_quarter[center, :, :],
               cmap="gist_gray_r",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[0].axis(False)
ax[0].set_title(r"$\textrm{Quarter Resolution}$")

# Line:
x = np.load("../data/LTE/x_regular_quarter.npy")
pix2Mm = (x.max() - x.min())*1e-6/len(x)
ax[0].hlines(y=4, xmin=6, xmax=6 + 1/pix2Mm, lw=2, color='w',
             path_effects=[pe.Stroke(linewidth=3, foreground="black"),pe.Normal()])

# Text:
ax[0].text(6, 6, r"\textbf{1 Mm}", color='w', fontsize=12,
           path_effects=[pe.Stroke(linewidth=1, foreground="black"),pe.Normal()])

ax[1].imshow(intensity_third[center, :, :],
               cmap="gist_gray_r",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[1].axis(False)
ax[1].set_title(r"$\textrm{Third Resolution}$")

im = ax[2].imshow(intensity_half[center, :, :],
               cmap="gist_gray_r",
               origin="lower",
               vmax=CMAX,
               vmin=CMIN)
ax[2].axis(False)
ax[2].set_title(r"$\textrm{Half Resolution}$")

fig.colorbar(im, fraction=0.05, pad=0.06)

plt.savefig("../img/compare_line/disk_centre_res.pdf")
"""

"""
fig, ax = plt.subplots(1, 2, figsize=(9, 3))

ax[0].plot(convergence_quarter, label=r"$\textrm{quarter~resolution}$")
ax[0].plot(convergence_third, label=r"$\textrm{third~resolution}$")
ax[0].plot(convergence_half, label=r"$\textrm{half~resolution}$")
ax[0].set_xlabel(r"$\textrm{Iteration}$")
ax[0].set_ylabel(r"$\textrm{Relative Change}$")
ax[0].set_yscale("log")
ax[0].legend()

ax[1].plot(convergence_5e5, label=r"$10^5~\textrm{sites}$")
ax[1].set_xlabel(r"$\textrm{Iteration}$")
ax[1].set_ylabel(r"$\textrm{Relative~Change}$")
ax[1].set_yscale("log")
ax[1].legend()

fig.tight_layout()
plt.savefig("../img/compare_line/convergence.pdf")
"""

"""
fig, ax = plt.subplots(1, 3, figsize=(9,3), constrained_layout=True)

ax[0].imshow(intensity_des_5e5[center, :, :],
          cmap="gist_gray_r",
          origin="lower",
          vmax=CMAX,
          vmin=CMIN)
ax[0].set_title(r"$\textrm{Destruction}$")
ax[0].axis(False)

# Line:
x = np.load("../data/LTE/x_regular_full.npy")
pix2Mm = (x.max() - x.min())*1e-6/len(x)
ax[0].hlines(y=4, xmin=6, xmax=6 + 1/pix2Mm, lw=2, color='w',
             path_effects=[pe.Stroke(linewidth=3, foreground="black"),pe.Normal()])

# Text:
ax[0].text(6, 6, r"\textbf{1 Mm}", color='w', fontsize=12,
           path_effects=[pe.Stroke(linewidth=1, foreground="black"),pe.Normal()])

ax[1].imshow(intensity_ext_5e5[center, :, :],
          cmap="gist_gray_r",
          origin="lower",
          vmax=CMAX,
          vmin=CMIN)

ax[1].set_title(r"$\textrm{Total~Extinction}$")
ax[1].axis(False)

im = ax[2].imshow(intensity_1e6[center, :, :],
          cmap="gist_gray_r",
          origin="lower",
          vmax=CMAX,
          vmin=CMIN)
ax[2].set_title(r"$\textrm{Continuum~Extinction}$")
ax[2].axis(False)

fig.colorbar(im, fraction=0.05, pad=0.06)
#plt.show()
plt.savefig("../img/compare_line/disk_centre_sites.pdf")
"""

fig, ax = plt.subplots(1, 2, figsize=(6,3), constrained_layout=True)

ax[0].imshow(intensity_1e6[right_wing, :, :],
          cmap="gist_gray_r",
          origin="lower",
          vmax=CMAX,
          vmin=CMIN)
ax[0].set_title(r"$\textrm{1~to~1}$")
ax[0].axis(False)

# Line:
x = np.load("../data/LTE/x_regular_full.npy")
pix2Mm = (x.max() - x.min())*1e-6/len(x)
ax[0].hlines(y=4, xmin=6, xmax=6 + 1/pix2Mm, lw=2, color='w',
             path_effects=[pe.Stroke(linewidth=3, foreground="black"),pe.Normal()])

# Text:
ax[0].text(6, 6, r"\textbf{1 Mm}", color='w', fontsize=12,
           path_effects=[pe.Stroke(linewidth=1, foreground="black"),pe.Normal()])

im=ax[1].imshow(intensity_1e6_1dot5[right_wing, :, :],
          cmap="gist_gray_r",
          origin="lower",
          vmax=CMAX,
          vmin=CMIN)

ax[1].set_title(r"$\textrm{Upsampled}$")
ax[1].axis(False)

fig.colorbar(im, fraction=0.05, pad=0.06)
plt.show()
# plt.savefig("../img/compare_line/disk_centre_sampling.png")
