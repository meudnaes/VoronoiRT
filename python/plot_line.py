import numpy as np
import matplotlib as mpl
# import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from brightness_temperature import *
from plot_searchlight import get_intensity, font_size, iunits

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


center = np.argmin(np.abs(wavelength - lambda0))
blue_wing = center - 11
red_wing = center + 11
continuum = np.argmax(wavelength)

plt.rcParams['text.usetex'] = True

PATH = "./linedata/"

CMAX = 100
CMIN = 0

CMAX_wing = 80
CMIN_wing = 0

CMAP = "gist_gray_r"
CMAP_CONT = "gist_gray_r"

lpad = 8

if __name__ == "__main__":

    intensity_half = get_intensity("half_res_ul7n12_disk_centre_1.npy", PATH)
    intensity_third = get_intensity("regular_third_disk_centre.npy", PATH)
    intensity_quarter = get_intensity("regular_quarter_disk_centre.npy", PATH)

    intensity_cont_ext_5e5 = get_intensity("voronoi_5e5_disk_centre.npy", PATH)
    intensity_cont_ext_5e5_1dot5 = get_intensity("voronoi_ul7n12_5e5_disk_centre_1dot5.npy", PATH)
    intensity_cont_ext_1e6 = get_intensity("voronoi_ul7n12_1e6_disk_centre_1.npy", PATH)
    intensity_cont_ext_1e6_1dot5 = get_intensity("voronoi_ul7n12_1e6_disk_centre_1dot5.npy", PATH)
    intensity_cont_ext_2e6 = get_intensity("voronoi_ul7n12_2e6_disk_centre_1.npy", PATH)
    intensity_cont_ext_3e6 = get_intensity("voronoi_ul7n12_3e6_disk_centre_1.npy", PATH)

    intensity_tot_ext_5e5 = get_intensity("total_ext_5e5_disk_centre_1.npy", PATH)
    intensity_tot_ext_1e6 = get_intensity("total_ext_1e6_disk_centre_1.npy", PATH)
    intensity_tot_ext_2e6 = get_intensity("total_ext_2e6_disk_centre_1.npy", PATH)
    intensity_tot_ext_3e6 = get_intensity("total_ext_3e6_disk_centre_1.npy", PATH)

    intensity_destruction_5e5 = get_intensity("destruction_5e5_disk_centre_1.npy", PATH)
    intensity_destruction_1e6 = get_intensity("destruction_1e6_disk_centre_1.npy", PATH)

    intensity_density_5e5 = get_intensity("density_5e5_disk_centre_1.npy", PATH)

    intensity_ionised_5e5 = get_intensity("ionised_hydrogen_5e5_disk_centre_1.npy", PATH)
    intensity_ionised_1e6 = get_intensity("ionised_hydrogen_1e6_disk_centre_1.npy", PATH)
    intensity_ionised_2e6 = get_intensity("ionised_hydrogen_2e6_disk_centre_1.npy", PATH)

    intensity_uniform_1e6 = get_intensity("uniform_1e6_disk_centre_1.npy", PATH)

    convergence_quarter = np.load(PATH+"regular_ul7n12_quarter.npy")
    convergence_half = np.load(PATH+"regular_ul7n12_half.npy")
    convergence_third = np.load(PATH+"regular_ul7n12_third.npy")

    convergence_cont_5e5 = np.load("./convergence/voronoi_ul7n12_5e5_convergence.npy")
    convergence_cont_1e6 = np.load("./convergence/voronoi_ul7n12_1e6_convergence.npy")
    convergence_cont_2e6 = np.load("./convergence/voronoi_ul7n12_2e6_convergence.npy")
    convergence_cont_3e6 = np.load("./convergence/voronoi_ul7n12_3e6_convergence.npy")

    convergence_ionised_5e5 = np.load("./convergence/ionised_hydrogen_5e5_convergence.npy")
    convergence_ionised_1e6 = np.load("./convergence/ionised_hydrogen_1e6_convergence.npy")
    convergence_ionised_2e6 = np.load("./convergence/ionised_hydrogen_2000000_convergence.npy")

    convergence_density_5e5 = np.load("./convergence/density_5e5_convergence.npy")

    convergence_destruction_5e5 = np.load("./convergence/destruction_5e5_convergence.npy")
    convergence_destruction_1e6 = np.load("./convergence/destruction_1e6_convergence.npy")

    convergence_tot_ext_5e5 = np.load("./convergence/total_ext_5e5_convergence.npy")
    convergence_tot_ext_1e6 = np.load("./convergence/total_ext_1e6_convergence.npy")
    convergence_tot_ext_2e6 = np.load("./convergence/total_ext_2e6_convergence.npy")
    convergence_tot_ext_3e6 = np.load("./convergence/total_ext_3e6_convergence.npy")

    convergence_uniform_1e6 = np.load("./convergence/uniform_1e6_convergence.npy")

    velocity = ((wavelength - lambda0)/lambda0*constants.c).to("km s-1")
    print("Velocity at blue wing: %.3f" %(velocity[blue_wing].value))
    print("Velocity at continuum: %.3f" %(velocity[continuum].value))

    CMAX_continuum = intensity_half[continuum, :, :].max()
    CMIN_continuum = intensity_half[continuum, :, :].min()

    font_size()

    # compare sampling methods
    fig, ax = plt.subplots(1, 2, figsize=(7.5,4), constrained_layout=True)

    # plot disk-centre intensity in wings and centre, and continuum
    ax[0].imshow(intensity_cont_ext_1e6[center, :, :],
                   cmap=CMAP,
                   origin="lower",
                   vmax=CMAX,
                   vmin=CMIN)
    ax[0].axis(False)
    ax[0].set_title(r"$\alpha^c~\textrm{sampling}$")

    im = ax[1].imshow(intensity_uniform_1e6[center, :, :],
                   cmap=CMAP,
                   origin="lower",
                   vmax=CMAX,
                   vmin=CMIN)
    ax[1].axis(False)
    ax[1].set_title(r"$U~\textrm{sampling}$")

    x = np.load("../data/LTE/x_regular_full.npy")
    pix2Mm = (x.max() - x.min())*1e-6/len(x)

    # Line:
    ax[0].hlines(y=16, xmin=20, xmax=20 + 1/pix2Mm, lw=3, color='w',
                 path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])

    # ax[0].vlines(x=(40+1/pix2Mm)/2, ymin=14, ymax=18, lw=1/pix2Mm-8.25, color='w',
                 # path_effects=[pe.Stroke(linewidth=1/pix2Mm-6.25, foreground="black"),pe.Normal()])

    # Text:
    ax[0].text(18, 20, r"\textbf{1 Mm}", color='w', fontsize=14,
               path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])

    # Line:
    ax[1].hlines(y=16, xmin=20, xmax=20 + 1/pix2Mm, lw=3, color='w',
                 path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])

    # Text:
    ax[1].text(18, 20, r"\textbf{1 Mm}", color='w', fontsize=14,
               path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])

    fig.colorbar(im, fraction=0.043, pad=0.04, label=iunits)

    fig.suptitle(r"$\textbf{Disk-centre intensity at line centre, irregular grid}$")
    plt.savefig("../img/compare_line/quick_compare.pdf")
    plt.close()

    ################################################################################
    ################################################################################
    ################################################################################
    # compare sampling methods
    fig, ax = plt.subplots(1, 4, figsize=(14.5,4), constrained_layout=True)

    # plot disk-centre intensity in wings and centre, and continuum
    ax[0].imshow(intensity_cont_ext_1e6[center, :, :],
                   cmap=CMAP,
                   origin="lower",
                   vmax=CMAX,
                   vmin=CMIN)
    ax[0].axis(False)
    ax[0].set_title(r"$\alpha^c~\textrm{sampling}$")

    ax[1].imshow(intensity_ionised_1e6[center, :, :],
                   cmap=CMAP,
                   origin="lower",
                   vmax=CMAX,
                   vmin=CMIN)
    ax[1].axis(False)
    ax[1].set_title(r"$N_\textrm{\small{H\,II}}^\textrm{\small{LTE}}~\textrm{sampling}$")

    ax[2].imshow(intensity_tot_ext_1e6[center, :, :],
                   cmap=CMAP,
                   origin="lower",
                   vmax=CMAX,
                   vmin=CMIN)
    ax[2].axis(False)
    ax[2].set_title(r"$\alpha^\textrm{tot}~\textrm{sampling}$")

    im = ax[3].imshow(intensity_destruction_1e6[center, :, :],
                   cmap=CMAP,
                   origin="lower",
                   vmax=CMAX,
                   vmin=CMIN)
    ax[3].axis(False)
    ax[3].set_title(r"$\varepsilon~\textrm{sampling}$")

    x = np.load("../data/LTE/x_regular_full.npy")
    pix2Mm = (x.max() - x.min())*1e-6/len(x)

    # Line:
    ax[0].hlines(y=16, xmin=20, xmax=20 + 1/pix2Mm, lw=3, color='w',
                 path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])

    # ax[0].vlines(x=(40+1/pix2Mm)/2, ymin=14, ymax=18, lw=1/pix2Mm-8.25, color='w',
                 # path_effects=[pe.Stroke(linewidth=1/pix2Mm-6.25, foreground="black"),pe.Normal()])

    # Text:
    ax[0].text(18, 20, r"\textbf{1 Mm}", color='w', fontsize=14,
               path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])

    # Line:
    ax[1].hlines(y=16, xmin=20, xmax=20 + 1/pix2Mm, lw=3, color='w',
                 path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])

    # Text:
    ax[1].text(18, 20, r"\textbf{1 Mm}", color='w', fontsize=14,
               path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])

    # Line:
    ax[2].hlines(y=16, xmin=20, xmax=20 + 1/pix2Mm, lw=3, color='w',
                 path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])


    # Text:
    ax[2].text(18, 20, r"\textbf{1 Mm}", color='w', fontsize=14,
               path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])

    # Line:
    ax[3].hlines(y=16, xmin=20, xmax=20 + 1/pix2Mm, lw=3, color='w',
                 path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])


    # Text:
    ax[3].text(18, 20, r"\textbf{1 Mm}", color='w', fontsize=14,
               path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])

    fig.colorbar(im, fraction=0.043, pad=0.04, label=iunits)

    fig.suptitle(r"$\textbf{Disk-centre intensity at line centre, irregular grid}$")
    # plt.show()
    plt.savefig("../img/compare_line/compare_sites.pdf")

    ################################################################################
    ################################################################################
    ################################################################################
    # Plot images over line wing, centre, and continuum, regular grid
    fig, ax = plt.subplots(1, 3, figsize=(13,4), constrained_layout=True)

    # plot disk-centre intensity in wings and centre, and continuum
    im = ax[0].imshow(intensity_half[blue_wing, :, :],
                   cmap=CMAP,
                   origin="lower",
                   vmax=CMAX_wing,
                   vmin=CMIN_wing)
    ax[0].axis(False)
    wl = wavelength[blue_wing]
    ax[0].set_title(r"$\textrm{Blue wing}~%.3f\,\textrm{nm}$" %wl)

    cbar = plt.colorbar(im, ax=ax[0], fraction=0.046, pad=0.04)
    cbar.set_label(iunits, rotation=90, labelpad=0)

    im = ax[1].imshow(intensity_half[center, :, :],
                   cmap=CMAP,
                   origin="lower",
                   vmax=CMAX,
                   vmin=CMIN)
    ax[1].axis(False)
    wl = wavelength[center]
    ax[1].set_title(r"$\textrm{Line centre}~%.3f\,\textrm{nm}$" %wl)

    cbar = plt.colorbar(im, ax=ax[1], fraction=0.046, pad=0.04)
    cbar.set_label(iunits, rotation=90, labelpad=0)

    im = ax[2].imshow(intensity_half[continuum, :, :],
                     cmap=CMAP_CONT,
                     origin="lower",
                     vmax=CMAX_continuum,
                     vmin=CMIN_continuum)
    ax[2].axis(False)
    wl = wavelength[continuum]
    ax[2].set_title(r"$\textrm{Continuum}~%.3f\,\textrm{nm}$" %wl)

    cbar = plt.colorbar(im, ax=ax[2], fraction=0.046, pad=0.04)
    cbar.set_label(iunits, rotation=90, labelpad=lpad)

    x = np.load("../data/LTE/x_regular_half.npy")
    pix2Mm = (x.max() - x.min())*1e-6/len(x)

    # Line:
    ax[0].hlines(y=8, xmin=10, xmax=10 + 1/pix2Mm, lw=3, color='w',
                 path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])

    # ax[0].vlines(x=(20+1/pix2Mm)/2, ymin=7, ymax=9, lw=1/pix2Mm, color='w',
                 # path_effects=[pe.Stroke(linewidth=1/pix2Mm, foreground="black"),pe.Normal()])

    # Text:
    ax[0].text(9, 10, r"\textbf{1 Mm}", color='w', fontsize=14,
               path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])

    # Line:
    ax[1].hlines(y=8, xmin=10, xmax=10 + 1/pix2Mm, lw=3, color='w',
                 path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])

    # ax[1].vlines(x=(20+1/pix2Mm)/2, ymin=7, ymax=9, lw=1/pix2Mm, color='w',
                 # path_effects=[pe.Stroke(linewidth=1/pix2Mm, foreground="black"),pe.Normal()])

    # Text:
    ax[1].text(9, 10, r"\textbf{1 Mm}", color='w', fontsize=14,
               path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])

    # Line:
    ax[2].hlines(y=8, xmin=10, xmax=10 + 1/pix2Mm, lw=3, color='w',
                 path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])

    # ax[2].vlines(x=(20+1/pix2Mm)/2, ymin=7, ymax=9, lw=1/pix2Mm, color='w',
                 # path_effects=[pe.Stroke(linewidth=1/pix2Mm, foreground="black"),pe.Normal()])

    # Text:
    ax[2].text(9, 10, r"\textbf{1 Mm}", color='w', fontsize=14,
               path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])

    fig.suptitle(r"$\textbf{Disk-centre intensity, regular Grid}$")
    plt.savefig("../img/compare_line/regular_disk_centre.pdf")


    ################################################################################
    ################################################################################
    ################################################################################
    # Plot images over line wing, centre, and continuum, irregular grid
    fig, ax = plt.subplots(1, 3, figsize=(13,4), constrained_layout=True)

    # plot disk-centre intensity in wings and centre, and continuum
    im = ax[0].imshow(intensity_cont_ext_3e6[blue_wing, :, :],
                   cmap=CMAP,
                   origin="lower",
                   vmax=CMAX_wing,
                   vmin=CMIN_wing)
    ax[0].axis(False)
    wl = wavelength[blue_wing]
    ax[0].set_title(r"$\textrm{Blue wing}~%.3f\,\textrm{nm}$" %wl)

    cbar = plt.colorbar(im, ax=ax[0], fraction=0.046, pad=0.04)
    cbar.set_label(iunits, rotation=90, labelpad=0)

    im = ax[1].imshow(intensity_cont_ext_3e6[center, :, :],
                   cmap=CMAP,
                   origin="lower",
                   vmax=CMAX,
                   vmin=CMIN)
    ax[1].axis(False)
    wl = wavelength[center]
    ax[1].set_title(r"$\textrm{Line centre}~%.3f\,\textrm{nm}$" %wl)

    cbar = plt.colorbar(im, ax=ax[1], fraction=0.046, pad=0.04)
    cbar.set_label(iunits, rotation=90, labelpad=0)

    im = ax[2].imshow(intensity_cont_ext_3e6[continuum, :, :],
                   cmap=CMAP_CONT,
                   origin="lower",
                   vmax=CMAX_continuum,
                   vmin=CMIN_continuum)
    ax[2].axis(False)
    wl = wavelength[continuum]
    ax[2].set_title(r"$\textrm{Continuum}~%.3f\,\textrm{nm}$" %wl)

    cbar = plt.colorbar(im, ax=ax[2], fraction=0.046, pad=0.04)
    cbar.set_label(iunits, rotation=90, labelpad=lpad)

    x = np.load("../data/LTE/x_regular_full.npy")
    pix2Mm = (x.max() - x.min())*1e-6/len(x)

    # Line:
    ax[0].hlines(y=16, xmin=20, xmax=20 + 1/pix2Mm, lw=3, color='w',
                 path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])

    # ax[0].vlines(x=(40+1/pix2Mm)/2, ymin=14, ymax=18, lw=1/pix2Mm-8.25, color='w',
                 # path_effects=[pe.Stroke(linewidth=1/pix2Mm-6.25, foreground="black"),pe.Normal()])

    # Text:
    ax[0].text(18, 20, r"\textbf{1 Mm}", color='w', fontsize=14,
               path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])

    # Line:
    ax[1].hlines(y=16, xmin=20, xmax=20 + 1/pix2Mm, lw=3, color='w',
                 path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])

    # ax[1].vlines(x=(40+1/pix2Mm)/2, ymin=14, ymax=18, lw=1/pix2Mm-8.25, color='w',
                 # path_effects=[pe.Stroke(linewidth=1/pix2Mm-6.25, foreground="black"),pe.Normal()])

    # Text:
    ax[1].text(18, 20, r"\textbf{1 Mm}", color='w', fontsize=14,
               path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])

    # Line:
    ax[2].hlines(y=16, xmin=20, xmax=20 + 1/pix2Mm, lw=3, color='w',
                 path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])

    # ax[2].vlines(x=(40+1/pix2Mm)/2, ymin=14, ymax=18, lw=1/pix2Mm-8.25, color='w',
                 # path_effects=[pe.Stroke(linewidth=1/pix2Mm-6.25, foreground="black"),pe.Normal()])

    # Text:
    ax[2].text(18, 20, r"\textbf{1 Mm}", color='w', fontsize=14,
               path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])

    fig.suptitle(r"$\textbf{Disk-centre intensity, irregular grid}$")
    plt.savefig("../img/compare_line/irregular_disk_centre.pdf")

    ################################################################################
    ################################################################################
    ################################################################################
    # compare regular resolutions
    fig, ax = plt.subplots(1, 3, figsize=(11.75,4), constrained_layout=True)

    ax[0].imshow(intensity_quarter[center, :, :],
                   cmap=CMAP,
                   origin="lower",
                   vmax=CMAX,
                   vmin=CMIN)
    ax[0].axis(False)
    ax[0].set_title(r"$\textrm{Quarter resolution}$")

    ax[1].imshow(intensity_third[center, :, :],
                   cmap=CMAP,
                   origin="lower",
                   vmax=CMAX,
                   vmin=CMIN)
    ax[1].axis(False)
    ax[1].set_title(r"$\textrm{One-third resolution}$")

    im = ax[2].imshow(intensity_half[center, :, :],
                   cmap=CMAP,
                   origin="lower",
                   vmax=CMAX,
                   vmin=CMIN)
    ax[2].axis(False)
    ax[2].set_title(r"$\textrm{Half resolution}$")

    # Line:
    x = np.load("../data/LTE/x_regular_quarter.npy")
    pix2Mm = (x.max() - x.min())*1e-6/len(x)
    ax[0].hlines(y=8/2, xmin=10/2, xmax=10/2 + 1/pix2Mm, lw=3, color='w',
                 path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])
    # Text:
    ax[0].text(9/2, 10/2, r"\textbf{1 Mm}", color='w', fontsize=14,
               path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])



    # Line:
    x = np.load("../data/LTE/x_regular_third.npy")
    pix2Mm = (x.max() - x.min())*1e-6/len(x)
    ax[1].hlines(y=8*2/3, xmin=10*2/3, xmax=10*2/3 + 1/pix2Mm, lw=3, color='w',
                 path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])

    # Text:
    ax[1].text(9*2/3, 10*2/3, r"\textbf{1 Mm}", color='w', fontsize=14,
               path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])

    # Line:
    x = np.load("../data/LTE/x_regular_half.npy")
    pix2Mm = (x.max() - x.min())*1e-6/len(x)
    ax[2].hlines(y=8, xmin=10, xmax=10 + 1/pix2Mm, lw=3, color='w',
                 path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])
    # Text:
    ax[2].text(9, 10, r"\textbf{1 Mm}", color='w', fontsize=14,
               path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])

    # fig.suptitle(r"$\textbf{Disk-centre Intensity line centre, Regular Grid}$")

    fig.colorbar(im, fraction=0.043, pad=0.04, label=iunits)

    plt.savefig("../img/compare_line/regular_resolutions.pdf")

    ################################################################################
    ################################################################################
    ################################################################################
    # plot convergence
    fig, ax = plt.subplots(1, 3, figsize=(14, 5.5), sharey=True)

    ax[0].plot(convergence_quarter, label=r"$\textrm{regular (1/4 res.)}$", color="k", ls="solid")
    ax[0].plot(convergence_ionised_5e5, label=r"$N_\textrm{H\,\small{II}}$", color="red", ls="solid")
    ax[0].plot(convergence_cont_5e5, label=r"$\alpha^\textrm{cont}$", color="blue", ls="dashed")
    ax[0].plot(convergence_tot_ext_5e5, label=r"$\alpha^\textrm{tot}$", color="gold", ls="solid")
    ax[0].plot(convergence_density_5e5, label=r"$\rho$", color="gray", ls="dashdot")
    ax[0].plot(convergence_destruction_5e5, label=r"$\varepsilon$", color="cyan", ls="solid")

    ax[1].plot(convergence_third, label=r"$\textrm{regular (1/3 res.)}$", color="k", ls="solid")
    ax[1].plot(convergence_destruction_1e6, label=r"$\varepsilon$", color="cyan", ls="solid")
    ax[1].plot(convergence_cont_1e6, label=r"$\alpha^\textrm{cont}$", color="blue", ls="dashed")
    ax[1].plot(convergence_ionised_1e6, label=r"$N_\textrm{H\,\small{II}}$", color="red", ls="solid")
    ax[1].plot(convergence_uniform_1e6, label=r"$U~\textrm{(uniform)}$", color="pink", ls="solid")

    ax[2].plot(convergence_half, label=r"$\textrm{regular (1/2 res.)}$", color="k", ls="solid")
    ax[2].plot(convergence_cont_3e6, label=r"$\alpha^\textrm{cont}$", color="blue", ls="dashed")

    # ax.plot(convergence_cont_2e6, label=r"$\alpha^\textrm{cont}~2\cdot 10^6~\textrm{sites}$", color="b", ls="dashdot")
    # ax.plot(convergence_tot_ext_2e6, label=r"$\alpha^\textrm{tot}~2\cdot 10^6~\textrm{sites}$", color="g", ls="dashdot")
    # ax.plot(convergence_tot_ext_1e6, label=r"$\alpha^\textrm{tot}~1\cdot 10^6~\textrm{sites}$", color="g", ls="dashed")

    ax[0].set_ylabel(r"$\textrm{Max rel. change,}~\max\left(1 - S_\textrm{new}/S_\textrm{old}\right)$")
    ax[0].set_yscale("log")

    ax[0].legend()
    ax[1].legend()
    ax[2].legend()

    ax[0].set_xlabel(r"$\textrm{Iterations}$")
    ax[1].set_xlabel(r"$\textrm{Iterations}$")
    ax[2].set_xlabel(r"$\textrm{Iterations}$")

    ax[0].set_title(r"$\sim 5\cdot 10^5~\textrm{points}$")
    ax[1].set_title(r"$\sim 10^6~\textrm{points}$")
    ax[2].set_title(r"$\sim 3\cdot10^6~\textrm{points}$")

    #ax.set_title(r"$\textrm{Convergence}$")
    fig.tight_layout()
    plt.savefig("../img/compare_line/convergence.pdf")

    ################################################################################
    ################################################################################
    ################################################################################
    # resolution irregular grid
    fig, ax = plt.subplots(1, 3, figsize=(11.75,4), constrained_layout=True)

    ax[0].imshow(intensity_cont_ext_5e5[center, :, :],
              cmap=CMAP,
              origin="lower",
              vmax=CMAX,
              vmin=CMIN)
    ax[0].set_title(r"$5\cdot 10^5~\textrm{sites}$")
    ax[0].axis(False)

    ax[1].imshow(intensity_cont_ext_2e6[center, :, :],
              cmap=CMAP,
              origin="lower",
              vmax=CMAX,
              vmin=CMIN)

    ax[1].set_title(r"$2 \cdot 10^6~\textrm{sites}$")
    ax[1].axis(False)

    im = ax[2].imshow(intensity_cont_ext_3e6[center, :, :],
              cmap=CMAP,
              origin="lower",
              vmax=CMAX,
              vmin=CMIN)
    ax[2].set_title(r"$3 \cdot 10^6~\textrm{sites}$")
    ax[2].axis(False)

    x = np.load("../data/LTE/x_regular_full.npy")
    pix2Mm = (x.max() - x.min())*1e-6/len(x)

    # Line:
    ax[0].hlines(y=16, xmin=20, xmax=20 + 1/pix2Mm, lw=3, color='w',
                 path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])

    # Text:
    ax[0].text(18, 20, r"\textbf{1 Mm}", color='w', fontsize=14,
               path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])

    # Line:
    ax[1].hlines(y=16, xmin=20, xmax=20 + 1/pix2Mm, lw=3, color='w',
                 path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])

    # Text:
    ax[1].text(18, 20, r"\textbf{1 Mm}", color='w', fontsize=14,
               path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])

    # Line:
    ax[2].hlines(y=16, xmin=20, xmax=20 + 1/pix2Mm, lw=3, color='w',
                 path_effects=[pe.Stroke(linewidth=5, foreground="black"),pe.Normal()])

    # Text:
    ax[2].text(18, 20, r"\textbf{1 Mm}", color='w', fontsize=14,
               path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])

    fig.colorbar(im, fraction=0.046, pad=0.04, label=iunits)
    # fig.suptitle(r"$\textbf{Disk-Centre~Intensity~\textit{I}}_{\lambda_0}$")
    # plt.show()
    plt.savefig("../img/compare_line/disk_centre_irregular_resolution.pdf")

    ################################################################################
    ################################################################################
    ################################################################################
    # plot all lines to highlight differences

    fig, ax = plt.subplots(1, 2, figsize=(10, 5.5), constrained_layout=True, sharey=True)

    I_regular = intensity_half.reshape(len(wavelength), -1)
    I_regular *= units.kW*units.m**(-2)*units.nm**(-1)

    I_irregular = intensity_cont_ext_3e6.reshape(len(wavelength), -1)
    I_irregular *= units.kW*units.m**(-2)*units.nm**(-1)


    Tb_regular = T_b(wavelength[:, np.newaxis]*units.nm, I_regular)
    Tb_irregular = T_b(wavelength[:, np.newaxis]*units.nm, I_irregular)

    ax[0].plot(wavelength[center-17:center+18],
               Tb_regular[center-17:center+18, ::4].value,
               color='k',
               lw=0.03,
               alpha=0.5,
               rasterized=True)
    ax[0].plot(wavelength[center-17:center+18],
               np.mean(Tb_regular[center-17:center+18], axis=1).value,
               color="crimson", label=r"$\textrm{spatial average}$")
    ax[0].axvline(lambda0, ls="dashed", color="royalblue", lw=0.75)
    ax[0].axvline(wavelength[blue_wing], ls="dashed", color="deepskyblue", lw=0.75)
    ax[0].set_xlabel(r"$\textrm{Wavelength [nm]}$")
    ax[0].set_ylabel(r"$\textrm{Brightness temperature [K]}$")
    ax[0].set_title(r"$\textrm{Regular grid}$")
    ax[0].legend(loc="upper right")
    ax[0].text(x=lambda0+0.001, y=6150, s=r"$\lambda_0$", color="royalblue")
    ax[0].text(x=wavelength[blue_wing]-0.006, y=6150,
               s=r"$\textrm{wing}$", color="deepskyblue", rotation="vertical")
    # ax[0].set_xticks(list(ax[0].get_xticks()) + [lambda0])
    # ax[0].set_xticklabels([r"$%.2f$" %x for x in list(ax[0].get_xticks())[:-1]] + [r"$\lambda_0$"])

    ax[1].plot(wavelength[center-17:center+18],
               Tb_irregular[center-17:center+18, ::16].value,
               color='k',
               lw=0.03,
               alpha=0.5,
               rasterized=True)
    ax[1].plot(wavelength[center-17:center+18],
               np.mean(Tb_irregular[center-17:center+18], axis=1).value,
               color="crimson", label=r"$\textrm{spatial average}$")
    ax[1].axvline(lambda0, ls="dashed", color="royalblue", lw=0.75)
    ax[1].axvline(wavelength[blue_wing], ls="dashed", color="deepskyblue", lw=0.75)
    ax[1].set_xlabel(r"$\textrm{Wavelength [nm]}$")
    ax[1].set_ylim(6000,12000)
    ax[1].set_title(r"$\textrm{Irregular grid}$")
    ax[1].legend(loc="upper right")
    ax[1].text(x=lambda0+0.001, y=6150, s=r"$\lambda_0$", color="royalblue")
    ax[1].text(x=wavelength[blue_wing]-0.006, y=6150,
               s=r"$\textrm{wing}$", color="deepskyblue", rotation="vertical")
    # ax[1].set_xticklabels([r"$%.2f$" %x for x in list(ax[0].get_xticks())[:-1]] + [r"$\lambda_0$"])
    # ax[1].set_xticks(list(ax[1].get_xticks()) + [lambda0])

    # fig.suptitle(r"$\textrm{Disk-Centre Intensity}$")
    plt.savefig("../img/compare_line/lines.pdf", dpi=300)
