"""
Extracts data from searchlighttest, stored in npy arrays. Filenames have the
format "I_(theta)_(phi)_(method).npy", where theta and phi are the angles that
gives the vector the ray is coming from, in degrees. Method is either "regular"
or "voronoi"
"""

import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

plt.rcParams['text.usetex'] = True

PATH = "../data/searchlight_data/"

iunits = r"$\textrm{kW}~\textrm{m}^{-2}~\textrm{nm}^{-1}~\textrm{sr}^{-1}$"

def get_intensity(fname, path=PATH):
    """
    fname::string
        name of datafile
    """

    intensity = np.load(path+fname)

    return intensity

def get_coordinates(method, path=PATH):
    """
    method::string
        voronoi or regular
    """

    x = np.load(path+"x_"+method+".npy")
    y = np.load(path+"y_"+method+".npy")

    return x, y

def font_size(SMALL_SIZE=14, MEDIUM_SIZE=16, BIGGER_SIZE=18):
    # Elin
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=BIGGER_SIZE)    # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

if __name__ == "__main__":
    font_size()

    intensity_i = np.roll(get_intensity("I_160_45_voronoi.npy"), 25, (0, 1))
    intensity_i = intensity_i[:305,:305].T
    x_i, y_i = get_coordinates("voronoi")
    x_i -= x_i[25]
    y_i -= y_i[25]
    i_extent = [x_i[0], x_i[305], y_i[0], y_i[305]]

    intensity_r = np.roll(get_intensity("I_160_45_regular.npy"), 4, (0, 1))
    intensity_r = intensity_r[:30,:30].T
    x_r, y_r = get_coordinates("regular")
    x_r -= x_r[2]
    y_r -= y_r[2]
    r_extent = [x_r[0], x_r[30], y_r[0], y_r[30]]

    """
    Create two heatmaps for both searchlight tests in one figure
    """
    vmax = intensity_r.max()
    fig, ax = plt.subplots(1, 2, figsize=(11,5), constrained_layout=True)
    im1 = ax[0].imshow(intensity_r,
                       cmap="magma",
                       origin="lower",
                       vmax=vmax,
                       extent=r_extent)
    ax[0].set_xlabel(r"$x~\textrm{[m]}$"); ax[0].set_ylabel(r"$y~\textrm{[m]}$")

    im2 = ax[1].imshow(intensity_i,
                       cmap="magma",
                       origin="lower",
                       vmax=vmax,
                       extent=i_extent)
    ax[1].set_xlabel(r"$x~\textrm{[m]}$"); ax[1].set_ylabel(r"$y~\textrm{[m]}$")

    plt.colorbar(im2, fraction=0.046, pad=0.04, label=iunits)
    plt.savefig("../img/searchlight_2D.pdf")
    plt.close()

    """
    Create two surface plots for both searchlight tests in one figure
    Following: https://matplotlib.org/stable/gallery/mplot3d/subplot3d.html
    """
    X_r, Y_r = np.meshgrid(x_r[0:30], y_r[0:30])


    font_size(SMALL_SIZE=12, MEDIUM_SIZE=14, BIGGER_SIZE=16)

    # set up a figure twice as wide as it is tall
    fig = plt.figure(figsize=(10, 5.1), constrained_layout=True)#figsize=plt.figaspect(0.5))

    # set up the axes for the first plot
    ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    ax1.grid(False)
    ax1.set_box_aspect((np.ptp(X_r), np.ptp(Y_r), np.ptp(intensity_i)))
    ax1.set_xlabel(r"$x~\textrm{[m]}$")
    ax1.set_ylabel(r"$y~\textrm{[m]}$")
    ax1.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax1.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax1.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax1.plot_surface(X_r, Y_r, intensity_r,
                           cmap="magma",
                           vmax = vmax,
                           rstride=1,
                           cstride=1)

    # ax.set_zlim(-1.01, 1.01)

    X_i, Y_i = np.meshgrid(x_i[0:305], y_i[0:305])

    # set up the axes for the second plot
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    ax2.grid(False)
    ax2.set_box_aspect((np.ptp(X_i), np.ptp(Y_i), np.ptp(intensity_i)))
    ax2.set_xlabel(r"$x~\textrm{[m]}$")
    ax2.set_ylabel(r"$y~\textrm{[m]}$")
    ax2.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax2.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax2.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    surf = ax2.plot_surface(X_i, Y_i, intensity_i,
                    cmap="magma",
                    vmax = vmax,
                    rstride=5,
                    cstride=5)
    #axdist = 10
    #ax1.dist = axdist
    #ax2.dist = axdist
    fig.colorbar(surf, fraction=0.046, pad=0.04, label=iunits)
    plt.savefig("../img/compare_searchlight/searchlight_3D.pdf", bbox_inches='tight', pad_inches=0.0)
    plt.close()
