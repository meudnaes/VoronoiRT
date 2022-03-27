"""
Extracts data from searchlighttest, stored in npy arrays. Filenames have the
format "I_(theta)_(phi)_(method).npy", where theta and phi are the angles that
gives the vector the ray is coming from, in degrees. Method is either "regular"
or "voronoi"
"""

import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

PATH = "../data/searchlight_data/"


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

if __name__ == "__main__":
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
    ax[0].set_xlabel("x"); ax[0].set_ylabel("y")

    im2 = ax[1].imshow(intensity_i,
                       cmap="magma",
                       origin="lower",
                       vmax=vmax,
                       extent=i_extent)
    ax[1].set_xlabel("x"); ax[1].set_ylabel("y")

    plt.colorbar(im2, fraction=0.05, pad=0.06)
    plt.savefig("../img/searchlight_2D.png", dpi=500)
    plt.close()

    """
    Create two surface plots for both searchlight tests in one figure
    Following: https://matplotlib.org/stable/gallery/mplot3d/subplot3d.html
    """
    X_r, Y_r = np.meshgrid(x_r[0:30], y_r[0:30])

    # set up a figure twice as wide as it is tall
    fig = plt.figure(figsize=(10, 6), constrained_layout=True)#figsize=plt.figaspect(0.5))

    # set up the axes for the first plot
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    ax.grid(False)
    ax.set_box_aspect((np.ptp(X_r), np.ptp(Y_r), np.ptp(intensity_i)))
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.plot_surface(X_r, Y_r, intensity_r,
                           cmap="magma",
                           vmax = vmax,
                           rstride=1,
                           cstride=1)

    # ax.set_zlim(-1.01, 1.01)

    X_i, Y_i = np.meshgrid(x_i[0:305], y_i[0:305])

    # set up the axes for the second plot
    ax = fig.add_subplot(1, 2, 2, projection='3d')
    ax.grid(False)
    ax.set_box_aspect((np.ptp(X_i), np.ptp(Y_i), np.ptp(intensity_i)))
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    surf = ax.plot_surface(X_i, Y_i, intensity_i,
                    cmap="magma",
                    vmax = vmax,
                    rstride=5,
                    cstride=5)
    fig.colorbar(surf, fraction=0.04, pad=0.06)
    plt.savefig("../img/searchlight_3D.png", dpi=500)
    plt.close()
