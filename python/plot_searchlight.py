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


def get_intensity(fname):
    """
    fname::string
        name of datafile
    """

    intensity = np.load(PATH+fname)

    return intensity

def get_coordinates(method):
    """
    method::string
        voronoi or regular
    """

    x = np.load(PATH+"x_"+method+".npy")
    y = np.load(PATH+"y_"+method+".npy")

    return x, y

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

vmax = intensity_r.max()
fig, ax = plt.subplots(1, 2, figsize=(10,5), constrained_layout=True)
im1 = ax[0].imshow(intensity_i,
                   cmap="magma",
                   origin="lower",
                   vmax=vmax,
                   extent=i_extent)

im2 = ax[1].imshow(intensity_r,
                   cmap="magma",
                   origin="lower",
                   vmax=vmax,
                   extent=r_extent)

plt.colorbar(im2, fraction=0.046, pad=0.04)#im2, cax=1)
# plt.savefig("searchlight_2D.png", dpi=500)
plt.close()

X_r, Y_r = np.meshgrid(x_r[0:30], y_r[0:30])

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.grid(False)
ax.plot_surface(X_r, Y_r, intensity_r,
                cmap="magma",
                vmax = vmax,
                rstride=1,
                cstride=1)

plt.show()
X_i, Y_i = np.meshgrid(x_i[0:305], y_i[0:305])

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.grid(False)
ax.plot_surface(X_i, Y_i, intensity_i,
                cmap="magma",
                vmax = vmax,
                rstride=5,
                cstride=5)
# show it
plt.show()
