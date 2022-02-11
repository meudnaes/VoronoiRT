"""
Extracts data from searchlighttest, stored in npy arrays. Filenames have the
format "I_(theta)_(phi)_(method).npy", where theta and phi are the angles that
gives the vector the ray is coming from, in degrees. Method is either "regular"
or "voronoi"
"""

import numpy as np
import matplotlib.pyplot as plt

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

intensity = get_intensity("I_20_15_voronoi.npy")
x, y = get_coordinates("voronoi")

plt.imshow(intensity.T, cmap="magma", origin="lower")
plt.show()
