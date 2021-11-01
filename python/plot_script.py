"""
Script to plot sites by layer. Works best when plotting interactively. Run
`voronoi.jl` first to generate sites.
"""

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.animation import ArtistAnimation
from IPython.display import HTML
from matplotlib_inline.backend_inline import set_matplotlib_formats

set_matplotlib_formats('svg')

layer = np.load("layer.npy") - 1
sites = np.load("p.npy")/1e6

colors = ["r",
          "orange",
          "y",
          "g",
          "b",
          "purple",
          "k"]

colored_sites = np.array(colors)[layer]

fig = plt.figure(figsize=(15,15))
ax = fig.add_subplot(projection='3d')

ax.scatter(sites[1,:],
           sites[2,:],
           sites[0,:],
           c=colored_sites,
           s=35)

ax.set_xlabel("x [Mm]")
ax.set_ylabel("y [Mm]")
ax.set_zlabel("z [Mm]")
ax.set_title("Sites by layer")

plt.show()
# plt.savefig("../img/sites_by_layer")
