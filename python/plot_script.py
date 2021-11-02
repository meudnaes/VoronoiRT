"""
Script for plotting and testing. Works best when plotting interactively. Run
`voronoi.jl` first to generate sites.
"""

import numpy as np
import matplotlib.pyplot as plt

from numpy import cos, sin
from matplotlib.animation import ArtistAnimation
from IPython.display import HTML
from matplotlib_inline.backend_inline import set_matplotlib_formats

set_matplotlib_formats('svg')

# Plot sites layer by layer
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

# Plot intersection site
cell_ID = np.load("cell_ID.npy") - 1
cell_neighbours = np.load("cell_neighbours.npy") - 1

p_r = sites[:, cell_ID]
positions = sites[:, cell_neighbours]

# Do an angle
theta = 30*180/np.pi
phi = 10*180/np.pi

# Unit vector towards upwind direction of the ray
k = -np.array([cos(theta), cos(phi)*sin(theta), sin(phi)*sin(theta)])

# Recipe from Camps et. al. (2013)
r = p_r[:]
s = np.zeros(len(cell_neighbours), np.float64)
for i in range(len(cell_neighbours)):
    # Natural neighbor position
    p_i = positions[:,i]

    # Calculate plane biseting site and neighbor
    # Normal vector
    n = p_i - p_r

    # Point on the plane
    p = (p_i + p_r)/2

    s[i] = n@(p - r)/(n@k)

q_index = np.where(s > 0, s, np.inf).argmin()
s_q = s[q_index]
q = r + s_q*k

indices = np.delete(np.arange(0, len(cell_neighbours)), q_index)

fig = plt.figure(figsize=(12,12))
ax = fig.add_subplot(projection='3d')

ax.scatter(positions[1,indices],
           positions[2,indices],
           positions[0,indices],
           c="blue",
           s=35)

ax.scatter(p_r[1],
           p_r[2],
           p_r[0],
           c="red",
           s=35)

ax.scatter(positions[1,q_index],
           positions[2,q_index],
           positions[0,q_index],
           c="green",
           s=35)

ax.quiver(q[1],
          q[2],
          q[0],
          -s_q*k[1],
          -s_q*k[2],
          -s_q*k[0],
          color="red",
          linewidth=3)

ax.scatter(q[1],
           q[2],
           q[0],
           c="red",
           s=50,
           marker="x",
           linewidth=3)

ax.set_xlabel("x [Mm]")
ax.set_ylabel("y [Mm]")
ax.set_zlabel("z [Mm]")
ax.set_title("SC Ray Intersection")

plt.show()
# plt.savefig("../img/ray_intersection")
