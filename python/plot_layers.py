import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from plot_searchlight import font_size

plt.rcParams['text.usetex'] = True

font_size()

sites = np.load("sites.npy")
layers_up = np.load("layers_up.npy") - 1
perm_up = np.load("perm_up.npy") - 1
layers_down = np.load("layers_down.npy") - 1
perm_down = np.load("perm_down.npy") - 1

if len(layers_up) >= len(layers_down):
  c_max = len(layers_up)
else:
  c_max = len(layers_down)

colors = cm.rainbow(np.linspace(0, 1, c_max))

colored_sites = np.zeros((len(perm_up), 4))
for i in range(len(perm_up)):
    idx = perm_up[i]
    layer = np.searchsorted(layers_up, i)
    colored_sites[idx,:] = colors[layer]

fig = plt.figure(figsize=(9.3,4), constrained_layout=True)
ax = fig.add_subplot(1, 2, 1, projection='3d')

ax.grid(False)

ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))


scatter = ax.scatter(sites[1,:],
                     sites[2,:],
                     sites[0,:],
                     c=colored_sites,
                     s=3)


ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.set_zlabel(r"$z$")
ax.set_title(r"$\textrm{Upward Layers}$")

colors = cm.rainbow(np.linspace(0, 1, c_max))

colored_sites = np.zeros((len(perm_down), 4))
for i in range(len(perm_down)):
    idx = perm_down[i]
    layer = np.searchsorted(layers_down, i)
    colored_sites[idx,:] = colors[layer]

ax = fig.add_subplot(1, 2, 2, projection='3d')

ax.grid(False)

ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

scatter = ax.scatter(sites[1,:],
                     sites[2,:],
                     sites[0,:],
                     c=colored_sites,
                     s=3)

ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.set_zlabel(r"$z$")
ax.set_title(r"$\textrm{Downward Layers}$")

cmap = mpl.cm.rainbow
bounds = np.linspace(0, c_max, c_max+1)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                  ax=ax,
                  fraction=0.046,
                  pad=0.2)
cb.set_label(r'$\rm{layer}$')

plt.savefig("../img/layers.pdf")
plt.close();
