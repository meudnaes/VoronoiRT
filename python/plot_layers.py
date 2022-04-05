import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True

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

fig = plt.figure(figsize=(13.25,6), constrained_layout=True)
ax = fig.add_subplot(1, 2, 1, projection='3d')

ax.grid(False)

ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))


scatter = ax.scatter(sites[1,:],
                     sites[2,:],
                     sites[0,:],
                     c=colored_sites,
                     s=5)


ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.set_zlabel(r"$z$")
ax.set_title(r"$\textrm{Upward Layers}$", fontsize=16)

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
                     s=5)

ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.set_zlabel(r"$z$")
ax.set_title(r"$\textrm{Downward Layers}$", fontsize=16)



cmap = mpl.cm.rainbow
norm = mpl.colors.Normalize(vmin=1, vmax=c_max)
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                  ax=ax,
                  fraction=0.04,
                  pad=0.1)
cb.set_label(r'$\rm{layer}$', fontsize=14)

plt.savefig("layers.pdf")
plt.close();
