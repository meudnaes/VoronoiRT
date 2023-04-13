import h5py
import numpy as np
#import cmasher as cmr
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import matplotlib.patches as patches

from matplotlib.colors import LogNorm
from plot_searchlight import font_size

def get_quantities(filename):
    with h5py.File(filename, "r") as f:
        # List all groups
        print("Keys: %s" % f.keys())

        temperature = f["temperature"][()].squeeze()
        Ne = f["electron_density"][()].squeeze()
        NH = f["hydrogen_populations"][()].squeeze()

    return temperature, Ne, NH

def get_coordinates(filename):
    with h5py.File(filename, "r") as f:
        # List all groups
        # print("Keys: %s" % f.keys())

        x = f["x"][()].flatten()
        y = f["y"][()].flatten()
        z = f["z"][()].flatten()

    return x, y, z

def B_field(filename):
    with h5py.File(filename, "r") as f:
        # List all groups
        # print("Keys: %s" % f.keys())

        B_x = f["B_x"][()].squeeze()
        B_y = f["B_y"][()].squeeze()
        B_z = f["B_z"][()].squeeze()

    return B_x, B_y, B_z

if __name__ == "__main__":
    plt.rcParams['text.usetex'] = True

    PATH = "/mn/stornext/u3/tiago/data/qs006023/"

    quarter = "bifrost_qs006023_s525_quarter.hdf5"
    half = "bifrost_qs006023_s525_half.hdf5"
    full = "bifrost_qs006023_s525.hdf5"

    font_size()

    x, y, z = get_coordinates(PATH+full)

    z_diff = z[:-1] - z[1:]
    print(np.average(z_diff))
    print(x[1]-x[0])
    print(y[1]-y[0])

    print(x[-1] - x[0])
    print(y[-1] - y[0])
    print(z[-1] - z[0])

    fig, ax = plt.subplots()
    ax.plot(z/1e6)
    ax.set_ylabel(r"\textrm{Height [Mm]}")
    ax.set_xlabel(r"\textrm{Index}")
    plt.close()

    surf = np.argmin(np.abs(z - 0))
    print("Surface", z[surf]/1e3)

    Bx, By, Bz = B_field(PATH+full)

    print("mean B: ", np.mean(np.sqrt(Bx[:,:,surf]**2 + By[:,:,surf]**2 + Bz[:,:,surf]**2))*1e3)

    temperature, Ne, NH = get_quantities(PATH+full)

    fig, ax = plt.subplots(2, 1, figsize=(6.75, 11), constrained_layout=True)

    im = ax[0].imshow(np.flip(temperature[:, :, surf].T, 1),
                      origin="lower",
                      norm=LogNorm(),
                      cmap="plasma")                     # "gist_heat"
    ax[0].set_xticks([])
    ax[0].set_yticks([])

    ax[0].set_title(r"$\textrm{Temperature at surface}$")
    cbar = plt.colorbar(im, ax=ax[0], fraction=0.046, pad=0.04)
    cbar.set_label(r"$T~\textrm{[K]}$", rotation=90)

    im = ax[1].imshow(np.flip((Bz[:, :, surf]*1e3).T, 1),
                      origin="lower",
                      cmap="gist_gray",
                      vmax=15,
                      vmin=-15)
                      # (Bz[:, :, surf]*1e3).max(),)
                      # (Bz[:, :, surf]*1e3).max())
    ax[1].set_xticks([])
    ax[1].set_yticks([])
    cbar = plt.colorbar(im, ax=ax[1], fraction=0.046, pad=0.04)
    cbar.set_label(r"$B_z~\textrm{[mT]}$", rotation=90)

    ax[1].set_title(r"$\textrm{Vertical magnetic field at surface}$")

    pix2Mm = (x.max() - x.min())*1e-6/len(x)

    # Scale:
    rect = patches.Rectangle(xy=[20, 14], width=1/pix2Mm, height=3, color='w',
                             path_effects=[pe.Stroke(linewidth=3, foreground="black"),pe.Normal()])
    ax[0].add_patch(rect)

    # Text:
    ax[0].text(27, 22, r"\textbf{1 Mm}", color='w', fontsize=14,
               path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])

    # Scale:
    rect = patches.Rectangle(xy=[20, 14], width=1/pix2Mm, height=3, color='w',
                             path_effects=[pe.Stroke(linewidth=3, foreground="black"),pe.Normal()])
    ax[1].add_patch(rect)


    # Text:
    ax[1].text(27, 22, r"\textbf{1 Mm}", color='w', fontsize=14,
               path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])

    plt.savefig("../img/atmos.pdf")
