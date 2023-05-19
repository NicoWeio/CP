# %%
import numpy as np
import matplotlib.pyplot as plt

# x, y = np.meshgrid(np.arange(0, 20), np.arange(0, 20))

# sum_elements = [
#     (2*(1-np.cos(n*np.pi))) / (n*np.pi * np.sinh(n*np.pi)) * np.sin(n*np.pi*x) * np.sinh(n*np.pi*y)
#     for n in range(1, 100)
# ]

# sum_elements_np = np.array(sum_elements)

# # phi = np.sum(sum_elements, axis=0)

phi = np.zeros((20, 20))
for x_index in range(20):
    for y_index in range(20):
        for n in range(1, 200):  # NOTE: higher n → overflows 😬
            x = x_index / 20
            y = y_index / 20
            phi[x_index, y_index] += (2*(1-np.cos(n*np.pi))) / (n*np.pi * np.sinh(n*np.pi)) * \
                np.sin(n*np.pi*x) * np.sinh(n*np.pi*y)

# fix x/y axes
phi = np.swapaxes(phi, 0, 1)

# %% plot


def plot_phi(my_phi):
    fig, ax = plt.subplots()
    absmax=np.max(np.abs(my_phi))
    im = ax.imshow(
        my_phi,
        # vmin=0, vmax=1,
        # vmin=np.min(my_phi), vmax=np.max(my_phi),
        vmin=-absmax, vmax=absmax,
        # cmap="hot",
        cmap="seismic",
        origin="lower",
        # extent=[0, 20, 0, 20],
    )

    plt.colorbar(im)
    return fig


plot_phi(phi).savefig("build/A1_c_analytical.png")

# %%
