"""
Using matplotlib, display an animated heatmap based on a CSV file, where each line is a point in time and contains nested arrays for x and y.

Example line:
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[0,0.5,0.625,0.65625,0.664062,0.666016,0.666504,0.666626,0.666656 …
"""
# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from pathlib import Path
import sys

if len(sys.argv) > 1:
    path_base = f"build/A1_{sys.argv[1]}_"
else:
    raise ValueError("Bitte Aufgabenteil spezifizieren: b/c/d/e")
# %% read data
txt = Path(path_base + "phi.txt").read_text()
data_phi = [
    [
        [
            float(q)
            for q in y_list.removesuffix(",").split(",")
        ]
        for y_list in row.removesuffix("|").split("|")
    ]
    for row in txt.strip().split("\n")
]
# convert to numpy array
data_phi = np.array(data_phi)

# %% read data
txt = Path(path_base + "E.txt").read_text()
data_E = [
    [
        [
            [
                float(q)
                for q in a.split(";")
            ]
            for a in y_list.removesuffix(",").split(",")
        ]
        for y_list in row.removesuffix("|").split("|")
    ]
    for row in txt.strip().split("\n")
]
# convert to numpy array
data_E = np.array(data_E)

# %% fix x/y axes
data_phi = np.swapaxes(data_phi, 1, 2)
data_E = np.swapaxes(data_E, 1, 2)
# NOTE: We also need to specify the origin="lower" parameter in imshow() to fix the axes

# %% plot
fig, ax = plt.subplots()
absmax = np.max(np.abs(data_phi))
im = ax.imshow(
    data_phi[0],
    # vmin=0, vmax=1,
    # vmin=np.min(data_phi), vmax=np.max(data_phi),
    vmin=-absmax, vmax=absmax,
    # cmap="hot",
    cmap="seismic",
    origin="lower",
)

# visualize E as vector field
# show every n-th vector
n = 1
x, y = np.meshgrid(np.arange(0, data_E.shape[1]), np.arange(0, data_E.shape[2]))
quiver = ax.quiver(x[::n, ::n], y[::n, ::n], data_E[0, ::n, ::n, 0], data_E[0, ::n, ::n, 1], color="black")
plt.colorbar(im)


def animate(i):
    fig.suptitle(f"t = {i+1}")  # NOTE: the initial state (t=0) is not actually included
    im.set_data(data_phi[i])
    quiver.set_UVC(data_E[i, ::n, ::n, 0], data_E[i, ::n, ::n, 1])
    return im


# %% plot final state
animate(len(data_phi) - 1)
plt.savefig(path_base + "final.pdf")

# %% animate
ani = FuncAnimation(fig, animate, frames=len(data_phi), interval=100, blit=False)
# plt.show()
ani.save(path_base + "anim.mp4", dpi=300)
