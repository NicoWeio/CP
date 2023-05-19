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

# path_base = "build/A1_b_"
# path_base = "build/A1_c_"
# path_base = "build/A1_d_"
path_base = "build/A1_e_"

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

# %% plot
fig, ax = plt.subplots()
im = ax.imshow(data_phi[0], vmin=0, vmax=1, cmap="hot")

# visualize E as vector field
# show every n-th vector
n = 1
x, y = np.meshgrid(np.arange(0, data_E.shape[1]), np.arange(0, data_E.shape[2]))
quiver = ax.quiver(x[::n, ::n], y[::n, ::n], data_E[0, ::n, ::n, 0], data_E[0, ::n, ::n, 1], color="white")


def animate(i):
    fig.suptitle(f"t = {i}")
    im.set_data(data_phi[i])
    quiver.set_UVC(data_E[i, ::n, ::n, 0], data_E[i, ::n, ::n, 1])
    return im


ani = FuncAnimation(fig, animate, frames=len(data_phi), interval=100, blit=False)
plt.show()
