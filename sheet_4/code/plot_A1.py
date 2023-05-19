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

# %% read data
txt = Path("build/A1_b.txt").read_text()
data = [
    [
        [
            float(q)
            for q in y_list.removesuffix(",").split(",")
        ]
        for y_list in row.removesuffix("|").split("|")
        # if y_list
    ]
    for row in txt.strip().split("\n")
]
# convert to numpy array
data = np.array(data)

# %% plot
fig, ax = plt.subplots()
im = ax.imshow(data[0], vmin=0, vmax=1, cmap="hot")

def animate(i):
    fig.suptitle(f"t = {i}")
    im.set_data(data[i])
    return im

ani = FuncAnimation(fig, animate, frames=len(data), interval=100, blit=False)
plt.show()
