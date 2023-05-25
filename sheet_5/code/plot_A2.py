import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.patches as patches
import numpy as np
import pandas as pd
from pathlib import Path

print("\nStart Animation of time evolution of u(x, y, t)...")

# PARAMETERS
path = "build/u_n.csv"
dx = 0.01
a = 1.5
b = 1.0

# read in data to numpy array with shape (t, x, y)
print(f'Reading data from "{path}" ...')
txt = Path(path).read_text()
data_u = [
    [
        [
            float(q)
            for q in y_list.removesuffix(",").split(",")
        ]
        for y_list in row.removesuffix("|").split("|")
    ]
    for row in txt.strip().split("\n")
]
data_u = np.array(data_u)
#data_u = np.transpose(data_u, axes=(0, 2, 1)) #(t, x, y)

print("Shape of data: ", data_u.shape)
box_height = data_u.shape[1]-1
box_width = data_u.shape[2]-1


# create 2d heatmap animation with shape of data_u = (t, x, y)
fig, ax = plt.subplots(figsize=(12, 8), dpi=200)
im = ax.imshow(data_u[0], origin="lower", cmap="RdGy", extent=[0, a, 0, b])
ax.set_title(f"Time step: 0", fontsize=10, loc="left")
#ax.axis("off")
fig.colorbar(im, ax=ax)

# set axis labels
ax.set_xlabel('x')
ax.set_ylabel('y')

# set axis limits
ax.set_xlim(0, a)
ax.set_ylim(0, b)

# create box around that shows boundaries
rect = patches.Rectangle((0, 0), a, b, linewidth=2, edgecolor='k', facecolor='none')
ax.add_patch(rect)

# Animation
def animate(i):
    print(f"Create frame number {i}", end="\r")
    im.set_data(data_u[i])
    ax.set_title(f"Time step: {i}", fontsize=10, loc="left")
    rect.set_width(a)
    rect.set_height(b)

    return im


anim = FuncAnimation(fig, animate, frames=np.arange(0, data_u.shape[0], 4), interval=100)
plt.show()
anim.save("build/A2_animation.mp4", writer="ffmpeg", fps=30)