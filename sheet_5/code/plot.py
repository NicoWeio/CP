import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import pandas as pd
from pathlib import Path

path = "build/u_n.csv"

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
# convert to numpy array
data_u = np.array(data_u)
#data_u = np.transpose(data_u, axes=(0, 2, 1)) #(t, x, y)
# take only every 10th time step
#data_u = data_u[::20]

print("Shape of data: ", data_u.shape)

# create 2d heatmap animation with shape of data_u = (t, x, y)
fig, ax = plt.subplots(figsize=(12, 8), dpi=300)
im = ax.imshow(data_u[0], origin="lower", cmap="hot")
ax.set_title(f"Time step: 0", fontsize=10, loc="left")
ax.axis("off")

def animate(i):
    print(f"Create frame number {i}", end="\r")
    im.set_data(data_u[i])
    ax.set_title(f"Time step: {i}", fontsize=10, loc="left")
    return im


anim = FuncAnimation(fig, animate, frames=data_u.shape[0], interval=100)
plt.show()
anim.save("build/A2_animation.mp4", writer="ffmpeg", fps=30)
