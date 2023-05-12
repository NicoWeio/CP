# %%
import matplotlib.animation as animation
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def load_merged_r(filename):
    # Load the data file
    with open(filename, 'r') as f:
        data = f.read().strip()

    # Split into lines
    lines = data.split('\n')

    # For each line, remove the trailing ' | '
    lines = [line.rstrip(' | ') for line in lines]

    # For each line, split into coordinate pairs
    pairs = [line.split(' | ') for line in lines]

    # Convert each pair into a numpy array
    pairs_arr = [
        [np.array(pair.split(','), dtype=float) for pair in line]
        for line in pairs
    ]

    # Join everything into a single numpy array
    return np.array(pairs_arr)


# Load data
r_merged = load_merged_r("build/r_merged.txt")
data = r_merged
r_merged

# %%

# Create a random NumPy array of shape (timesteps, particles, xy)
# timesteps = 100
# particles = 10
# xy = 2
# data = np.random.rand(timesteps, particles, xy)

# Create a figure and axis object
fig, ax = plt.subplots()

# Define a function to update the plot for each frame

x_max = np.max(data[:, :, 0])
y_max = np.max(data[:, :, 1])

def update(frame):
    ax.clear()
    ax.set_xlim(0, x_max)
    ax.set_ylim(0, y_max)
    ax.set_title('Frame {}'.format(frame))
    ax.scatter(data[frame, :, 0], data[frame, :, 1])


# Create the animation using FuncAnimation
ani = animation.FuncAnimation(fig, update, frames=len(data), interval=100)

# Show the animation
plt.show()

# Save the animation (mp4)
ani.save('build/animation.mp4', writer='ffmpeg', fps=(100), dpi=600)

# %%
