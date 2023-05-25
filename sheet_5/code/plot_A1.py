import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# PARAMETERS
xmin = -10
xmax = 10
dt = 0.02

# Read data from csv file, last column drop
df = pd.read_csv('build/A1_psi.csv', sep=',', header=None)
df = df.drop(df.columns[[-1,]], axis=1)
data = df.to_numpy() #each line contains the data of one time step

# create 1D animation
fig, ax = plt.subplots(dpi = 200)

# set axis limits
ax.set_xlim(xmin, xmax)
ax.set_ylim(data.min(), data.max())

# set axis labels
ax.set_xlabel('x')
ax.set_ylabel(fr'$\psi$')

# set title
ax.set_title(r'Time evolution of $\psi$(x, t = 0)')

# create empty line
line, = ax.plot([], [], lw=2, color='green')

# initialization function
def init():
    line.set_data([], [])
    return line,

# animation function
def animate(i):
    x = np.linspace(xmin, xmax, len(data[i]))
    y = data[i]
    line.set_data(x, y)

    #update title
    ax.set_title(fr'Time evolution of $\psi$(x, t = {i*dt:.2f})')

    return line

# call the animator
# show only every 10th frame
anim = FuncAnimation(fig, animate, init_func=init, frames=np.arange(0, len(data), 1), interval=1)


plt.show()

# save animation as mp4 but only hold every 10th frame
anim.save('build/A1_psi.mp4', writer='ffmpeg', fps=100)
#anim.save('A1_psi.mp4', fps=30, extra_args=['-vcodec', 'libx264'], savefig_kwargs={'pad_inches':0.05}, dpi=200)

