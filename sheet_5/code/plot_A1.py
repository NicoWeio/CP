import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

print("\nStart Animation of time evolution of probability density...")
# PARAMETERS
xmin = -10
xmax = 10
dt = 0.02

# Read data_rho from csv file, last column drop
df_rho = pd.read_csv('build/A1_rho.csv', sep=',', header=None)
df_psiRe = pd.read_csv('build/A1_psiRe.csv', sep=',', header=None)
df_psiIm = pd.read_csv('build/A1_psiIm.csv', sep=',', header=None)

df_rho = df_rho.drop(df_rho.columns[[-1,]], axis=1)
df_psiRe = df_psiRe.drop(df_psiRe.columns[[-1,]], axis=1)
df_psiIm = df_psiIm.drop(df_psiIm.columns[[-1,]], axis=1)

data_rho = df_rho.to_numpy() #each line contains the data of one time step
data_psiRe = df_psiRe.to_numpy()
data_psiIm = df_psiIm.to_numpy()

# create 1D animation
fig, ax = plt.subplots(dpi = 200)

# set axis limits
ax.set_xlim(xmin, xmax)
ax.set_ylim(data_rho.min(), data_rho.max())

# set axis labels
ax.set_xlabel('x')
ax.set_ylabel(fr'$\psi$')

# set title
ax.set_title(r'Time evolution of $|\psi$(x, t = 0)$|^2$')

# create empty line
line, = ax.plot([], [], lw=3, color='black', label=r'$|\psi|^2$')
line2, = ax.plot([], [], lw=2, color='blue', label=r'$Re(\psi)$')
line3, = ax.plot([], [], lw=2, color='orange', label=r'$Im(\psi)$')

plt.legend()

# initialization function
def init():
    line.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    return line,

# animation function
def animate(i):
    print(f"Create frame number {i}", end="\r")
    x = np.linspace(xmin, xmax, len(data_rho[i]))

    line.set_data(x, data_rho[i])
    line2.set_data(x, data_psiRe[i])
    line3.set_data(x, data_psiIm[i])

    #update title
    ax.set_title(fr'Time evolution of $|\psi$(x, t = {i*dt:.2f})$|^2$')
    
    return line

# call the animator
# show only every 10th frame
anim = FuncAnimation(fig, animate, init_func=init, frames=np.arange(1, len(data_rho), 1), interval=1)
plt.show()
# save animation as mp4 but only hold every 10th frame
anim.save('build/A1_psi.mp4', writer='ffmpeg', fps=50)
#anim.save('A1_psi.mp4', fps=30, extra_args=['-vcodec', 'libx264'], savefig_kwargs={'pad_inches':0.05}, dpi=200)

# plot last frame
plt.plot(np.linspace(xmin, xmax, len(data_rho[-1])), data_rho[-1], color='green')
plt.xlabel('x')
plt.ylabel(fr'$|\psi|^2$')
plt.title(r'$|\psi$(x, t = 10)$|^2$')
plt.savefig('build/A1_t10.pdf', dpi=200)