import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
plt.style.use('ggplot')

#e.g. line 0.916164,0.11087 | 0.339836,4.88557 | 0.819589,5.93698 |
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

# load data
data_local = load_merged_r('build/A2_localbest.txt')
data_full = load_merged_r('build/A2_full.txt')
data_global = load_merged_r('build/A2_globalbest.txt')

#print shape
print(f'local: {data_local.shape}')
print(f'full: {data_full.shape}')
print(f'global: {data_global.shape}')
NUM_PARTICLES = data_local.shape[1]
NUM_ITER = data_local.shape[0]

# b) plot global best and local best
plt.rcParams.update({'font.size': 14}) #set fontsize
plt.figure(figsize=(10, 10), dpi=200)
colors = plt.cm.get_cmap('tab10', NUM_PARTICLES)

iter_lin = np.arange(0, NUM_ITER)
for i in range(NUM_PARTICLES):
    plt.plot(iter_lin, data_local[:,i,2], label=f'local best {i+1}', color=colors(i), linewidth=1.5)
plt.plot(iter_lin, data_global[:,0,2], '--', label='global best', color='red', linewidth=2)

plt.xlim(0, 20)
plt.xlabel('iteration')
plt.ylabel('f(x,y)')
plt.legend()
plt.tight_layout()
plt.savefig("build/A2_b.pdf")

# c) 3D animation of the moving particles
def f(x):
    return (x[0] - 1.9)**2 + (x[1] - 2.1)**2 + 2*np.cos(4*x[0] + 2.8) + 3*np.sin(2*x[1] + 0.6)
plt.rcParams.update({'font.size': 10}) #set fontsize
fig = plt.figure(figsize=(10, 10), dpi=200)
ax = plt.axes(projection='3d', azim=30, elev=30)

# set limits
ax.set_xlim3d(-5, 5)
ax.set_ylim3d(-5, 5)
ax.set_zlim3d(0, 10)

# set labels
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('f(x,y)')

# animation function
def animate(i):
    print(f'Animate frame {i}', end='\r')
    # clear the axes
    ax.clear()

    # set limits
    ax.set_xlim3d(-5, 5)
    ax.set_ylim3d(-5, 5)
    ax.set_zlim3d(0, 100)

    # set labels
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('f(x,y)')

    # plot the function
    x = np.linspace(-5, 5, 100)
    y = np.linspace(-5, 5, 100)
    X, Y = np.meshgrid(x, y)
    Z = f([X, Y])
    ax.plot_surface(X, Y, Z, cmap='twilight', edgecolor='grey', alpha=0.8, facecolors=plt.cm.twilight(Z/np.max(Z)))
    #viridis, seismic

    # animate the trace of the particles
    for j in range(NUM_PARTICLES):
        ax.plot(data_local[:i,j,0], data_local[:i,j,1], data_local[:i,j,2], label=f'local best {j+1}', color=colors(j), linewidth=1.5, alpha=0.7)
    ax.plot(data_global[:i,0,0], data_global[:i,0,1], data_global[:i,0,2], '--', label='global best', color='red', linewidth=2)

    # set legend outside
    ax.legend(loc='upper left', fontsize=10)#bbox_to_anchor=(1.05, 1)

# animation show first 50 frames
anim = animation.FuncAnimation(fig, animate, frames=NUM_ITER, interval=100)
anim.save('build/A2_c.mp4', writer='ffmpeg', fps=20)

# snapshot of the last frame
animate(100)
plt.savefig("build/A2_c.pdf")
#plt.show()