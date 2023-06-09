import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def map_log(N: int, x0: float, r: float):
    for _ in range(N):
        x0 = r*x0*(1-x0)

    return x0

def map_cube(N: int, x0: float, r: float):
    for _ in range(N):
        x0 = r*x0 - x0**3 

    return x0

if __name__ == "__main__":
    # set values
    # r space
    spacing = 1E-3
    min_r = 0.0
    max_r = 4.0
    
    dx = 0.01 # step size for initial values x0 in range [0,1]
    N = 1000 # number of iterations

    print(f'Number of iterations: {N}')
    print(f'Spacing: {spacing}')
    print(f'Min r: {min_r}')
    print(f'Max r: {max_r}')
    print(f'Step size for initial values x0: {dx} in range [0,1]')

    # create array with r values
    r_space_log = np.arange(min_r, max_r, spacing)
    r_space_cube = np.arange(0, 3, spacing)
    print(f'Shape of r: {r_space_log.shape}')
    
    x0_space = np.arange(0.0, 1.0+dx, dx) #create array with initial values x0 including 1

    # create empty array for x
    x_space_log = np.zeros((len(r_space_log), len(x0_space)))
    x_space_cube = np.zeros((len(r_space_cube), len(x0_space)))
    print(f'Shape of full x array: {x_space_log.shape}')

    # calc logistic map for each x0 and r
    for ind, x0 in enumerate(x0_space):
        x_space_log[:, ind] = map_log(N, x0, r_space_log)  # vectorized calculation
        x_space_cube[:, ind] = map_cube(N, x0, r_space_cube)  # vectorized calculation

    # remove data points: x < 0 or x > 1
    x_space_log[(x_space_log < 1E-8) | (x_space_log > 1.0)] = np.nan
    x_space_cube[(x_space_cube < 1E-8) | (x_space_cube > 1.0)] = np.nan
    

    # TODO: calc Feigenbaum constant
    

    # plot logistic map
    plt.figure(figsize=(10,8), dpi=200, rasterized=True)

    for ind, x0 in enumerate(x0_space):
        label_name = f'x0 = {x0:.1f}' if ind % 10 == 0 else None
        plt.plot(r_space_log, x_space_log[:, ind], '.', linewidth=0.2, markersize=0.2, label=label_name, color=f'C{ind}')

    plt.xlabel('r')
    plt.ylabel('x')
    plt.title('Bifurcation diagram (logistic map)')
    #plt.legend(fontsize=10, loc='upper left')
    plt.savefig("build/A1_bifurkation_log.pdf")

    # plot logistic map
    plt.figure(figsize=(10,8), dpi=200, rasterized=True)

    for ind, x0 in enumerate(x0_space):
        label_name = f'x0 = {x0:.1f}' if ind % 10 == 0 else None
        plt.plot(r_space_cube, x_space_cube[:, ind], '.', linewidth=0.2, markersize=0.2, label=label_name, color=f'C{ind}')

    plt.xlabel('r')
    plt.ylabel('x')
    plt.title('Bifurcation diagram (cubic map)')
    #plt.legend(fontsize=10, loc='upper left')
    plt.savefig("build/A1_bifurkation_cube.pdf")
    
    plt.show()