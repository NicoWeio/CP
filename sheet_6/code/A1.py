import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def map_log(N: int, x0: float, r: float):
    for _ in range(N):
        x0 = r*x0*(1-x0)

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
    r_space = np.arange(min_r, max_r, spacing)
    print(f'Shape of r: {r_space.shape}')

    # create array with initial values x0 including 1
    x0_space = np.arange(0.0, 1.0+dx, dx)

    # create empty array for x
    x_space = np.zeros((len(r_space), len(x0_space)))
    print(f'Shape of full x array: {x_space.shape}')
    
    # create empty array for x
    x_space = np.zeros((len(r_space), len(x0_space)))
    print(f'Shape of full x array: {x_space.shape}')

    # calc logistic map for each x0 and r
    for ind, x0 in enumerate(x0_space):
        x_space[:, ind] = map_log(N, x0, r_space)  # vectorized calculation

    # remove data points: x < 0 or x > 1
    x_space[(x_space < 1E-8) | (x_space > 1.0)] = np.nan
    

    # TODO: calc Feigenbaum constant
    

    # plot
    plt.figure(figsize=(10,8), dpi=200)

    for ind, x0 in enumerate(x0_space):
        label_name = f'x0 = {x0:.1f}' if ind % 10 == 0 else None
        plt.plot(r_space, x_space[:, ind], '.', linewidth=0.2, markersize=0.2, label=label_name, color=f'C{ind}')

    #plt.legend(fontsize=5, loc='upper left')
    plt.xlabel('r')
    plt.ylabel('x')
    plt.title('Bifurcation diagram')
    plt.legend(fontsize=10, loc='upper left')
    plt.savefig("build/bifurkation_log.pdf")
    plt.show()
