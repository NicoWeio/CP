import matplotlib.pyplot as plt
import numpy as np

def map(N: int, x0: float, r: float):
    for n in N:
        x0 = r*x0*(1-x0)

    return x0


if __name__ == "__main__":
    # set values
    spacing = 1E-3
    min_r = 0.0
    max_r = 4.0

    r = np.arange(0, 4, spacing)
    print(r.shape)