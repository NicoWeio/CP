import numpy as np
from scipy.optimize import minimize

def f(x):
    #f(x,y) = (x − 1.9)2 + (y − 2.1)2 + 2 cos(4x + 2.8) + 3 sin(2y + 0.6)
    return (x[0] - 1.9)**2 + (x[1] - 2.1)**2 + 2*np.cos(4*x[0] + 2.8) + 3*np.sin(2*x[1] + 0.6)

# minimize
res = minimize(f, [1.6, 2.0])
print(f"Minimum: {res.x} at {res.fun}")