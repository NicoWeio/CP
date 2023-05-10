import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load data
r1 = pd.read_csv('build/r1.csv').to_numpy()
r2 = pd.read_csv('build/r2.csv').to_numpy()
r3 = pd.read_csv('build/r3.csv').to_numpy()
r4 = pd.read_csv('build/r4.csv').to_numpy()


plt.figure(figsize=(8,8), dpi=200)
plt.plot(r1 [:,0], r1 [:,1], 'ro', label='r1')
plt.show()

plt.figure(figsize=(8,8), dpi=200)
plt.plot(r2 [:,0], r2 [:,1], 'ro', label='r2')
plt.show()

plt.figure(figsize=(8,8), dpi=200)
plt.plot(r3 [:,0], r3 [:,1], 'ro', label='r3')
plt.show()

plt.figure(figsize=(8,8), dpi=200)
plt.plot(r4 [:,0], r4 [:,1], 'ro', label='r4')
plt.show()