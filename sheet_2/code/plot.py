import matplotlib.pyplot as plt
import numpy as np

S, r, z, phi =np.genfromtxt("build/A1.csv", unpack=True, delimiter=",")

plt.plot(r, z, label='Kurve')
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.legend(loc='best')

# in matplotlibrc leider (noch) nicht möglich
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot.pdf')


