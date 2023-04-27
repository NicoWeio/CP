import matplotlib.pyplot as plt
import numpy as np

S, r, z, phi =np.genfromtxt("build/A1.csv", unpack=True, delimiter=",")

plt.plot(S, r, label='Kurve')
plt.xlabel(r"$S'= \frac{S}{a}$")
plt.ylabel(r"$r'= \frac{r}{a}$")
plt.legend(loc='best')

# in matplotlibrc leider (noch) nicht möglich
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot.pdf')
plt.clf()

