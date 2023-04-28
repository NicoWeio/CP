import matplotlib.pyplot as plt
import numpy as np

S, r, z, phi = np.genfromtxt("build/A1.csv", unpack=True, delimiter=",")

s_ = S[22:]
S_end = s_[np.where(np.amin(r[22:]) == r[22:])[0][0]]

plt.plot(S, r)
plt.vlines(S_end, ymax=np.amax(r), ymin=0, ls="--", colors="r")
plt.xlabel(r"$S'= \frac{S}{a}$")
plt.ylabel(r"$r'= \frac{r}{a}$")

# in matplotlibrc leider (noch) nicht möglich
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/A1_r.pdf')
plt.clf()


plt.plot(S, z)
plt.vlines(S_end, ymax=np.amax(z), ymin=0, ls="--", colors="r")
plt.xlabel(r"$S'= \frac{S}{a}$")
plt.ylabel(r"$z'= \frac{z}{a}$")


# in matplotlibrc leider (noch) nicht möglich
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/A1_z.pdf')
plt.clf()

plt.plot(S, phi)
plt.vlines(S_end, ymax=np.amax(phi), ymin=0, ls="--", colors="r")
plt.xlabel(r"$S'= \frac{S}{a}$")
plt.ylabel(r"$\psi$")

# in matplotlibrc leider (noch) nicht möglich
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/A1_psi.pdf')
plt.clf()

#
print(S_end)
print(np.amax(r[:np.where(np.amin(r[22:]) == r[22:])[0][0]]))
print(np.amax(z[:np.where(np.amin(r[22:]) == r[22:])[0][0]]))
print(np.amax(phi[:np.where(np.amin(r[22:]) == r[22:])[0][0]]))
