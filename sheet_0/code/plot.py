import matplotlib.pyplot as plt
import numpy as np

#A3
#a)
t_euler1, y_euler1 = np.genfromtxt("build/A3_euler1.csv", unpack= True, delimiter = ",")
t_sym1, y_sym1 = np.genfromtxt("build/A3_symeuler1.csv", unpack= True, delimiter = ",")


fig = plt.figure()

# Create zoom-in plot
ax = plt.plot(t_euler1, y_euler1, "--", label='Euler')
ax = plt.plot(t_sym1, y_sym1, "--", label='symmetrischer Euler')
ax = plt.plot(t_sym1, np.exp(-t_euler1), "--", label='analytisch')
plt.xlabel('x', labelpad = 15)
plt.ylabel('y', labelpad = 15)
plt.xlabel(r'$t$')
plt.ylabel(r'$y$')
plt.title(r"$\Delta t = 0.1$")
plt.legend(loc='best')

# Create zoom-out plot
ax_new = fig.add_axes([0.4, 0.3, 0.47, 0.4]) # the position of zoom-out plot compare to the ratio of zoom-in plot 
plt.plot(t_euler1, y_euler1, ",")
plt.plot(t_sym1, y_sym1, ",")
plt.plot(t_sym1, np.exp(-t_euler1), ",")
plt.xlim(8, 10)
plt.ylim(-0.002, 0.002)
plt.savefig('build/a3_a.pdf')
plt.clf()

#differenz
plt.plot(t_euler1, np.exp(-t_euler1)- y_euler1, ",", label='Euler')
plt.plot(t_sym1, np.exp(-t_euler1)- y_sym1, ",", label='symmetrischer Euler')
plt.xlabel(r'$t$')
plt.ylabel(r'$y_{analytisch}-y$')
plt.title(r"$\Delta t = 0.1$, Differenz")
plt.legend(loc='best')

plt.tight_layout()
plt.savefig('build/a3_a_diff.pdf')
plt.clf()

#b)
t_euler2, y_euler2 = np.genfromtxt("build/A3_euler2.csv", unpack= True, delimiter = ",")
t_sym2, y_sym2 = np.genfromtxt("build/A3_symeuler2.csv", unpack= True, delimiter = ",")


fig = plt.figure()

# Create zoom-in plot
ax = plt.plot(t_euler2, y_euler2, "--", label='Euler')
ax = plt.plot(t_sym2, y_sym2, "--", label='symmetrischer Euler')
ax = plt.plot(t_sym2, np.exp(-t_euler2), "--", label='analytisch')
plt.xlabel('x', labelpad = 15)
plt.ylabel('y', labelpad = 15)
plt.xlabel(r'$t$')
plt.ylabel(r'$y$')
plt.title(r"$\Delta t = 0.1$")
plt.legend(loc='best')

## Create zoom-out plot
#ax_new = fig.add_axes([0.4, 0.3, 0.47, 0.4]) # the position of zoom-out plot compare to the ratio of zoom-in plot 
#plt.plot(t_euler1, y_euler1, ",")
#plt.plot(t_sym1, y_sym1, ",")
#plt.plot(t_sym1, np.exp(-t_euler1), ",")
#plt.xlim(8, 10)
#plt.ylim(-0.002, 0.002)
plt.savefig('build/a3_b.pdf')
plt.clf()

#differenz
plt.plot(t_euler2, np.exp(-t_euler2)- y_euler2, ",", label='Euler')
plt.plot(t_sym2, np.exp(-t_euler2)- y_sym2, ",", label='symmetrischer Euler')
plt.xlabel(r'$t$')
plt.ylabel(r'$y_{analytisch}-y$')
plt.title(r"$\Delta t = 0.1$, Differenz")
plt.legend(loc='best')
plt.savefig('build/a3_b_diff.pdf')
plt.clf()