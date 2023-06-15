import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
plt.style.use('ggplot')

# constants
NUM_BINS_EX1 = 10
NUM_BINS_EX2 = 20
ALPHA = 0.5
N = 10**5

# TASK 1
# a)
# read in csv, last line is empty, one column
r_a = np.genfromtxt('build/A1_a.csv')
r_b = np.genfromtxt('build/A1_b.csv')
r_c = np.genfromtxt('build/A1_c.csv')
r_d = np.genfromtxt('build/A1_d.csv')

plt.figure(figsize=(6, 4), dpi=200)
plt.hist(r_a, bins=NUM_BINS_EX1, label=r'$a) r0 = 1234, a = 20, c = 120, m = 6075$', alpha=ALPHA, histtype='step', linewidth=2.0)
plt.hist(r_b, bins=NUM_BINS_EX1, label=r'$b) r0 = 1234, a = 137, c = 187, m = 256$', alpha=ALPHA, histtype='step', linewidth=2.0)
plt.hist(r_c, bins=NUM_BINS_EX1, label=r'$c) r0 = 123456789, a = 65539, c = 0, m = 2^{31}$', alpha=ALPHA, histtype='step', linewidth=2.0)
plt.hist(r_d, bins=NUM_BINS_EX1, label=r'$d) r0 = 1234, a = 7^5, c = 0, m = 2^{31} − 1$', alpha=ALPHA, histtype='step', linewidth=2.0)

#plt.title(r'Histogram of LCG RNG with $10^5$ numbers')
plt.xlabel('Random numbers')
plt.ylabel('# of occurences')
plt.legend(fontsize='medium', loc='lower center')
plt.tight_layout()
plt.savefig('build/A1_hist_lcg.pdf')

# b)
# pairs of consecutive values
pairs_a = np.column_stack((r_a[1:], r_a[:-1]))  # (rn, rn-1) pairs
pairs_b = np.column_stack((r_b[1:], r_b[:-1]))
pairs_c = np.column_stack((r_c[1:], r_c[:-1]))
pairs_d = np.column_stack((r_d[1:], r_d[:-1]))

# plot the pairs as points in a scatter plot
fig, ax = plt.subplots(2, 2, figsize=(7, 7), dpi=200)
ax[0,0].scatter(pairs_a[:6075, 0], pairs_a[:6075, 1], s=0.2, color='C0')
ax[0,1].scatter(pairs_b[:256, 0], pairs_b[:256, 1], s=0.2, color='C1')
ax[1,0].scatter(pairs_c[:N//2, 0], pairs_c[:N//2, 1], s=0.2, color='C2', rasterized=True)
ax[1,1].scatter(pairs_d[:N//2, 0], pairs_d[:N//2, 1], s=0.2, color='C3', rasterized=True)

ax[0,0].set_title(r'$a) r0 = 1234, a = 20, c = 120, m = 6075$', fontsize='small')
ax[0,1].set_title(r'$b) r0 = 1234, a = 137, c = 187, m = 256$', fontsize='small')
ax[1,0].set_title(r'$c) r0 = 123456789, a = 65539, c = 0, m = 2^{31}$', fontsize='small')
ax[1,1].set_title(r'$d) r0 = 1234, a = 7^5, c = 0, m = 2^{31} − 1$', fontsize='small')

ax[0,0].set_xlabel(r'$r_n$')
ax[0,1].set_xlabel(r'$r_n$')
ax[1,0].set_xlabel(r'$r_n$')
ax[1,1].set_xlabel(r'$r_n$')

ax[0,0].set_ylabel(r'$r_{n-1}$')
ax[0,1].set_ylabel(r'$r_{n-1}$')
ax[1,0].set_ylabel(r'$r_{n-1}$')
ax[1,1].set_ylabel(r'$r_{n-1}$')

ax[0,0].grid()
ax[0,1].grid()
ax[1,0].grid()
ax[1,1].grid()

# title of full figure
#fig.suptitle(r'Correlation between two immediately consecutively generated numbers (LCG):', fontsize='medium')
fig.tight_layout()
fig.savefig('build/A1_corr_lcg.pdf')



# TASK 2
# numpy normal distribution
np.random.seed(1234)
rnorm0 = np.random.normal(size=N)

# a) Box-Mueller
rnorm1 = np.genfromtxt('build/A2_boxmueller.csv')
mean_a = np.mean(rnorm1)
var_a = np.var(rnorm1)
plt.figure(figsize=(6, 4), dpi=200)
plt.hist(rnorm1, bins=NUM_BINS_EX2, density=True, label='Box-Mueller method', alpha=ALPHA, histtype='step', linewidth=4.0)
plt.hist(rnorm0, bins=NUM_BINS_EX2, density=True, label='numpy.random', alpha=ALPHA, histtype='step', linewidth=4.0)
plt.plot(np.linspace(-4,4,1000), norm.pdf(np.linspace(-4,4,1000)), linestyle='--', label=r'$f(x) = \frac{1}{\sqrt{2\pi}} \exp(-\frac{x^2}{2})$')
#plt.title(f'Box-Mueller RNG with mean {mean_a:.3f} and variance {var_a:.3f}')
plt.xlabel('Random numbers')
plt.ylabel('Probability density')
plt.legend(fontsize='medium')
plt.savefig('build/A2_hist_boxmueller.pdf') 

# b) Central limit theorem
rnorm2 = np.genfromtxt('build/A2_central.csv')
mean_b = np.mean(rnorm2)
var_b = np.var(rnorm2)
plt.figure(figsize=(6, 4), dpi=200)
plt.hist(rnorm2, bins=NUM_BINS_EX2, density=True, label='Central limit', alpha=ALPHA, histtype='step', linewidth=4.0)
plt.hist(rnorm0, bins=NUM_BINS_EX2, density=True, label='numpy.random', alpha=ALPHA, histtype='step', linewidth=4.0)
plt.plot(np.linspace(-4,4,1000), norm.pdf(np.linspace(-4,4,1000)), linestyle='--', label=r'$f(x) = \frac{1}{\sqrt{2\pi}} \exp(-\frac{x^2}{2})$')
#plt.title(f'Central limit with mean {mean_b:.3f} and variance {var_b:.3f}')
plt.xlabel('Random numbers')
plt.ylabel('Probability density')
plt.legend(fontsize='medium')
plt.savefig('build/A2_hist_central.pdf')

# c) Von Neumann rejection: sin(x)/2
r_neumann = np.genfromtxt('build/A2_neumann.csv')
plt.figure(figsize=(6, 4), dpi=200)
plt.hist(r_neumann, bins=NUM_BINS_EX2, density=True, label='Von Neumann rejection', alpha=ALPHA, histtype='step', linewidth=4.0)
plt.plot(np.linspace(0,np.pi,1000), 0.5*np.sin(np.linspace(0,np.pi,1000)), linestyle='--', label=r'$f(x) = \frac{1}{2} \sin(x)$')
#plt.title(r'Von Neumann rejection for $f(x) = \frac{1}{2} \sin(x)$')
plt.xlabel('Random numbers')
plt.ylabel('Probability density')
plt.xticks([0, np.pi/2, np.pi], [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$'])
plt.legend()
plt.tight_layout()
plt.savefig('build/A2_hist_neumann.pdf')


# d) Inversion method: 3x^2
r_inv = np.genfromtxt('build/A2_inverse.csv')
plt.figure(figsize=(6, 4), dpi=200)
plt.hist(r_inv, bins=NUM_BINS_EX2, label='Inversion method', density=True, alpha=ALPHA, histtype='step', linewidth=4.0)
plt.plot(np.linspace(0,1,1000), 3*np.linspace(0,1,1000)**2, linestyle='--', label=r'$f(x) = 3x^2$')
#plt.title(r'Inversion method for $f(x) = 3x^2$')
plt.xlabel('Random numbers')
plt.ylabel('Probability density')
plt.legend()

plt.savefig('build/A2_hist_inversion.pdf')

#plt.show()