import matplotlib.pyplot as plt
import numpy as np

# a)
# read in csv, last line is empty, one column
r_a = np.genfromtxt('build/A1_a.csv', delimiter=',', skip_footer=1)
r_b = np.genfromtxt('build/A1_b.csv', delimiter=',', skip_footer=1)
r_c = np.genfromtxt('build/A1_c.csv', delimiter=',', skip_footer=1)
r_d = np.genfromtxt('build/A1_d.csv', delimiter=',', skip_footer=1)

NUM_BINS = 10
ALPHA = 0.5
N = 10**5

plt.figure(figsize=(6, 4), dpi=200)
plt.hist(r_a, bins=NUM_BINS, label=r'$a) r0 = 1234, a = 20, c = 120, m = 6075$', alpha=ALPHA, histtype='step', linewidth=4.0)
plt.hist(r_b, bins=NUM_BINS, label=r'$b) r0 = 1234, a = 137, c = 187, m = 256$', alpha=ALPHA, histtype='step', linewidth=4.0)
plt.hist(r_c, bins=NUM_BINS, label=r'$c) r0 = 123456789, a = 65539, c = 0, m = 2^{31}$', alpha=ALPHA, histtype='step', linewidth=4.0)
plt.hist(r_d, bins=NUM_BINS, label=r'$d) r0 = 1234, a = 7^5, c = 0, m = 2^{31} − 1$', alpha=ALPHA, histtype='step', linewidth=4.0)

plt.title(r'Histogram of LCG RNG with $10^5$ numbers')
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

# title of full figure
fig.suptitle(r'Correlation between two immediately consecutively generated numbers (LCG):', fontsize='medium')
fig.tight_layout()
fig.savefig('build/A1_corr_lcg.pdf')


# TASK 2
# a) Box-Mueller
rnorm = np.genfromtxt('build/A2_boxmueller.csv')
mean_a = np.mean(rnorm)
var_a = np.var(rnorm)
plt.figure(figsize=(6, 4), dpi=200)
plt.hist(rnorm, bins=NUM_BINS, label=r'$\mu = 0, \sigma = 1$', alpha=ALPHA, histtype='step', linewidth=4.0)
plt.title(f'Box-Mueller RNG with mean {mean_a:.3f} and variance {var_a:.3f}')
plt.xlabel('Random numbers')
plt.ylabel('# of occurences')
plt.savefig('build/A2_hist_boxmueller.pdf') 

# b) Central limit theorem
rnorm = np.genfromtxt('build/A2_central.csv')
mean_b = np.mean(rnorm)
var_b = np.var(rnorm)
plt.figure(figsize=(6, 4), dpi=200)
plt.hist(rnorm, bins=NUM_BINS, label=r'$\mu = 0, \sigma = 1$', alpha=ALPHA, histtype='step', linewidth=4.0)
plt.title(f'Central limit theorem with mean {mean_b:.3f} and variance {var_b:.3f}')
plt.xlabel('Random numbers')
plt.ylabel('# of occurences')
plt.savefig('build/A2_hist_central.pdf')



plt.show()