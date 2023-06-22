import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('ggplot')
plt.rcParams.update({'font.size': 14}) #set fontsize

print("\nPlotting...\n")

# data
df_bisection = pd.read_csv('build/A1_bisection.csv')
df_newton = pd.read_csv('build/A1_newton.csv')

# plot values x,y,z against step
# increase pointer size
plt.figure(figsize=(10, 6), dpi=200)
plt.plot(df_bisection['step'], df_bisection['x'], '.-', label=r'Bisection $x$', linewidth=3.0, markersize=10)
plt.plot(df_bisection['step'], df_bisection['y'], '.-', label=r'Bisection $y$', linewidth=3.0, markersize=10)
plt.plot(df_bisection['step'], df_bisection['z'], '.-', label=r'Bisection $z$', linewidth=3.0, markersize=10)
plt.plot(df_newton['step'], df_newton['x'], '.-', label=r'Newton $x$', linewidth=3.0, markersize=10)

plt.xlabel('Step')
plt.ylabel('Value')
plt.legend()
plt.savefig('build/A1_values.pdf')

# plot error against step
plt.figure(figsize=(10, 6), dpi=200)
plt.plot(df_bisection['step'][1:], df_bisection['error'][1:], '.-', label=r'Bisection', linewidth=3.0, markersize=10)
plt.plot(df_newton['step'][1:], df_newton['error'][1:], '.-', label=r'Newton', linewidth=3.0, markersize=10)

#plt.yscale('log')
plt.xlabel('Step')
plt.ylabel('Error')
plt.legend()
plt.savefig('build/A1_error.pdf')

#plt.show()