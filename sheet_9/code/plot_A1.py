import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('ggplot')

# data
df_bisection = pd.read_csv('build/A1_bisection.csv')
df_newton = pd.read_csv('build/A1_newton.csv')

# plot values x,y,z against step
plt.figure(figsize=(10, 6), dpi=200)
plt.plot(df_bisection['step'], df_bisection['x'], '.-', label=r'Bisection $x$')
plt.plot(df_bisection['step'], df_bisection['y'], '.-', label=r'Bisection $y$')
plt.plot(df_bisection['step'], df_bisection['z'], '.-', label=r'Bisection $z$')
plt.plot(df_newton['step'], df_newton['x'], '.-', label=r'Newton $x$')

plt.xlabel('Step')
plt.ylabel('Value')
plt.legend()
plt.savefig('build/A1_values.pdf')

# plot error against step
plt.figure(figsize=(10, 6), dpi=200)
plt.plot(df_bisection['step'][1:], df_bisection['error'][1:], '.-', label=r'Bisection')
plt.plot(df_newton['step'][1:], df_newton['error'][1:], '.-', label=r'Newton')

#plt.yscale('log')
plt.xlabel('Step')
plt.ylabel('Error')
plt.legend()
plt.savefig('build/A1_error.pdf')

#plt.show()