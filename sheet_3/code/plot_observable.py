import pandas as pd
import matplotlib.pyplot as plt

# read in b)set.tsv
df = pd.read_csv("build/d)set.tsv", sep="\t")
# 't', 'T', 'Ekin', 'Epot', 'vSx', 'vSy'

# create subplots for (t,T) and (t, E)
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 6))
#plt.rcParams.update({'font.size': 18})


# plot T(t)
ax1.plot(df["t"], df["T"])
ax2.set_xlabel("$t$")
ax1.set_ylabel("$T$")
ax1.grid()


# plot E(t)
ax2.plot(df["t"], df["Ekin"] + df["Epot"], label="$E_{total}$")
ax2.plot(df["t"], df["Ekin"], label="$E_{kin}$")
ax2.plot(df["t"], df["Epot"], label="$E_{pot}$")

ax2.legend()
ax2.set_xlabel("$t$")
ax2.set_ylabel("$E$")
ax2.grid()

plt.tight_layout()
plt.savefig("build/observable.pdf")