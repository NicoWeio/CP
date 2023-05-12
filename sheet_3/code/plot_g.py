import pandas as pd
import matplotlib.pyplot as plt

# read in b)set.tsv
df = pd.read_csv("build/d)g.tsv", sep="\t")
# 'rBin', 'g'

# plot these as a histogram (lineplot with step)
plt.figure(figsize=(8, 6), dpi=300)
plt.step(df["rBin"], df["g"], where='post')
plt.xlabel("$r$")
plt.ylabel("$g(r)$")
plt.grid()
plt.savefig("build/g.pdf")