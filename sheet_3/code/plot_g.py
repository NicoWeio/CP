# %%
import pandas as pd
# import seaborn as sns
import matplotlib.pyplot as plt

# read in b)set.tsv
df = pd.read_csv("build/b)g.tsv", sep="\t")
# 'rBin', 'g'

# plot these as a histogram (lineplot with step)
plt.step(df["rBin"], df["g"], where='post')
plt.xlabel("$r$")
plt.ylabel("$g(r)$")
plt.grid()
plt.show()
