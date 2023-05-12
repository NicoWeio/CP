# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# read in b)set.tsv
df = pd.read_csv("build/b)set.tsv", sep="\t")
# 't', 'T', 'Ekin', 'Epot', 'vSx', 'vSy'

# plot the energies
sns.lineplot(data=df, x="t", y="Ekin", label="Ekin")
# |Epot|
abs_Epot = df["Epot"].abs()
sns.lineplot(data=df, x="t", y=abs_Epot, label="|Epot|")
plt.xlabel("$t$")
plt.ylabel("$E$")
plt.grid()
plt.show()
# %%
