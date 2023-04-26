import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

df_r1 = pd.read_csv("build/euler_r1.csv")
df_r2 = pd.read_csv("build/euler_r2.csv")


plt.figure(figsize=(6,6), dpi=200)
plt.plot(df_r1.iloc[:,0].to_numpy(), df_r1.iloc[:,1].to_numpy(), 'b.', label="Mass 1")
plt.plot(df_r2.iloc[:,0].to_numpy(), df_r2.iloc[:,1].to_numpy(), 'r.', label="Mass 2")

plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.savefig("build/A2a_euler.pdf")
plt.show()