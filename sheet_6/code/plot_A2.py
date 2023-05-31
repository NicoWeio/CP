import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# read in data
path1 = f"build/A2_r20.csv"
path2 = f"build/A2_r28.csv"

df1 = pd.read_csv(path1)
df2 = pd.read_csv(path2)

# plot x,y projection
# 1st initial condition
plt.figure(figsize=(8,8), dpi=200)
plt.scatter(df1["x"], df1["y"], s=0.1, c="black")
plt.xlabel("x")
plt.ylabel("y")
plt.title("x,y projection for r=20")
plt.savefig("build/A2_xyprojection1.pdf")
plt.show()

# 2nd initial condition
plt.figure(figsize=(8,8), dpi=300)
plt.scatter(df2["x"], df2["y"], s=0.1, c="black")
plt.xlabel("x")
plt.ylabel("y")
plt.title("x,y projection for r=28")
plt.savefig("build/A2_xyprojection2.pdf")
plt.show()

