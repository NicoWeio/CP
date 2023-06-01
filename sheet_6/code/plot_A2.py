import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# read in data
path1 = f"build/A2_r20.csv"
path2 = f"build/A2_r28.csv"

df1 = pd.read_csv(path1)
df2 = pd.read_csv(path2)

# plot x,y projection <<<<<<<<<<
# r=20
plt.figure(figsize=(8,8), dpi=200)
plt.scatter(df1["x"], df1["y"], s=0.1, c=df1["z"])
plt.xlabel("x")
plt.ylabel("y")
plt.title("x,y projection for r=20")
cbar = plt.colorbar()
cbar.ax.set_ylabel('z-value')
plt.savefig("build/A2_xy_r20.pdf")
plt.clf()

# r=28
plt.figure(figsize=(8,8), dpi=200)
plt.scatter(df2["x"], df2["y"], s=0.1, c=df2["z"])
plt.xlabel("x")
plt.ylabel("y")
plt.title("x,y projection for r=28")
cbar = plt.colorbar()
cbar.ax.set_ylabel('z-value')
plt.savefig("build/A2_xy_r28.pdf")
plt.clf()

# poincare section Z=20=const <<<<<<<<<<
# r=20
# calc linear interpolation between point before and after Z=20
# and find intersection with Z=20
# find first point before Z=20

# 3d plot <<<<<<<<<<
# r=20
fig = plt.figure(figsize=(8,8), dpi=200)    
ax = fig.add_subplot(111, projection='3d')
ax.scatter(df1["x"], df1["y"], df1["z"], s=0.1, c="black")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_title("3d plot for r=20")
plt.savefig("build/A2_3d_r20.pdf")
plt.clf()

# r=28
fig = plt.figure(figsize=(8,8), dpi=200)    
ax = fig.add_subplot(111, projection='3d')
ax.scatter(df1["x"], df1["y"], df1["z"], s=0.1, c="black")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_title("3d plot for r=28")
plt.savefig("build/A2_3d_r28.pdf")
plt.clf()