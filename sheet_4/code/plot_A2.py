import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

#Parameters
n = 40 # hold every n-th frame

# Read data
df = pd.read_csv("build/A1.csv", header=None)

# hold every n-th frame
df = df.iloc[::n, :]

# remove last column (line ends with ,)
t = df.iloc[:, 0].to_numpy()
df.drop(df.columns[[0]], axis=1, inplace=True)

# print total frames
print("Total frames: ", len(df))

# create FuncAnimation to display the one dimensional diffusion equation
fig, ax = plt.subplots()

# set axis limits
ax.set_xlim(0, len(df.columns)-1)
ax.set_ylim(df.min().min(), df.max().max()+0.01)

ax.set_xlabel("x")
ax.set_ylabel("u(x,t)")

# plot initial condition
line, = ax.plot(df.iloc[0, :], color="black")

# update function
def update(frame):
    line.set_ydata(df.iloc[frame, :])
    
    if frame % 10 == 0:
        ax.set_title("Time = " + str(t[frame].round(2)))
    return line, ax

# create animation
anim = FuncAnimation(fig, update, frames=len(df), interval=1)
plt.show()
# save animation
anim.save("build/A2.mp4", dpi=300, fps=100)