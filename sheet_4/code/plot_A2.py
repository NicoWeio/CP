import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def animate(df_data, n):
    # hold every n-th frame
    df = df_data.iloc[::n, :].copy()

    # remove last column (line ends with ,)
    t = df.iloc[:, 0].to_numpy()
    df.drop(df.columns[[0]], axis=1, inplace=True)

    # print total frames
    print("Total frames: ", len(df))
    print("Number bins: ", len(df.columns))

    # create FuncAnimation to display the one dimensional diffusion equation
    fig, ax = plt.subplots()

    # set ticks constant for all animations
    #ax.set_xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90])

    # set axis limits
    ax.set_xlim(0, len(df.columns)-1)
    ax.set_ylim(df.min().min()-0.01, df.max().max()+0.01)

    ax.set_xlabel(r"$x$", fontsize=14)
    ax.set_ylabel(r"$u (x,t)$", fontsize=14)

    # plot initial condition
    line, = ax.plot(np.arange(len(df.columns)), df.iloc[0, :], color="black", linewidth=2)

    # update function
    def update(frame):
        line.set_ydata(df.iloc[frame, :])
        
        if frame % 10 == 0:
            ax.set_title("t = " + str(t[frame].round(2)))
        return line, ax

    # create animation
    anim = FuncAnimation(fig, update, frames=len(df), interval=1)
    return anim

if __name__=="__main__":
    # time of animation = len(df)/(n*fps)
    
    # a)
    print("\nPlot a)")
    df_a = pd.read_csv("build/A2_a).csv", header=None)
    animation = animate(df_a, n=50)
    plt.show()  # Get the current axes
    #animation.save("build/A2_a).mp4", dpi=300, fps=100)

    # b)
    # stable
    print("\nPlot b) small dt")
    df_b_slow = pd.read_csv("build/A2_b)slow.csv", header=None)
    animation = animate(df_b_slow, n=3)
    plt.show()  # Get the current axes
    #animation.save("build/A2_b)slow.mp4", dpi=300, fps=100)

    #unstable
    print("\nPlot b) large dt")
    df_b_fast = pd.read_csv("build/A2_b)fast.csv", header=None)
    animation = animate(df_b_fast, n=1)
    plt.show()  # Get the current axes
    #animation.save("build/A2_b)fast.mp4", dpi=300, fps=10)

    # c)
    # u1
    print("\nPlot c) u1")
    df_c_u1 = pd.read_csv("build/A2_c)u1.csv")
    animation = animate(df_c_u1, n=10)
    #animation.save("build/A2_c)u1.mp4", dpi=300, fps=100)
    plt.show()

    # u2
    print("\nPlot c) u2")
    df_c_u2 = pd.read_csv("build/A2_c)u2.csv")
    animation = animate(df_c_u2, n=200)
    ax = plt.show()  # Get the current axes
    #animation.save("build/A2_c)u2.mp4", dpi=300, fps=100)

