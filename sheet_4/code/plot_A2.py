import pandas as pd
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
    return anim

if __name__=="__main__":
    # time of animation = len(df)/(n*fps)

    # a)
    print("\nPlot a)")
    df_a = pd.read_csv("build/A2_a).csv", header=None)
    animation = animate(df_a, n=40)
    plt.show()
    #animation.save("build/A2_a).mp4", dpi=300, fps=100)

    # b)
    # stable
    print("\nPlot b) small dt")
    df_b_slow = pd.read_csv("build/A2_b)slow.csv", header=None)
    animation = animate(df_b_slow, n=1)
    plt.show()
    #animation.save("build/A2_b)slow.mp4", dpi=300, fps=30)

    #unstable
    print("\nPlot b) large dt")
    df_b_fast = pd.read_csv("build/A2_b)fast.csv", header=None)
    animation = animate(df_b_fast, n=1)
    plt.show()
    #animation.save("build/A2_b)fast.mp4", dpi=300, fps=30)