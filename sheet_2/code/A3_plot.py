import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def create_plot(h, df_euler_r1, df_euler_r2, df_verlet_r1, df_verlet_r2):
    # plot results of euler algorithm

    plt.figure(figsize=(6,6), dpi=200)
    plt.plot(df_euler_r1.iloc[:,0].to_numpy(), df_euler_r1.iloc[:,1].to_numpy(), 'b.', label="mass 1")
    plt.plot(df_euler_r2.iloc[:,0].to_numpy(), df_euler_r2.iloc[:,1].to_numpy(), 'r.', label="mass 2")

    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.savefig(f"build/A2_euler_h{h}.pdf")
    plt.clf()

    # plot results of verlet algorithm

    plt.figure(figsize=(6,6), dpi=200)
    plt.plot(df_verlet_r1.iloc[:,0].to_numpy(), df_verlet_r1.iloc[:,1].to_numpy(), 'b.', label="mass 1")
    plt.plot(df_verlet_r2.iloc[:,0].to_numpy(), df_verlet_r2.iloc[:,1].to_numpy(), 'r.', label="mass 2")

    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.savefig(f"build/A2_verlet_h{h}.pdf")
    plt.clf()

# a)
h_list = ["1_000000", "0_100000", "0_010000"]
for h in h_list:
    df_euler_r1 = pd.read_csv(f"build/euler_r1_h{h}.csv")
    df_euler_r2 = pd.read_csv(f"build/euler_r2_h{h}.csv")
    df_verlet_r1 = pd.read_csv(f"build/verlet_r1_h{h}.csv")
    df_verlet_r2 = pd.read_csv(f"build/verlet_r2_h{h}.csv")

    create_plot(h, df_euler_r1, df_euler_r2, df_verlet_r1, df_verlet_r2)


# c)
#df_euler_r1_pos = pd.read_csv(f"build/euler_r1_h0_050000.csv")
#df_euler_r1_neg = pd.read_csv(f"build/euler_r1_h-0_050000.csv")
#df_euler_r1 = df_euler_r1_pos.append(df_euler_r1_neg, ignore_index=True)
#
#df_verlet_r1_pos = pd.read_csv(f"build/verlet_r1_h0_050000.csv")
#df_verlet_r1_neg = pd.read_csv(f"build/verlet_r1_h-0_050000.csv")
#df_verlet_r1 = df_verlet_r1_pos.append(df_verlet_r1_neg, ignore_index=True)

#create_plot("_nr_c", df_euler_r1, df_euler_r2, df_verlet_r1, df_verlet_r2)

#plt.xlabel("x")
#plt.ylabel("y")
#plt.legend()
#plt.savefig(f"build/A2_euler_h{h}.pdf")
#plt.clf()