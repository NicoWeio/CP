import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# 1. Schreiben Sie eine Routine, die zu gegebenem r und x0 die jeweilige Abbildung N mal ausführt.
# Dieses “Warmlaufen” soll garantieren, dass man zu einem Fixpunkt/Orbit konvergiert ist, sofern
# dieser existiert. Überlegen Sie sich, welche Werte r > 0 maximal annehmen darf, damit xn in
# jedem Schritt innerhalb der jeweiligen vorgegebenen Intervalle bleibt.

# 2. Schreiben Sie nun eine Routine, die nach dem Warmlaufen weiteriteriert. Zählen sie dabei die
# Orbits auf welche sich die Punkte für die verschiedenen r verteilen.

# 3. Variieren Sie r in Schritten ∆r = 1 · 10−3 und erzeugen Sie so ein Bifurkationsdiagramm für die
# jeweilige Abbildung. Testen Sie dabei verschiedene Startwerte x0 aus, da bestimmte Fixpunk-
# te/Orbits gegebenenfalls nur von bestimmten Startwerten aus erreicht werden können. Speichern
# sie diejenigen Werte von r, bei denen sich die Anzahl der Orbits verdoppelt. Bestimmen sie dar-
# aus die Feigenbaumkonstante. Bestimmen Sie zudem ein r∞ , ab welchem sich keine “echten”
# Orbits mehr einstellen und das Bifurkationsdiagramm chaotisch wird. Da die Orbitgröße ab
# einem gewissen Punkt sehr rapide mit r wächst, ist es eine gute Näherung bereits ab einer
# Orbitgröße von G > 64 Chaos anzunehmen um r∞ zu bestimmen.

def map_log(N: int, x0: float, r: float):
    for _ in range(N):
        x0 = r*x0*(1-x0)

    return x0

if __name__ == "__main__":
    # set values
    spacing = 1E-2
    min_r = 0.0
    max_r = 5.0

    # create array with r values
    r_space = np.arange(min_r, max_r, spacing)
    print(f'Shape of r: {r_space.shape}')

    # create array with initial values x0
    x0_space = np.arange(0.0, 1.0, 0.1)

    # create empty array for x
    x_space = np.zeros((len(r_space), len(x0_space)))
    print(f'Shape of full x array: {x_space.shape}')

    # calc logistic map for each x0 and r
    # iterate over x0
    for ind,x0 in enumerate(x0_space):
        x_space[:,ind] = np.array([map_log(1000, x0, r) for r in r_space]) # iterate over r
    
    #df = pd.DataFrame(x_space, columns=[f'x0={x0}' for x0 in x0_space])
    #print(df)

    # remove data points: x < 0 or x > 1
    x_space[x_space < 1E-8] = np.nan
    x_space[x_space > 1.0] = np.nan

    # TODO: calc Feigenbaum constant
    

    # plot
    plt.figure(figsize=(10,8), dpi=200)

    #plt.plot(r_space, x_space, 'k.', markersize=0.1)
    for ind, x0 in enumerate(x0_space):
        plt.plot(r_space, x_space[:, ind], '.', markersize=0.2, label=f'x0={x0:.1f}')

    plt.legend(fontsize=6)
    plt.xlabel('r')
    plt.ylabel('x')
    plt.savefig("build/bifurkation_log.pdf")
    plt.show()
