import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 1. Read in data (x,y,z) <<<<<<<<<<
path1 = f'build/A2_r20.csv' # r=20
path2 = f'build/A2_r28.csv' # r=28

df1 = pd.read_csv(path1)
df2 = pd.read_csv(path2)

# 2. Plot x,y projection <<<<<<<<<<
# r=20
plt.figure(figsize=(8,8), dpi=200)
plt.scatter(df1['x'], df1['y'], s=0.1, c=df1['z'])
plt.xlabel('x')
plt.ylabel('y')
plt.title('x,y projection for r=20')
cbar = plt.colorbar()
cbar.ax.set_ylabel('z-value')
plt.savefig('build/A2_xy_r20.pdf')

# r=28
plt.figure(figsize=(8,8), dpi=200)
plt.scatter(df2['x'], df2['y'], s=0.1, c=df2['z'])
plt.xlabel('x')
plt.ylabel('y')
plt.title('x,y projection for r=28')
cbar = plt.colorbar()
cbar.ax.set_ylabel('z-value')
plt.savefig('build/A2_xy_r28.pdf')

# 3. 3d plot <<<<<<<<<<
# r=20
fig = plt.figure(figsize=(8,8), dpi=200)    
ax = fig.add_subplot(111, projection='3d')
ax.scatter(df1['x'], df1['y'], df1['z'], s=0.1, c='black')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title('Trajectories for r=20')
plt.savefig('build/A2_3d_r20.pdf')

# r=28
fig = plt.figure(figsize=(8,8), dpi=200)    
ax = fig.add_subplot(111, projection='3d')
ax.scatter(df2['x'], df2['y'], df2['z'], s=0.1, c='black')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title('Trajectories for r=28')
plt.savefig('build/A2_3d_r28.pdf')

# 4. Poincaré-Schnitt at z = const = 20 <<<<<<<<<<<
z_const = 20
print(f'Poincare cut at z = {z_const}')

def poincare_cut(x, y, z, z_const):

    X_poincare = []							# Habe hier mal etwas Code ergänzt, mit dem man auf die richtige Lösung kommt
    Y_poincare = []
    Z_poincare = []

    for i in range(1, len(z)):
        if z[i] <= z_const and z[i-1] > z_const:
            t = (z_const - z[i-1]) / (z[i] - z[i-1])
            X_poincare.append(x[i-1] + t * (x[i] - x[i-1]))
            Y_poincare.append(y[i-1] + t * (y[i] - y[i-1]))
            Z_poincare.append(z_const)
    
    return X_poincare, Y_poincare

    # last value over z_const									
    #over_indices = np.where(z > z_const)[0]

    # first value under z_const
    #under_indices = np.where(z < z_const)[0]

    # linear interpolation between last value over z_const and first value under z_const
    #interp_index = np.interp(z_const, z[over_indices], over_indices)				# Verstehe nicht, wie hier interpoliert wird. Fehlt hier nicht etwas, damit zwischen 
                                                    # z_i < z_const < z_{i+1} interpoliert werden kann?   												 
    #print(interp_index)

    # x and y values of the interpolated point
    #interp_x = np.interp(z_const, z[over_indices], x[over_indices])				# siehe oben
    #interp_y = np.interp(z_const, z[over_indices], y[over_indices])

    # integer indices for the poincare cut
    #poincare_indices = np.concatenate((under_indices, [interp_index])).astype(int)

    # coordinates of the poincare cut
    #poincare_x = x[poincare_indices]
    #poincare_y = y[poincare_indices]

    #return poincare_x, poincare_y


# r=20
poincare_x1, poincare_y1 = poincare_cut(df1['x'].to_numpy(), df1['y'].to_numpy(), df1['z'].to_numpy(), z_const) # arguments: x, y, z, z_const
print('\nr=30')
print(f'Interpolated x: {poincare_x1}')
print(f'Interpolated y: {poincare_y1}')

plt.figure(figsize=(8,8), dpi=200)
plt.scatter(poincare_x1, poincare_y1, s=0.1, c='black')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Poincaré-cut at z = 20 for r=20')
plt.savefig('build/A2_poincare_r20.pdf')
# r=28
poincare_x2, poincare_y2 = poincare_cut(df2['x'].to_numpy(), df2['y'].to_numpy(), df2['z'].to_numpy(), z_const) # arguments: x, y, z, z_const
print('\nr=28')
print(f'Interpolated x: {poincare_x2}')
print(f'Interpolated y: {poincare_y2}')

plt.figure(figsize=(8,8), dpi=200)
plt.scatter(poincare_x2, poincare_y2, s=0.1, c='black')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Poincaré-cut at z = 20 for r=28')
plt.savefig('build/A2_poincare_r28.pdf')

# 5. show all figures <<<<<<<<<<
plt.show()
