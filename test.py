import matplotlib.pyplot as plt
from coupled_bimolecular import *


N = 100
J = 4
level = 0
w = 1.
T = 1.
X0,X1,Y0,Y1,t_grid = path_coupled(N,J,level,w,T)
print("#################################################################")
print("t_grid = "+str(t_grid))

N0 = int(N/(pow(J,level)))     # fine grid
N1 = int(N/(pow(J,level+1))) # coarse grid

plt.plot(range(len(t_grid)),t_grid,'k-')
#plt.plot(t_grid,X0[:,10],'k-')
#plt.plot(t_grid,X1[:,10],'r-')
plt.show()
