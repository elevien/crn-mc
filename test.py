import matplotlib.pyplot as plt
from coupled_bimolecular import *


N = 50
J = 4
level = 0
w = 1.
T = 1.
X,Y,t_grid = path_uncoupled(N,w,T)
#for i in events:
#    print(i.wait_absolute)

plt.plot(t_grid,X[:,5])
plt.plot(t_grid,X[:,1])
plt.plot(t_grid,X[:,2])
print(t_grid)
ax = plt.gca()

#ax.imshow(X, extent=[0,100,0,1], aspect=100,interpolation="none")
plt.show()
