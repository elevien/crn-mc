import matplotlib.pyplot as plt
from coupled_bimolecular import *


N = 100
J = 4
level = 0
w = 1.
T = 1.
X,Y,t_grid = path_uncoupled(N,w,T)
#for i in events:
#    print(i.wait_absolute)

plt.plot(t_grid,X[:,10])
plt.show()
