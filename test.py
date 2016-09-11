import matplotlib.pyplot as plt
from coupled_bimolecular import *


N = 100
J = 4
level = 0
w = 1.
T = 1.
mesh_fine = int(N/(pow(J,level)))     # fine grid
mesh_coarse  = int(N/(pow(J,level+1))) # coarse grid
X0,Y0,X1,Y1,clock = path_coupled(N,J,level,w,T)
#for i in events:
#    print(i.wait_absolute)

plt.plot(range(mesh_fine),X0[-1],'k-')

#x0 =0
#for i in range(J):
#    x0 = x0 + X0[:,1*J+i]
#plt.plot(clock,x0,'k-')
#ax = plt.gca()

#ax.imshow(X1, extent=[0,100,0,1], aspect=100,interpolation="none")
plt.show()
