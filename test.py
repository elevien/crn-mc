import matplotlib.pyplot as plt
from coupled_bimolecular import *


N = 50
J = 10
level = 0
w = 1.
T = 1.
mesh_fine = int(N/(pow(J,level)))     # fine grid
mesh_coarse  = int(N/(pow(J,level+1))) # coarse grid
X0,Y0,X1,Y1,clock = path_coupled(N,J,level,w,T)
#for i in events:
#    print(i.wait_absolute)

plt.plot(range(mesh_fine),X0[0],'k--')
plt.plot(range(mesh_fine),Y0[0],'r--')

plt.plot(range(mesh_fine),X0[-1],'k-')
plt.plot(range(mesh_fine),Y0[-1],'r-')
#plt.plot(range(mesh_fine),Y0[-1],'g-')

X1_onfine = zeros(mesh_fine)
for i in range(mesh_coarse):
    X1_onfine[J*i:J*(i+1)] = X1[-1,i]/J
plt.step(range(mesh_fine),X1_onfine,'g-')
#plt.plot(clock,X0[:,0],'k-')
#plt.plot(clock,X0[:,1],'r+')
#plt.plot(clock,X1[:,0]/J,'g--')
#x0 =0
#for i in range(J):
#    x0 = x0 + X0[:,1*J+i]
#plt.plot(clock,x0,'k-')
#ax = plt.gca()

#ax.imshow(X1, extent=[0,100,0,1], aspect=100,interpolation="none")
plt.show()
