import sys
sys.path.insert(0, './../')
from mesh import *
from model import *
from simulation import *
from pylab import *

Nx = 100
Np = 50
L = 1.
D0 = pow(10,-2.)/pow((L/Nx),2)
D1 = pow(10,-5.)/pow((L/Nx),2)
Nspecies = 2 #(U,V)
mesh = make_lattice1d(Nx,L)
model = Model(Nspecies,mesh)
ic = Np*ones((Nspecies,Nx))
model.system_state = ic

model.add_diffusions(0,D0)
model.add_diffusions(1,D1)

r = array([0,0])
p = array([1,0])
model.add_reaction(r,p,2*pow(10,1.))

r = array([1,0])
p = array([0,0])
model.add_reaction(r,p,2.)

r = array([0,0])
p = array([0,1])
model.add_reaction(r,p,6.25*pow(10,1.))

r = array([2,1])
p = array([3,0])
model.add_reaction(r,p,6.25*pow(10,-8.))

path,clock = gillespie(model,3)

#plt.imshow(path[:,0], aspect='auto', interpolation="none")
plt.plot(range(Nx),path[-1,0],'k-')
#plt.plot(range(Nx),path[-1,1],'r--')
#plt.plot(range(Nx),path[-1,2],'k-')
#plt.plot(range(Nx),path[-1,Nspecies+2],'k+')
ax = plt.gca()
#savetxt('./../../output/schnakenberg1d-path.csv',path,delimiter=',')
#savetxt('./../../output/schnakenberg1d-clock.csv',clock,delimiter=',')
#savetxt('./../../output/schnakenberg1d-mesh.csv',mesh,delimiter=',')
plt.show()
