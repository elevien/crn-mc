import sys
from pylab import *
from crn_mc.mesh import *
from crn_mc.model import *
from crn_mc.simulation.paths import *
from crn_mc.simulation.montecarlo import *


Nx = 1
Np = 100
L = 1.
T = 10
Nspecies = 3 #(U,V)
mesh = make_lattice1d(Nx,L)
print(mesh)
model = Model(Nspecies,mesh)
model_hybrid = ModelHybrid(Nspecies,mesh)
model_hybrid.systemState = np.zeros((Nspecies,1))
model_hybrid.systemState[0,0] = Np
model_hybrid.systemState[1,0] = 100.
model.systemState = np.zeros((Nspecies,1))
model.systemState[0,0] = Np
model.systemState[1,0] = 100.



r = array([0,1,0])
p = array([1,1,0])
model_hybrid.addreactionFast(r,p,1.*Np)
model.addreaction(r,p,1.*Np)

r = array([1,0,0])
p = array([0,0,0])
model_hybrid.addreactionFast(r,p,1.)
model.addreaction(r,p,1.)

r = array([0,1,0])
p = array([0,0,1])
# works if this is fast reaction channel
model_hybrid.addreactionSlow(r,p,1.)
model.addreaction(r,p,1.)

r = array([0,0,1])
p = array([0,1,0])
model_hybrid.addreactionSlow(r,p,1.)
model.addreaction(r,p,1.)



path_hybrid,clock_hybrid = chv(model_hybrid,T,pow(Np,-2.),'lsoda',1.)
path,clock = gillespie(model,T)
plt.plot(clock_hybrid,path_hybrid[:,1],'k-')
#plt.plot(clock_hybrid,Np*path_hybrid[:,1],'r-')
plt.plot(clock,path[:,1],'r-')


ax = plt.gca()
plt.show()
