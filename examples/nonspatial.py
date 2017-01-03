import sys
from pylab import *
from crn_mc.mesh import *
from crn_mc.model import *
from crn_mc.simulation.paths import *


Nx = 1
Np = 100
L = 1.
T = 100
Nspecies = 3 #(U,V)
mesh = make_lattice1d(Nx,L)
print(mesh)
model = Model(Nspecies,mesh)
model_hybrid = ModelHybrid(Nspecies,mesh)
model_hybrid.system_state = np.zeros((Nspecies,1))
model_hybrid.system_state[0,0] = Np
model_hybrid.system_state[1,0] = 20.
model.system_state = np.zeros((Nspecies,1))
model.system_state[0,0] = Np
model.system_state[1,0] = 20.



r = array([0,1,0])
p = array([1,1,0])
model_hybrid.add_reaction_fast(r,p,1.*Np)
model.add_reaction(r,p,1.*Np)

r = array([1,0,0])
p = array([0,0,0])
model_hybrid.add_reaction_fast(r,p,1.)
model.add_reaction(r,p,1.)

r = array([0,1,0])
p = array([0,0,1])
# works if this is fast reaction channel
model_hybrid.add_reaction_slow(r,p,1.)
model.add_reaction(r,p,1.)

r = array([0,0,1])
p = array([0,1,0])
model_hybrid.add_reaction_slow(r,p,1.)
model.add_reaction(r,p,1.)



path_hybrid,clock_hybrid = chv(model_hybrid,T,pow(Np,-2.),'lsoda',100.)
path,clock = gillespie(model,T)
plt.plot(clock_hybrid,path_hybrid[:,1],'k-')
#plt.plot(clock_hybrid,Np*path_hybrid[:,1],'r-')
plt.plot(clock,path[:,1],'r-')


ax = plt.gca()
plt.show()
