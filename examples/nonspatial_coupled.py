import sys
sys.path.insert(0, './../')
from mesh import *
from model import *
from simulation import *
from pylab import *


Nx = 1
Np = 50
L = 1.
T = 7
Nspecies = 4 #(U,V)
mesh = make_lattice1d(Nx,L)
model = ModelHybridSplitCoupled(Nspecies,mesh)
model.system_state[0,0] = Np
model.system_state[1,0] = Np
model.system_state[2,0] = 2.
model.system_state[Nspecies,0] = Np
model.system_state[Nspecies+1,0] = Np
model.system_state[Nspecies+2,0] = 2.



r = array([0,0,1,0])
p = array([1,0,1,0])
model.add_reaction_fast(r,p,Np)

r = array([2,0,0,0])
p = array([1,0,0,0])
model.add_reaction_fast(r,p,1./Np)

r = array([0,2,0,0])
p = array([1,1,0,0])
model.add_reaction_fast(r,p,1./Np)

r = array([1,0,0,1])
p = array([1,1,0,1])
model.add_reaction_fast(r,p,1.)


r = array([1,0,1,0])
p = array([1,0,0,1])
model.add_reaction_slow(r,p,1./Np)

r = array([0,1,0,1])
p = array([0,1,1,0])
model.add_reaction_slow(r,p,4./Np)



path,clock= gillespie_hybrid(model,T,0.001)
plt.plot(clock,path[:,1],'k--')
plt.plot(clock,path[:,Nspecies+1],'k-')


ax = plt.gca()
plt.show()
