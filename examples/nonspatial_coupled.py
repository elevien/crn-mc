import sys
sys.path.insert(0, './../')
from mesh import *
from model import *
from simulation import *
from pylab import *


Nx = 1
Np = 50
L = 1.
T = 10
Nspecies = 3 #(U,V)
mesh = make_lattice1d(Nx,L)
model = ModelHybridSplitCoupled(Nspecies,mesh)
model.system_state[0,0] = Np
model.system_state[1,0] = 2.
model.system_state[3,0] = Np
model.system_state[4,0] = 2.



r = array([1,1,0])
p = array([2,1,0])
model.add_reaction_fast(r,p,1.)

r = array([2,0,0])
p = array([1,0,0])
model.add_reaction_fast(r,p,1./Np)

r = array([0,1,0])
p = array([0,0,1])
model.add_reaction_slow(r,p,1.)

r = array([0,0,1])
p = array([0,1,0])
model.add_reaction_slow(r,p,4.)



path,clock= gillespie_hybrid(model,T,0.01)
plt.plot(clock,path[:,0],'r+')
plt.plot(clock,path[:,3],'k-')


ax = plt.gca()
plt.show()
