import sys
sys.path.insert(0, './../')
from mesh import *
from model import *
from simulation import *
from pylab import *

Nx = 1
Np = 50
L = 1.
Nspecies = 3 #(U,V)
mesh = make_lattice1d(Nx,L)
model = ModelHybrid(Nspecies,mesh)
model.system_state = np.zeros((Nspecies,1))
model.system_state[0,0] = Np
model.system_state[1,0] = 1.



r = array([0,1,0])
p = array([1,1,0])
model.add_reaction_fast(r,p,500.)

r = array([1,0,0])
p = array([0,0,0])
model.add_reaction_fast(r,p,20.)

r = array([1,1,0])
p = array([0,0,1])
model.add_reaction_slow(r,p,2.4)

r = array([0,0,1])
p = array([0,1,0])
model.add_reaction_slow(r,p,2.3)


path,clock = gillespie_hyrbid(model,10.,0.01)

plt.plot(clock,path[:,0],'k-')
plt.plot(clock,path[:,1],'r-')

ax = plt.gca()
plt.show()
