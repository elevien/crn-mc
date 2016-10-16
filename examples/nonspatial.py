import sys
sys.path.insert(0, './../')
from mesh import *
from model import *
from simulation import *
from pylab import *

Nx = 1
Np = 500
L = 1.
T = 10
Nspecies = 3 #(U,V)
mesh = make_lattice1d(Nx,L)
model = Model(Nspecies,mesh)
model_hybrid = ModelHybrid(Nspecies,mesh)
model_hybrid.system_state = np.zeros((Nspecies,1))
model_hybrid.system_state[0,0] = Np
model_hybrid.system_state[1,0] = 2.
model.system_state = np.zeros((Nspecies,1))
model.system_state[0,0] = Np
model.system_state[1,0] = 2.



r = array([1,1,0])
p = array([2,1,0])
model_hybrid.add_reaction_fast(r,p,1.)
model.add_reaction(r,p,1.)

r = array([2,0,0])
p = array([1,0,0])
model_hybrid.add_reaction_fast(r,p,1./Np)
model.add_reaction(r,p,1./Np)

r = array([0,1,0])
p = array([0,0,1])
model_hybrid.add_reaction_slow(r,p,1.)
model.add_reaction(r,p,7.)

r = array([0,0,1])
p = array([0,1,0])
model_hybrid.add_reaction_slow(r,p,4.)
model.add_reaction(r,p,4.)



path_hybrid,clock_hybrid = gillespie_hybrid(model_hybrid,T,0.0001)
path,clock = gillespie(model,T)
plt.plot(clock_hybrid,path_hybrid[:,0],'k-')
#plt.plot(clock_hybrid,Np*path_hybrid[:,1],'r-')
plt.plot(clock,path[:,0],'r-')


ax = plt.gca()
plt.show()
