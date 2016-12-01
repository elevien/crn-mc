import sys
sys.path.insert(0, './../')
from mesh import *
from model import *
from simulation import *
from pylab import *
from timer import *

Nx = 1
L = 1.
T = 5.
Np = 100.
Nspecies = 3 #(U,V)
mesh = make_lattice1d(Nx,L)
m1 = ModelHybridSplitCoupled(Nspecies,mesh)
m2 = ModelHybrid(Nspecies,mesh)
m3 = Model(Nspecies,mesh)

m1.system_state[0,0] = Np
m1.system_state[1,0] = 1.
m1.system_state[Nspecies,0] = Np
m1.system_state[Nspecies+1,0] = 1.
m2.system_state[0,0] = Np
m2.system_state[1,0] = 1.
m3.system_state[0,0] = Np
m3.system_state[1,0] = 1.

r = array([0,1,0])
p = array([1,1,0])
m1.add_reaction_fast(r,p,1.*Np)
m2.add_reaction_fast(r,p,1.*Np)
m3.add_reaction(r,p,1.)

r = array([1,0,0])
p = array([0,0,0])
m1.add_reaction_fast(r,p,0.6)
m2.add_reaction_fast(r,p,0.6)
m3.add_reaction(r,p,1.)

r = array([0,1,0])
p = array([0,0,1])
m1.add_reaction_slow(r,p,1.)
m2.add_reaction_slow(r,p,1.)
m3.add_reaction(r,p,1.)

r = array([0,0,1])
p = array([0,1,0])
m1.add_reaction_slow(r,p,1.)
m2.add_reaction_slow(r,p,1.)
m3.add_reaction(r,p,1.)


with timer() as t:
    path2,clock2 = strang_split(m2,T,0.1,0.0001,'lsoda')
    path3,clock3 = gillespie_hybrid(m2,T,0.1,0.0001,'lsoda')
tcpu = t.secs
print(tcpu)
print(path2)
plt.plot(clock2,path2[:,0],'k-')
plt.plot(clock3,path3[:,0],'r-')
plt.show()
