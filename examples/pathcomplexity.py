import sys
sys.path.insert(0, './../')
from mesh import *
from model import *
from simulation import *
from pylab import *
from timer import *

Nx = 1
L = 1.
T = 1.
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

r = array([0,0,0])
p = array([1,0,0])
m1.add_reaction_fast(r,p,1.)
m2.add_reaction_fast(r,p,1.)
m3.add_reaction(r,p,1.)

r = array([0,1,0])
p = array([0,0,1])
m1.add_reaction_slow(r,p,1.)
m2.add_reaction_slow(r,p,1.)
m3.add_reaction(r,p,0.1)

r = array([0,0,1])
p = array([0,1,0])
m1.add_reaction_slow(r,p,0.8)
m2.add_reaction_slow(r,p,0.8)
m3.add_reaction(r,p,0.8)
