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

def model_setup(Np):
    # set initial conditions
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
    m3.add_reaction(r,p,1.*Np)

    r = array([1,0,0])
    p = array([0,0,0])
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
    m3.add_reaction(r,p,0.1)



h = 0.01
delta = 1.5
Np_range = array([100,200,300])
T_coupled = zeros(len(Np_range))
T_crude = zeros(len(Np_range))
Q_coupled = zeros(len(Np_range))
Q_crude = zeros(len(Np_range))

for i in range(len(Np_range)):
    print(i)
    Np = Np_range[i]
    model_setup(Np)
    with timer() as t:
        Q_coupled[i] = mc_hyrbidCoupled(m1,m2,1,Np,delta,h)
    T_coupled[i] = t.secs
    with timer() as t:
        Q_crude[i] = mc_crude(m3,1,Np,delta)
    T_crude[i] = t.secs

fname1 = '/Users/E/Dropbox/RESEARCH/coupled_mc/paperdata/Q_coupled.csv'
fname2 = '/Users/E/Dropbox/RESEARCH/coupled_mc/paperdata/Q_crude.csv'
fname3 = '/Users/E/Dropbox/RESEARCH/coupled_mc/paperdata/T_coupled.csv'
fname4 = '/Users/E/Dropbox/RESEARCH/coupled_mc/paperdata/T_crude.csv'
fname5 = '/Users/E/Dropbox/RESEARCH/coupled_mc/paperdata/Np_range.csv'
savetxt(fname1,Q_coupled,delimiter=',')
savetxt(fname2,Q_crude,delimiter=',')
savetxt(fname3,T_coupled,delimiter=',')
savetxt(fname4,T_crude,delimiter=',')
savetxt(fname5,Np_range,delimiter=',')
