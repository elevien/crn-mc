import sys
sys.path.insert(0, './../')
from mesh import *
from model import *
from simulation import *
from pylab import *
from timer import *


Nx = 1
Np = 5
L = 1.
T = 100.
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
m3.add_reaction(r,p,1.)

r = array([0,0,1])
p = array([0,1,0])
m1.add_reaction_slow(r,p,0.8)
m2.add_reaction_slow(r,p,0.8)
m3.add_reaction(r,p,0.8)



h = 0.1
delta = 1.5
q1 = mc_hyrbidCoupled(m1,m2,1,Np,delta,h)
q2 = mc_crude(m3,1,Np,delta)
print(q1)
print(q2)
# with timer() as t;
#     path1,clock1= gillespie_hybrid(m1,T,1./pow(Np,2),0.1,'lsoda')
# print(t.secs)
#
# with timer() as t:
#     path2,clock2= gillespie_hybrid(m2,T,1./pow(Np,2),0.1,'lsoda')
# print(t.secs)


#plt.plot(clock1,path1[:,0],'g-')
#plt.plot(clock2,path2[:,0],'g--')
##plt.plot(clock1,path1[:,Nspecies],'k-',alpha=0.4)


#ax = plt.gca()
#plt.show()
