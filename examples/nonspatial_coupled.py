import sys
sys.path.insert(0, './../')
from mesh import *
from model import *
from simulation import *
from pylab import *
from timer import *


Nx = 1
Np = 100
L = 1.
T = 100.
Nspecies = 4 #(U,V)
mesh = make_lattice1d(Nx,L)
m1 = ModelHybridSplitCoupled(Nspecies,mesh)
m2 = ModelHybrid(Nspecies,mesh)
m1.system_state[0,0] = Np
m1.system_state[1,0] = Np
m1.system_state[2,0] = 2.
m1.system_state[Nspecies,0] = Np
m1.system_state[Nspecies+1,0] = Np
m1.system_state[Nspecies+2,0] = 2.
m2.system_state[0,0] = Np
m2.system_state[1,0] = Np
m2.system_state[2,0] = 2.




r = array([1,0,1,0])
p = array([0,1,1,0])
m1.add_reaction_fast(r,p,1.)
m2.add_reaction_fast(r,p,1.)

r = array([0,1,0,1])
p = array([1,0,0,1])
m1.add_reaction_fast(r,p,1.)
m2.add_reaction_fast(r,p,1.)

r = array([0,0,1,0])
p = array([0,0,0,1])
m1.add_reaction_slow(r,p,1.)
m2.add_reaction_slow(r,p,1.)

r = array([0,0,0,1])
p = array([0,0,1,0])
m1.add_reaction_slow(r,p,0.8)
m2.add_reaction_slow(r,p,0.8)


#Np_range = [10,20,30,40]
#for n in Np_range:

    # s = zeros(10)
    # for k in range(10):
    #     path,clock= gillespie_hybrid(model,T,1./pow(n,2))
    #     s[k] = path[-1,1]-path[-1,Nspecies+1]
    # print(var(s/10.))
with timer() as t:
    path1,clock1= gillespie_hybrid(m1,T,1./pow(Np,2))
print(t.secs)

with timer() as t:
    path2,clock2= gillespie_hybrid(m2,T,1./pow(Np,2))
print(t.secs)


plt.plot(clock1,path1[:,1],'g-')
plt.plot(clock2,path2[:,1],'g--')
plt.plot(clock1,path1[:,Nspecies+1],'k-',alpha=0.4)


ax = plt.gca()
plt.show()
