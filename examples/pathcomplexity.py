import sys
sys.path.insert(0, './../')
from mesh import *
from model import *
from simulation import *
from pylab import *
from timer import *

Nx = 1
L = 1.
T = 10.
Np =100.
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
m1.add_reaction_fast(r,p,0.6)
m2.add_reaction_fast(r,p,0.6)
m3.add_reaction(r,p,0.6)

r = array([2,1,0])
p = array([2,0,1])
m1.add_reaction_slow(r,p,3./pow(Np,2))
m2.add_reaction_slow(r,p,3./pow(Np,2))
m3.add_reaction(r,p,3./pow(Np,2))

r = array([0,0,1])
p = array([0,1,0])
m1.add_reaction_slow(r,p,1.)
m2.add_reaction_slow(r,p,1.)
m3.add_reaction(r,p,1.)

#with timer() as t:
#    path2,clock2 = strang_split(m2,T,0.001,0.001,'lsoda')
#tcpu = t.secs
#print(tcpu)

m2.system_state[0,0] = Np
m2.system_state[1,0] = 1.


with timer() as t:
    path3,clock3 = chv(m1,T,pow(Np,-1.),'lsoda',20.)
tcpu = t.secs
print(tcpu)

#with timer() as t:
#    path4,clock4 = gillespie(m3,T)
#tcpu = t.secs
#print(tcpu)


#plt.plot(clock2,path2[:,0],'k-')
plt.plot(clock3,path3[:,0],'r-')
plt.plot(clock3,path3[:,3],'k-')
#plt.plot(clock4,path4[:,0],'b--')
plt.show()
