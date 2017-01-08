import sys
from pylab import *
from crn_mc.mesh import *
from crn_mc.model import *
from crn_mc.simulation.paths import *
from crn_mc.simulation.montecarlo import *


Nx = 1
L = 1
T = 5.
mesh = make_lattice1d(Nx,L)
systemSize = 10.
m = Model(mesh,systemSize)
X1 = m.addspecies("M",1.,[10.])
X2 = m.addspecies("D",1.,[10.])
X3 = m.addspecies("RNA",0.,[10.])
X4 = m.addspecies("DNA",0.,[10.])
X5 = m.addspecies("DNA.D",1.,[10.])
X6 = m.addspecies("DNA.2D",1.,[10.])
OPEN = m.addspecies("0",0.,[10.])

m.addreaction([["RNA",1]],[["RNA",1],["M",1]],2.,0.,SLOW)
m.addreaction([["M",1]],[["0",1]],2.,0.,SLOW)
m.addreaction([["DNA.D",1]],[["RNA",1],["DNA.D",1]],2.,0.,SLOW)
m.addreaction([["RNA",1]],[["0",1]],2.,0.,SLOW)
m.addreaction([["DNA",1],["D",1]],[["DNA.D",1]],2.,0.,SLOW)
m.addreaction([["DNA.D",1]],[["DNA",1],["D",1]],2.,0.,SLOW)
m.addreaction([["DNA.D",1],["D",1]],[["DNA.2D",1]],2.,0.,SLOW)
m.addreaction([["M",2]],[["D",1]],2.,0.,FAST)
m.addreaction([["D",1]],[["M",2]],2.,0.,FAST)


# path,clock = makepath_coupled(m,T,pow(systemSize,-2.),'lsoda',1.,0)
# plt.plot(clock,path[:,0],'r-+')
# plt.plot(clock,path[:,1],'r--+')
# plt.plot(clock,path[:,0+m.dimension],'k-')
# plt.plot(clock,path[:,1+m.dimension],'k--')
delta = 1.5
Q,err=montecarlo_coupled(m,T,delta,0,0,'lsoda',10.)
Q2,err2=montecarlo_crude(m,T,delta,0,0,'lsoda',10.)
print(Q)
print(Q2)
plt.plot(err,'k-')
plt.plot(err2,'r-')
plt.show()

plt.show()
