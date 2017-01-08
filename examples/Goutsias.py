import sys
from pylab import *
from crn_mc.mesh import *
from crn_mc.model import *
from crn_mc.simulation.paths import *
from crn_mc.simulation.montecarlo import *


Nx = 1
L = 1
T = 10.
mesh = make_lattice1d(Nx,L)
systemSize = 40.
m = Model(mesh,systemSize)
X1 = m.addspecies("M",1.,[2.])
X2 = m.addspecies("D",1.,[2.])
X3 = m.addspecies("RNA",0.,[2.])
X4 = m.addspecies("DNA",0.,[2.])
X5 = m.addspecies("DNA.D",1.,[2.])
X6 = m.addspecies("DNA.2D",1.,[2.])
OPEN = m.addspecies("0",0.,[2.])

ic = [2.,2.,2.,2.,2.,2.,2.]

m.addreaction([["RNA",1]],[["RNA",1],["M",1]],2.,0.,SLOW)
m.addreaction([["M",1]],[["0",1]],2.,0.,SLOW)
m.addreaction([["DNA.D",1]],[["RNA",1],["DNA.D",1]],2.,0.,SLOW)
m.addreaction([["RNA",1]],[["0",1]],2.,0.,SLOW)
m.addreaction([["DNA",1],["D",1]],[["DNA.D",1]],2.,0.,SLOW)
m.addreaction([["DNA.D",1]],[["DNA",1],["D",1]],2.,0.,SLOW)
m.addreaction([["DNA.D",1],["D",1]],[["DNA.2D",1]],2.,0.,SLOW)
m.addreaction([["M",2]],[["D",1]],2.,1.,FAST)
m.addreaction([["D",1]],[["M",2]],2.,1.,FAST)

# path_exact,clock_exact = makepath(m,T,pow(systemSize,-2.),'lsoda',1.,0)
# for i in range(m.dimension):
#     m.systemState[i].value[0]= ic[i]
# path,clock = makepath_coupled(m,T,pow(systemSize,-2.),'lsoda',1.,0)
# plt.plot(clock,path[:,1],'r-')
# plt.plot(clock,path[:,1+m.dimension],'k-',alpha=0.4)
# plt.plot(clock_exact,path_exact[:,1],'g--')

delta = 1.
Q,err=montecarlo_coupled(m,T,delta,0,0,'lsoda',10.)
for i in range(m.dimension):
    m.systemState[i].value[0]= ic[i]
Q2,err2=montecarlo_crude(m,T,delta,0,0,'lsoda',10.)
print(Q)
print(Q2)
plt.plot(err,'k-')
plt.plot(err2,'r-')
plt.show()

plt.show()
