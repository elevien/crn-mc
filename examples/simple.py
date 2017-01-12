import sys
from pylab import *
from crn_mc.mesh import *
from crn_mc.model import *
from crn_mc.simulation.paths import *
from crn_mc.simulation.montecarlo import *


Nx = 1
L = 1
T = 20.
mesh = make_lattice1d(Nx,L)
systemSize = 50.
m = Model(mesh,systemSize)

#
m.addspecies("A",exponent=1.)
m.addspecies("B",exponent=1.)
m.addspecies("C",exponent=0.)
m.addspecies("D",exponent=0.)
m.addreaction([["A",1],["C",1]],[["B",1],["C",1]],1.,exponent=1.)
m.addreaction([["B",1],["D",1]],[["A",1],["D",1]],1.,exponent=1.)
m.addreaction([["C",1]],[["D",1]],0.5,exponent=0.)
m.addreaction([["D",1]],[["C",1]],0.3,exponent=0.)
# set initial data
ic = [1.,0.,1.,0.]
for i in range(m.dimension):
    m.systemState[i].value[0]=ic[i]
for e in m.events:
    print(e)
#
delta = 2.
# WHY DO THESE NEED TO BE IN THIS ORDER?
Q2,standdev2 = montecarlo(m,T,delta,method='lsoda',sample_rate = 4.,
                                    estimator = 'coupled')
Q1,standdev1 = montecarlo(m,T,delta,method='lsoda',sample_rate = 4.,
                                    estimator = 'crude',path_type='exact')

plt.plot(standdev1)
plt.plot(standdev2)
print(Q1)
print(Q2)
print(standdev1)
print(standdev2)

# path,clock= makepath(m,T,pow(systemSize,-2.),sample_rate = 1.,path_type='coupled')
# plt.plot(clock,path[:,0],'k-')
# plt.plot(clock,path[:,4],'b+')
# print(path.dtype)

plt.show()
