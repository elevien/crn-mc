import sys
from pylab import *
from crn_mc.mesh import *
from crn_mc.model import *
from crn_mc.simulation.paths import *
from crn_mc.simulation.montecarlo import *


Nx = 1
L = 1
T = 100.
mesh = make_lattice1d(Nx,L)
systemSize = 100.
m = Model(mesh,systemSize)

#
X1 = m.addspecies("A",exponent=1.)
X2 = m.addspecies("B",exponent=1.)
X3 = m.addspecies("C",exponent=0.)
X4 = m.addspecies("D",exponent=0.)

m.addreaction([["A",1],["C",1]],[["B",1],["C",1]],1.,exponent=1.)
m.addreaction([["B",1],["D",1]],[["A",1],["D",1]],1.,exponent=0.)
m.addreaction([["C",1]],[["D",1]],0.5,exponent=0.)
m.addreaction([["D",1]],[["C",1]],0.3,exponent=0.)
# set initial data
ic = [1.,0.,1.,0.]
for i in range(m.dimension):
    m.systemState[i].value[0]=ic[i]
for e in m.events:
    print(e)

delta = 1.8
Q,standdev = montecarlo(m,T,delta,method='lsoda',sample_rate = 0.,estimator = 'crude')
plt.plot(standdev)
print(Q)

#path_exact,clock_exact = makepath(m,T,pow(systemSize,-2.),sample_rate = 1.,path_type='coupled')
#plt.plot(clock_exact,path_exact[:,0],'k-')
#plt.plot(clock_exact,path_exact[:,4],'b+')

plt.show()
