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
systemSize = 20.
m = Model(mesh,systemSize)
# from HYE-WON KANG AND THOMAS G. KURTZ 2013
X1 = m.addspecies("M",exponent=1.)
X2 = m.addspecies("D",exponent=1.)
X3 = m.addspecies("RNA",exponent=0.)
X4 = m.addspecies("DNA",exponent=0.)
X5 = m.addspecies("DNA.D",exponent=0.)
X6 = m.addspecies("DNA.2D",exponent=0.)
OPEN = m.addspecies("0",exponent=0.)
m.addreaction([["RNA",1]],[["RNA",1],["M",1]],4.30,exponent=-1.)
m.addreaction([["M",1]],[["0",1]],7.0,exponent=-1.)
m.addreaction([["DNA.D",1]],[["RNA",1],["DNA.D",1]],7.15,exponent=-1.)
m.addreaction([["RNA",1]],[["0",1]],0.39,exponent=-1.)
m.addreaction([["DNA",1],["D",1]],[["DNA.D",1]],1.99,exponent=0.)
m.addreaction([["DNA.D",1]],[["DNA",1],["D",1]],0.479,exponent=0.)
m.addreaction([["DNA.D",1],["D",1]],[["DNA.2D",1]],199.,exponent=-2.)
m.addreaction([["DNA.2D",1]],[["DNA.D",1],["D",1]],8.77e10-8,exponent=-2.)
m.addreaction([["M",2]],[["D",1]],8.30,exponent=1.)
m.addreaction([["D",1]],[["M",2]],0.55,exponent=1.)

# set initial data
ic = [1.,1.,0.,2.,0.,0.,0.]
for i in range(m.dimension):
    m.systemState[i].value[0]= ic[i]
for e in m.events:
    print(e)
delta = 2.
Q2,standdev2 = montecarlo(m,T,delta,method='lsoda',sample_rate = 4.,
                                    estimator = 'coupled')
Q1,standdev1 = montecarlo(m,T,delta,method='lsoda',sample_rate = 4.,
                                    estimator = 'crude',path_type='exact')

plt.plot(standdev1)
plt.plot(standdev2)
plt.show()
print(Q1)
print(Q2)
print(standdev1)
print(standdev2)
