import sys
from pylab import *
from crn_mc.mesh import *
from crn_mc.model import *
from crn_mc.simulation.paths import *
from crn_mc.simulation.montecarlo import *


Nx = 1
L = 1
T = 0.05
mesh = make_lattice1d(Nx,L)
systemSize = 100.
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
ic = [1.,1.,0.,50.,0.,0.,0.]
for i in range(m.dimension):
    m.systemState[i].value[0]= ic[i]

for e in m.events:
    print(e)


path_exact,clock_exact = makepath(m,T,pow(systemSize,-2.),sample_rate = 20.,treatment='exact')
for i in range(m.dimension):
    m.systemState[i].value[0]= ic[i]
#path,clock = makepath_coupled(m,T,pow(systemSize,-2.),sample_rate = 10.)
#plt.plot(clock,path[:,2],'r-')
#plt.plot(clock,path[:,1+m.dimension],'k-',alpha=0.4)

plt.plot(clock_exact,path_exact[:,3],'k-+')
plt.plot(clock_exact,path_exact[:,4],'k-+')
#plt.plot(clock_exact,path_exact[:,2],'g--')

# delta = 1.
# Q,err=montecarlo_coupled(m,T,delta,0,0,'lsoda',10.)
# for i in range(m.dimension):
#     m.systemState[i].value[0]= ic[i]
# Q2,err2=montecarlo_crude(m,T,delta,0,0,'lsoda',10.)
# print(Q)
# print(Q2)
# plt.plot(err,'k-')
# plt.plot(err2,'r-')
plt.show()
