import sys
from pylab import *
from crn_mc.mesh import *
from crn_mc.model import *
from crn_mc.simulation.paths import *
from crn_mc.simulation.montecarlo import *
from timer import *


Nx = 1
L = 1
T = 10.
mesh = make_lattice1d(Nx,L)
systemSize = 100.
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
#WHY DO THESE NEED TO BE IN THIS ORDER? NEED TO RESET EVENT TYPES?
print('running coupled monte carlo ...')
with timer(verbose=False) as t:
    Q2,standdev2,event_count2 = montecarlo(m,T,delta,ode_method='lsoda',sample_rate = 4.,
                                    estimator = 'coupled',path_type='hybrid')
print("   time         = "+ str(t.secs))
print("   samples      = "+ str(len(standdev2)))
print("   event_count  = "+ str(event_count2))
for e in m.events:
    print(e)
print('running crude monte carlo ...')
with timer(verbose=False) as t:
    Q1,standdev1,event_count1 = montecarlo(m,T,delta,ode_method='lsoda',sample_rate = 4.,
                                    estimator = 'crude',path_type='exact')
print("   time         = "+ str(t.secs))
print("   samples      = "+ str(len(standdev1)))
print("   event_count  = "+ str(event_count1))

plt.plot(standdev1,'k-')
plt.plot(standdev2,'g-')
print(Q1)
print(Q2)
print(standdev1)
print(standdev2)
# with timer(verbose=False) as t:
#     pt = 'hybrid'
#     path,clock= makepath(m,T,pow(systemSize,-2.),sample_rate = 5.,path_type=pt,ode_method='vode')
# print("   time   ("+pt+") = "+ str(t.secs))
# print("   length ("+pt+") = "+ str(len(clock)))
# plt.plot(clock,path[:,0],'k-')
# plt.show()
