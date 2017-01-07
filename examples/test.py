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

systemSize = 100.
m = Model(mesh,systemSize)
m.addSpecies("A",1.,array([1.]))
m.addSpecies("B",1.,array([0.]))
m.addSpecies("C",0.,array([1.]))
m.addSpecies("D",0.,array([0.]))


r = array([["A",1],["C",1]])
p = array([["B",1],["C",1]])
m.addReaction(r,p,2.,1.,FAST)

r = array([["B",1],["D",1]])
p = array([["A",1],["D",1]])
m.addReaction(r,p,2.,1.,FAST)

r = array([["D",1]])
p = array([["C",1]])
m.addReaction(r,p,1.,0.,SLOW)

r = array([["C",1]])
p = array([["D",1]])
m.addReaction(r,p,1.,0.,SLOW)

m2 = ModelHybridSplitCoupled(m)

    #print(e)

print("model-----------------------------")
for e in m.events:
    e.updateRate()
    print(e)

print("coupling-----------------------------")
for e in m2.events:
    e.updateRate()
    print(e)
##for s in m2.systemState:
    #print(s.name)

#path,clock = gillespie(m,T,0)
path,clock = chv1d(m,T,pow(systemSize,-2.),'lsoda',10.,0)
print('done')

plt.plot(clock,path[:,0],'k-')

#ax = plt.gca()
plt.show()
