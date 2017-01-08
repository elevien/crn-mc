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
m.addReaction(r,p,1.,0.,FAST)

r = array([["C",1]])
p = array([["D",1]])
m.addReaction(r,p,1.,0.,FAST)

#m2 = ModelHybridSplitCoupled(m)


#path,clock = gillespie(m,T,0)
#path,clock = pathHybrid(m,T,pow(systemSize,-2.),'lsoda',10.,0)
path,clock = pathHybridSplitCoupled(m,T,pow(systemSize,-2.),'lsoda',10.,0)
#print(path)

plt.plot(clock,path[:,0],'k-')
plt.plot(clock,path[:,1],'k--')
#print(path)
plt.show()
