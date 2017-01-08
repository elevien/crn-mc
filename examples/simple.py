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
m.addspecies("A",1.,array([1.]))
m.addspecies("B",1.,array([0.]))
m.addspecies("C",0.,array([1.]))
m.addspecies("D",0.,array([0.]))


r = array([["A",1],["C",1]])
p = array([["B",1],["C",1]])
m.addreaction(r,p,2.,1.,FAST)

r = array([["B",1],["D",1]])
p = array([["A",1],["D",1]])
m.addreaction(r,p,2.,1.,FAST)

r = array([["D",1]])
p = array([["C",1]])
m.addreaction(r,p,1.,0.,SLOW)

r = array([["C",1]])
p = array([["D",1]])
m.addreaction(r,p,1.,0.,SLOW)

delta = 2
#sol2,Er2 = mc_crude(m,T,systemSize,delta,0,lambda: makepath(m,T,pow(systemSize,-delta),'lsoda',5.,0))
#print(Er2)
#path,clock = gillespie(m,T,0)
path,clock = makepath(m,T,pow(systemSize,-2.),'lsoda',1.,0)
path2,clock2 = makepath_coupled(m,T,pow(systemSize,-2.),'lsoda',1.,0)

plt.plot(clock,path[:,0],'k-')
plt.plot(clock2,path2[:,0],'r-')
plt.plot(clock2,path2[:,4],'r--')


#plt.plot(clock,path[:,2],'r-+')
#plt.plot(clock,path[:,4+2],'r--+')
#print(path)
plt.show()
