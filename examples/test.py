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

systemSize = 10.
m = Model(mesh,systemSize)
m.addSpecies("A",1.,mesh,array([1.]))
m.addSpecies("B",0.,mesh,array([1.]))
m.addSpecies("C",0.,mesh,array([1.]))
m.addSpecies("D",0.,mesh,array([1.]))

r = array([["B",1.]])
p = array([["A",1.],["B",1.]])
m.addReaction(r,p,2.,1.,'FAST')

r = array([["A",1.]])
p = array([["D",1.]])
m.addReaction(r,p,2.,1.,'SLOW')

r = array([["B",1.]])
p = array([["C",1.]])
m.addReaction(r,p,2.,0.,'SLOW')

r = array([["C",1.]])
p = array([["B",1.]])
m.addReaction(r,p,2.,0.,'SLOW')

#path,clock = gillespie(m,T,0)
path,clock = chv1d(m,T,pow(systemSize,-2.),'lsoda',50.,0)

plt.plot(clock,path[:,0],'k-')

ax = plt.gca()
plt.show()
