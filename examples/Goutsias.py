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
m.addspecies("M",1.,array([1.]))
m.addspecies("D",1.,array([0.]))
m.addspecies("RNA",0.,array([1.]))
m.addspecies("DNA",0.,array([0.]))
m.addspecies("DNA-D",0.,array([0.]))
m.addspecies("DNA-2D",0.,array([0.]))
