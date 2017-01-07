import sys
from pylab import *
from crn_mc.mesh import *
from crn_mc.model import *
from crn_mc.simulation.paths import *
from crn_mc.simulation.montecarlo import *
from Goutsias import *

Np = 100
model,model_hybrid,model_coupled = Goutsias(Np)

T = 10000.
path_coupled,clock_coupled =chv(model_coupled,T,pow(Np,-2.),'lsoda',0.)
#path,clock = gillespie(model,T)

species = 3
plt.plot(clock_coupled,path_coupled[:,0],'k--')
plt.plot(clock_coupled,path_coupled[:,0+6],'k-')



ax = plt.gca()
plt.show()
