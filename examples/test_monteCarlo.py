import sys
from pylab import *
from crn_mc.mesh import *
from crn_mc.model import *
from crn_mc.simulation.paths import *
from crn_mc.simulation.montecarlo import *
from Goutsias import *

Np = 100
model,model_hybrid,model_coupled = Goutsias(Np)

T = 0.1
delta = 1.
species = 0
sol,Er = mc_crude(model_hybrid,T,Np,delta,species,lambda: chv(model_hybrid,T,pow(Np,-delta),'lsoda',0.))
#sol,Er = mc_crude(model,T,Np,delta,lambda:gillespie(model,T))
sol2,Er2 = mc_hyrbidCoupled(model_coupled,model_hybrid,T,Np,delta,0.,species)
#print(sol)
plt.plot(Er,'ko')
plt.plot(Er2,'r-+')
plt.show()
print(Er)
print(Er2)
plt.show()
