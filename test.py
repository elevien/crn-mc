from mesh import *
from model import *
from simulation import *
from pylab import *


Nx = 4
J =2
Nspecies = 2
mesh,coupling = make_lattice1d_coupled(Nx,J)
model = SplitCoupled(Nspecies,mesh,coupling)
model_uncoupled = Model(Nspecies,mesh)


ic = make_coupledSS(5*ones((Nspecies,Nx)),coupling)
model.system_state = ic
model_uncoupled.system_state = 5*ones((Nspecies,Nx))
model.add_diffusions(0)
model.add_diffusions(1)
model_uncoupled.add_diffusions(0)
model_uncoupled.add_diffusions(1)

r = array([-1,-1])
p = array([0,0])
model.add_reaction(r,p)
model_uncoupled.add_reaction(r,p)




path,clock = next_reaction(model,10)
path_uncoupled,clock_uncoupled = next_reaction(model_uncoupled,10)


plt.plot(clock,path[:,0,J],'k-')
plt.plot(clock,path[:,2,J]/J,'r+')
plt.plot(clock,path_uncoupled[:,0,J],'g--')
ax = plt.gca()
ax.set_ylim([1, 20])
plt.show()
