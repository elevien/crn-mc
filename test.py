from mesh import *
from model import *
from simulation import *
from pylab import *


Nx = 8
Np = 30
L = 1.
J =2
D0 = 10.
D1 = 3.

Nspecies = 3
mesh,coupling = make_lattice1d_coupled(Nx,L,J)
model = SplitCoupled(Nspecies,mesh,coupling)
model_uncoupled = Model(Nspecies,mesh)


ic = make_coupledSS(Np*ones((Nspecies,Nx)),coupling)
model.system_state = ic
model_uncoupled.system_state = Np*ones((Nspecies,Nx))
print("Diffusion coupled ---------------------------------------------------")
model.add_diffusions(0,D0)
#model.add_diffusions(1,D1)
print("Diffusion uncoupled ---------------------------------------------------")
model_uncoupled.add_diffusions(0,D0)
#model_uncoupled.add_diffusions(1,D1)

r = array([1,1,0])
p = array([0,0,1])
model.add_reaction(r,p)
model_uncoupled.add_reaction(r,p)

r = array([0,0,0])
p = array([1,0,0])
model.add_reaction(r,p)
model_uncoupled.add_reaction(r,p)

r = array([0,0,0])
p = array([0,1,0])
model.add_reaction(r,p)
model_uncoupled.add_reaction(r,p)

print("Running coupled ---------------------------------------------------")
path,clock = next_reaction(model,10)
print("Running uncoupled ---------------------------------------------------")
path_uncoupled,clock_uncoupled = next_reaction(model_uncoupled,10)



plt.plot(clock,path[:,2,J],'k-')
plt.plot(clock,path[:,Nspecies+2,J]/J,'r-')
plt.plot(clock,path_uncoupled[:,2,J],'g--')

#plt.plot(range(Nx),path[-1,2],'k-')
#plt.plot(range(Nx),path[-1,Nspecies+2],'k+')
ax = plt.gca()
#ax.set_ylim([0,4*Np])
plt.show()
