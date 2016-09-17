from mesh import *
from model import *
from simulation import *
from pylab import *


Nx = 10
J =2
Nspecies = 2
mesh,coupling = make_lattice1d_coupled(Nx,J)
print(coupling)
model = SplitCoupled(Nspecies,mesh,coupling)

model.system_state = make_coupledSS(1000*ones((Nspecies,Nx)),coupling)
model.add_diffusions(0)

r = array([1,1])
p = array([0,0])
model.add_reaction(r,p)



path,clock = next_reaction(model,10)


for species in path[-1]:
    print(species)
