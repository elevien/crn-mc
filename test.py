from mesh import *
from model import *
from simulation import *
from pylab import *


Nx = 4
J = 2
Nspecies = 2
mesh,coupling = make_lattice1d_coupled(Nx,J)

model = SplitCoupled(Nspecies,mesh,coupling)

model.system_state = 100*ones((2*Nspecies,Nx))
model.add_diffusions(0)

r = array([1,1])
p = array([0,0])
model.add_reaction(r,p)



path,clock = next_reaction(model,10)


for species in path[-1]:
    print(species)
