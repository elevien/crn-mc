from mesh import *
from model import *
from simulation import *
from pylab import *


Nx = 20
J =5
Nspecies = 2
mesh,coupling = make_lattice1d_coupled(Nx,J)
model = SplitCoupled(Nspecies,mesh,coupling)

ic = make_coupledSS(10*ones((Nspecies,Nx)),coupling)
model.system_state = ic
model.add_diffusions(0)

r = array([-1,-1])
p = array([0,0])
model.add_reaction(r,p)



path,clock = next_reaction(model,10)

print("initial:")
for species in ic:
    print("species ="+str(species))

print("final:")
for species in path[-1]:
    print("species ="+str(species))
