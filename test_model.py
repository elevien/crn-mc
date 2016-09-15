from mesh import *
from model import *
from simulation import *
from pylab import *


Nx = 50
Nspecies = 2
mesh_fine = make_lattice1d(Nx)
mesh_coarse = make_lattice1d(Nx/2)

model = Model(Nspecies,mesh_fine)
model.system_state = 50*ones((Nspecies,Nx))
model.add_diffusions(0)
model.add_diffusions(1)

r = array([1,1])
p = array([0,0])
model.add_reaction(r,p)



path,clock = next_reaction(model,10)
#print(clock)
#print(path)

plt.plot(range(Nx),path[-1,0],'k-')
plt.plot(range(Nx),path[-1,1],'k--')


plt.show()
