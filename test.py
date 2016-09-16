from mesh import *
from model import *
from simulation import *
from pylab import *


Nx = 50
J = 2;
Nspecies = 2
mesh = make_lattice1d(Nx)
coupling = []*(Nx/J)
for i in coupling:
    i = range()


#model = Model(Nspecies,mesh_fine)
#model.system_state = 50*ones((Nspecies,Nx))
#model.add_diffusions(0)
#model.add_diffusions(1)

#r = array([1,1])
#p = array([0,0])
#model.add_reaction(r,p)



path,clock = next_reaction(model,10)
#print(clock)
#print(path)

plt.plot(range(Nx),path[-1,0],'k-')
plt.plot(range(Nx),path[-1,1],'k--')


plt.show()
