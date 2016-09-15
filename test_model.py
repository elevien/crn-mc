from mesh import *
from model import *
from simulation import *
from pylab import *


Nx = 5
mesh = make_lattice1d(Nx)
model = Model(3,mesh)
model.system_state = 30*ones((3,Nx))
#model.init_diffusions(1)
#model.init_diffusions(2)


r = array([1,1,0])
p = array([0,0,3])
model.add_reaction(r,p)


path,clock = next_reaction(model,10)
#print(clock)
#print(path)

plt.plot(clock,path[:,0,3],'k-')
plt.plot(clock,path[:,2,3],'k--')
plt.show()
