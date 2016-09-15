from mesh import *
from model import *
from simulation import *
from pylab import *


Nx = 100
mesh = make_lattice1d(Nx)
model = Model(5,mesh)
model.system_state = 30*ones((5,Nx))
model.init_diffusions(1)
model.init_diffusions(2)


r = array([1,1,0,0,0])
p = array([0,0,3,0,0])
model.add_reaction(r,p)

r = array([0,0,1,0,0])
p = array([0,0,0,1,1])
model.add_reaction(r,p)

r = array([0,0,0,1,1])
p = array([2,0,0,0,0])
model.add_reaction(r,p)

r = array([0,0,0,0,1])
p = array([0,0,0,0,0])
model.add_reaction(r,p)


path,clock = next_reaction(model,10)
#print(clock)
#print(path)

plt.plot(range(Nx),path[-1,0],'k-')
plt.plot(range(Nx),path[-1,1],'k--')
plt.plot(range(Nx),path[-1,2],'g-')
plt.plot(range(Nx),path[-1,3],'r-')
plt.plot(range(Nx),path[-1,4],'b-')

plt.show()
