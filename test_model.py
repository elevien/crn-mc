from mesh import *
from model import *
from simulation import *
from pylab import *


Nx = 5
mesh = make_lattice1d(Nx)
model = Model(3,mesh)
model.system_state = ones((3,Nx))
model.init_diffusions(1)
model.init_diffusions(2)

r = array([1,1,0])
p = array([0,0,1])
model.add_reaction(r,p,"produce C")
for i in model.events:
    print(i)
