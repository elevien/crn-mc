import sys
sys.path.insert(0, './../')
from mesh import *
from model import *
from simulation import *
from pylab import *

Nx = 1
Np = 50
L = 1.
Nspecies = 2 #(U,V)
mesh = make_lattice1d(Nx,L)
model = Model(Nspecies,mesh)
model.system_state = Np*np.ones((Nspecies,1))


r = array([1,1])
p = array([0,0])
model.add_reaction(r,p,1.)

r = array([0,1])
p = array([1,0])
model.add_reaction(r,p,2.)

r = array([0,0])
p = array([1,0])
model.add_reaction(r,p,2.4)

r = array([0,0])
p = array([0,1])
model.add_reaction(r,p,2.3)


path,clock = gillespie(model,10.)

plt.plot(clock,path[:,0],'k-')
plt.plot(clock,path[:,1],'r-')

ax = plt.gca()
plt.show()
