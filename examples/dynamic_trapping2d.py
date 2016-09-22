import sys
sys.path.insert(0, './../')
from mesh import *
from model import *
from simulation import *
from pylab import *

Nx = 100
Np = 50
L = 1.
D0 = pow(10,-3.)/pow((L/Nx),2)
D1 = 0.
Nspecies = 3
mesh = make_lattice1d(Nx,L)
model = Model(Nspecies,mesh)
ic = Np*ones((Nspecies,Nx))
model.system_state = ic

model.add_diffusions(0,D0)
model.add_diffusions(1,D1)

# absorption by trap
r = array([1,1,0])
p = array([0,0,1])
model.add_reaction(r,p,4.)

# trap opening and closing
r = array([0,1,0])
p = array([0,0,1])
model.add_reaction(r,p,2.)

r = array([0,0,1])
p = array([0,1,0])
model.add_reaction(r,p,5.)

path,clock = gillespie(model,2)

plt.imshow(path[:,1], aspect='auto')
#plt.plot(range(Nx),path[-1,0],'k-')
#plt.plot(range(Nx),path[-1,1],'r--')
#plt.plot(range(Nx),path[-1,2],'k-')
#plt.plot(range(Nx),path[-1,Nspecies+2],'k+')
ax = plt.gca()
#savetxt('./../../output/schnakenberg1d-path.csv',path,delimiter=',')
#savetxt('./../../output/schnakenberg1d-clock.csv',clock,delimiter=',')
#savetxt('./../../output/schnakenberg1d-mesh.csv',mesh,delimiter=',')
plt.show()
