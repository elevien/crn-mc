import sys
sys.path.insert(0, './../')
from mesh import *
from model import *
from simulation import *
from pylab import *

# setup model
Nx = 50
Np = 50
L = 1.
D0 = 1./pow((L/Nx),2)
D1 = 0.
J = 1.
T = 0.00001
Nspecies = 3
mesh = make_lattice1d(Nx,L)
m_u = Model(Nspecies,mesh)


icx = zeros(Nx)
icx[Nx/2] = Np
icy1 = zeros(Nx)
icy2 = ones(Nx)

m_u.system_state[0] = icx
m_u.system_state[1] = icy1
m_u.system_state[2] = icy2


m_u.add_diffusions(0,D0)

# absorption by trap
r = array([1,1,0])
p = array([0,0,1])
m_u.add_reaction(r,p,1.)

# trap opening and closing
r = array([0,1,0])
p = array([0,0,1])
m_u.add_reaction(r,p,2.)

r = array([0,0,1])
p = array([0,1,0])
m_u.add_reaction(r,p,5.)

path = mc_crude(m_u,m_u.system_state,T,1000,4)

plt.plot(range(Nx),path[0])
#plt.plot(range(Nx),path[-1,0],'k-')
#plt.plot(range(Nx),path[-1,1],'r--')
#plt.plot(range(Nx),path[-1,2],'k-')
#plt.plot(range(Nx),path[-1,Nspecies+2],'k+')
ax = plt.gca()
#savetxt('./../../output/schnakenberg1d-path.csv',path,delimiter=',')
#savetxt('./../../output/schnakenberg1d-clock.csv',clock,delimiter=',')
#savetxt('./../../output/schnakenberg1d-mesh.csv',mesh,delimiter=',')
plt.show()
