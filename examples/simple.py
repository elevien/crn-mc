import sys
from pylab import *
from crn_mc.mesh import *
from crn_mc.model import *
from crn_mc.simulation.paths import *
from crn_mc.simulation.montecarlo import *


Nx = 1
L = 1.
T = 10.
mesh = make_lattice1d(Nx,L)

systemSize = 50.
m = Model(mesh,systemSize)
A = m.addspecies("A",1.,[1.])
B = m.addspecies("B",1.,[0.])
C = m.addspecies("C",0.,[1.])
D = m.addspecies("D",0.,[0.])


r = [["C",1]]
p = [["B",1],["C",1]]
m.addreaction(r,p,2.,1.,FAST)

r = [["B",1],["D",1]]
p = [["A",1],["D",1]]
m.addreaction(r,p,2.,1.,FAST)

r = [["D",1]]
p = [["C",1]]
m.addreaction(r,p,1.,0.,SLOW)

r = [["C",1]]
p = [["D",1]]
m.addreaction(r,p,1.,0.,SLOW)

#path,clock = makepath(m,T,pow(systemSize,-2.),'lsoda',1.,0)

delta = 1.5
Q,err=montecarlo_coupled(m,T,delta,0,0,'lsoda',10.)
Q2,err2=montecarlo_crude(m,T,delta,0,0,'lsoda',10.)
print(Q)
print(Q2)
plt.plot(err,'k-')
plt.plot(err2,'r-')
plt.show()
#path2,clock2 = makepath_coupled(m,T,pow(systemSize,-2.),'lsoda',1.,0)

#plt.plot(clock,path[:,0],'k-')
#plt.plot(clock2,path2[:,0],'r-')
#plt.plot(clock2,path2[:,4],'r--')

#print(path)
#plt.show()
