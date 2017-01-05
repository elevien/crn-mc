import sys
from pylab import *
from crn_mc.mesh import *
from crn_mc.model import *
from crn_mc.simulation.paths import *
from crn_mc.simulation.montecarlo import *

# there are 10 reaction
# X1 = M
# X2 = D
# X3 = RNA
# X4 = DNA unbound
# X5 = DNA bound at 1 site
# X6 = DNA bound at 2 stites

# stoichiometric vectors

r1 = array([0,0,1,0,0,0])
p1 = array([1,0,1,0,0,0])

r2 = array([1,0,0,0,0,0])
p2 = array([0,0,0,0,0,0])

r3 = array([0,0,0,0,1,0])
p3 = array([0,0,1,0,1,0])

r4 = array([0,0,1,0,0,0])
p4 = array([0,0,0,0,0,0])

r5 = array([0,1,0,1,0,0])
p5 = array([0,0,0,0,1,0])

r6 = array([0,0,0,0,1,0])
p6 = array([0,1,0,1,0,0])

r7 = array([0,1,0,0,1,0])
p7 = array([0,0,0,0,0,1])

r8 = array([0,0,0,0,0,1])
p8 = array([0,1,0,0,1,0])

r9 = array([2,0,0,0,0,0])
p9 = array([1,0,0,0,0,0])

r10 = array([1,0,0,0,0,0])
p10 = array([2,0,0,0,0,0])


# scaling exponents


g1 = 1.
g2 = 1.
g3 = 0.
g4 = 0.
g5 = 0.
g6 = 0.



z1 = 4.3
z2 = 0.07
z3 = 7.15
z4 = 0.39
z5 = 1.99
z6 = 47.9
z7 = 0.0199
z8 = 8.7
z9 = 0.083
z10 = 0.5
#
# b1 = -1.
# b2 = -2.
# b3 = -1.
# b4 = -1.
# b5 = -1.
# b6 = 0.
# b7 = -3.
# b8 = -2.
# b9 = -1
# b10 = 0.

a1 = -1.
a2 = -1.
a3 = -1.
a4 = -1.
a5 = 0.
a6 = 0.
a7 = -2.
a8 = -2.
a9 = 1.
a10 = 1.

# setup model
Nx = 1
Np = 100.
L = 1.
T = 1
Nspecies = 6 #(U,V)
mesh = make_lattice1d(Nx,L)
model = Model(Nspecies,mesh)
model_coupled = ModelHybridSplitCoupled(Nspecies,mesh)

# add reactions
model.addReaction(r1,p2,z1*pow(Np,a1))
model.addReaction(r2,p2,z2*pow(Np,a2))
model.addReaction(r3,p3,z3*pow(Np,a3))
model.addReaction(r4,p4,z4*pow(Np,a4))
model.addReaction(r5,p5,z5*pow(Np,a5))
model.addReaction(r6,p6,z6*pow(Np,a6))
model.addReaction(r7,p7,z7*pow(Np,a7))
model.addReaction(r8,p8,z8*pow(Np,a8))
model.addReaction(r9,p9,z9*pow(Np,a9))
model.addReaction(r10,p10,z10*pow(Np,a10))

model_coupled.addReactionSlow(r1,p2,z1*pow(Np,a1))
model_coupled.addReactionSlow(r2,p2,z2*pow(Np,a2))
model_coupled.addReactionSlow(r3,p3,z3*pow(Np,a3))
model_coupled.addReactionSlow(r4,p4,z4*pow(Np,a4))
model_coupled.addReactionSlow(r5,p5,z5*pow(Np,a5))
model_coupled.addReactionSlow(r6,p6,z6*pow(Np,a6))
model_coupled.addReactionSlow(r7,p7,z7*pow(Np,a7))
model_coupled.addReactionSlow(r8,p8,z8*pow(Np,a8))
model_coupled.addReactionFast(r9,p9,z9*pow(Np,a9))
model_coupled.addReactionFast(r10,p10,z10*pow(Np,a10))

# set intial conditions
model.systemState[0,0] = 10*pow(Np,g1)
model.systemState[1,0] = 10*pow(Np,g2)
model.systemState[2,0] = 10.
model.systemState[3,0] = 10.
model.systemState[4,0] = 10.
model.systemState[5,0] = 10.

# set intial conditions
model_coupled.systemState[0,0] = 10*pow(Np,g1)
model_coupled.systemState[1,0] = 10*pow(Np,g2)
model_coupled.systemState[2,0] = 10.
model_coupled.systemState[3,0] = 10.
model_coupled.systemState[4,0] = 10.
model_coupled.systemState[5,0] = 10.
model_coupled.systemState[0+6,0] = 10*pow(Np,g1)
model_coupled.systemState[1+6,0] = 10*pow(Np,g2)
model_coupled.systemState[2+6,0] = 10.
model_coupled.systemState[3+6,0] = 10.
model_coupled.systemState[4+6,0] = 10.
model_coupled.systemState[5+6,0] = 10.



# simulate
path_coupled,clock_coupled =chv(model_coupled,T,pow(Np,-2.),'lsoda',10.)
path,clock = gillespie(model,T)

plt.plot(clock_coupled,path_coupled[:,0]*pow(Np,-g1),'k+')
plt.plot(clock_coupled,path_coupled[:,1]*pow(Np,-g2),'k+')

plt.plot(clock,path[:,0]*pow(Np,-g1),'r-')
plt.plot(clock,path[:,1]*pow(Np,-g2),'k-')
# plt.plot(clock,path[:,0+6]*pow(Np,-g1),'r-')
# plt.plot(clock,path[:,1+6]*pow(Np,-g2),'k-')

ax = plt.gca()
plt.show()
