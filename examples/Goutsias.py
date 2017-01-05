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
T = 10
Nspecies = 6 #(U,V)
mesh = make_lattice1d(Nx,L)
model = Model(Nspecies,mesh)
model_coupled = ModelHybrid(Nspecies,mesh)

# add reactions
model.add_reaction(r1,p2,z1*pow(Np,a1))
model.add_reaction(r2,p2,z2*pow(Np,a2))
model.add_reaction(r3,p3,z3*pow(Np,a3))
model.add_reaction(r4,p4,z4*pow(Np,a4))
model.add_reaction(r5,p5,z5*pow(Np,a5))
model.add_reaction(r6,p6,z6*pow(Np,a6))
model.add_reaction(r7,p7,z7*pow(Np,a7))
model.add_reaction(r8,p8,z8*pow(Np,a8))
model.add_reaction(r9,p9,z9*pow(Np,a9))
model.add_reaction(r10,p10,z10*pow(Np,a10))

model_coupled.add_reaction_slow(r1,p2,z1*pow(Np,a1))
model_coupled.add_reaction_slow(r2,p2,z2*pow(Np,a2))
model_coupled.add_reaction_slow(r3,p3,z3*pow(Np,a3))
model_coupled.add_reaction_slow(r4,p4,z4*pow(Np,a4))
model_coupled.add_reaction_slow(r5,p5,z5*pow(Np,a5))
model_coupled.add_reaction_slow(r6,p6,z6*pow(Np,a6))
model_coupled.add_reaction_slow(r7,p7,z7*pow(Np,a7))
model_coupled.add_reaction_slow(r8,p8,z8*pow(Np,a8))
model_coupled.add_reaction_fast(r9,p9,z9*pow(Np,a9))
model_coupled.add_reaction_fast(r10,p10,z10*pow(Np,a10))

# set intial conditions
model.system_state = np.zeros((Nspecies,1))
model.system_state[0,0] = 10*pow(Np,g1)
model.system_state[1,0] = 10*pow(Np,g2)
model.system_state[2,0] = 10.
model.system_state[3,0] = 10.
model.system_state[4,0] = 10.
model.system_state[5,0] = 10.

# set intial conditions
model_coupled.system_state = np.zeros((Nspecies,1))
model_coupled.system_state[0,0] = 10*pow(Np,g1)
model_coupled.system_state[1,0] = 10*pow(Np,g2)
model_coupled.system_state[2,0] = 10.
model_coupled.system_state[3,0] = 10.
model_coupled.system_state[4,0] = 10.
model_coupled.system_state[5,0] = 10.
# model_coupled.system_state[0+6,0] = 10*pow(Np,g1)
# model_coupled.system_state[1+6,0] = 10*pow(Np,g2)
# model_coupled.system_state[2+6,0] = 10.
# model_coupled.system_state[3+6,0] = 10.
# model_coupled.system_state[4+6,0] = 10.
# model_coupled.system_state[5+6,0] = 10.



# simulate
path_coupled,clock_coupled =chv(model_coupled,T,pow(Np,-2.),'lsoda',10.)
path,clock = gillespie(model,T)

plt.plot(clock_coupled,path_coupled[:,4]*pow(Np,-g1),'r--')
plt.plot(clock_coupled,path_coupled[:,5]*pow(Np,-g2),'k--')

plt.plot(clock,path[:,4]*pow(Np,-g1),'r-')
plt.plot(clock,path[:,5]*pow(Np,-g2),'k-')
# plt.plot(clock,path[:,0+6]*pow(Np,-g1),'r-')
# plt.plot(clock,path[:,1+6]*pow(Np,-g2),'k-')

ax = plt.gca()
plt.show()
