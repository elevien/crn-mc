import sys
sys.path.insert(0, './../')
from mesh import *
from model import *
from simulation import *
from pylab import *


# make models
Nx = 10
Np = 1000
L = 1.
J = 2
T = 10.
Nspecies = 2
mesh,coupling = make_lattice1d_coupled(Nx,L,J)
initial_condition = (Np/Nx)*ones((Nspecies,Nx))

Ntests = 4
D_range = linspace(0.01,2.,Ntests)
var_u = np.zeros(Ntests)
var_c = np.zeros(Ntests)

for m in range(Ntests):

    model_c = ModelSplitCoupled(Nspecies,mesh,coupling)
    model_u = Model(Nspecies,mesh)

    model_c.system_state = make_coupledSS(initial_condition,coupling)
    model_u.system_state = initial_condition

    model_u.add_diffusions(0,D_range[m])
    model_c.add_diffusions(0,D_range[m])

    r = array([1,0])
    p = array([0,1])
    model_u.add_reaction(r,p,1.)
    model_c.add_reaction(r,p,0.1)

    r = array([0,1])
    p = array([1,0])
    model_u.add_reaction(r,p,0.3)
    model_c.add_reaction(r,p,0.3)

    Ntests_2 = 10

    print("D "+str(D_range[m])+"-----------------------")

    err_c = np.zeros(Ntests_2)
    err_u = np.zeros(Ntests_2)
    for k in range(Ntests_2):
        print("test  # "+str(k))
        path_c,clock_c = gillespie(model_c,T)
        path_u,clock_u = gillespie(model_u,T)
        path_cc_proj = np.zeros(Nx)
        for i in range(int(Nx/J)):
            path_cc_proj[i*J] = path_c[-1,Nspecies,J*i]/J
            path_cc_proj[i*J+1] = path_c[-1,Nspecies,J*i]/J
        err_c[k] = abs(np.max(path_c[-1,0] - path_cc_proj))
        err_u[k] = abs(np.max(path_c[-1,0] - path_u[-1,0]))
    var_u[m] = var(err_u)
    var_c[m] = var(err_c)


print("uncoupled var = "+str(var_u))
print("coupled var = "+str(var_c))

plt.plot(D_range,var_u,'k-',label='uncoupled')
plt.plot(D_range,var_c,'r-',label='coupled')
plt.legend(bbox_to_anchor=(0.9, 0.9), borderaxespad=0.)

plt.show()
