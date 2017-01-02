import sys
sys.path.insert(0, './../')
from mesh import *
from model import *
from simulation import *
from pylab import *
from timer import *


Nx = 1
L = 1.
T = 5.
Nspecies = 3 #(U,V)
mesh = make_lattice1d(Nx,L)


def model_setup(Np):

    m1 = ModelHybridSplitCoupled(Nspecies,mesh)
    m2 = ModelHybrid(Nspecies,mesh)
    m3 = Model(Nspecies,mesh)

    # set initial conditions
    m1.system_state[0,0] = Np
    m1.system_state[1,0] = 1.
    m1.system_state[2,0] = 0.
    m1.system_state[Nspecies,0] = Np
    m1.system_state[Nspecies+1,0] = 1.
    m1.system_state[Nspecies+2,0] = 0.
    m2.system_state[0,0] = Np
    m2.system_state[1,0] = 1.
    m2.system_state[2,0] = 0.
    m3.system_state[0,0] = Np
    m3.system_state[1,0] = 1.
    m3.system_state[2,0] = 0.

    r = array([0,1,0])
    p = array([1,1,0])
    m1.add_reaction_fast(r,p,1.*Np)
    m2.add_reaction_fast(r,p,1.*Np)
    m3.add_reaction(r,p,1.)

    r = array([1,0,0])
    p = array([0,0,0])
    m1.add_reaction_fast(r,p,1.)
    m2.add_reaction_fast(r,p,1.)
    m3.add_reaction(r,p,1.)

    r = array([1,1,0])
    p = array([1,0,1])
    m1.add_reaction_slow(r,p,4./Np)
    m2.add_reaction_slow(r,p,4./Np)
    m3.add_reaction(r,p,4./Np)

    r = array([0,0,1])
    p = array([0,1,0])
    m1.add_reaction_slow(r,p,5.)
    m2.add_reaction_slow(r,p,5.)
    m3.add_reaction(r,p,5.)

    return m1,m2,m3



delta = 1.1
Np_range = array([100,150,200,250])
T_coupled = zeros(len(Np_range))
T_crude = zeros(len(Np_range))
Q_coupled = zeros(len(Np_range))
Q_crude = zeros(len(Np_range))
V_coupled = zeros(len(Np_range))

fname1 = '/Users/E/Dropbox/RESEARCH/coupled_mc/code/output/Q_coupled_gene_var.csv'
fname2 = '/Users/E/Dropbox/RESEARCH/coupled_mc/code/output/Q_crude_gene_var.csv'
fname3 = '/Users/E/Dropbox/RESEARCH/coupled_mc/code/output/T_coupled-gene_fb.csv'
fname4 = '/Users/E/Dropbox/RESEARCH/coupled_mc/code/output/T_crude_gene_var.csv'

fname5 = '/Users/E/Dropbox/RESEARCH/coupled_mc/code/output/Np_range-T.csv'
fname6 = '/Users/E/Dropbox/RESEARCH/coupled_mc/code/output/V_coupled-genefb_s500.csv'


#fname5 = 'Np_range_gene_var.csv'
#fname6 = 'V_coupled_gene_d1p1.csv'

for i in range(len(Np_range)):
    print(i)
    Np = Np_range[i]
    m1,m2,m3= model_setup(Np)
    #q = zeros(500)
    #for j in range(500):
    #    m1,m2,m3= model_setup(Np)
    #    path,clock = chv(m1,T,pow(Np,-delta),'lsoda',0.)
    #    q[j] = abs(path[-1,0]-path[-1,3])


    #V_coupled[i] = var(q/Np)
    with timer() as t:
        Q_coupled[i] = mc_hyrbidCoupled(m1,m2,T,Np,delta,0.)
    T_coupled[i] = t.secs

    #    Q_crude[i] = mc_crude(m3,T,Np,delta)
    #T_crude[i] = t.secs
    #q[j] = Q_coupled[i]
    #V_coupled[i] = var(q/Np)

savetxt(fname5,Np_range,delimiter=',')
savetxt(fname3,T_coupled,delimiter=',')
#plt.plot(Np_range,V_coupled,'r-')
#plt.plot(Np_range,T_crude,'k-')
#plt.show()

#savetxt(fname1,Q_coupled,delimiter=',')
#savetxt(fname2,Q_crude,delimiter=',')
#savetxt(fname3,T_coupled,delimiter=',')
#savetxt(fname4,T_crude,delimiter=',')
