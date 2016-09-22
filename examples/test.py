import sys
sys.path.insert(0, './../')
from mesh import *
from model import *
from simulation import *
from pylab import *


def test2d():
    Nx = 3
    Ny = 2
    Lx = 1.
    Ly = 1.
    Np = 10
    Nspecies = 3
    D0 = 10.
    D1 = 3.

    mesh = make_lattice2d(Nx,Nx,Lx,Ly)
    print(mesh.topology)
    model= Model(Nspecies,mesh)
    model.add_diffusions(0,D0)
    model.add_diffusions(1,D1)
    model.system_state = Np*ones((Nspecies,Nx*Ny))
    path,clock = next_reaction(model,10)
    print(path)

    return None

def test1d_uncoupled():
    Nx = 100
    Np = 50
    L = 1.
    D0 = pow(10,-3.)/pow((L/Nx),2)
    D1 = pow(10,-1.)/pow((L/Nx),2)
    Nspecies = 2 #(U,V)
    mesh = make_lattice1d(Nx,L)
    model = Model(Nspecies,mesh)
    ic = Np*ones((Nspecies,Nx))
    model.system_state = ic

    model.add_diffusions(0,D0)
    model.add_diffusions(1,D1)

    r = array([0,0])
    p = array([1,0])
    model.add_reaction(r,p,4*pow(10,3.))

    r = array([1,0])
    p = array([0,0])
    model.add_reaction(r,p,2.)

    r = array([0,0])
    p = array([0,1])
    model.add_reaction(r,p,1.2*pow(10,3.))

    r = array([2,1])
    p = array([3,0])
    model.add_reaction(r,p,6.25*pow(10,-8.))

    path,clock = gillespie(model,10)
    plt.plot(range(Nx),path[-1,0],'k-')
    plt.plot(range(Nx),path[-1,1],'r--')
    print(clock[-1])
    #plt.plot(range(Nx),path[-1,2],'k-')
    #plt.plot(range(Nx),path[-1,Nspecies+2],'k+')
    ax = plt.gca()
    #ax.set_ylim([0,2*Np])
    #name1 = './../../paperdata/catalyst/Q_coupled_d2.csv'
    #fname2 = './../../paperdata/catalyst/Q_crude_d2.csv'

    #savetxt(fname1,Q_coupled,delimiter=',')
    #savetxt(fname2,Q_crude,delimiter=',')
    plt.show()


def test_var_coupled():
    Nx = 20
    Np = 80
    L = 1.
    D0 = pow(10,-3.)/pow((L/Nx),2)
    D1 = pow(10,-1.)/pow((L/Nx),2)
    J = 2

    Nspecies = 2
    mesh,coupling = make_lattice1d_coupled(Nx,L,J)
    model = SplitCoupled(Nspecies,mesh,coupling)
    model_uncoupled = Model(Nspecies,mesh)


    ic = make_coupledSS(Np*ones((Nspecies,Nx)),coupling)
    model.system_state = ic
    model_uncoupled.system_state = Np*ones((Nspecies,Nx))
    model.add_diffusions(0,D0)
    model.add_diffusions(1,D1)
    model_uncoupled.add_diffusions(0,D0)
    model_uncoupled.add_diffusions(1,D1)

    r = array([0,0])
    p = array([1,0])
    model.add_reaction(r,p,4*pow(10,3.))
    model_uncoupled.add_reaction(r,p,4*pow(10,3.))

    r = array([1,0])
    p = array([0,0])
    model.add_reaction(r,p,2.)
    model_uncoupled.add_reaction(r,p,2.)

    r = array([0,0])
    p = array([0,1])
    model.add_reaction(r,p,1.2*pow(10,3.))
    model_uncoupled.add_reaction(r,p,1.2*pow(10,3.))

    r = array([2,1])
    p = array([3,0])
    model.add_reaction(r,p,6.25*pow(10,-8.))
    model_uncoupled.add_reaction(r,p,6.25*pow(10,-8.))

    Ntests = 10
    err = np.zeros(Ntests)
    err_uncoupled = np.zeros(Ntests)
    for k in range(Ntests):
        print("test  # "+str(k))
        path,clock = gillespie(model,10)
        path_uncoupled,clock_uncoupled = gillespie(model_uncoupled,10)
        path_mod = np.zeros(Nx)
        for i in range(int(Nx/J)):
            path_mod[i*J] = path[-1,Nspecies,J*i]/J
            path_mod[i*J+1] = path[-1,Nspecies,J*i]/J
        err_uncoupled[k] = abs(np.max(path[-1,0] - path_uncoupled[-1,0]))
        err[k] = abs(np.max(path[-1,0] - path_mod))
    plt.plot(range(Ntests),err_uncoupled,'k-',label='uncoupled')
    plt.plot(range(Ntests),err,'r-',label='coupled')

    print("uncoupled var = "+str(var(err_uncoupled)))
    print("coupled var = "+str(var(err)))

    plt.legend(bbox_to_anchor=(0.9, 0.9), borderaxespad=0.)

    plt.show()




def test1d_coupled():
    Nx = 90
    Np = 75
    L = 1.
    D0 = pow(10,-3.)/pow((L/Nx),2)
    D1 = pow(10,-1.)/pow((L/Nx),2)
    J = 2

    Nspecies = 3
    mesh,coupling = make_lattice1d_coupled(Nx,L,J)
    model = SplitCoupled(Nspecies,mesh,coupling)
    model_uncoupled = Model(Nspecies,mesh)


    ic = make_coupledSS(Np*ones((Nspecies,Nx)),coupling)
    model.system_state = ic
    model_uncoupled.system_state = Np*ones((Nspecies,Nx))
    model.add_diffusions(0,D0)
    model.add_diffusions(1,D1)
    model_uncoupled.add_diffusions(0,D0)
    model_uncoupled.add_diffusions(1,D1)

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

    path,clock = gillespie(model,10)
    path_uncoupled,clock_uncoupled = gillespie(model_uncoupled,10)
    print(clock[-1])

    path_mod = np.zeros(Nx)
    for i in range(int(Nx/J)):
        path_mod[i*J] = path[-1,Nspecies,J*i]/J
        path_mod[i*J+1] = path[-1,Nspecies,J*i]/J
    plt.plot(range(Nx),path[-1,0],'k-',label='coupled - fine grid')
    plt.plot(range(Nx),path_mod,'r-',label='coupled - coarse grid')
    plt.plot(range(Nx),path_uncoupled[-1,0],'k--',label='uncoupled')

    #plt.plot(clock,path[:,0,J],'k-',label='coupled - fine grid')
    #plt.plot(clock,path[:,Nspecies,J]/J,'r-',label='coupled - coarse grid')
    #plt.plot(clock,path_uncoupled[:,0,J],'k--',label='uncoupled')
    plt.xlabel('Voxel', fontsize=20)
    plt.ylabel('U', fontsize=20)
    plt.legend(bbox_to_anchor=(0.9, 0.9), borderaxespad=0.)
    ax = plt.gca()
    #ax.set_ylim([1.,4*Np])
    savefig('./../../output/DynamicTrap.pdf', bbox_inches='tight')
    plt.show()
    return None

test1d_coupled()
