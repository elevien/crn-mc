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
    Nx = 20
    Np = 10
    L = 1.
    D0 = pow(10,-3.)
    D1 = pow(10,-1.)
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

    path,clock = next_reaction(model,10)
    plt.plot(range(Nx),path[-1,0],'k-')
    plt.plot(range(Nx),path[-1,1],'k--')

    #plt.plot(range(Nx),path[-1,2],'k-')
    #plt.plot(range(Nx),path[-1,Nspecies+2],'k+')
    ax = plt.gca()
    #ax.set_ylim([0,2*Np])
    plt.show()


def test1d_coupled():
    Nx = 20
    Np = 10
    L = 1.
    D0 = pow(10,-3.)
    D1 = pow(10,-1.)
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

    path,clock = next_reaction(model,10)
    path_uncoupled,clock_uncoupled = next_reaction(model_uncoupled,10)

    #plt.plot(range(Nx),path[-1,0],'k-')
    #plt.plot(range(Nx),path[-1,Nspecies]/J,'r+')
    #plt.plot(range(Nx),path_uncoupled[-1,0],'k--')

    plt.plot(clock,path[:,0,J],'k-')
    plt.plot(clock,path[:,Nspecies,J],'r-')
    plt.plot(clock,path_uncoupled[:,0,J],'k--')
    ax = plt.gca()
    #ax.set_ylim([1.,2*Np])
    plt.show()
    return None

test1d_coupled()
