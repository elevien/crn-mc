import numpy as np

class Mesh:
    def __init__(self,dimension,topology,geometry):
        self.size = len(topology)
        self.dimension = dimension
        self.topology = topology  # adjaceny matrix
        self.geometry = geometry  # list of (x,y,z)

def make_lattice1d(Nx):
    topology = np.zeros((Nx,Nx))
    d = np.ones(Nx-1)
    topology = np.diag(d,1)+np.diag(d,-1)
    mesh  = Mesh(1,topology,[])
    return mesh

def make_lattice1d_coupled(Nx,J):
    topology = np.zeros((Nx,Nx))
    coupling = np.zeros((Nx,Nx))
    d = np.ones(Nx-1)
    topology = np.diag(d,1)+np.diag(d,-1)
    mesh  = Mesh(1,topology,[])

    for i in range(int(Nx/J)):
        coupling[i*J:(i+1)*J,i*J:(i+1)*J] = np.ones((J,J))
    return mesh,coupling

# need to implement
def make_lattice2d(Nx,Ny):
    return None

def make_lattice3d(Nx,Ny):
    return None
