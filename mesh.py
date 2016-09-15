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


# need to implement
def make_lattice2d(Nx,Ny):
    return None

def make_lattice3d(Nx,Ny):
    return None
