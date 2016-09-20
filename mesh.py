import numpy as np

class Mesh:
    def __init__(self,dimension,topology,geometry):
        self.Nvoxels = len(topology)
        self.dimension = dimension
        self.topology = topology  # adjaceny matrix (numpy array), 0 along main diagonal, 1 elsewhere
        # only really works for regular grids
        self.geometry = geometry  #  numpy array of Nvoxels pairs (volume,x,(y,(z)))

def get_coarseMesh_voxel(voxel,coupling):
    # returns the coarse mesh voxel associated with
    # voxel by the coupling

    # by convention I take the coarse mesh voxel to by the smallest
    # index coupled to voxel according to coupling
    i = 0
    while coupling[voxel,i]<1:
        i = i+1
    return i


def make_lattice1d(Nx,L):
    # generates uniform 1d lattice on [0,L]
    topology = np.zeros((Nx,Nx))
    d = np.ones(Nx-1)
    topology = np.diag(d,1)+np.diag(d,-1)
    geometry = np.zeros((Nx,2))
    h = L/Nx
    print(h)
    geometry[:,0] = h*np.ones(Nx)
    geometry[:,1] = np.linspace(0,L-h,Nx)
    print(geometry)
    mesh  = Mesh(1,topology,geometry)
    return mesh

def make_lattice1d_coupled(Nx,L,J):
    mesh = make_lattice1d(Nx,L)
    coupling = np.zeros((Nx,Nx))
    for i in range(int(Nx/J)):
        coupling[i*J:(i+1)*J,i*J:(i+1)*J] = np.ones((J,J))
    return mesh,coupling

# need to implement
def make_lattice2d(Nx,Ny):
    topology = np.zeros((Nx*Ny,Nx*Ny))
    d = np.ones(Nx-1)
    topology = np.diag(d,1)+np.diag(d,-1)
    geometry = np.zeros((Nx,2))
    h = Nx/L
    geometry[:,0] = h*np.ones(Nx)
    geometry[:,1] = linspace(0,L-h,Nx)
    mesh  = Mesh(1,topology,geometry)

    return None

def make_lattice3d(Nx,Ny):
    return None
