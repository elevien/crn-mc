import numpy as np
from mesh import *
from events import *

class Model:
    def __init__(self,Nspecies,mesh):
        self.mesh = mesh
        self.Nspecies = Nspecies
        self.system_state = np.zeros((self.Nspecies,self.mesh.Nvoxels))
        self.events = []

    def add_reaction(self,reactants,products,intensity):
        for i in range(self.mesh.Nvoxels):
            reaction = Reaction(self,i,reactants,products,intensity)
            self.events.append(reaction)
        return None

    def add_diffusions(self,species,diffusivity):
        for i in range(self.mesh.Nvoxels):
            for j in range(self.mesh.Nvoxels):
                if self.mesh.topology[i,j]>0:
                    diffusion =Diffusion(self,i,j,species,diffusivity)
                    self.events.append(diffusion)
        return None

class ModelSplitCoupled(Model):
    def __init__(self,Nspecies,mesh,coupling):
        # coupling is a reducible matrix encoding which voxel
        # are grouped to make coarse mesh. See make_lattice1d_coupled(Nx,J)
        # for example

        self.coupling = coupling
        self.mesh = mesh # the fine mesh
        self.Nspecies = Nspecies
        # the arrays for the coarse and fine mesh are the same size
        # but only certian voxels store information for
        # coarse mesh. See get_coarseMesh_voxel(voxel,coupling)

        self.system_state = np.zeros((2*self.Nspecies,self.mesh.Nvoxels))
        self.events = []

    def add_reaction(self,reactants,products,intensity):
        for i in range(self.mesh.Nvoxels):
            reaction = Reaction_SplitCommon(self,i,reactants,products,intensity)
            self.events.append(reaction)
            reaction = Reaction_SplitFine(self,i,reactants,products,intensity)
            self.events.append(reaction)
            if i == get_coarseMesh_voxel(i,self.coupling): # only do this once coarse mesh voxel
                reaction = Reaction_SplitCoarse(self,i,reactants,products,intensity)
                self.events.append(reaction)



    def add_diffusions(self,species,diffusivity):
        for i in range(self.mesh.Nvoxels):
            for j in range(self.mesh.Nvoxels):
                if self.mesh.topology[j,i]>0:
                    if self.coupling[j,i]>0: # if i and j are in the same "subvoxel"
                        # add regular diffusion on fine mesh
                        diffusion = Diffusion(self,i,j,species,diffusivity)
                        self.events.append(diffusion)
                    else:
                        # add coupled diffusion
                        # various diffusion classes handle getting the "stoichiometry" right
                        diffusion = Diffusion_SplitCommon(self,i,j,species,diffusivity)
                        self.events.append(diffusion)
                        diffusion = Diffusion_SplitCoarse(self,i,j,species,diffusivity)
                        self.events.append(diffusion)
                        diffusion = Diffusion_SplitFine(self,i,j,species,diffusivity)
                        self.events.append(diffusion)
        return None

def make_coupledSS(system_state,coupling):
    # convert a system_state into a system_state for the coupled model
    # by redistributing the particles
    Nspecies = len(system_state)
    Nx = len(coupling[0])
    system_state_coupled = [None]*2*Nspecies
    system_state_coupled[0:Nspecies] = system_state
    system_state_coupled[Nspecies:2*Nspecies] = np.zeros((Nspecies,Nx))
    # this could be cleaned up a bit
    for i in range(Nx):
        if i == get_coarseMesh_voxel(i,coupling):
            for k in range(Nspecies):
                n = 0
                for j in range(Nx):
                    if coupling[i,j]>0:
                        n = n + system_state[k,j]
                system_state_coupled[Nspecies+k][i] = n
    return system_state_coupled
