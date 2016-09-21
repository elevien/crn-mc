
import numpy as np
from mesh import *

global exp_max
exp_max =  1000000000.
def rho(u,v):
    return u - min(u,v)

def exponential0(rate):
    if (rate <= 0):
        print("Warning: next reaction at t = infinity")
        return exp_max
    else:
        return np.random.exponential(1./rate)

class Event:
    def __init__(self,model):
        self.model = model
        self.time_internal = 0.
        self.wait_internal = exponential0(1.)
        self.update_rate()
        if self.rate>0:
            self.wait_absolute = (self.wait_internal-self.time_internal)/self.rate
        else:
            self.wait_absolute = exp_max

    def fire(self,delta):
        self.update_rate()
        self.wait_internal = exponential0(1.)
        self.time_internal = self.time_internal + self.rate*delta
        self.update_wait_absolute()
        return None

    def no_fire(self,delta):
        self.update_rate()
        # wait_internal remains unchanged
        self.time_internal = self.time_internal + self.rate*delta
        self.update_wait_absolute()
        return None

    def update_rate(self):
        return None

    def update_wait_absolute(self):
        if self.rate>0:
            self.wait_absolute = self.wait_internal/self.rate
        else:
            self.wait_absolute = exp_max
        return None

class Diffusion(Event):
    def __init__(self,model,voxel,voxel_in,species,diffusivity):
        self.voxel = voxel
        self.voxel_in = voxel_in
        self.species = species
        self.diffusivity = diffusivity
        self.stoichiometric_coeffs = np.zeros((len(model.system_state),model.mesh.Nvoxels))
        self.stoichiometric_coeffs[species] = -np.identity(model.mesh.Nvoxels)[self.voxel]+np.identity(model.mesh.Nvoxels)[self.voxel_in]
        super().__init__(model)

    def __str__(self):
        return "Diffusion of species "+str(self.species)+" from voxel "+str(self.voxel)+" to "+str(self.voxel_in)

    def update_rate(self):
        self.rate = self.diffusivity*self.model.system_state[self.species][self.voxel]
        return None


class Reaction(Event):
    def __init__(self,model,voxel,reactants,products,intensity):
        self.voxel = voxel
        self.reactants = reactants
        self.products = products
        self.intensity = intensity
        self.stoichiometric_coeffs = np.zeros((len(model.system_state),model.mesh.Nvoxels))
        self.stoichiometric_coeffs[:,self.voxel] = products-reactants
        super().__init__(model)
    def __str__(self):
        return "Reaction in voxel "+str(self.voxel)

    def update_rate(self):
        # works for order =1,2
        a = self.intensity
        for i in range(self.model.Nspecies):
            if self.reactants[i]>0:
                a = a*pow(self.model.system_state[i,self.voxel],self.reactants[i])
        self.rate = a
        return None

#-----------------------------------------------------------------------------------------
# split coupling



class Diffusion_SplitCommon(Event):
    def __init__(self,model,voxel,voxel_in,species,diffusivity):
        self.voxel = voxel
        self.voxel_in = voxel_in
        self.voxel_coarse = get_coarseMesh_voxel(self.voxel,model.coupling)
        self.voxel_in_coarse = get_coarseMesh_voxel(self.voxel_in,model.coupling)
        self.species = species
        self.diffusivity = diffusivity
        self.stoichiometric_coeffs = np.zeros((len(model.system_state),model.mesh.Nvoxels))
        self.stoichiometric_coeffs[species] = -np.identity(model.mesh.Nvoxels)[voxel]+np.identity(model.mesh.Nvoxels)[voxel_in]

        self.stoichiometric_coeffs[model.Nspecies+species] = -np.identity(model.mesh.Nvoxels)[self.voxel_coarse]+np.identity(model.mesh.Nvoxels)[self.voxel_in_coarse]
        super().__init__(model)

    def update_rate(self):
        a1 = self.model.system_state[self.species][self.voxel]
        a2 = self.model.system_state[self.model.Nspecies+self.species][self.voxel_coarse]
        self.rate = self.diffusivity*min(a1,a2)
        return None


class Diffusion_SplitFine(Event):
    def __init__(self,model,voxel,voxel_in,species,diffusivity):
        self.voxel = voxel
        self.voxel_in = voxel_in
        self.voxel_coarse = get_coarseMesh_voxel(self.voxel,model.coupling)
        self.voxel_in_coarse = get_coarseMesh_voxel(self.voxel_in,model.coupling)
        self.species = species
        self.diffusivity = diffusivity
        self.stoichiometric_coeffs = np.zeros((len(model.system_state),model.mesh.Nvoxels))
        self.stoichiometric_coeffs[species] = -np.identity(model.mesh.Nvoxels)[self.voxel]+np.identity(model.mesh.Nvoxels)[self.voxel_in]
        super().__init__(model)

    def update_rate(self):
        a1 = self.model.system_state[self.species][self.voxel]
        a2 = self.model.system_state[self.model.Nspecies+self.species][self.voxel_coarse]
        self.rate = self.diffusivity*rho(a1,a2)
        return None

class Diffusion_SplitCoarse(Event):
    def __init__(self,model,voxel,voxel_in,species,diffusivity):
        self.voxel = voxel
        self.voxel_in = voxel_in
        self.voxel_coarse = get_coarseMesh_voxel(self.voxel,model.coupling)
        self.voxel_in_coarse = get_coarseMesh_voxel(self.voxel_in,model.coupling)
        self.species = species
        self.diffusivity = diffusivity
        self.stoichiometric_coeffs = np.zeros((len(model.system_state),model.mesh.Nvoxels))
        self.stoichiometric_coeffs[model.Nspecies+species] = -np.identity(model.mesh.Nvoxels)[self.voxel_coarse]+np.identity(model.mesh.Nvoxels)[self.voxel_in_coarse]
        super().__init__(model)

    def update_rate(self):
        a1 = self.model.system_state[self.species][self.voxel]
        a2 = self.model.system_state[self.model.Nspecies+self.species][self.voxel_coarse]
        self.rate = self.diffusivity*rho(a2,a1)
        return None


class Reaction_SplitCommon(Event):
    def __init__(self,model,voxel,reactants,products,intensity):
        self.voxel = voxel
        self.voxel_coarse = get_coarseMesh_voxel(voxel,model.coupling)
        self.reactants = reactants
        self.products = products
        self.intensity = intensity
        self.stoichiometric_coeffs = np.zeros((len(model.system_state),model.mesh.Nvoxels))
        self.stoichiometric_coeffs[0:model.Nspecies,self.voxel] = products-reactants
        self.stoichiometric_coeffs[model.Nspecies:2*model.Nspecies,self.voxel_coarse] = products-reactants
        super().__init__(model)

    def update_rate(self):
        # works for order =1,2
        a1 = self.intensity
        a2 = self.intensity
        for i in range(self.model.Nspecies):
            if self.reactants[i]>0:
                a1 = a1*pow(self.model.system_state[i][self.voxel],self.reactants[i])
                a2 = a2*pow(self.model.system_state[self.model.Nspecies+i][self.voxel_coarse],self.reactants[i])
        self.rate = min(a1,a2/np.count_nonzero(self.model.coupling[self.voxel]))
        return None

class Reaction_SplitFine(Event):
    def __init__(self,model,voxel,reactants,products,intensity):
        self.voxel = voxel
        self.voxel_coarse = get_coarseMesh_voxel(voxel,model.coupling)
        self.reactants = reactants
        self.products = products
        self.intensity = intensity
        self.stoichiometric_coeffs = np.zeros((len(model.system_state),model.mesh.Nvoxels))
        self.stoichiometric_coeffs[0:model.Nspecies,self.voxel] = products-reactants
        super().__init__(model)

    def update_rate(self):
        # works for order =1,2
        a1 = self.intensity
        a2 = self.intensity
        for i in range(self.model.Nspecies):
            if self.reactants[i]>0:
                a1 = a1*pow(self.model.system_state[i][self.voxel],self.reactants[i])
                a2 = a2*pow(self.model.system_state[self.model.Nspecies+i][self.voxel_coarse],self.reactants[i])
        self.rate = rho(a1,a2/len(self.model.coupling[self.voxel]))
        return None



class Reaction_SplitCoarse(Event):
    def __init__(self,model,voxel,reactants,products,intensity):
        self.voxel = voxel
        self.voxel_coarse = get_coarseMesh_voxel(voxel,model.coupling)
        self.reactants = reactants
        self.products = products
        self.intensity = intensity
        self.stoichiometric_coeffs = np.zeros((len(model.system_state),model.mesh.Nvoxels))
        self.stoichiometric_coeffs[model.Nspecies:2*model.Nspecies,self.voxel_coarse] = products-reactants
        super().__init__(model)

    def update_rate(self):
        # works for order =1,2

        #
        a1 = 0.
        for j in range(self.model.mesh.Nvoxels):
            aa = self.intensity
            for i in range(self.model.Nspecies):
                if self.reactants[i]>0 and self.model.coupling[self.voxel,j]>0:
                    aa = aa*pow(self.model.system_state[i][j],self.reactants[i])
            a1 = a1 + aa

        a2 = self.intensity

        for i in range(self.model.Nspecies):
            if self.reactants[i]>0:
                a2 = a2*pow(self.model.system_state[self.model.Nspecies+i][self.voxel_coarse],self.reactants[i])
        self.rate = rho(a2,a1)
        return None
