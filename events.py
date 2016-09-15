
import numpy as np
global exp_max
exp_max =  1000000000.
def rho(u,v):
    return u - min(u,v)

def exponential0(rate):
    if (rate <= 0):
        return exp_max
    else:
        return np.random.exponential(1./rate)

class Event:
    def __init__(self,model):
        self.model = model
        self.time_internal = 0.
        self.wait_internal = exponential0(1.)
        self.update_rate()
        self.wait_absolute = (self.wait_internal-self.time_internal)/self.rate

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
        self.rate = 0.
        return None

    def update_wait_absolute(self):
        if self.rate>0:
            self.wait_absolute = self.wait_internal/self.rate
        else:
            self.wait_absolute = exp_max
        return None

class Diffusion(Event):
    def __init__(self,model,voxel_in,voxel_out,species):
        self.voxel_out = voxel_out
        self.voxel_in = voxel_in
        self.species = species
        self.stoichiometric_coeffs = np.zeros((model.species,model.mesh.size))
        self.stoichiometric_coeffs[species] = -np.identity(model.mesh.size)[voxel_out]+np.identity(model.mesh.size)[voxel_in]
        super().__init__(model)

    def __str__(self):
        return "Diffusion of species "+str(self.species)+" from voxel "+str(self.voxel_in)+" to "+str(self.voxel_out)

    def update_rate(self):
        self.rate = self.model.system_state[self.species][self.voxel_out]
        return None


class Reaction(Event):
    def __init__(self,model,voxel,reactants,products):
        self.voxel = voxel
        self.reactants = reactants
        self.products = products
        self.stoichiometric_coeffs = np.zeros((model.species,model.mesh.size))
        self.stoichiometric_coeffs[:,self.voxel] = products-reactants
        super().__init__(model)
    def __str__(self):
        return "Reaction in voxel "+str(self.voxel)

    def update_rate(self):
        # works for order =1,2
        a = 1.
        for i in range(self.model.species):
            if self.reactants[i]>0:
                a = a*self.model.system_state[i,self.voxel]
        self.rate = a
        return None

#-----------------------------------------------------------------------------------------
# split coupling

class Diffusion_SplitCommon(Event):
    def __init__(self,model,voxel_in,voxel_out,species):
        self.voxel_out = voxel_out
        self.voxel_in = voxel_in
        self.species = species
        self.stoichiometric_coeffs = np.zeros((model.species*2,model.mesh.size))
        self.stoichiometric_coeffs[species] = -np.identity(model.mesh.size)[voxel_out]+np.identity(model.mesh.size)[voxel_in]
        self.stoichiometric_coeffs[2*species] = -np.identity(model.mesh.size)[voxel_out]+np.identity(model.mesh.size)[voxel_in]
        super().__init__(model)

    def update_rate(self,system_state):
        a1 = self.model.system_state[self.species][self.voxel_out]
        a2 = self.model.system_state[2*self.species][self.voxel_out]
        self.rate = min(a1,a2)
        return None


class Diffusion_SplitFine(Event):
    def __init__(self,model,voxel_in,voxel_out,species):
        self.voxel_out = voxel_out
        self.voxel_in = voxel_in
        self.species = species
        self.stoichiometric_coeffs = np.zeros((model.species*2,model.mesh.size))
        self.stoichiometric_coeffs[species] = -np.identity(model.mesh.size)[voxel_out]+np.identity(model.mesh.size)[voxel_in]
        super().__init__(model)

    def update_rate(self,system_state):
        a1 = system_state[self.species][self.voxel_out_fine]
        a2 = system_state[2*self.species][self.voxel_out_coarse]
        self.rate = rho(a1,a2)
        return None

class Diffusion_SplitCoarse(Event):
    def __init__(self,model,voxel_in,voxel_out,species):
        self.voxel_out = voxel_out
        self.voxel_in = voxel_in
        self.species = species
        self.stoichiometric_coeffs = np.zeros((model.species*2,model.mesh.size))
        self.stoichiometric_coeffs[2*species] = -np.identity(model.mesh.size)[voxel_out]+np.identity(model.mesh.size)[voxel_in]
        super().__init__(model)

    def update_rate(self,system_state):
        a1 = system_state[self.species][self.voxel_out_fine]
        a2 = system_state[2*self.species][self.voxel_out_coarse]
        self.rate = rho(a2,a1)
        return None


class Reaction_SplitCommon(Event):
    def __init__(self,i,reactants,products):
        self.voxel = voxel
        self.reactants = reactants
        self.products = products
        self.stoichiometric_coeffs = np.zeros((2*model.species,model.mesh.size))
        self.stoichiometric_coeffs[0:self.species,self.voxel] = products-reactants
        self.stoichiometric_coeffs[self.species*2:self.species,self.voxel] = products-reactants
        super().__init__(model)

    def update_rate(self,system_state):
        # works for order =1,2
        a1 = 1.
        a2 = 1.
        for i in range(self.model.species):
            if self.reactants[i]>0:
                a1 = a1*self.model.system_state[i,self.voxel]
                a2 = a2*self.model.system_state[i,self.voxel]
        self.rate = min(a1,a2/len(self.coupling[0]))
        return None

class Reaction_SplitFine(Event):
    def __init__(self,i,reactants,products):
        self.voxel = voxel
        self.reactants = reactants
        self.products = products
        self.stoichiometric_coeffs = np.zeros((2*model.species,model.mesh.size))
        self.stoichiometric_coeffs[0:self.species,self.voxel] = products-reactants
        self.stoichiometric_coeffs[self.species*2:self.species,self.voxel] = products-reactants
        super().__init__(model)

    def update_rate(self,system_state):
        # works for order =1,2
        a1 = 1.
        a2 = 1.
        for i in range(self.model.species):
            if self.reactants[i]>0:
                a1 = a1*self.model.system_state[i,self.voxel]
                a2 = a2*self.model.system_state[i,self.voxel]
        self.rho = min(a1,a2/len(self.coupling[0]))
        return None



class Reaction_SplitCoarse(Event):
    def __init__(self,mesh_fine,mesh_coarse,refinement,voxels_fine,voxel_coarse,order,stoichiometric_coeffs,system_state):
        self.voxel = voxel # this is the coarse "voxel"
        self.reactants = reactants
        self.products = products
        self.stoichiometric_coeffs = np.zeros((2*model.species,model.mesh.size))
        self.stoichiometric_coeffs[0:self.species,self.voxel] = products-reactants
        self.stoichiometric_coeffs[self.species*2:self.species,self.voxel] = products-reactants
        super().__init__(model)

    def update_rate(self,system_state):
        # works for order =1,2
        a1 = 0.
        for j in self.model.coupling[self.voxel]:
            aa = 1.
            for i in range(self.model.species):
                if self.reactants[i]>0:
                    aa = aa*self.model.system_state[i,j]
            a1 = a1 + aa

        a2 = 1.
        for i in range(self.order):
            for i in range(self.model.species):
                if self.reactants[i]>0:
                    a2 = a2*self.model.system_state[i,self.voxel[0]]
        self.rate = rho(a1,a2)
        return None
