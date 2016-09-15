
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
        self.update_rate(model.system_state)
        self.wait_absolute = (self.wait_internal-self.time_internal)/self.rate

    def fire(self,delta):
        self.update_rate(model.system_state)
        self.wait_internal = exponential0(1.)
        self.time_internal = self.time_internal + self.rate*delta
        self.update_wait_absolute()
        return None

    def no_fire(self,delta):
        self.update_rate(model.system_state)
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

    def update_rate(self,system_state):
        self.rate = self.model.system_state[self.species][self.voxel_out]
        return None


class Reaction(Event):
    def __init__(self,model,voxel,reactants,products,name):
        self.voxel = voxel
        self.name  = name
        self.reactants = reactants
        self.products = products
        self.stoichiometric_coeffs = np.zeros((model.species,model.mesh.size))
        self.stoichiometric_coeffs[:,self.voxel] = products-reactants
        super().__init__(model)
    def __str__(self):
        return "Reaction {"+self.name+"} in voxel "+str(self.voxel)

    def update_rate(self,system_state):
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
    def __init__(self,mesh_fine,mesh_coarse,voxel_out_fine,voxel_in_fine,voxel_out_coarse,voxel_in_coarse,species,system_state):
        self.voxel_out_fine = voxel_out_fine
        self.voxel_in_fine = voxel_in_fine
        self.voxel_out_coarse = voxel_out_coarse
        self.voxel_in_coarse = voxel_in_coarse
        self.species = species
        self.stoichiometric_coeffs = [\
         -np.identity(mesh_fine)[voxel_out_fine]+np.identity(mesh_fine)[voxel_in_fine],0,\
         -np.identity(mesh_coarse)[voxel_out_coarse]+np.identity(mesh_coarse)[voxel_in_coarse],0\
         ]
        super().__init__(system_state)

    def update_rate(self,system_state):
        a1 = system_state[self.species][self.voxel_out_fine]
        a2 = system_state[1+self.species][self.voxel_out_coarse]
        self.rate = 20*min(a1,a2)
        return None


class Diffusion_SplitFine(Event):
    def __init__(self,mesh_fine,mesh_coarse,voxel_out_fine,voxel_in_fine,voxel_out_coarse,voxel_in_coarse,species,system_state):
        self.voxel_out_fine = voxel_out_fine
        self.voxel_in_fine = voxel_in_fine
        self.voxel_out_coarse = voxel_out_coarse
        self.voxel_in_coarse = voxel_in_coarse
        self.species = species
        self.stoichiometric_coeffs = [\
         -np.identity(mesh_fine)[voxel_out_fine]+np.identity(mesh_fine)[voxel_in_fine],0,0,0\
         ]
        super().__init__(system_state)

    def update_rate(self,system_state):
        a1 = system_state[self.species][self.voxel_out_fine]
        a2 = system_state[1+self.species][self.voxel_out_coarse]
        self.rate = 20*rho(a1,a2)
        return None

class Diffusion_SplitCoarse(Event):
    def __init__(self,mesh_fine,mesh_coarse,voxel_out_fine,voxel_in_fine,voxel_out_coarse,voxel_in_coarse,species,system_state):
        self.voxel_out_fine = voxel_out_fine
        self.voxel_in_fine = voxel_in_fine
        self.voxel_out_coarse = voxel_out_coarse
        self.voxel_in_coarse = voxel_in_coarse
        self.species = species
        self.stoichiometric_coeffs= [\
         0,0,-np.identity(mesh_coarse)[voxel_out_coarse]+np.identity(mesh_coarse)[voxel_in_coarse],0\
         ]
        super().__init__(system_state)

    def update_rate(self,system_state):
        a1 = system_state[self.species][self.voxel_out_fine]
        a2 = system_state[2+self.species][self.voxel_out_coarse]
        self.rate = 20*rho(a2,a1)
        return None


class Reaction_SplitCommon(Event):
    def __init__(self,mesh_fine,mesh_coarse,refinement,voxel_fine,voxel_coarse,order,stoichiometric_coeffs,system_state):
        self.voxel_fine = voxel_fine
        self.voxel_coarse = voxel_coarse
        self.order = order
        self.stoichiometric_coeffs = stoichiometric_coeffs
        self.refinement = refinement
        super().__init__(system_state)

    def update_rate(self,system_state):
        # works for order =1,2
        a1 = 1.
        for i in range(self.order):
            a1 = a1*system_state[i][self.voxel_fine]

        a2 = 1.
        for i in range(self.order):
            a2 = a1*system_state[2+i][self.voxel_coarse]

        self.rate = min(a1,a2/self.refinement)
        return None

class Reaction_SplitFine(Event):
    def __init__(self,mesh_fine,mesh_coarse,refinement,voxel_fine,voxel_coarse,order,stoichiometric_coeffs,system_state):
        self.voxel_fine = voxel_fine
        self.voxel_coarse = voxel_coarse
        self.order = order
        self.stoichiometric_coeffs = stoichiometric_coeffs
        self.refinement = refinement
        super().__init__(system_state)

    def update_rate(self,system_state):
        # works for order =1,2
        a1 = 1.
        for i in range(self.order):
            a1 = a1*system_state[i][self.voxel_fine]

        a2 = 1.
        for i in range(self.order):
            a2 = a1*system_state[2+i][self.voxel_coarse]

        self.rate = rho(a1,a2/self.refinement)
        return None


class Reaction_SplitCoarse(Event):
    def __init__(self,mesh_fine,mesh_coarse,refinement,voxels_fine,voxel_coarse,order,stoichiometric_coeffs,system_state):
        self.voxels_fine = voxels_fine
        self.voxel_coarse = voxel_coarse
        self.order = order
        self.stoichiometric_coeffs = stoichiometric_coeffs
        self.refinement = refinement
        super().__init__(system_state)

    def update_rate(self,system_state):
        # works for order =1,2
        a1 = 0.
        for j in self.voxels_fine:
            aa = 1.
            for i in range(self.order):
                aa = aa*system_state[i][j]
            a1 = a1 + aa

        a2 = 1.
        for i in range(self.order):
            a2 = a1*system_state[2+i][self.voxel_coarse]

        self.rate = rho(a1,a2)
        return None
