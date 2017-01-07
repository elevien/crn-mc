import numpy as np
#from .mesh import *


"""

"""

 class Species:
     def __init__(self,name,exponent,mesh,value):
         self.name = name
         self.exponent = exponent
         self.value = value
         self.mesh = mesh


class Model:

    """
    Class wrapping all the static information of a biohchemical model
    """

    def __init__(self,Nspecies,mesh):
        self.mesh = mesh
        self.Nspecies = Nspecies
        self.ss_d1 = Nspecies
        self.ss_d2 = mesh.Nvoxels
        # self.systemState = []
        self.systemState = np.zeros((self.Nspecies,self.mesh.Nvoxels))
        self.events = []

    # def addSpecies(self,name,exponent,mesh,value):
    #     species = Species(name,exponent,mesh,value)
    #     self.systemState.append(species)
    #     return None

    def addReaction(self,reactants,products,intensity,scale):
        """ Add new reaction
        Input:
            - reactants [numpy array]
            - products [numpy array]
            - intensity [float]
        """
        for i in range(self.mesh.Nvoxels):
            reaction = Reaction(self,i,reactants,products,intensity,scale)
            self.events.append(reaction)
        return None

    def addDiffusions(self,species,diffusivity):

        """ Add diffusive channel
        Input:
            - reactants [numpy array]
            - products [numpy array]
            - intensity [float]
        """
        for i in range(self.mesh.Nvoxels):
            for j in range(self.mesh.Nvoxels):
                if self.mesh.topology[i,j]>0:
                    diffusion = Diffusion(self,i,j,species,diffusivity)
                    self.events.append(diffusion)
        return None

    def makeStateMatrix():
        #
        state = np.zeros((self.Nspecies,self.mesh.Nvoxels))
        for i in



class ModelHybrid(Model):

    """
    Clas wrapping all the static information of a biohchemical model
    """
    def __init__(self,Nspecies,mesh):
        self.mesh = mesh

        self.Nspecies = Nspecies
        self.ss_d1 = Nspecies
        self.ss_d2 = mesh.Nvoxels
        self.systemState = np.zeros((self.Nspecies,self.mesh.Nvoxels))

        self.eventsFast = []
        self.eventsSlow = []

    def addReactionSlow(self,reactants,products,intensity,scale):
        for i in range(self.mesh.Nvoxels):
            reaction = Reaction(self,i,reactants,products,intensity,scale)
            self.eventsSlow.append(reaction)
        return None

    def addReactionFast(self,reactants,products,intensity):
        for i in range(self.mesh.Nvoxels):
            reaction = Reaction(self,i,reactants,products,intensity,scale)
            self.eventsFast.append(reaction)
        return None

class ModelHybridSplitCoupled(Model):
    # not yet implemented for spatial model
    def __init__(self,Nspecies,mesh):
        self.mesh = mesh
        #self.Nspecies_fast = Nspecies_fast
        #self.Nspecies_slow = Nspecies_slow
        self.Nspecies = Nspecies
        self.ss_d1 = 2*Nspecies
        self.ss_d2 = mesh.Nvoxels
        self.systemState = np.zeros((2*self.Nspecies,self.mesh.Nvoxels))
        self.eventsFast = []
        self.eventsSlow = []
        self.eventsUncoupled = []

    def addReactionSlowUncoupled(self,reactants,products,intensity,scale):
        for i in range(self.mesh.Nvoxels):
            reaction = ReactionHybridFast_Exact(self,i,reactants,products,intensity)
            self.eventsSlow.append(reaction)
        return None

    def addReactionSlowCoupled(self,reactants,products,intensity,scale):
        for i in range(self.mesh.Nvoxels):
            reaction_common = ReactionHybridSlow_SplitCommon(self,i,reactants,products,intensity)
            self.eventsSlow.append(reaction_common)
            reaction_fast = ReactionHybridSlow_SplitHybrid(self,i,reactants,products,intensity)
            self.eventsSlow.append(reaction_fast)
            reaction_slow = ReactionHybridSlow_SplitExact(self,i,reactants,products,intensity)
            self.eventsSlow.append(reaction_slow)
        return None

    def addReactionFast(self,reactants,products,intensity,scale):
        for i in range(self.mesh.Nvoxels):
            reaction = ReactionHybridFast_Hybrid(self,i,reactants,products,intensity)
            self.eventsFast.append(reaction)
            reaction = ReactionHybridFast_Exact(self,i,reactants,products,intensity)
            self.eventsSlow.append(reaction)
        return None

####################################################################################


global exp_max
exp_max =  10e20

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
        self.updateRate()

    def updateRate(self):
        return None


class Diffusion(Event):
    def __init__(self,model,voxel,voxel_in,species,diffusivity):
        self.voxel = voxel
        self.voxel_in = voxel_in
        self.species = species
        self.diffusivity = diffusivity
        self.direction = np.zeros((len(model.systemState),model.mesh.Nvoxels))
        self.direction[species] = -np.identity(model.mesh.Nvoxels)[self.voxel]+np.identity(model.mesh.Nvoxels)[self.voxel_in]
        super().__init__(model)

    def __str__(self):
        return "Diffusion of species "+str(self.species)+" from voxel "+str(self.voxel)+" to "+str(self.voxel_in)

    def updateRate(self):
        self.rate = self.diffusivity*self.model.systemState[self.species][self.voxel]
        return None


class Reaction(Event):
    def __init__(self,model,voxel,reactants,products,intensity,scale):
        self.voxel = voxel
        self.reactants = reactants
        self.products = products
        self.intensity = intensity
        self.scale = scale
        self.direction = np.zeros((len(model.systemState),model.mesh.Nvoxels))
        self.direction[:,self.voxel] = products-reactants
        super().__init__(model)
    def __str__(self):
        return "Reaction in voxel "+str(self.voxel)

    def updateRate(self):
        # works for order =1,2
        a = self.intensity
        for i in range(self.model.Nspecies):
            if self.reactants[i]>0:
                a = a*pow(self.model.systemState[i,self.voxel],self.reactants[i])
        self.rate = a
        return None

#-----------------------------------------------------------------------------------------
# hybrid split coupling


class ReactionHybridFast_Hybrid(Event):
    def __init__(self,model,voxel,reactants,products,intensity,scale):
        self.voxel = voxel
        self.reactants = reactants
        self.products = products
        self.intensity = intensity
        self.scale = scale
        self.direction = np.zeros((len(model.systemState),model.mesh.Nvoxels))
        self.direction[0:model.Nspecies,self.voxel] = products-reactants
        super().__init__(model)
    def __str__(self):
        return "Reaction in voxel "+str(self.voxel)

    def updateRate(self):
        a1 = self.intensity
        for i in range(self.model.Nspecies):
            if self.reactants[i]>0:
                a1 = a1*pow(self.model.systemState[i,self.voxel],self.reactants[i])
        self.rate = a1
        return None


class ReactionHybridFast_Exact(Event):
    def __init__(self,model,voxel,reactants,products,intensity,scale):
        self.voxel = voxel
        self.reactants = reactants
        self.products = products
        self.intensity = intensity
        self.scale = scale
        self.direction = np.zeros((len(model.systemState),model.mesh.Nvoxels))
        self.direction[model.Nspecies:2*model.Nspecies,self.voxel] = products-reactants
        super().__init__(model)
    def __str__(self):
        return "Reaction in voxel "+str(self.voxel)

    def updateRate(self):
        a2 = self.intensity
        for i in range(self.model.Nspecies):
            if self.reactants[i]>0:
                a2 = a2*pow(self.model.systemState[self.model.Nspecies+i,self.voxel],self.reactants[i])
        self.rate = a2
        return None


class ReactionHybridSlow_SplitCommon(Event):
    def __init__(self,model,voxel,reactants,products,intensity,scale):
        self.voxel = voxel
        self.reactants = reactants
        self.products = products
        self.intensity = intensity
        self.scale = scale
        self.direction = np.zeros((len(model.systemState),model.mesh.Nvoxels))
        self.direction[0:model.Nspecies,self.voxel] = products-reactants
        self.direction[model.Nspecies:2*model.Nspecies,self.voxel] = products-reactants
        super().__init__(model)
    def __str__(self):
        return "Reaction in voxel "+str(self.voxel)

    def updateRate(self):
        a1 = self.intensity
        a2 = self.intensity
        for i in range(self.model.Nspecies):
            if self.reactants[i]>0:
                a1 = a1*pow(self.model.systemState[i,self.voxel],self.reactants[i])
                a2 = a2*pow(self.model.systemState[self.model.Nspecies+i,self.voxel],self.reactants[i])
        self.rate = min(a1,a2)
        return None

class ReactionHybridSlow_SplitHybrid(Event):
    def __init__(self,model,voxel,reactants,products,intensity,scale):
        self.voxel = voxel
        self.reactants = reactants
        self.products = products
        self.intensity = intensity
        self.scale = scale
        self.direction = np.zeros((len(model.systemState),model.mesh.Nvoxels))
        self.direction[0:model.Nspecies,self.voxel] = products-reactants
        super().__init__(model)


    def updateRate(self):
        a1 = self.intensity
        a2 = self.intensity
        for i in range(self.model.Nspecies):
            if self.reactants[i]>0:
                a1 = a1*pow(self.model.systemState[i,self.voxel],self.reactants[i])
                a2 = a2*pow(self.model.systemState[self.model.Nspecies+i,self.voxel],self.reactants[i])
        self.rate = rho(a1,a2)
        return None

class ReactionHybridSlow_SplitExact(Event):
    def __init__(self,model,voxel,reactants,products,intensity,scale):
        self.voxel = voxel
        self.reactants = reactants
        self.products = products
        self.intensity = intensity
        self.scale = scale
        self.direction = np.zeros((len(model.systemState),model.mesh.Nvoxels))
        self.direction[model.Nspecies:2*model.Nspecies,self.voxel] = products-reactants
        super().__init__(model)

    def updateRate(self):
        a1 = self.intensity
        a2 = self.intensity
        for i in range(self.model.Nspecies):
            if self.reactants[i]>0:
                a1 = a1*pow(self.model.systemState[i,self.voxel],self.reactants[i])
                a2 = a2*pow(self.model.systemState[self.model.Nspecies+i,self.voxel],self.reactants[i])
        self.rate = rho(a2,a1)
        return None
