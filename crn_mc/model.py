import numpy as np
from .mesh import *
from .events import *
from .species import *


class Model:

    """
    Class wrapping all the static information of a biohchemical model
    """

    def __init__(self,mesh,systemSize):
        self.mesh = mesh
        self.systemSize = systemSize
        self.Nspecies = 0
        self.systemState = []
        #self.systemState = np.zeros((self.Nspecies,self.mesh.Nvoxels))
        self.events = []

    def addSpecies(self,name,exponent,value):
        scale = pow(self.systemSize,-exponent)
        species = Species(name,scale,self.mesh,value)
        self.systemState.append(species)
        self.Nspecies = self.Nspecies+1
        return None

    def addReaction(self,reactants_vect,products_vect,intensity,scale,speed):
        """ Add new reaction
        Input:
            - reactants [numpy array]
            - products [numpy array]
            - intensity [float]
        """
        reactants = []
        products = []
        for r in reactants_vect:
            name = r[0]
            coeff = r[1]
            species = list(filter(lambda s: s.name == r[0], self.systemState))[0]
            reactants.append([species,coeff])
        for p in products_vect:
            name = p[0]
            coeff = p[1]
            species = list(filter(lambda s: s.name == p[0], self.systemState))[0]
            products.append([species,coeff])

        for i in range(self.mesh.Nvoxels):
            scale = pow(self.systemSize,-scale)
            reaction = Reaction(i,reactants,products,intensity,scale,speed)
            self.events.append(reaction)
        return None

    def getStateInVoxel(self,voxel):
        #
        state = np.zeros(len(self.systemState))
        for i in range(len(self.systemState)):
            state[i] = self.systemState[i].value[voxel]
        return state

    def react(self,reaction):
        for r in reaction.reactants:
            r[0].value[reaction.voxel] = r[0].value[reaction.voxel]-r[0].scale*float(r[1])
        for p in reaction.products:
            p[0].value[reaction.voxel] = p[0].value[reaction.voxel]+p[0].scale*float(p[1])
        return None

    # class ModelHybridSplitCoupled(Model):
    #
    #     def __init__(self,mesh,systemSize):
    #         elf.mesh = mesh
    #         self.systemSize = systemSize
    #         self.Nspecies = 0
    #         self.systemState = []
    #         #self.systemState = np.zeros((self.Nspecies,self.mesh.Nvoxels))
    #         self.events = []
    #
    #     def addSpecies(self,name,exponent,value):
    #         scale = pow(self.systemSize,-exponent)
    #         species = Species(name,scale,self.mesh,value)
    #         self.systemState.append(species)
    #         self.Nspecies = self.Nspecies+1
    #         return None
