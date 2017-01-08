import numpy as np
import copy
from .mesh import *
from .events import *
from .species import *


class Model:
    """ Contains all the static information about a biochemical model """

    def __init__(self,mesh,systemSize):
        self.mesh = mesh
        self.systemSize = systemSize
        self.Nspecies = 0
        self.systemState = []
        self.events = []

    def addspecies(self,name,exponent,value):
        """ Adds new species to the model """

        scale = pow(self.systemSize,-exponent)
        species = Species(name,scale,self.mesh,value)
        self.systemState.append(species)
        self.Nspecies = self.Nspecies+1
        return None

    def addreaction(self,reactants_vect,products_vect,intensity,exponent,speed):
        """ Add new reaction to the model

        Input:
            - reactants -- array of tuples (Species name,Integer>0)
            - products -- array of tuples (Species name,Integer>0)
            - intensity -- float
            - exponent -- float
            - speed -- "FAST","SLOW"...

        Output:
            None
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
            scale = pow(self.systemSize,-exponent)
            reaction = Reaction(i,reactants,products,intensity,scale,speed)
            self.events.append(reaction)
        return None

    def getstate(self,voxel):
        """ Return state of each species in voxel """

        state = np.zeros(len(self.systemState))
        for i in range(len(self.systemState)):
            state[i] = self.systemState[i].value[voxel]
        return state
