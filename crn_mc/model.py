import numpy as np
import copy
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

class ModelHybridSplitCoupled(Model):


    def __init__(self,model):
        self.mesh = copy.deepcopy(model.mesh)
        self.systemSize = copy.deepcopy(model.systemSize)
        self.Nspecies = 0
        self.systemState = []
        self.events = []
        #self.systemState = copy.deepcopy(model.systemState)
        #self.events = copy.deepcopy(model.events)
        self.setupCoupling(model)


    def setupCoupling(self,model):
            #  """
            # Set up the coupling between model and it's exact
            # version. For each fast reaction we make a slow reaction in
            # the coupled model, while each slow reaction is brown up into
            # 3 reaction channels via the split coupling method.
            #
        # Keyword arguments:
        # model -- the model to couple
        #  """
        #for s in model.systemState:
        #    self.addSpecies(s.name,s.scale,s.value)

        for s in model.systemState:
            self.addSpecies(s.name,s.scale,s.value)
            self.addSpecies(s.name+'-coupled',s.scale,s.value)

        for e in model.events:
            if e.speed == FAST:
                # make SLOW version for coupled model
                e_old = copy.deepcopy(e)
                e_new = copy.deepcopy(e)
                for s in e_new.reactants:
                    newname = s[0].name+'-coupled'
                    species = list(filter(lambda s: s.name == newname, self.systemState))[0]
                    s[0] = species
                for s in e_new.products:
                    newname = s[0].name+'-coupled'
                    species = list(filter(lambda s: s.name == newname, self.systemState))[0]
                    s[0] = species
                e_new.speed = SLOW
                self.events.append(e_old)
                self.events.append(e_new)
            elif e.speed == SLOW:
                # make only slow reactions
                e_slow = copy.deepcopy(e)
                for s in e_slow.reactants:
                    newname = s[0].name+'-coupled'
                    species = list(filter(lambda s: s.name == newname, self.systemState))[0]
                    s[0] = species
                for s in e_slow.products:
                    newname = s[0].name+'-coupled'
                    species = list(filter(lambda s: s.name == newname, self.systemState))[0]
                    s[0] = species
                e_fast = copy.deepcopy(e)
                for s in e_fast.reactants:
                    newname = s[0].name
                    species = list(filter(lambda s: s.name == newname, self.systemState))[0]
                    s[0] = species
                for s in e_fast.products:
                    newname = s[0].name
                    species = list(filter(lambda s: s.name == newname, self.systemState))[0]
                    s[0] = species
                # make common part
                e_common = copy.deepcopy(e)
                for s in e_common.reactants:
                    newname = s[0].name+'-coupled'
                    species = list(filter(lambda s: s.name == newname, self.systemState))[0]
                    coeff =  s[1]
                    e_common.reactants = np.append(e_common.reactants,[[species,coeff]],axis=0)
                for s in e_common.products:
                    newname = s[0].name+'-coupled'
                    species = list(filter(lambda s: s.name == newname, self.systemState))[0]
                    coeff =  s[1]
                    e_common.products = np.append(e_common.products,[[species,coeff]],axis=0)
                print(e_common)
                e_slow.speed = COUPLED_SLOW
                e_common.speed = COUPLED_COMMON
                e_fast.speed = COUPLED_FAST
                self.events.append(e_slow)
                self.events.append(e_fast)
                self.events.append(e_common)
        return None
        #
        #     elif e.speed == SMALL:
