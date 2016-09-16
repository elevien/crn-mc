import numpy as np
from mesh import *
from events import *

class Model:
    def __init__(self,species,mesh):
        self.mesh = mesh
        self.species = species
        self.system_state = [None]*self.species
        self.events = []

    def add_reaction(self,reactants,products):
        for i in range(self.mesh.size):
            reaction = Reaction(self,i,reactants,products)
            self.events.append(reaction)
        return None

    def add_diffusions(self,species):
        for i in range(self.mesh.size):
            for j in range(self.mesh.size):
                if self.mesh.topology[j,i]>0:
                    diffusion =Diffusion(self,i,j,species)
                    self.events.append(diffusion)
        return None



class SplitCoupled(Model):
    def __init__(self,species,mesh,coupling):
        self.coupling = coupling
        self.mesh = mesh # the fine mesh
        self.species = species
        self.system_state = [None]*self.species
        self.events = []

    def add_reaction(self,reactants,products):
        for i in range(self.mesh.size):
            reaction = Reaction_SplitCommon(self,i,reactants,products)
            self.events.append(reaction)
            reaction = Reaction_SplitFine(self,i,reactants,products)
            self.events.append(reaction)
        for i in range(self.coupling):
            reaction = Reaction_SplitCoarse(self,i,reactants,products)
            self.events.append(reaction)



    def add_diffusions(self,species):
        for i in range(self.mesh.size):
            for j in range(self.mesh.size):
                if self.mesh.topology[j,i]>0:
                    if self.coupling[j,i]>0:
                        # add regular diffusion on fine mesh
                        diffusion = Diffusion(self,i,j,species)
                        self.events.append(diffusion)
                    else:
                        # add coupled diffusion
                        # particiles on the coarse mesh live on a fine mesh but never
                        # touch the fine mesh voxels
                        diffusion = Diffusion_SplitCommon(self,i,j,species)
                        self.events.append(diffusion)
                        diffusion = Diffusion_SplitCoarse(self,i,j,species)
                        self.events.append(diffusion)
                        diffusion = Diffusion_SplitFine(self,i,j,species)
                        self.events.append(diffusion)
        return None
