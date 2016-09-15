import numpy as np
from mesh import *
from events import *

class Model:
    def __init__(self,species,mesh):
        self.mesh = mesh
        self.species = species
        self.system_state = np.zeros((self.species,self.mesh.size))
        self.events = []

    def add_reaction(self,reactants,products,name):
        for i in range(self.mesh.size):
            reaction = Reaction(self,i,reactants,products,name)
            self.events.append(reaction)
        return None

    def init_diffusions(self,species):
        for i in range(self.mesh.size):
            for j in range(self.mesh.size):
                if self.mesh.topology[j,i]>0:
                    diffusion =Diffusion(self,i,j,species)
                    self.events.append(diffusion)
        return None
