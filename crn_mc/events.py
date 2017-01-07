import numpy as np
from .mesh import *


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
    def __init__(self):
        self.updateRate()

    def updateRate(self):
        return None


class Reaction(Event):
    """

    A single stochastic reaction channel

    """
    def __init__(self,voxel,reactants,products,intensity,scale,speed):
        """
        Input:
            - voxel [float]
            - reactants [list of tuples (Species,integer>0)]
            - product [list of tuples (Species,integer>0)]
        """
        self.voxel = voxel
        self.reactants = reactants
        self.products = products
        self.intensity = intensity
        self.scale = scale
        self.speed = speed
        super().__init__()
    def __str__(self):
        return "Reaction in voxel "+str(self.voxel)

    def updateRate(self):
        # works for order =1,2
        if self.speed=="SLOW":
            a = self.intensity
            for s in self.reactants:
                    base = s[0].value[self.voxel]
                    exp = float(s[1])
                    a = a*pow(base,exp)
            self.rate = a/self.scale
        elif self.speed=="FAST":
            a = self.intensity
            for s in self.reactants:
                    base = s[0].value[self.voxel]
                    exp = float(s[1])
                    a = a*pow(base,exp)
            self.rate = a
        return None
