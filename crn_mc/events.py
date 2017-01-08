import numpy as np
from .mesh import *
from .species import *



EXP_MAX =  10e20
FAST = 0
SLOW = 1
SMALL = 2
COUPLED_FAST = 3
COUPLED_SLOW = 4
COUPLED_COMMON = 5


class Event:
    """ Base class for all dynamics events in the model. """
    def __init__(self):
        self.updaterate()
    def updaterate(self):
        return None

class NullEvent:
    """ Null event represents nothing happening. Useful in path generation. """
    def __init__(self):
        super().__init__()
    def __str__(self):
        return "Null Reaction"
    def updaterate(self):
        return None


class Reaction(Event):
    """ An event representing the reaction of some number of molecules. """

    def __init__(self,voxel,reactants,products,intensity,scale,speed):
        """ Make a reaction.

        Input:
            - voxel -- integer
            - reactants -- list of tuples of the form (species,integer)
            - product  -- list of tuples of the form (species,integer)

        """

        self.voxel = voxel
        self.reactants = reactants
        self.products = products
        self.intensity = intensity
        self.scale = scale
        self.speed = speed
        self.rate = 0.
        super().__init__()

    def __str__(self):
        s = "In voxel "+str(self.voxel)+" with rate = "+str(self.rate)+": "
        for r in self.reactants:
            s = s+str(r[1])+"*"+r[0].name + " +"
        s = s[:-1]
        s = s + " -> "
        for p in self.products:
            s = s+str(p[1])+"*"+p[0].name + " +"
        s = s[:-1]
        return s

    def updaterate(self):
        """ Update the reaction rate based on the values of species involves. """

        if self.speed == SLOW:
            self.computerate_slow()
        elif self.speed == FAST:
            self.computerate_fast()

        return None

    def computerate_slow(self):
        """ Compute the reaction rate for slow channels. """

        a = self.intensity
        for s in self.reactants:
                base = s[0].value[self.voxel]
                exp = float(s[1])
                a = a*pow(base,exp)
        self.rate = a/self.scale
        return None

    def computerate_fast(self):
        """ Compute fast rates based on leading order term in ma rates. """

        # note that this expression should NOT involve any scales
        a = self.intensity
        for s in self.reactants:
            base = s[0].value[self.voxel]
            exp = float(s[1])
            a = a*pow(base,exp)
        self.rate = a
        return None

    def react(self):
        """ update species involved in reaction accoding to stoichiometry. """

        for r in self.reactants:
            r[0].value[self.voxel] = r[0].value[self.voxel]-r[0].scale*float(r[1])
        for p in self.products:
            p[0].value[self.voxel] = p[0].value[self.voxel]+p[0].scale*float(p[1])
        return None
