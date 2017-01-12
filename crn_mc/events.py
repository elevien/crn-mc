import numpy as np
from .mesh import *
from .species import *



EXP_MAX =  10e20
FAST = "FAST"
SLOW = "SLOW"
NULL = "NULL"


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

    def __init__(self,voxel,reactants,products,intensity,scale):
        """ Reaction constructor.

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

        # determine how to treat the reaction in the thermodynamic limit
        self.hybridType = 0.
        if self.scale == 1.:
            self.hybridType = SLOW
        else:
            self.hybridType = FAST

        # check if reaction vanishes is thermodynamic limit

        # doesn't handle events with \gamma_i  = rho_j = 0 for some (S_0) species
        # but \gamma_i > \rho_j for some (S_1) species
        # in such instances reaction only appears in slow species evolution

        for p in products:
            if p[0].scale>self.scale:
                self.hybridType = NULL

        self.rate = 0.
        super().__init__()

    def __str__(self):
        s = "With scale "+str(self.scale)+" of type "+self.hybridType+": "
        for r in self.reactants:
            s = s+str(r[1])+""+r[0].name + " +"
        s = s[:-1]
        s = s + " -> "
        for p in self.products:
            s = s+str(p[1])+""+p[0].name + " +"
        s = s[:-1]
        return s

    def updaterate(self):
        """ Update the reaction rate based on the values of species involves. """

        if self.hybridType == SLOW:
            self.computerate_slow()
        else:
            self.computerate_fast()

        return None

    def computerate_slow(self):
        """ Compute the reaction rate for slow channels. """

        a = self.intensity
        for s in self.reactants:
                base = s[0].value[self.voxel]
                exp = float(s[1])
                a = a*pow(base,exp)
        self.rate = a*self.scale
        return None

    def computerate_fast(self):
        """ Compute fast rates based on leading order term in ma rates. """

        # this is the right hand side of the reaction rate equations
        # note that the expression should NOT involve any scales
        # todo: add high order terms in rate
        a = self.intensity
        for s in self.reactants:
            base = s[0].value[self.voxel]
            exp = float(s[1])
            a = a*pow(base,exp)
        self.rate = a
        return None

    def react(self):
        """ update species involved in reaction accoding to stoichiometry. """
        if self.hybridType == NULL:
            # NULL reactions need special path, since they don't alter all
            # their products and reactants
            for r in self.reactants:
                if r[0].scale == self.scale:
                    r[0].value[self.voxel] = r[0].value[self.voxel]-(1./r[0].scale)*float(r[1])
            for p in self.products:
                if p[0].scale == self.scale:
                    p[0].value[self.voxel] = p[0].value[self.voxel]+(1./p[0].scale)*float(p[1])
        else:

            for r in self.reactants:
                r[0].value[self.voxel] = r[0].value[self.voxel]-(1./r[0].scale)*float(r[1])
            for p in self.products:
                p[0].value[self.voxel] = p[0].value[self.voxel]+(1./p[0].scale)*float(p[1])
        return None
