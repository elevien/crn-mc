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

    def updateRate(self):

        # works for order =1,2
        if self.speed == SLOW:
            self.computeRateSlow()
        if self.speed == SMALL:
            self.computeRateSlow()
        elif self.speed == FAST:

            self.computeRateFast()
        elif self.speed == COUPLED_SLOW:

            self.computeRateCoupledSlow()
        elif self.speed == COUPLED_FAST:
            self.computeRateCoupledFast()
        elif self.speed == COUPLED_COMMON:
            self.computeRateCoupledCommon()
        return None

    def computeRateSlow(self):
        a = self.intensity
        for s in self.reactants:
                base = s[0].value[self.voxel]
                exp = float(s[1])
                a = a*pow(base,exp)
        self.rate = a/self.scale
        return None

    def computeRateFast(self):
        a = self.intensity
        for s in self.reactants:
            base = s[0].value[self.voxel]
            exp = float(s[1])
            a = a*pow(base,exp)
        self.rate = a
        return None


    def computeRateCoupledSlow(self):
        a_slow = self.intensity
        a_fast = self.intensity

        for s in self.reactants[0:len(self.reactants)]:
            base = s[0].value[self.voxel]
            exp = float(s[1])
            a_slow = a_slow*pow(base,exp)

        for s in self.reactants[len(self.reactants):2*len(self.reactants)]:
            base = s[0].value[self.voxel]
            exp = float(s[1])
            a_fast = a_fast*pow(base,exp)

        self.rate = a_slow - min(a_slow,a_fast)
        return None

    def computeRateCoupledFast(self):
        a_slow = self.intensity
        a_fast = self.intensity

        for s in self.reactants[0:len(self.reactants)]:
            base = s[0].value[self.voxel]
            exp = float(s[1])
            a_slow = a_slow*pow(base,exp)
        for s in self.reactants[len(self.reactants):2*len(self.reactants)]:
            base = s[0].value[self.voxel]
            exp = float(s[1])
            a_fast = a_fast*pow(base,exp)

        self.rate = a_fast - min(a_slow,a_fast)
        return None

    def computeRateCoupledCommon(self):
        a_slow = self.intensity
        a_fast = self.intensity

        for s in self.reactants[0:len(self.reactants)]:
            base = s[0].value[self.voxel]
            exp = float(s[1])
            a_slow = a_slow*pow(base,exp)
        for s in self.reactants[len(self.reactants):2*len(self.reactants)]:
            base = s[0].value[self.voxel]
            exp = float(s[1])
            a_fast = a_fast*pow(base,exp)

        self.rate = min(a_slow,a_fast)
        return None
