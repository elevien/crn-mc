import numpy as np
from .mesh import *


class Species:

    def __init__(self,name,scale,mesh,value):
        self.name = name
        self.scale = scale # {system size}^{exponent}
        self.mesh = mesh
        self.value = value

    def __str__(self):
        return self.name
