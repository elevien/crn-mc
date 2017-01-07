import numpy as np
from .mesh import *


class Species:

    def __init__(self,name,scale,mesh,value):
        self.name = name
        self.scale = scale
        self.value = value
        self.mesh = mesh

    def __str__(self):
        return self.name
