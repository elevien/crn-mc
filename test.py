import matplotlib.pyplot as plt
from coupled_bimolecular import *


N = 100
J = 4
level = 0
w = 1.
T = 1.
X0,X1,Y0,Y1 = path_coupled(N,J,level,w,T)
