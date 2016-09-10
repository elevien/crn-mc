import matplotlib.pyplot as plt
from coupled_bimolecular import *


N = 100
J = 4
level = 0
w = 1.
T = 1.
events = path_uncoupled(N,w,T)
for i in events:
    print(i.wait_absolute)
