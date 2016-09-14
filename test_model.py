from mesh import *
from pylab import *

topology = array([[0,1],[1,0]])
M = Mesh(topology)
for i in M.voxels:
    print(i.neighbors)
