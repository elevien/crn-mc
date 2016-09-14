from pylab import *
class Mesh:
    def __init__(self,topology):
        self.size = len(topology)
        self.topology = topology
        self.voxels = []

        for i in range(self.size):
            self.voxels.append(Voxel(1.,[]))

        for i in range(self.size):
            neighbors = []
            for j in range(self.size):
                if self.topology[j,i]>0:
                    neighbors.append(self.voxels[j])
            self.voxels[i].add_neighbor(Voxel(1.,neighbors))

class Voxel:
    def __init__(self,volume,neighbors):
        self.volume = volume
        self.neighbors = []

    def add_neighbor(self,voxel):
        self.neighbors.append(voxel)
