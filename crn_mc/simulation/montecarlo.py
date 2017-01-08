from ..mesh import *
from ..model import *
from .paths import *
import numpy as np
from scipy.integrate import ode
import copy



def montecarlo_crude(model,T,delta,species,voxel,method,sample_rate):
    """ Obtains statistics of model using a crude monte carlo esimator. """
    samples = []
    standdev = []
    eps = pow(model.systemSize,-delta)
    h = eps

    M0 = 10
    Mmax = 10e5

    # get the intial conditions
    ic = np.zeros(model.dimension)
    for j in range(model.dimension):
        ic[j] = model.systemState[j].value[0]

    # do a fixed number of prelimnary simulations
    for i in range(M0):
        path,clock= makepath(model,T,h,method,sample_rate,voxel)
        samples = np.append(samples,path[-1,species])
        # reset initial conditions
        for j in range(model.dimension):
            model.systemState[j].value[0] = ic[j]
    i = 0
    standdev.append(np.std(samples))
    while standdev[i-1]>eps and i<Mmax:
        path,clock= makepath(model,T,h,method,sample_rate,voxel)
        samples = np.append(samples,path[-1,species])
        # reset initial conditions
        for j in range(model.dimension):
            model.systemState[j].value[0] = ic[j]
        standdev.append(np.std(samples/(i+M0)))
        i = i+1
    return sum(samples/(i+M0)),standdev[1:-1]

def montecarlo_coupled(model,T,delta,species,voxel,method,sample_rate):
    """ Obtains statistics of model using a coupled monte carlo esimator. """
    samples = []
    standdev = []
    eps = pow(model.systemSize,-delta)
    h = eps

    M0 = 10
    Mmax = 10e5

    # get the intial conditions
    ic = np.zeros(model.dimension)
    for j in range(model.dimension):
        ic[j] = model.systemState[j].value[0]

    # do a fixed number of prelimnary simulations
    for i in range(M0):
        path,clock= makepath_coupled(model,T,h,method,sample_rate,voxel)
        sample = path[-1,species+model.dimension]-path[-1,species]
        samples = np.append(samples,sample)
        # reset initial conditions
        for j in range(model.dimension):
            model.systemState[j].value[0] = ic[j]
    i = 0
    standdev.append(np.std(samples))
    while standdev[i-1]>eps and i<Mmax:
        path,clock= makepath_coupled(model,T,h,method,sample_rate,voxel)
        sample = path[-1,species+model.dimension]-path[-1,species]
        samples = np.append(samples,sample)
        # reset initial conditions
        for j in range(model.dimension):
            model.systemState[j].value[0] = ic[j]
        standdev.append(np.std(samples/(i+M0)))
        i = i+1
    Q1 = sum(samples/(i+M0))
    Q2,a = montecarlo_crude(model,T,delta,species,voxel,method,sample_rate)
    return Q1+Q2,standdev[1:-1]
