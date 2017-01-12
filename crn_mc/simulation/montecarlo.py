from ..mesh import *
from ..model import *
from .paths import *
import numpy as np
from scipy.integrate import ode
import copy


def identity(arg):
    return arg

def montecarlo(model,T,delta,method='lsoda',sample_rate =0.,estimator='crude',
        path_type='hybrid',func = identity,*args,**kwargs):

    voxel = 0.
    if estimator == 'crude':
        return montecarlo_crude(model,T,func,delta,voxel,method,sample_rate,path_type)
    elif estimator == 'coupled':
        return montecarlo_coupled(model,T,func,delta,voxel,method,sample_rate)

def montecarlo_crude(model,T,func,delta,voxel,method,sample_rate,path_type):
    """ Obtains statistics of model using a crude monte carlo esimator. """
    M0 = 10
    Mmax = 10e5
    samples = np.zeros((Mmax,model.dimension))
    standdev = np.zeros(Mmax)
    # for storing list new standard deviations
    new_standdevs = np.zeros(model.dimension)
    eps = pow(model.systemSize,-delta)
    h = eps
    # get the intial conditions
    ic = np.zeros(model.dimension)
    for j in range(model.dimension):
        ic[j] = model.systemState[j].value[0]
    i = 0.
    while (standdev[i-1]>eps or i<M0) and i<Mmax:
        path,clock= makepath(model,T,h,method=method,sample_rate = sample_rate,
            path_type=path_type)
        for j in range(model.dimension):
            # evalute f on each species to obtain samples
            samples[i,j] = func(path[-1,j])
            new_standdevs[j] = np.std(samples[:,j]/i)
            # reset initial conditions
            model.systemState[j].value[0] = ic[j]
        standdev[i] = max(new_standdevs)
        i = i+1
    return sum(samples/(i+M0)),standdev[1:i]

def montecarlo_coupled(model,T,func,delta,voxel,method,sample_rate):
    """ Obtains statistics of model using a coupled monte carlo esimator. """
    M0 = 10
    Mmax = 10e5
    samples = np.zeros((Mmax,model.dimension))
    standdev = np.zeros(Mmax)
    # for storing list new standard deviations
    new_standdevs = np.zeros(model.dimension)
    eps = pow(model.systemSize,-delta)
    h = eps
    # get the intial conditions
    ic = np.zeros(model.dimension)
    for j in range(model.dimension):
        ic[j] = model.systemState[j].value[0]
    i = 0.
    while (standdev[i-1]>eps or i<M0) and i<Mmax:
        path,clock=makepath(model,T,h,method = method,
                    sample_rate = sample_rate,
                    path_type='coupled')

        for j in range(model.dimension):
            # evalute f on each species to obtain samples
            samples[i,j] = path[-1,j]-path[-1,j+model.dimension]
            new_standdevs[j] = np.std(samples[:,j]/i)
            # reset initial conditions (remember this model is copied)
            model.systemState[j].value[0] = ic[j]
        standdev[i] = max(new_standdevs)
        i = i+1
    Q1 = sum(samples/i)
    Q2,standdev2 = montecarlo_crude(model,T,func,delta,voxel,method,sample_rate,'hybrid')
    return Q1+Q2,standdev[1:i]
