from ..mesh import *
from ..model import *
from .paths import *
import numpy as np
from scipy.integrate import ode
import copy


def identity(arg):
    return arg

def montecarlo(model,T,delta,ode_method='lsoda',sample_rate =0.,estimator='crude',
        path_type='hybrid',func = identity,min_samples=10,max_samples=10e5,*args,**kwargs):
    voxel = 0.

    print('######################################################')
    print('                 running monte carlo                  ')
    print('######################################################')
    print('')
    print('params ')
    print('------------------------------------------------------')
    print(' estimator   = '+estimator)
    print(' sample_rate = '+str(sample_rate))
    print(' ode_method  = '+ode_method)
    print(' T           = '+str(T))
    print(' delta       = '+str(delta))
    print(' max_samples = '+str(max_samples))
    print(' min_samples = '+str(min_samples))
    print('')

    if estimator == 'crude':
        estimate,standdev,event_count= montecarlo_crude(model,T,func,delta,voxel,ode_method,
            sample_rate,path_type,min_samples,max_samples)
    elif estimator == 'coupled':
        estimate,standdev,event_count=  montecarlo_coupled(model,T,func,delta,voxel,ode_method,
            sample_rate,min_samples,max_samples)

    print('ouput')
    print('------------------------------------------------------')
    print(' estimate     = '+str(estimate))
    print(' event_count  = '+str(event_count))
    print('')
    return estimate,standdev,event_count

def montecarlo_crude(model,T,func,delta,voxel,ode_method,sample_rate,path_type,
        min_samples,max_samples):
    """ Obtains statistics of model using a crude monte carlo esimator. """
    M0 = 10
    Mmax = 10e5
    samples = np.zeros((Mmax,model.dimension))
    event_count = 0.
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

    print('model setup ')
    print('------------------------------------------------------')
    print(' system size = '+str(eps))
    print(' eps         = '+str(eps))
    print(' events:')
    for e in model.events:
        print(' '+e.__str__())
    print('')

    print('generating samples...')
    while (standdev[i-1]>eps or i<M0) and i<Mmax:
        path,clock= makepath(model,T,h,ode_method=ode_method,sample_rate = sample_rate,
            path_type=path_type)
        event_count = event_count+len(clock)
        for j in range(model.dimension):
            # evalute f on each species to obtain samples
            samples[i,j] = func(path[-1,j])
            new_standdevs[j] = np.std(samples[:,j]/i)
            # reset initial conditions
            model.systemState[j].value[0] = ic[j]
        standdev[i] = max(new_standdevs)
        i = i+1
    if standdev[i-1]>eps:
        print(' estimator failed to converge')
    else:
        print(' success!')
    print(' ')
    return sum(samples/(i+M0)),standdev[1:i],event_count

def montecarlo_coupled(model,T,func,delta,voxel,ode_method,sample_rate,
        min_samples,max_samples):
    """ Obtains statistics of model using a coupled monte carlo esimator. """
    M0 = 10
    Mmax = 10e5
    samples = np.zeros((Mmax,model.dimension))
    event_count = 0.
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

    print('model setup ')
    print('------------------------------------------------------')
    print(' system size = '+str(eps))
    print(' eps         = '+str(eps))
    print(' events:')
    for e in model.events:
        print(' '+e.__str__())
    print('')

    print('generating samples...')
    while (standdev[i-1]>eps or i<M0) and i<Mmax:
        path,clock=makepath(model,T,h,ode_method = ode_method,
                    sample_rate = sample_rate,
                    path_type='coupled')
        event_count = event_count+len(clock)
        for j in range(model.dimension):
            # evalute f on each species to obtain samples
            samples[i,j] = func(path[-1,j])-func(path[-1,j+model.dimension])
            new_standdevs[j] = np.std(samples[:,j]/i)
            # reset initial conditions (remember this model is copied)
            model.systemState[j].value[0] = ic[j]
        standdev[i] = max(new_standdevs)
        i = i+1
    if standdev[i-1]>eps:
        print(' estimator failed to converge')
    else:
        print(' success!')
    print(' ')
    Q1 = sum(samples/i)
    Q2,standdev2,event_count2 = montecarlo_crude(model,T,func,delta,voxel,ode_method,
            sample_rate,'hybrid',min_samples,max_samples)
    return Q1+Q2,standdev[1:i],event_count+event_count2
