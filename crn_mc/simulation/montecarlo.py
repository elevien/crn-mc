from ..mesh import *
from ..model import *
from .paths import *
import sys
import numpy as np
from scipy.integrate import ode
import copy


def identity(arg):
    return arg

def montecarlo(model,T,delta,ode_method='lsoda',sample_rate =0.,estimator='crude',
        path_type='hybrid',func = identity,min_samples=10,max_samples=10e5,
        output=sys.stdout,*args,**kwargs):
    voxel = 0.

    print('######################################################',file=output)
    print('                 running monte carlo                  ',file=output)
    print('######################################################',file=output)
    print('',file=output)
    print('params ',file=output)
    print('------------------------------------------------------',file=output)
    print(' estimator   =',estimator,file=output)
    print(' sample_rate =',str(sample_rate),file=output)
    print(' ode_method  =',ode_method,file=output)
    print(' T           =',str(T),file=output)
    print(' delta       =',str(delta),file=output)
    print(' max_samples =',str(max_samples),file=output)
    print(' min_samples =',str(min_samples),file=output)
    print('',file=output)

    if estimator == 'crude':
        estimate,standdev,event_count= montecarlo_crude(model,T,func,delta,voxel,ode_method,
            sample_rate,path_type,min_samples,max_samples,output)
    elif estimator == 'coupled':
        estimate,standdev,event_count=  montecarlo_coupled(model,T,func,delta,voxel,ode_method,
            sample_rate,min_samples,max_samples,output)

    print('output',file=output)
    print('------------------------------------------------------',file=output)
    print(' estimate     = '+str(estimate),file=output)
    print(' event_count  = '+str(event_count),file=output)
    print('',file=output)
    return estimate,standdev,event_count

def montecarlo_crude(model,T,func,delta,voxel,ode_method,sample_rate,path_type,
        min_samples,max_samples,output):
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

    print('model setup ',file=output)
    print('------------------------------------------------------',file=output)
    print(' system size = '+str(eps),file=output)
    print(' eps         = '+str(eps),file=output)
    print(' events:',file=output)
    for e in model.events:
        print(' '+e.__str__(),file=output)
    print('',file=output)

    print('generating samples...',file=output)
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
        print(' estimator failed to converge',file=output)
    else:
        print(' success!',file=output)
    print(' ',file=output)
    return sum(samples/(i+M0)),standdev[1:i],event_count

def montecarlo_coupled(model,T,func,delta,voxel,ode_method,sample_rate,
        min_samples,max_samples,output):
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

    print('model setup ',file=output)
    print('------------------------------------------------------',file=output)
    print(' system size = '+str(eps),file=output)
    print(' eps         = '+str(eps),file=output)
    print(' events:',file=output)
    for e in model.events:
        print(' '+e.__str__(),file=output)
    print('',file=output)

    print('generating samples...',file=output)
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
        print(' estimator failed to converge',file=output)
    else:
        print(' success!',file=output)
    print(' ',file=output)
    Q1 = sum(samples/i)
    Q2,standdev2,event_count2 = montecarlo_crude(model,T,func,delta,voxel,ode_method,
            sample_rate,'hybrid',min_samples,max_samples,output)
    return Q1+Q2,standdev[1:i],event_count+event_count2
