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
        output_file=sys.stdout,*args,**kwargs):
    voxel = 0.

    print('######################################################',file=output_file)
    print('                 running monte carlo                  ',file=output_file)
    print('######################################################',file=output_file)
    print('',file=output_file)
    print('params ',file=output_file)
    print('------------------------------------------------------',file=output_file)
    print(' estimator   =',estimator,file=output_file)
    print(' sample_rate =',str(sample_rate),file=output_file)
    print(' ode_method  =',ode_method,file=output_file)
    print(' T           =',str(T),file=output_file)
    print(' delta       =',str(delta),file=output_file)
    print(' max_samples =',str(max_samples),file=output_file)
    print(' min_samples =',str(min_samples),file=output_file)
    print('',file=output_file)

    if estimator == 'crude':
        estimate,standdev,event_count= montecarlo_crude(model,T,func,delta,voxel,ode_method,
            sample_rate,path_type,min_samples,max_samples,output_file)
    elif estimator == 'coupled':
        estimate,standdev,event_count=  montecarlo_coupled(model,T,func,delta,voxel,ode_method,
            sample_rate,min_samples,max_samples,output_file)

    print('results',file=output_file)
    print('------------------------------------------------------',file=output_file)
    print(' estimate     = '+str(estimate),file=output_file)
    print(' event_count  = '+str(event_count),file=output_file)
    print('',file=output_file)
    return estimate,standdev,event_count

def montecarlo_crude(model,T,func,delta,voxel,ode_method,sample_rate,path_type,
        min_samples,max_samples,output_file):
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

    print('model setup ',file=output_file)
    print('------------------------------------------------------',file=output_file)
    print(' system size = '+str(model.systemSize),file=output_file)
    print(' eps         = '+str(eps),file=output_file)
    print(' events:',file=output_file)
    for e in model.events:
        print(' '+e.__str__(),file=output_file)
    print('',file=output_file)

    print('generating samples...',file=output_file)
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
        print(' estimator failed to converge',file=output_file)
    else:
        print(' success!',file=output_file)
    print(' ',file=output_file)
    return sum(samples/(i+M0)),standdev[1:i],event_count

def montecarlo_coupled(model,T,func,delta,voxel,ode_method,sample_rate,
        min_samples,max_samples,output_file):
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

    print('model setup ',file=output_file)
    print('------------------------------------------------------',file=output_file)
    print(' system size = '+str(model.systemSize),file=output_file)
    print(' eps         = '+str(eps),file=output_file)
    print(' events:',file=output_file)
    for e in model.events:
        print(' '+e.__str__(),file=output_file)
    print('',file=output_file)

    print('generating samples...',file=output_file)
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
        print(' estimator failed to converge',file=output_file)
    else:
        print(' success!',file=output_file)
    print(' ',file=output_file)
    Q1 = sum(samples/i)
    Q2,standdev2,event_count2 = montecarlo_crude(model,T,func,delta,voxel,ode_method,
            sample_rate,'hybrid',min_samples,max_samples,output_file)
    return Q1+Q2,standdev[1:i],event_count+event_count2
