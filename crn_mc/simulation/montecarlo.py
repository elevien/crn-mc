from ..mesh import *
from ..model import *
from .paths import *
import sys
import numpy as np
from scipy.integrate import ode
import copy
import json


def identity(arg):
    return arg

def montecarlo(model,initial_data,T,delta,ode_method='lsoda',sample_rate =0.,estimator='crude',
        path_type='hybrid',func = identity,min_samples=10,max_samples=10e5,
        output_file=sys.stdout,*args,**kwargs):
    voxel = 0.

    if estimator == 'crude':
        estimate,standdev,event_count= montecarlo_crude(model,initial_data,T,func,delta,voxel,ode_method,
            sample_rate,path_type,min_samples,max_samples,output_file)
    elif estimator == 'coupled':
        estimate,standdev,event_count=  montecarlo_coupled(model,initial_data,T,func,delta,voxel,ode_method,
            sample_rate,min_samples,max_samples,output_file)
    # this is all for output
    params_dict = {'estimator':estimator,'sample_rate':sample_rate,
        'ode_method':ode_method,'T':T,'delta':delta,
        'max_samples':max_samples,'min_samples':min_samples}
    model_dict = {'system_size':model.systemSize,'events':list([e.__str__() for e in model.events]),
        'initial_data':list(initial_data)}
    results_dict = {'estimate':{model.systemState[i].name:estimate[i] for i in range(model.dimension)},
        'event_count':event_count,'standdev':list(standdev)}
    output = {'params':params_dict,'model':model_dict,'results':results_dict}
    #print('\n',file = output_file)
    print(json.dumps(output),file = output_file)


def montecarlo_crude(model,initial_data,T,func,delta,voxel,ode_method,sample_rate,path_type,
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
    i = 0.
    for j in range(model.dimension):
        model.systemState[j].value[0] = initial_data[j]


    while (standdev[i-1]>eps or i<M0) and i<Mmax:
        path,clock= makepath(model,T,h,ode_method=ode_method,sample_rate = sample_rate,
            path_type=path_type)
        event_count = event_count+len(clock)
        for j in range(model.dimension):
            # evalute f on each species to obtain samples
            samples[i,j] = func(path[-1,j])
            new_standdevs[j] = np.std(samples[:,j]/i)
            # reset initial conditions
            model.systemState[j].value[0] = initial_data[j]
        standdev[i] = max(new_standdevs)
        i = i+1
    return sum(samples/(i+M0)),standdev[1:i],event_count

def montecarlo_coupled(model,initial_data,T,func,delta,voxel,ode_method,sample_rate,
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
    i = 0.

    for j in range(model.dimension):
        model.systemState[j].value[0] = initial_data[j]
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
            model.systemState[j].value[0] = initial_data[j]
        standdev[i] = max(new_standdevs)
        i = i+1

    Q1 = sum(samples/i)
    Q2,standdev2,event_count2 = montecarlo_crude(model,initial_data,T,func,delta,voxel,ode_method,
            sample_rate,'hybrid',min_samples,max_samples,output_file)
    return Q1+Q2,standdev[1:i],event_count+event_count2
