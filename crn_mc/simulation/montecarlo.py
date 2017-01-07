from ..mesh import *
from ..model import *
from .paths import *
import numpy as np
from scipy.integrate import ode
import copy



def mc_crude(model,T,Np,delta,species,path_method):
    Q = []
    eps = pow(Np,-delta)
    Err = []
    M0 = 10
    Mmax = 100000
    # do a fixed number of prelimnary simulations
    for i in range(M0):
        path,clock= path_method()
        Q = np.append(Q,path[-1,species])
    i = 0
    Err.append(np.std(Q))
    while Err[i-1]>eps and i<Mmax:
        path,clock= path_method()
        Q = np.append(Q,path[-1,species])
        Err.append(np.std(Q))
        i = i+1
        #print(Err[i-1])
    return sum(Q),Err

def mc_hyrbidCoupled(model_coupled,model_hybrid,T,Np,delta,sample_rate,species):
    Q = []
    eps = pow(Np,-delta)
    #eps = 1.0
    Err = []
    M0 = 10
    Mmax = 100000
    # do a fixed number of prelimnary simulations
    for i in range(M0):
        path,clock = chv(model_coupled,T,eps,'lsoda',sample_rate)
        Q = np.append(Q,path[-1,species]-path[-1,species+model_coupled.Nspecies])
    i = 0
    Err.append(np.std(Q))
    while Err[i-1]>eps and i<Mmax:
        path,clock = chv(model_coupled,T,eps,'lsoda',sample_rate)
        Q = np.append(Q,path[-1,species]-path[-1,species+model_coupled.Nspecies])
        Err.append(np.std(Q))
        i = i+1
    q2,Err2 = mc_crude(model_hybrid,T,Np,delta,species,lambda: chv(model_hybrid,T,eps,'lsoda',sample_rate))

    return sum(Q)+q2,Err

def mc_crudeDiffusions(model,T,eps,delta):

    clock_quantized = np.linspace(0,T,resolution)
    Nt  = len(clock_quantized)
    average_quantized = np.zeros((len(model.systemState),model.mesh.Nvoxels))

    for i in range(Nruns):
        print("run "+str(i)+"/"+str(Nruns))
        model.syste_state = initial_conditions
        path,clock = gillespie(model,T)
        average_quantized = average_quantized + path[-1]

    average_quantized = average_quantized/Nruns

    return average_quantized

def mc_splitCoupledDiffusions(models,initial_conditions,T,runs,resolution):

    level = len(models)
    clock_quantized = np.linspace(0,T,resolution)
    Nt  = len(clock_quantized)
    path_quantized = np.zeros((Nt,len(model.systemState),model.mesh.Nvoxels))

    for i in range(levels):
        for j in range(runs[i]):
            path_average_o,clock_o = mc_crude(model,Nruns,resolution)
            path_average = path_average + path_average_o
            path_clock = path_clock + path_clock_o

    return path_quantized,clock_quantized

def quantize_path(path,clock,clock_quantized):
    path_quantized = np.zeros((len(clock_quantized),len(path[0,:,0]),len(path[0,0])))
    path_quantized[0] = path[0]
    i = 1.
    for k in range(len(clock_quantized)):
        while clock[i] < clock_quantized[k] and i+1<len(clock):
            path_quantized[k] = path[i]
            i = i+1

    return path_quantized
