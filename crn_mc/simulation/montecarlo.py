from ..mesh import *
from ..model import *
from .paths import *
import numpy as np
from scipy.integrate import ode
import copy



def mc_crude(model,T,Np,delta):
    Q = 0
    eps = pow(Np,-delta)
    Nruns = int(pow(eps,-2))
    for k in range(Nruns):
        path,clock= gillespie(model,T)
        Q = Q+path[-1,0]
    Q = Q/float(Nruns)
    return Q

def mc_hyrbidCoupled(model_coupled,model_hybrid,T,Np,delta,sample_rate):
    eps = pow(Np,-delta)
    h = eps
    Nruns_coupled = int(pow(eps,1/(2*delta)-2))
    Nruns_hybrid = int(pow(eps,-2))
    Q_coupled = 0
    Q_hybrid = 0
    for k in range(Nruns_coupled):
        path,clock = chv(model_coupled,T,eps,'lsoda',sample_rate)
        Q_coupled = Q_coupled-(path[-1,0]-path[-1,model_coupled.Nspecies])
    Q_coupled = Q_coupled/float(Nruns_coupled)

    #for k in range(Nruns_hybrid):
    #    path,clock = chv(model_hybrid,T,eps,'lsoda',sample_rate)
    #    Q_hybrid = Q_hybrid+path[-1,0]
    #Q_hybrid = Q_hybrid/float(Nruns_hybrid)
    return Q_coupled #+Q_hybrid

def mc_crudeDiffusions(model,T,eps,delta):

    clock_quantized = np.linspace(0,T,resolution)
    Nt  = len(clock_quantized)
    average_quantized = np.zeros((len(model.system_state),model.mesh.Nvoxels))

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
    path_quantized = np.zeros((Nt,len(model.system_state),model.mesh.Nvoxels))

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
