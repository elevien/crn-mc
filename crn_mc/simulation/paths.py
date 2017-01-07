from ..mesh import *
from ..model import *
import numpy as np
from scipy.integrate import ode
import copy


global Nt
Nt =  10e5


def gillespie1d(model,T):
    path = np.zeros((Nt,len(model.systemState),model.mesh.Nvoxels))
    clock = np.zeros(Nt)
    path[0,:] = model.systemState
    k = 1
    for e in model.events:
        e.updateRate()
    while (k<Nt) and (clock[k-1]<T):
        # compute aggregate rate
        agg_rate = sum((e.rate for e in model.events))
        delta = exponential0(agg_rate)
        # find next reaction
        r =  np.random.rand()
        firing_event = find_reaction(model.events,agg_rate,r)
        direction = firing_event.direction
        # update system state
        clock[k] = clock[k-1]+delta
        model.systemState =  model.systemState + direction
        path[k][:] = model.systemState

        # update rates
        for e in model.events:
            e.updateRate()
        k = k+1
    #print("k = "+str(k))
    return path[0:k-1],clock[0:k-1]

def chvRHS(t,y,m,sample_rate):
    m.systemState = y[0:len(m.systemState)].reshape(m.ss_d1,m.mesh.Nvoxels)
    for e in m.eventsFast:
        e.updateRate()
    agg_rate = sum((e.rate for e in m.eventsSlow))
    rhs = np.zeros(len(m.systemState)+1)
    for e in m.eventsFast:
        rhs[0:len(m.systemState)] = rhs[0:len(m.systemState)]\
         + e.direction[:,0].reshape(len(m.systemState),)*e.rate
    rhs[len(m.systemState)] = 1.
    rhs = rhs/(agg_rate+sample_rate)
    return rhs


def chv1d(model,T,h,method,sample_rate):

    # there is a bug here. Making sample rate large has problems

    path = np.zeros((Nt,len(model.systemState),model.mesh.Nvoxels))
    clock = np.zeros(Nt)
    path[0,:] = model.systemState
    k = 0
    tj = ode(chvRHS).set_integrator(method,atol = h,rtol = h)
    tj.set_f_params(model,sample_rate)

    while (k+1<Nt) and (clock[k]<T):
        k = k+1
        s1 = exponential0(1)
        # solve
        y0 = np.append(model.systemState.reshape(model.ss_d1*model.ss_d2,),0)
        tj.set_initial_value(y0,0)
        tj.integrate(s1)
        ys1 = tj.y

        model.systemState = ys1[0:len(model.systemState)].reshape(model.ss_d1,model.mesh.Nvoxels)
        t_next = tj.y[len(model.systemState)]

        for e in model.eventsSlow:
            e.updateRate()
        for e in model.eventsFast:
            e.updateRate()

        # update slow species
        r = np.random.rand()
        agg_rate = sum((e.rate for e in model.eventsSlow))
        if r>sample_rate/(agg_rate+sample_rate):
            firing_event = find_reaction(model.eventsSlow,agg_rate,r)
            direction = firing_event.direction
            model.systemState = model.systemState + direction
        clock[k] = clock[k-1] + t_next
        path[k][:] = model.systemState

    # now find the value of the continous part at exactly T
    rre = ode(rre_f).set_integrator(method,atol = h,rtol = h)
    rre.set_f_params(model)
    rre.set_initial_value(path[k-1][:].reshape(model.ss_d1*model.ss_d2,),0)
    s1 = T-clock[k-1]
    rre.integrate(s1)
    model.systemState = rre.y.reshape(model.ss_d1,model.mesh.Nvoxels)
    clock[k] = T
    path[k][:] = model.systemState

    return path[0:k+1],clock[0:k+1]

def find_reaction(events,agg_rate,r):
    s = 0.
    for e in events:
        s = s+e.rate
        if r<s/agg_rate:
            return e
