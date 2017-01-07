from ..mesh import *
from ..model import *
import numpy as np
from scipy.integrate import ode
import copy


global Nt
Nt =  10e5


def gillespie(model,T,voxel):
    path = np.zeros((Nt,len(model.systemState)))
    path[0][:] = model.getStateInVoxel(0)
    clock = np.zeros(Nt)

    k = 1
    for e in model.events:
        e.updateRate()
    agg_rate = sum((e.rate for e in model.events))
    while (k<Nt) and (clock[k-1]<T) and (agg_rate >0):
        # compute aggregate rate
        delta = exponential0(agg_rate)

        # find next reaction
        r =  np.random.rand()
        firing_event = find_reaction(model.events,agg_rate,r)
        # update system state

        clock[k] = clock[k-1]+delta
        model.react(firing_event)
        path[k][:] = model.getStateInVoxel(0)

        # update rates
        for e in model.events:
            e.updateRate()
        k = k+1
        agg_rate = sum((e.rate for e in model.events))
    #print("k = "+str(k))
    return path[0:k-1],clock[0:k-1]

def chvRHS(t,y,m,sample_rate):
    for i in range(m.Nspecies):
        m.systemState[i].value[0] = y[i]
    for e in m.events:
        e.updateRate()
    slow = filter(lambda e: e.speed == SLOW, m.events)
    agg_rate = 0.
    for s in slow:
        agg_rate = agg_rate + s.rate

    rhs = np.zeros(m.Nspecies+1)
    fast = filter(lambda e: e.speed == FAST, m.events)
    for e in fast:
        for i in range(m.Nspecies):
            name = m.systemState[i].name
            r = list(filter(lambda e: e[0].name == name, e.reactants))
            p = list(filter(lambda e: e[0].name == name, e.products))
            direction = 0.
            if r:
                direction = direction - float(r[0][1])
            if p:
                direction = direction + float(p[0][1])
            rhs[i] = rhs[i]+ direction*e.rate
    rhs[len(m.systemState)] = 1.
    rhs = rhs/(agg_rate+sample_rate)
    return rhs


def chv1d(model,T,h,method,sample_rate,voxel):

    # there is a bug here. Making sample rate large has problems

    path = np.zeros((Nt,len(model.systemState)))
    path[0][:] = model.getStateInVoxel(0)
    clock = np.zeros(Nt)

    k = 0
    tj = ode(chvRHS).set_integrator(method,atol = h,rtol = h)
    tj.set_f_params(model,sample_rate)

    while (k+1<Nt) and (clock[k]<T):
        k = k+1
        s1 = exponential0(1)
        # solve
        y0 = np.append(model.getStateInVoxel(0),0)
        tj.set_initial_value(y0,0)
        tj.integrate(s1)
        ys1 = tj.y

        for i in range(model.Nspecies):
            model.systemState[i].value[0] = ys1[i]
        t_next = tj.y[model.Nspecies]

        for e in model.events:
            e.updateRate()
            print(e)
            print(e.rate)

        # update slow species
        r = np.random.rand()
        slow = filter(lambda e: e.speed == SLOW, model.events)
        agg_rate = 0.
        for s in slow:
            agg_rate = agg_rate + s.rate
        if r>sample_rate/(agg_rate+sample_rate):
            firing_event = find_reaction(model.events,agg_rate,r)
            #print(firing_event)
            model.react(firing_event)
        clock[k] = clock[k-1] + t_next
        path[k][:] = model.getStateInVoxel(0)
        #print(path[k][:])
        #print(clock[k])

    # need to find the value of the continous part at exactly T

    #clock[k] = T
    #path[k][:] = model.getStateInVoxel(0)

    return path[0:k+1],clock[0:k+1]

def find_reaction(events,agg_rate,r):
    s = 0.
    for e in events:
        if e.speed == SLOW:
            s = s+e.rate
            if r<s/agg_rate:
                return e
