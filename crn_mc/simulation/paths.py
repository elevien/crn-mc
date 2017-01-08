from ..mesh import *
from ..model import *
import inspect
import numpy as np
from scipy.integrate import ode
import copy


global Nt
Nt =  10e5


def gillespie(model,T,voxel):
    path = np.zeros((Nt,len(model.systemState)))
    path[0][:] = model.getstate(0)
    clock = np.zeros(Nt)

    k = 1
    for e in model.events:
        e.updaterate()
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
        path[k][:] = model.getstate(0)

        # update rates
        for e in model.events:
            e.updaterate()
        k = k+1
        agg_rate = sum((e.rate for e in model.events))
    #print("k = "+str(k))
    return path[0:k-1],clock[0:k-1]

def chvrhs(t,y,model,sample_rate):
    for i in range(model.Nspecies):
        model.systemState[i].value[0] = y[i]
    for e in model.events:
        e.updaterate()
    slow = filter(lambda e: e.speed == SLOW, model.events)
    agg_rate = 0.
    for s in slow:
        agg_rate = agg_rate + s.rate

    rhs = np.zeros(model.Nspecies+1)
    fast = filter(lambda e: e.speed == FAST, model.events)
    for e in fast:
        for i in range(model.Nspecies):
            name = model.systemState[i].name
            r = list(filter(lambda e: e[0].name == name, e.reactants))
            p = list(filter(lambda e: e[0].name == name, e.products))
            direction = 0.
            if r:
                direction = direction - float(r[0][1])
            if p:
                direction = direction + float(p[0][1])
            rhs[i] = rhs[i]+ direction*e.rate
    rhs[len(model.systemState)] = 1.

    rhs = rhs/(agg_rate+sample_rate)
    return rhs



def makepath(model,T,h,method,sample_rate,voxel):

    # there is a bug here. Making sample rate large has problems

    path = np.zeros((Nt,len(model.systemState)))
    path[0][:] = model.getstate(0)
    clock = np.zeros(Nt)

    k = 0
    tj = ode(chvrhs).set_integrator(method,atol = h,rtol = h)
    tj.set_f_params(model,sample_rate)

    while (k+1<Nt) and (clock[k]<T):
        k = k+1
        s1 = exponential0(1)
        # solve
        y0 = np.append(model.getstate(0),0)
        tj.set_initial_value(y0,0)
        tj.integrate(s1)
        ys1 = tj.y

        for i in range(model.Nspecies):
            model.systemState[i].value[0] = ys1[i]
        t_next = tj.y[model.Nspecies]

        for e in model.events:
            e.updaterate()
        # update slow species
        r = np.random.rand()
        slow_events = filter(lambda e: e.speed == SLOW, model.events)
        agg_rate = 0.
        for s in slow_events:
            agg_rate = agg_rate + s.rate
        if r>sample_rate/(agg_rate+sample_rate):
            firing_event = findreaction(model.events,agg_rate,r)
            print(firing_event)
            model.react(firing_event)
        clock[k] = clock[k-1] + t_next
        path[k][:] = model.getstate(0)
    return path[0:k+1],clock[0:k+1]

def chvrhs_coupled(t,y,model_hybrid,model_exact,sample_rate):
    for i in range(model_exact.Nspecies):
        model_hybrid.systemState[i].value[0] = y[i]
    for i in range(model_hybrid.Nspecies):
        model_exact.systemState[i].value[0] = y[i+model_exact.Nspecies]
    for e in model_exact.events:
        e.updaterate()
    for e in model_hybrid.events:
        e.updaterate()
    agg_rate = 0.
    for i in range(len(model_exact.events)):
        rate_hybrid = model_hybrid.events[i].rate
        rate_exact = model_exact.events[i].rate
        agg_rate = agg_rate + rate_hybrid + rate_exact - min(rate_hybrid,rate_exact)

    rhs = np.zeros(2*model_exact.Nspecies+1)
    fast = filter(lambda e: e.speed == FAST, model_hybrid.events)

    for e in fast:
        for i in range(model_exact.Nspecies):
            name = model_exact.systemState[i].name
            r = list(filter(lambda e: e[0].name == name, e.reactants))
            p = list(filter(lambda e: e[0].name == name, e.products))
            direction = 0.
            if r:
                direction = direction - float(r[0][1])
            if p:
                direction = direction + float(p[0][1])
            rhs[i] = rhs[i]+ direction*e.rate
    rhs[2*model_exact.Nspecies] = 1.

    rhs = rhs/(agg_rate+sample_rate)
    return rhs

def res(x,y):
    return x - min(x,y)


def makepath_coupled(model_hybrid,T,h,method,sample_rate,voxel):

    # make copy of model with exact dynamics
    model_exact = copy.deepcopy(model_hybrid)
    for e in model_exact.events:
        e.speed = SLOW

    # setup integrator
    path = np.zeros((Nt,2*model_hybrid.Nspecies))
    path[0][0:model_hybrid.Nspecies] = model_hybrid.getstate(0)
    path[0][model_hybrid.Nspecies:2*model_hybrid.Nspecies] = model_exact.getstate(0)
    clock = np.zeros(Nt)

    k = 0
    tj = ode(chvrhs_coupled).set_integrator(method,atol = h,rtol = h)
    tj.set_f_params(model_hybrid,model_exact,sample_rate)
    y0 = np.zeros(2*model_hybrid.Nspecies+1)

    while (k+1<Nt) and (clock[k]<T):
        k = k+1
        s1 = exponential0(1)
        # solve
        y0[0:model_hybrid.Nspecies] = model_hybrid.getstate(0)
        y0[model_hybrid.Nspecies:2*model_hybrid.Nspecies] = model_exact.getstate(0)
        y0[2*model_hybrid.Nspecies] = 0.
        tj.set_initial_value(y0,0)
        tj.integrate(s1)
        ys1 = tj.y

        for i in range(model_hybrid.Nspecies):
            model_hybrid.systemState[i].value[0] = ys1[i]
        for i in range(model_hybrid.Nspecies):
            model_exact.systemState[i].value[0] = ys1[i+model_hybrid.Nspecies]
        t_next = tj.y[2*model_hybrid.Nspecies]


        for e in model_hybrid.events:
            e.updaterate()
        for e in model_exact.events:
            e.updaterate()

        # update slow species
        r = np.random.rand()
        agg_rate = 0.
        for i in range(len(model_hybrid.events)):
            if model_hybrid.events[i].speed == SLOW:
                hybrid_rate = model_hybrid.events[i].rate
                exact_rate = model_exact.events[i].rate
                agg_rate = agg_rate + res(hybrid_rate,exact_rate )
                agg_rate = agg_rate + res(exact_rate,hybrid_rate )
                agg_rate = agg_rate + min(hybrid_rate,exact_rate )
            elif model_hybrid.events[i].speed == FAST:
                agg_rate = agg_rate + model_exact.events[i].rate
            else:
                print("PROBLEM")


        # find reaction
        if r>sample_rate/(agg_rate+sample_rate):
            firing_event_hybrid,firing_event_exact  = findreaction_coupled(model_hybrid.events,model_exact.events,agg_rate,r)

            if isinstance(firing_event_hybrid,Reaction):
                model_hybrid.react(firing_event_hybrid)
            if isinstance(firing_event_exact,Reaction):
                model_exact.react(firing_event_exact)
        clock[k] = clock[k-1] + t_next
        path[k][0:model_hybrid.Nspecies] = model_hybrid.getstate(0)
        path[k][model_hybrid.Nspecies:2*model_hybrid.Nspecies] = model_exact.getstate(0)
    return path[0:k+1],clock[0:k+1]

def findreaction(events,agg_rate,r):
    rate_sum = 0.
    for e in events:
        if e.speed == SLOW:
            rate_sum = rate_sum +e.rate
            if r<rate_sum/agg_rate:
                return e

null = NullReaction()
def findreaction_coupled(events_hybrid,events_exact,agg_rate,r):
    rate_sum = 0.
    for i in range(len(events_hybrid)):
        if events_hybrid[i].speed == SLOW:
            exact_rate = events_exact[i].rate
            hybrid_rate  = events_hybrid[i].rate
            rate_sum = rate_sum + res(hybrid_rate,exact_rate)
            if r<rate_sum/agg_rate:
                return events_hybrid[i],null
            rate_sum = rate_sum + res(exact_rate,exact_rate)
            if r<rate_sum/agg_rate:
                return null,events_exact[i]
            rate_sum = rate_sum + min(hybrid_rate,exact_rate)
            if r<rate_sum/agg_rate:
                return events_hybrid[i],events_exact[i]
        elif events_hybrid[i].speed == FAST:
            exact_rate = events_exact[i].rate
            rate_sum = rate_sum + exact_rate
            if r<rate_sum/agg_rate:
                return null,events_exact[i]
        else:
            print("PROBLEM")

    print(rate_sum/agg_rate)
    return null,null
