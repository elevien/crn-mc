from ..mesh import *
from ..model import *
import inspect
import numpy as np
from scipy.integrate import ode
import copy


global Nt
Nt =  10e5


def tryexponential(rate):
    """ Trys to compute exponential. """
    try:
        return np.random.exponential(1./rate)
    except ValueError:
        print("next jump time is at infinity")

def chvrhs(t,y,model,sample_rate):
    for i in range(model.dimension):
        model.systemState[i].value[0] = y[i]
    for e in model.events:
        e.updaterate()
    slow = filter(lambda e: e.speed == SLOW, model.events)
    agg_rate = 0.
    for s in slow:
        agg_rate = agg_rate + s.rate

    rhs = np.zeros(model.dimension+1)
    fast = filter(lambda e: e.speed == FAST, model.events)
    for e in fast:
        for i in range(model.dimension):
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


def makepath(model,T,h,method='lsoda',sample_rate = 0.,*args,**kwargs):
    """ Compute paths of hybrid model using CHV method. """
    voxel = 0.
    path = np.zeros((Nt,len(model.systemState)))
    path[0][:] = model.getstate(0)
    clock = np.zeros(Nt)

    k = 0
    tj = ode(chvrhs).set_integrator(method,atol = h,rtol = h)
    tj.set_f_params(model,sample_rate)

    while (k+1<Nt) and (clock[k]<T):
        k = k+1
        s1 = tryexponential(1)
        # solve
        y0 = np.append(model.getstate(0),0)
        tj.set_initial_value(y0,0)
        tj.integrate(s1)
        ys1 = tj.y

        for i in range(model.dimension):
            model.systemState[i].value[0] = ys1[i]
        t_next = tj.y[model.dimension]

        for e in model.events:
            e.updaterate()
        # update slow species
        r = np.random.rand()
        stochastic_events = getstochasticevents(model)
        agg_rate = 0.
        for e in stochastic_events:
            agg_rate = agg_rate + e.rate
        if r>sample_rate/(agg_rate+sample_rate):
            firing_event = findreaction(model.events,agg_rate,r)
            firing_event.react()
        clock[k] = clock[k-1] + t_next
        path[k][:] = model.getstate(0)
    return path[0:k+1],clock[0:k+1]

def getstochasticevents(model):
    stochastic_events = []
    for e in model.events:
        if e.scale == 0.:
            stochastic_events.append(e)
    return stochastic_events



def chvrhs_coupled(t,y,model_hybrid,model_exact,sample_rate):
    for i in range(model_exact.dimension):
        model_hybrid.systemState[i].value[0] = y[i]
    for i in range(model_hybrid.dimension):
        model_exact.systemState[i].value[0] = y[i+model_exact.dimension]
    for e in model_exact.events:
        e.updaterate()
    for e in model_hybrid.events:
        e.updaterate()
    agg_rate = 0.
    for i in range(len(model_exact.events)):
        rate_hybrid = model_hybrid.events[i].rate
        rate_exact = model_exact.events[i].rate
        agg_rate = agg_rate + rate_hybrid + rate_exact - min(rate_hybrid,rate_exact)

    rhs = np.zeros(2*model_exact.dimension+1)
    fast = filter(lambda e: e.speed == FAST, model_hybrid.events)

    for e in fast:
        for i in range(model_exact.dimension):
            name = model_exact.systemState[i].name
            r = list(filter(lambda e: e[0].name == name, e.reactants))
            p = list(filter(lambda e: e[0].name == name, e.products))
            direction = 0.
            if r:
                direction = direction - float(r[0][1])
            if p:
                direction = direction + float(p[0][1])
            rhs[i] = rhs[i] + direction*e.rate
    rhs[2*model_exact.dimension] = 1.

    rhs = rhs/(agg_rate+sample_rate)
    return rhs

def res(x,y):
    return x - min(x,y)


def makepath_coupled(model_hybrid,T,h,method='lsoda',sample_rate = 0.,*args,**kwargs):
    """ Compute paths of coupled exact-hybrid model using CHV method. """
    voxel = 0
    # make copy of model with exact dynamics
    model_exact = copy.deepcopy(model_hybrid)
    for e in model_exact.events:
        e.speed = SLOW

    # setup integrator
    path = np.zeros((Nt,2*model_hybrid.dimension))
    path[0][0:model_hybrid.dimension] = model_hybrid.getstate(0)
    path[0][model_hybrid.dimension:2*model_hybrid.dimension] = model_exact.getstate(0)
    clock = np.zeros(Nt)

    k = 0
    tj = ode(chvrhs_coupled).set_integrator(method,atol = h,rtol = h)
    tj.set_f_params(model_hybrid,model_exact,sample_rate)
    y0 = np.zeros(2*model_hybrid.dimension+1)

    while (k+1<Nt) and (clock[k]<T):
        k = k+1
        s1 = tryexponential(1)
        # solve
        y0[0:model_hybrid.dimension] = model_hybrid.getstate(0)
        y0[model_hybrid.dimension:2*model_hybrid.dimension] = model_exact.getstate(0)
        y0[2*model_hybrid.dimension] = 0.
        tj.set_initial_value(y0,0)
        tj.integrate(s1)
        ys1 = tj.y

        for i in range(model_hybrid.dimension):
            model_hybrid.systemState[i].value[0] = ys1[i]
        for i in range(model_hybrid.dimension):
            model_exact.systemState[i].value[0] = ys1[i+model_hybrid.dimension]
        t_next = tj.y[2*model_hybrid.dimension]


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
                firing_event_hybrid.react()
            if isinstance(firing_event_exact,Reaction):
                firing_event_exact.react()
        clock[k] = clock[k-1] + t_next
        path[k][0:model_hybrid.dimension] = model_hybrid.getstate(0)
        path[k][model_hybrid.dimension:2*model_hybrid.dimension] = model_exact.getstate(0)
    return path[0:k+1],clock[0:k+1]

def findreaction(events,agg_rate,r):
    rate_sum = 0.
    for e in events:
        if e.speed == SLOW:
            rate_sum = rate_sum +e.rate
            if r<rate_sum/agg_rate:
                return e

null = NullEvent()
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
    return null,null
