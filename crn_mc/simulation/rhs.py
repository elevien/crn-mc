from ..mesh import *
from ..model import *
from .timer import *
import copy,json
import numpy as np
from scipy.integrate import ode


# Right hand sides --------------------------------------------------------
# curretly spending too much time inside this function. perhaps don't
# use filter?

def chvrhs_hybrid(t,y,model,sample_rate):
    for i in range(model.dimension):
        model.systemState[i].value[0] = y[i]
    for e in model.events:
        e.updaterate()

    #MIXED = filter(lambda e: e.hybridType == MIXED, model.events)
    agg_rate = 0.
    for i in range(model.dimension):
        if model.events[i].hybridType == SLOW or model.events[i].hybridType == MIXED:
            agg_rate = agg_rate + model.events[i].rate
    #for s in MIXED:
    #    agg_rate = agg_rate + s.rate

    rhs = np.zeros(model.dimension+1)
    fast = filter(lambda e: e.hybridType == FAST, model.events)
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
    for i in range(len(model_hybrid.events)):
        if model_hybrid.events[i].hybridType == SLOW or model_hybrid.events[i].hybridType == MIXED:
            hybrid_rate = model_hybrid.events[i].rate
            exact_rate = model_exact.events[i].rate
            agg_rate = agg_rate + res(hybrid_rate,exact_rate )
            agg_rate = agg_rate + res(exact_rate,hybrid_rate )
            agg_rate = agg_rate + min(hybrid_rate,exact_rate )
        elif model_hybrid.events[i].hybridType == FAST or model_hybrid.events[i].hybridType == VITL:
            agg_rate = agg_rate + model_exact.events[i].rate


    rhs = np.zeros(2*model_exact.dimension+1)
    fast = filter(lambda e: e.hybridType == FAST, model_hybrid.events)
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

def rrerhs(t,y,model,sample_rate):
    """rhs of determistic part of equations, i.e. the rhs of reaction rate equations"""
    for i in range(model.dimension):
        model.systemState[i].value[0] = y[i]
    for e in model.events:
        e.updaterate()
    rhs = np.zeros(model.dimension)
    fast = filter(lambda e: e.hybridType == FAST, model.events)
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
    return rhs
