from model import *
import numpy as np


global Nt
Nt =  5000.


def next_reaction(model,T):
    path = np.zeros((Nt,len(model.system_state),model.mesh.Nvoxels))
    clock = np.zeros(Nt)
    path[0,:] = model.system_state
    k = 1
    while (k<Nt) and (clock[k-1]<T):
        firing_event = min(model.events, key=lambda e: e.wait_absolute)
        badrate = firing_event.wait_absolute
        m = model.events.index(firing_event)
        delta = firing_event.wait_absolute
        stoichiometric_coeffs = firing_event.stoichiometric_coeffs

        # update system
        clock[k] = clock[k-1]+delta
        model.system_state =  model.system_state + stoichiometric_coeffs

        # fire events
        firing_event.fire(delta)
        model.events.pop(m)
        for e in model.events:
            e.no_fire(delta)
        model.events.append(firing_event)
        if len(model.system_state[model.system_state  <0]) >0:  ## DB
            print("Warning: negative species count from event = " + str(firing_event))  ## DB
            break;
        path[k][:] = model.system_state
        k = k+1
    return path[0:k-1],clock[0:k-1]

def gillespie(model,T):
    path = np.zeros((Nt,len(model.system_state),model.mesh.Nvoxels))
    clock = np.zeros(Nt)
    path[0,:] = model.system_state
    k = 1
    while (k<Nt) and (clock[k-1]<T):
        # compute aggregate rate
        agg_rate = sum((e.rate for e in model.events))
        delta = exponential0(agg_rate)

        # find next reaction
        r =  np.random.rand()
        firing_event = binary_search(model.events,agg_rate,r)
        #print(firing_event)
        stoichiometric_coeffs = firing_event.stoichiometric_coeffs
        # update system state
        clock[k] = clock[k-1]+delta
        model.system_state =  model.system_state + stoichiometric_coeffs
        path[k][:] = model.system_state

        # update rates
        for e in model.events:
            e.update_rate()
        k = k+1
    return path[0:k-1],clock[0:k-1]


def binary_search(events,agg_rate,r):
    s = 0.
    for e in events:
        s = s+e.rate
        if r<s/agg_rate:
            return e

def hybrid(model,T):
    return None

def tau_leaping(model,T):
    return None
