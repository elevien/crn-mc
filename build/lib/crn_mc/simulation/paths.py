from ..mesh import *
from ..model import *
import numpy as np
from scipy.integrate import ode
import copy


global Nt
Nt =  500000.


def next_reaction(model,T):
    path = np.zeros((Nt,len(model.systemState),model.mesh.Nvoxels))
    clock = np.zeros(Nt)
    path[0,:] = model.systemState
    k = 1
    while (k<Nt) and (clock[k-1]<T):
        firing_event = min(model.events, key=lambda e: e.wait_absolute)
        badrate = firing_event.wait_absolute
        m = model.events.index(firing_event)
        delta = firing_event.wait_absolute
        direction = firing_event.direction

        # update system
        clock[k] = clock[k-1]+delta
        model.systemState =  model.systemState + direction

        # fire events
        firing_event.fire(delta)
        model.events.pop(m)
        for e in model.events:
            e.no_fire(delta)
        model.events.append(firing_event)
        if len(model.systemState[model.systemState  <0]) >0:  ## DB
            print("Warning: negative species count from event = " + str(firing_event))  ## DB
            break;
        path[k][:] = model.systemState
        k = k+1
    return path[0:k-1],clock[0:k-1]

def gillespie(model,T):
    path = np.zeros((Nt,len(model.systemState),model.mesh.Nvoxels))
    clock = np.zeros(Nt)
    path[0,:] = model.systemState
    k = 1
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

def rre_f(t,y,m):
    # make copy of model
    m.systemState = y.reshape(m.ss_d1,m.mesh.Nvoxels)
    for e in m.eventsFast:
        e.updateRate()
    rates = np.zeros(len(m.systemState))
    for e in m.eventsFast:
        rates = rates + e.direction[:,0].reshape(len(m.systemState),)*e.rate
    #print(rates.tolist())
    return rates

def chv_f(t,y,m,sample_rate):
    m.systemState = y[0:len(m.systemState)].reshape(m.ss_d1,m.mesh.Nvoxels)
    for e in m.eventsFast:
        e.updateRate()
    agg_rate = sum((e.rate for e in m.eventsSlow))+sample_rate
    rhs = np.zeros(len(m.systemState)+1)
    for e in m.eventsFast:
        rhs[0:len(m.systemState)] = rhs[0:len(m.systemState)]\
         + e.direction[:,0].reshape(len(m.systemState),)*e.rate
    rhs[len(m.systemState)] = 1.
    rhs = rhs/agg_rate
    return rhs


def chv(model,T,h,method,sample_rate):
    print('here')
    path = np.zeros((Nt,len(model.systemState),model.mesh.Nvoxels))
    clock = np.zeros(Nt)
    path[0,:] = model.systemState
    k = 0
    tj = ode(chv_f).set_integrator(method,atol = h,rtol = h)
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


def strang_split(model,T,h0,h,method):
    clock = np.arange(0,T,h0)
    path = np.zeros((len(clock),len(model.systemState),model.mesh.Nvoxels))
    path[0,:] = model.systemState

    # setup ODE integrator
    rre = ode(rre_f).set_integrator(method,atol = h,rtol = h)
    rre.set_f_params(model)

    for k in range(len(clock)):
        tY = clock[k]
        # gillespie 1/2 step
        while tY<clock[k]+h0/2.:
            agg_rate = sum((e.rate for e in model.eventsSlow))
            delta = exponential0(agg_rate)
            if delta<h0/2.:
                tY = tY+delta
                # find next reaction
                r = np.random.rand()
                firing_event = find_reaction(model.eventsSlow,agg_rate,r)
                direction = firing_event.direction
                # fire slow reaction and update system state
                model.systemState = model.systemState + direction
                for e in model.eventsFast:
                    e.updateRate()
                for e in model.eventsSlow:
                    e.updateRate()

        # integrate 1 step
        rre.set_initial_value(model.systemState,0)
        rre.integrate(h0)
        model.systemState = rre.y
        for e in model.eventsFast:
            e.updateRate()
        for e in model.eventsSlow:
            e.updateRate()

        # gillespie 1/2 step
        tY = clock[k]+h0/2.
        while tY<clock[k]+h0:
            agg_rate = sum((e.rate for e in model.eventsSlow))
            delta = exponential0(agg_rate)

            if delta<h0/2.:
                tY = tY+delta
                # find next reaction
                r =  np.random.rand()
                firing_event = find_reaction(model.eventsSlow,agg_rate,r)
                direction = firing_event.direction
                # fire slow reaction and update system state

                model.systemState =  model.systemState + direction
                for e in model.eventsFast:
                    e.updateRate()
                for e in model.eventsSlow:
                    e.updateRate()

        # store path
        path[k][:] = model.systemState
    return path,clock

def gillespie_hybrid(model,T,h1,h2,method):
    path = np.zeros((Nt,len(model.systemState),model.mesh.Nvoxels))
    clock = np.zeros(Nt)
    path[0,:] = model.systemState
    k = 1
    rre = ode(rre_f).set_integrator(method,atol = h1,rtol = h1)
    rre.set_f_params(model)
    while (k<Nt) and (clock[k-1]<T):
        # compute aggregate rate
        agg_rate = sum((e.rate for e in model.eventsSlow))
        delta = exponential0(agg_rate)
        if delta<h2:
            # find next reaction
            r =  np.random.rand()
            firing_event = find_reaction(model.eventsSlow,agg_rate,r)
            direction = firing_event.direction

            # fire slow reaction and update system state
            clock[k] = clock[k-1]+delta
            model.systemState =  model.systemState + direction
            path[k][:] = model.systemState

            # integrate
            rre.set_initial_value(model.systemState,clock[k])
            rre.integrate(rre.t+delta)
            model.systemState = rre.y
            path[k][:] = model.systemState

        else:
            # integrate
            rre.set_initial_value(model.systemState,clock[k])
            rre.integrate(rre.t+h2)
            clock[k] = clock[k-1]+h2
            model.systemState = rre.y
            path[k][:] = model.systemState

        # update rates
        for e in model.eventsFast:
            e.updateRate()
        for e in model.eventsSlow:
            e.updateRate()
        k = k+1
    #print("k = "+str(k))
    return path[0:k-1],clock[0:k-1]

def find_reaction(events,agg_rate,r):
    s = 0.
    for e in events:
        s = s+e.rate
        if r<s/agg_rate:
            return e



def tau_leaping(model,T):
    return None
