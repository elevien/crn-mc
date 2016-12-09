from model import *
import numpy as np
from scipy.integrate import ode
import copy


global Nt
Nt =  500000.


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
        stoichiometric_coeffs = firing_event.stoichiometric_coeffs
        # update system state
        clock[k] = clock[k-1]+delta
        model.system_state =  model.system_state + stoichiometric_coeffs
        path[k][:] = model.system_state

        # update rates
        for e in model.events:
            e.update_rate()
        k = k+1
    #print("k = "+str(k))
    return path[0:k-1],clock[0:k-1]

def rre_f(t,y,m):
    # make copy of model
    m.system_state = y
    rates = np.zeros(len(m.system_state))
    for e in m.events_fast:
        rates = rates + e.stoichiometric_coeffs[:,0]*e.rate
    #print(rates.tolist())
    return rates.tolist()

def chv_f(t,y,m,sample_rate):
    m.system_state = y[0:len(m.system_state)].reshape(m.ss_d1,m.mesh.Nvoxels)
    for e in m.events_fast:
        e.update_rate()
    agg_rate = sum((e.rate for e in m.events_slow))+sample_rate
    rhs = np.zeros(len(m.system_state)+1)
    for e in m.events_fast:
        rhs[0:len(m.system_state)] = rhs[0:len(m.system_state)]\
         + e.stoichiometric_coeffs[:,0].reshape(len(m.system_state),)*e.rate
    rhs[len(m.system_state)] = 1.
    rhs = rhs/agg_rate
    return rhs


def chv(model,T,h,method,sample_rate):
    path = np.zeros((Nt,len(model.system_state),model.mesh.Nvoxels))
    clock = np.zeros(Nt)
    path[0,:] = model.system_state
    k = 0
    tj = ode(chv_f).set_integrator(method,atol = h,rtol = h)
    tj.set_f_params(model,sample_rate)

    while (k<Nt) and (clock[k-1]<T):
        k = k+1
        s1 = exponential0(1)
        # solve
        y0 = np.append(model.system_state.reshape(model.ss_d1*model.ss_d2,),0)
        tj.set_initial_value(y0,0)
        tj.integrate(s1)
        ys1 = tj.y

        model.system_state = ys1[0:len(model.system_state)].reshape(model.ss_d1,model.mesh.Nvoxels)
        t_next = tj.y[len(model.system_state)]

        for e in model.events_slow:
            e.update_rate()
        for e in model.events_fast:
            e.update_rate()

        # update slow species
        r = np.random.rand()
        agg_rate = sum((e.rate for e in model.events_slow))
        if r>sample_rate/(agg_rate+sample_rate):
            firing_event = binary_search(model.events_slow,agg_rate,r)
            stoichiometric_coeffs = firing_event.stoichiometric_coeffs
            model.system_state = model.system_state + stoichiometric_coeffs
        clock[k] = clock[k-1] + t_next
        path[k][:] = model.system_state

    return path[0:k-1],clock[0:k-1]


def strang_split(model,T,h0,h,method):
    clock = np.arange(0,T,h0)
    path = np.zeros((len(clock),len(model.system_state),model.mesh.Nvoxels))
    path[0,:] = model.system_state

    # setup ODE integrator
    rre = ode(rre_f).set_integrator(method,atol = h,rtol = h)
    rre.set_f_params(model)

    for k in range(len(clock)):
        tY = clock[k]
        # gillespie 1/2 step
        while tY<clock[k]+h0/2.:
            agg_rate = sum((e.rate for e in model.events_slow))
            delta = exponential0(agg_rate)
            if delta<h0/2.:
                tY = tY+delta
                # find next reaction
                r = np.random.rand()
                firing_event = binary_search(model.events_slow,agg_rate,r)
                stoichiometric_coeffs = firing_event.stoichiometric_coeffs
                # fire slow reaction and update system state
                model.system_state = model.system_state + stoichiometric_coeffs
                for e in model.events_fast:
                    e.update_rate()
                for e in model.events_slow:
                    e.update_rate()

        # integrate 1 step
        rre.set_initial_value(model.system_state,0)
        rre.integrate(h0)
        model.system_state = rre.y
        for e in model.events_fast:
            e.update_rate()
        for e in model.events_slow:
            e.update_rate()

        # gillespie 1/2 step
        tY = clock[k]+h0/2.
        while tY<clock[k]+h0:
            agg_rate = sum((e.rate for e in model.events_slow))
            delta = exponential0(agg_rate)

            if delta<h0/2.:
                tY = tY+delta
                # find next reaction
                r =  np.random.rand()
                firing_event = binary_search(model.events_slow,agg_rate,r)
                stoichiometric_coeffs = firing_event.stoichiometric_coeffs
                # fire slow reaction and update system state

                model.system_state =  model.system_state + stoichiometric_coeffs
                for e in model.events_fast:
                    e.update_rate()
                for e in model.events_slow:
                    e.update_rate()

        # store path
        path[k][:] = model.system_state
    return path,clock

def gillespie_hybrid(model,T,h1,h2,method):
    path = np.zeros((Nt,len(model.system_state),model.mesh.Nvoxels))
    clock = np.zeros(Nt)
    path[0,:] = model.system_state
    k = 1
    rre = ode(rre_f).set_integrator(method,atol = h1,rtol = h1)
    rre.set_f_params(model)
    while (k<Nt) and (clock[k-1]<T):
        # compute aggregate rate
        agg_rate = sum((e.rate for e in model.events_slow))
        delta = exponential0(agg_rate)
        if delta<h2:
            # find next reaction
            r =  np.random.rand()
            firing_event = binary_search(model.events_slow,agg_rate,r)
            stoichiometric_coeffs = firing_event.stoichiometric_coeffs

            # fire slow reaction and update system state
            clock[k] = clock[k-1]+delta
            model.system_state =  model.system_state + stoichiometric_coeffs
            path[k][:] = model.system_state

            # integrate
            rre.set_initial_value(model.system_state,clock[k])
            rre.integrate(rre.t+delta)
            model.system_state = rre.y
            path[k][:] = model.system_state

        else:
            # integrate
            rre.set_initial_value(model.system_state,clock[k])
            rre.integrate(rre.t+h2)
            clock[k] = clock[k-1]+h2
            model.system_state = rre.y
            path[k][:] = model.system_state

        # update rates
        for e in model.events_fast:
            e.update_rate()
        for e in model.events_slow:
            e.update_rate()
        k = k+1
    #print("k = "+str(k))
    return path[0:k-1],clock[0:k-1]

def binary_search(events,agg_rate,r):
    s = 0.
    for e in events:
        s = s+e.rate
        if r<s/agg_rate:
            return e



def tau_leaping(model,T):
    return None



def mc_crude(model,T,Np,delta):
    Q = 0
    eps = pow(Np,-delta)
    Nruns = int(pow(eps,-2))
    for k in range(Nruns):
        path,clock= gillespie(model,T)
        Q = Q+path[-1,0]
    Q = Q/float(Nruns)
    return Q

def mc_hyrbidCoupled(model_coupled,Qz,T,Np,delta,h):
    eps = pow(Np,-delta)
    Nruns_coupled = int(pow(eps,1/delta-2))
    Nruns_hybrid = int(pow(eps,-2))
    Q_coupled = 0
    Q_hybrid = 0
    for k in range(Nruns_coupled):
        path,clock= gillespie_hybrid(model_coupled,T,eps,h,'lsoda')
        Q_coupled = Q_coupled-(path[-1,0]-path[-1,model_coupled.Nspecies])
    Q_coupled = Q_coupled/float(Nruns_coupled)

    # for k in range(Nruns_hybrid):
    #     path,clock= gillespie_hybrid(model_hybrid,T,eps,h,'dop853')
    #     Q_hybrid = Q_hybrid+path[-1,0]
    Q_hybrid = Q_hybrid/float(Nruns_hybrid)
    return Q_coupled+Qz*Np

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
