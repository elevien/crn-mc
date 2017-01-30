from ..mesh import *
from ..model import *
from .timer import *
from .rhs import *
import copy,json
import numpy as np
from scipy.integrate import ode



# should move this to arguments
global Nt
Nt =  10e5


# Helper functions --------------------------------------------------------
def tryexponential(rate):
    """ Trys to compute exponential. """
    try:
        return np.random.exponential(1./rate)
    except ValueError:
        print("next jump time is at infinity")

def getstochasticevents(model):
    stochastic_events = []
    for e in model.events:
        if e.hybridType != FAST:
            stochastic_events.append(e)
    return stochastic_events

def getstochasticevents_hybrid(model):
    stochastic_events = []
    for e in model.events:
        if e.hybridType == SLOW or e.hybridType == MIXED:
            stochastic_events.append(e)
    return stochastic_events

def findreaction_gillespie(events,agg_rate,r):
    rate_sum = 0.
    for e in events:
        rate_sum = rate_sum + e.rate
        if r<rate_sum/agg_rate:
            return e

def findreaction_hybrid(events,agg_rate,r):
    rate_sum = 0.
    for e in events:
        if e.hybridType == SLOW or e.hybridType == MIXED:
            rate_sum = rate_sum +e.rate
            if r<rate_sum/agg_rate:
                return e

null = NullEvent()
def findreaction_coupled(events_hybrid,events_exact,agg_rate,r):
    rate_sum = 0.
    for i in range(len(events_hybrid)):
        if events_hybrid[i].hybridType == SLOW or events_hybrid[i].hybridType == MIXED:
            exact_rate = events_exact[i].rate
            hybrid_rate  = events_hybrid[i].rate
            rate_sum = rate_sum + res(hybrid_rate,exact_rate)
            if r<rate_sum/agg_rate:
                return events_hybrid[i],null
            rate_sum = rate_sum + res(exact_rate,hybrid_rate)
            if r<rate_sum/agg_rate:
                return null,events_exact[i]
            rate_sum = rate_sum + min(hybrid_rate,exact_rate)
            if r<rate_sum/agg_rate:
                return events_hybrid[i],events_exact[i]
        elif events_hybrid[i].hybridType == FAST or events_hybrid[i].hybridType == VITL:
            exact_rate = events_exact[i].rate
            rate_sum = rate_sum + exact_rate
            if r<rate_sum/agg_rate:
                return null,events_exact[i]
        #else:
        #    print("PROBLEM")
    return null,null


# path generation ---------------------------------------------------------

def makepath(model,T,h = None,ode_method='lsoda',sample_rate = 0.,
                        path_type = 'hybrid',output_file=None,*args,**kwargs):
    if h == None:
        h = 1./model.systemSize
    cpu_time = 0.
    if path_type == 'hybrid':
        with timer() as t:
            path, clock = makepath_hybrid(model,T,h,ode_method,sample_rate)
        cpu_time = t.secs
        # the reason for computing path entry here is that it differs for esimators
        path_entry = {model.systemState[i].name:list(path[:,i])
            for i in range(model.dimension)}
    elif path_type == 'exact':
        with timer() as t:
            path, clock = makepath_exact(model,T)
        cpu_time = t.secs
        path_entry = {model.systemState[i].name:list(path[:,i])
            for i in range(model.dimension)}
    elif path_type == 'coupled':
        with timer() as t:
            path, clock = makepath_coupled(model,T,h,ode_method,sample_rate)
        cpu_time = t.secs
        # here this is a but different
        path_entry_hybrid = {model.systemState[i].name:list(path[:,i])
            for i in range(model.dimension)}
        path_entry_exact = {model.systemState[i].name+'*':
            list(path[:,i+model.dimension])for i in range(model.dimension)}
        path_entry = dict(path_entry_hybrid, **path_entry_exact)
    if output_file != None:
        params_dict = {'path_type':path_type,'sample_rate':sample_rate,'ode_method':ode_method,
                    'T':T,'h':h}
        model_info = {'system_size':model.systemSize,
            'events':list([e.__str__() for e in model.events])}
        results_dict = {
            'path':path_entry,'clock':list(clock),'cpu_time':cpu_time}
        output = {'params':params_dict,'model':model_info,'results':results_dict}
        print(json.dumps(output),file = output_file)
    return path,clock


def makepath_exact(model,T):
    """ Compute exact path using Gillespie algorithm. """
    voxel = 0.
    for e in model.events:
        e.hybridType = SLOW
        e.updaterate()
    path = np.zeros((Nt,len(model.systemState)))
    path[0][:] = model.getstate(0)
    clock = np.zeros(Nt)
    k = 0
    while (k+1<Nt) and (clock[k]<T):
        k = k+1
        for e in model.events:
            e.updaterate()
        r = np.random.rand()
        agg_rate = 0.
        for e in model.events:
            agg_rate = agg_rate + e.rate
        t_next = tryexponential(agg_rate)
        firing_event = findreaction_gillespie(model.events,agg_rate,r)
        firing_event.react()
        clock[k] = clock[k-1] + t_next
        path[k][:] = model.getstate(0)
    return path[0:k+1],clock[0:k+1]

def makepath_hybrid(model,T,h,ode_method,sample_rate):
    """ Compute paths of model. """
    voxel = 0.
    path = np.zeros((Nt,len(model.systemState)))
    path[0][:] = model.getstate(0)
    clock = np.zeros(Nt)
    for e in model.events:
        e.sethybridtype()
        e.updaterate()

    # for hybrid paths use chv ode_method
    k = 0
    tj = ode(chvrhs_hybrid).set_integrator(ode_method,atol = h,rtol = h)
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
        stochastic_events = getstochasticevents_hybrid(model)
        agg_rate = 0.
        for e in stochastic_events:
            agg_rate = agg_rate + e.rate
        if r>sample_rate/(agg_rate+sample_rate):
            firing_event = findreaction_hybrid(model.events,agg_rate,r)
            firing_event.react()
        clock[k] = clock[k-1] + t_next
        path[k][:] = model.getstate(0)
    # now compute value at exactly T
    tj = ode(rrerhs).set_integrator(ode_method,atol = h,rtol = h)
    tj.set_f_params(model,sample_rate)
    y0 = path[k-1][:]
    tj.set_initial_value(y0,clock[k-1])
    tj.integrate(T)
    yT = tj.y
    for i in range(model.dimension):
        model.systemState[i].value[0] = yT[i]
    clock[k] = T
    path[k][:] = model.getstate(0)
    return path[0:k+1],clock[0:k+1]


def makepath_coupled(model_hybrid,T,h,ode_method,sample_rate):
    """ Compute paths of coupled exact-hybrid model using CHV ode_method. """
    voxel = 0
    # make copy of model with exact dynamics
    model_exact = copy.deepcopy(model_hybrid)
    for e in model_hybrid.events:
        e.sethybridtype()
        e.updaterate()
    for e in model_exact.events:
        e.hybridType = SLOW
        e.updaterate()

    # setup integrator
    path = np.zeros((Nt,2*model_hybrid.dimension))
    path[0][0:model_hybrid.dimension] = model_hybrid.getstate(0)
    path[0][model_hybrid.dimension:2*model_hybrid.dimension] = model_exact.getstate(0)
    clock = np.zeros(Nt)

    k = 0
    tj = ode(chvrhs_coupled).set_integrator(ode_method,atol = h,rtol = h)
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
            if model_hybrid.events[i].hybridType == SLOW or model_hybrid.events[i].hybridType == MIXED:
                hybrid_rate = model_hybrid.events[i].rate
                exact_rate = model_exact.events[i].rate
                agg_rate = agg_rate + res(hybrid_rate,exact_rate )
                agg_rate = agg_rate + res(exact_rate,hybrid_rate )
                agg_rate = agg_rate + min(hybrid_rate,exact_rate )
            elif model_hybrid.events[i].hybridType == FAST or model_hybrid.events[i].hybridType == VITL:
                agg_rate = agg_rate + model_exact.events[i].rate

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
