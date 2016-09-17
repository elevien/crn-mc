from model import *
import numpy as np

def next_reaction(model,T):
    Nt = 500;
    path = np.zeros((Nt,len(model.system_state),model.mesh.size))
    clock = np.zeros(Nt)
    path[0,:] = model.system_state
    k = 1
    while (k<Nt):
        firing_event = min(model.events, key=lambda e: e.wait_absolute)
        m = model.events.index(firing_event)
        delta = firing_event.wait_absolute
        stoichiometric_coeffs = firing_event.stoichiometric_coeffs
        firing_event.fire(delta)
        model.events.pop(m)
        for e in model.events:
            e.no_fire(delta)
        model.events.append(firing_event)
        # update system
        clock[k] = clock[k-1]+delta

        model.system_state =  model.system_state + stoichiometric_coeffs
        if np.mean(model.system_state[model.system_state  <= 0]) <0:
            print("Warning: negative species count from event = " + str(firing_event))
            print(model.system_state[firing_event.species][firing_event.voxel_out])
            print(np.argmin(model.system_state[firing_event.species]))
        path[k][:] = model.system_state
        k = k+1
    return path,clock
