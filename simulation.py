from model import *

def next_reaction(model,T):
    Nt = 100;
    while (Nt<T):
        ss[:]= model.system_state
        firing_event = min(model.events, key=lambda e: e.wait_absolute)
        m = model.events.index(firing_event)
        delta = firing_event.wait_absolute
        stoichiometric_coeffs = firing_event.stoichiometric_coeffs
        firing_event.fire(system_state,delta)
        model.events.pop(m)
        for e in events:
            e.no_fire(model.system_state,delta)
        model.events.append(firing_event)
        # update system
        clock[k] = clock[k-1]+delta
       system_state[k][:] = ss + stoichiometric_coeffs
       k = k+1
   return system_state,clock
