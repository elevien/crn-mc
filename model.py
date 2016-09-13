class Model:
    def __init__(self,species,mesh,system_state):
        self.mesh = mesh
        self.species = species
        self.system_state
        self.events = []

    def add_reaction(self,stoichiometric_coeffs,voxel_in,voxel_out):
        reaction = Reaction_SplitFine(self,stoichiometric_coeffs)
        self.events.append(reaction)
        Return None

    def init_reactions(self,stoichiometric_coeffs,voxel_in,voxel_out):
        Return None

    def init_diffusions(self,stoichiometric_coeffs,voxel_in,voxel_out):
        for voxel in mesh:
            for neighbor in voxel.neighbor
                diffusion =Diffusion(self,mesh,voxel,neighbor,system_state)
                self.events.append(reaction)
        Return None

    def next_reaction(self,T):
        Nt = 100;
        while (Nt<T):
            ss[:]= system_state
            firing_event = min(events, key=lambda e: e.wait_absolute)
            m = events.index(firing_event)
            delta = firing_event.wait_absolute
            stoichiometric_coeffs = firing_event.stoichiometric_coeffs
            firing_event.fire(system_state,delta)
            events.pop(m)
            for e in events:
                e.no_fire(system_state,delta)
            events.append(firing_event)
            # update system
            clock[k] = clock[k-1]+delta
            system_state[k][:] = ss + stoichiometric_coeffs
            k = k+1
        return system_state,clock
