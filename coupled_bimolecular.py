from pylab import *

global exp_max
exp_max =  1000000000.
def rho(u,v):
    return u - min(u,v)

def exponential0(rate):
    if (rate <= 0):
        return exp_max
    else:
        return exponential(1./rate)

def get_reaction(N0,N1,l,J,r_next):
    if r_next == 0:
        # reaction: common
        v0 = -identity(N0)[l]
        v1 = -identity(N1)[l/J]
        u0 = -identity(N0)[l]
        u1 = -identity(N1)[l/J]
    elif r_next == 1:
        # reaction: only fine grid
        v0 = -identity(N0)[l]
        v1 = 0
        u0 = -identity(N0)[l]
        u1 = 0
    elif r_next == 2:
        # reaction: only coarse grid
        v0 = 0.
        v1 = -identity(N1)[l]
        u0 = 0.
        u1 = -identity(N1)[l]
    elif r_next == 3:
        # left diffusion: common
        v0 = identity(N0)[l]-identity(N0)[l+1]
        v1 = identity(N1)[l/J]-identity(N1)[l/J+1]
        u0 = 0
        u1 = 0
    elif r_next == 4:
        # left diffusion: only fine grid
        v0 = identity(N0)[l]-identity(N0)[l+1]
        v1 = 0
        u0 = 0
        u1 = 0
    elif r_next == 5:
        # left diffusion: only coarse grid
        v0 = 0
        v1 = identity(N1)[l]-identity(N1)[l+1]
        u0 = 0
        u1 = 0
    elif r_next == 6:
        # left diffusion: only fine grid
        v0 = identity(N0)[l]-identity(N0)[l+1]
        v1 = 0
        u0 = 0
        u1 = 0
    elif r_next == 7:
        # right diffusion: common
        v0 = -identity(N0)[l]+identity(N0)[l+1]
        v1 = -identity(N1)[l/J]+identity(N1)[l/J+1]
        u0 = 0
        u1 = 0
    elif r_next == 8:
        # right diffusion: only fine grid
        v0 = -identity(N0)[l]+identity(N0)[l+1]
        v1 = 0
        u0 = 0
        u1 = 0
    elif r_next == 9:
        # right diffusion: only coarse grid
        v0 = 0
        v1 = -identity(N1)[l]+identity(N1)[l+1]
        u0 = 0
        u1 = 0
    elif r_next == 10:
        # right diffusion: only fine grid
        v0 = -identity(N0)[l]+identity(N0)[l+1]
        v1 = 0
        u0 = 0
        u1 = 0
    return v0,v1,u0,u1

def path_coupled(N,J,level,w,t_max):
    Np = 100.
    Nt = 1000.
    t_grid = zeros(Nt)

    N0 = int(N/(pow(J,level)))     # fine grid
    N1 = int(N/(pow(J,level+1))) # coarse grid
    print('N0 = ' + str(N0))
    print('N1 = ' + str(N1))

    X0 = zeros((Nt,N0))
    X1 = zeros((Nt,N1))

    Y0 = zeros((Nt,N0))
    Y1 = zeros((Nt,N1))

    x0 = Np*ones(N0)*(pow(J,level))
    x1 = Np*ones(N1)*(pow(J,level+1))
    y0 = Np*ones(N0)*(pow(J,level))
    y1 = Np*ones(N1)*(pow(J,level+1))

    # only one reaction: X+Y->0
    t_r_0c1 = zeros(N0)
    t_r_0m1 = zeros(N0)
    t_r_1m0 = zeros(N1)

    a_r_0c1 = zeros(N0)
    a_r_0m1 = zeros(N0)
    a_r_1m0 = zeros(N1)

    #note: no diffusion for Y moleculers
    a_dx_0c1_l = zeros(N1-1) # left diffusion
    a_dx_0m1_l = zeros(N1-1)
    a_dx_1m0_l = zeros(N1-1)

    t_dx_0c1_l = zeros(N1-1) # left diffusion
    t_dx_0m1_l = zeros(N1-1)
    t_dx_1m0_l = zeros(N1-1)


    a_dx_0c1_r = zeros(N1-1) # right diffusions
    a_dx_0m1_r = zeros(N1-1)
    a_dx_1m0_r = zeros(N1-1)

    t_dx_0c1_r = zeros(N1-1) # right diffusions
    t_dx_0m1_r = zeros(N1-1)
    t_dx_1m0_r = zeros(N1-1)

    # these are the uncoupled diffusion channels on the fine grid
    a_dx_0_l = zeros((N1-1)*(J-1))
    a_dx_0_r = zeros((N1-1)*(J-1))

    t_dx_0_l = zeros((N1-1)*(J-1))
    t_dx_0_r = zeros((N1-1)*(J-1))


    # compute initial reaction rates and times
    # reaction kinetics

    for i in range(N1):
        a_sum = 0.
        for m in range(J):
            a_r_0c1[i*J+m] = min(x0[i*J+m]*y0[i*J+m],x1[i]*y1[i]*pow(J,-1.))
            a_r_0m1[i*J+m] = rho(x0[i*J+m]*y0[i*J+m],x1[i]*y1[i]*pow(J,-1.))
            t_r_0c1[i*J+m] = exponential0(a_r_0c1[i*J+m])
            t_r_0m1[i*J+m] = exponential0(a_r_0m1[i*J+m])
            a_sum = a_sum +x0[i*J+m]*y0[i*J+m]
        a_r_0m1[i] = rho(a_sum,x1[i]*y1[i])
        t_r_1m0[i] = exponential0(a_r_0m1[i])

    # diffusion
    for i in range(N1-1):
        a_dx_0c1_r[i] = min(w*x0[i*J],w*x1[i])
        a_dx_0m1_r[i] = rho(w*x0[i*J],w*x1[i])
        a_dx_1m0_r[i] = rho(w*x1[i],w*x0[i*J])
        a_dx_0c1_l[i] = min(w*x0[(i+1)*J],w*x1[i+1])
        a_dx_0m1_l[i] = rho(w*x0[(i+1)*J],w*x1[i+1])
        a_dx_1m0_l[i] = rho(w*x1[i+1],w*x0[(i+1)*J])

        t_dx_0c1_r[i] = exponential0(a_dx_0c1_r[i])
        t_dx_0m1_r[i] = exponential0(a_dx_0m1_r[i])
        t_dx_1m0_r[i] = exponential0(a_dx_1m0_r[i])
        t_dx_0c1_l[i] = exponential0(a_dx_0c1_l[i])
        t_dx_0m1_l[i] = exponential0(a_dx_0m1_l[i])
        t_dx_1m0_l[i] = exponential0(a_dx_1m0_l[i])

        for m in range(J-1):
            a_dx_0_r[i*(J-1)+m] = w*x0[i*J+m]
            a_dx_0_l[i*(J-1)+m] = w*x0[i*J+m+1]

            t_dx_0_r[i*(J-1)+m] = exponential0(a_dx_0_r[i*(J-1)+m])
            t_dx_0_l[i*(J-1)+m] = exponential0(a_dx_0_l[i*(J-1)+m])


    # find min time and it's index

    # generate 11 reactions
    T = array([min(t_r_0c1),min(t_r_0m1),min(t_r_1m0),
     min(t_dx_0c1_l),min(t_dx_0m1_l),min(t_dx_1m0_l),min(t_dx_0_l),
     min(t_dx_0c1_r),min(t_dx_0m1_r),min(t_dx_1m0_r),min(t_dx_0_r)])

    R = array([argmin(t_r_0c1),argmin(t_r_0m1),argmin(t_r_1m0),\
     argmin(t_dx_0c1_l),argmin(t_dx_0m1_l),argmin(t_dx_1m0_l),argmin(t_dx_0_l),\
     argmin(t_dx_0c1_r),argmin(t_dx_0m1_r),argmin(t_dx_1m0_r),argmin(t_dx_0_r)])

    #print('THE PROBLEM ='+ str(t_dx_0_l))
    print('T = '+str(T))
    print('R = '+str(R))

    t_next = min(T)
    r_next = argmin(T) # which type of reaction
    k = R[argmin(T)] # location of reaction
    v0,v1,u0,u1 = get_reaction(N0,N1,k,J,r_next)

    print('type = '+str(r_next))
    #print('v0 = '+str(v0))
    #print('v1 = '+str(v1))
    #print('u0 = '+str(u0))
    #print('u1 = '+str(v1))
    # recaluate rates

    # recaluate rates
    X0[1][:] = x0 + v0
    X1[1][:] = x1 + v1
    Y0[1][:] = y0 + u0
    Y1[1][:] = y1 + u1
    t_grid[1] = t_next

    k = 2.
    while (k<Nt):
        x0[:] = X0[k-1]
        x1[:] = X1[k-1]

        y0[:] = Y0[k-1]
        y1[:] = Y1[k-1]


        for i in range(N1):
            a_sum = 0.
            for m in range(J):
                a_r_0c1_new = min(x0[i*J+m]*y0[i*J+m],x1[i]*y1[i]*pow(J,-1.))
                a_r_0m1_new = rho(x0[i*J+m]*y0[i*J+m],x1[i]*y1[i]*pow(J,-1.))

                if a_r_0c1_new == 0:
                    t_r_0c1[i*J+m] = exp_max
                else:
                    if r_next == 0 or a_r_0c1[i*J+m]==0:
                        t_r_0c1[i*J+m] = exponential0(a_r_0c1_new)
                    else:
                        t_r_0c1[i*J+m] = (a_r_0c1[i*J+m]/a_r_0c1_new)*(t_r_0c1[i*J+m]-t_next)+t_next

                if a_r_0m1_new==0:
                    t_r_0m1[i*J+m] = exp_max
                else:
                    if r_next == 1 or a_r_0m1[i*J+m]==0:
                        t_r_0m1[i*J+m] = exponential0(a_r_0m1_new)
                    else:
                        t_r_0m1[i*J+m] = (a_r_0m1[i*J+m]/a_r_0m1_new)*(t_r_0m1[i*J+m]-t_next)+t_next

                a_r_0c1[i*J+m] = a_r_0c1_new
                a_r_0m1[i*J+m] = a_r_0m1_new

                a_sum = a_sum +x0[i*J+m]*y0[i*J+m]

            a_r_1m0_new = rho(a_sum,x1[i]*y1[i])
            if a_r_1m0_new==0:
                    t_r_1m0[i] = exp_max
            else:
                if r_next == 2 or a_r_1m0[i]==0:
                    t_r_1m0[i]  = exponential0(a_r_1m0_new)
                else:
                    t_r_1m0[i] = (a_r_1m0[i]/a_r_1m0_new)*(t_r_1m0[i]-t_next)+t_next

            a_r_0m1[i] = a_r_0m1_new

            # diffusion
            for i in range(N1-1):

                a_dx_0c1_l_new = min(w*x0[(i+1)*J],w*x1[i+1])
                a_dx_0m1_l_new = rho(w*x0[(i+1)*J],w*x1[i+1])
                a_dx_1m0_l_new = rho(w*x1[i+1],w*x0[(i+1)*J])


                if a_dx_0c1_l_new ==0:
                    t_dx_0c1_l[i] = exp_max
                else:
                    if r_next == 3 or a_dx_0c1_l[i]==0:
                        t_dx_0c1_l[i] = exponential0(a_dx_0c1_l_new)
                    else:
                        t_dx_0c1_l[i] = (a_dx_0c1_l[i]/a_dx_0c1_l_new)*(t_dx_0c1_l[i]-t_next)+t_next

                if a_dx_0m1_l_new==0:
                    t_dx_0m1_l[i] = exp_max
                else:
                    if r_next == 4 or a_dx_0m1_l[i]==0:
                        t_dx_0m1_l[i] = exponential0(a_dx_0m1_l_new)
                    else:
                        t_dx_0m1_l[i] = (a_dx_0m1_l[i]/a_dx_0m1_l_new)*(t_dx_0m1_l[i]-t_next)+t_next

                if a_dx_1m0_l_new==0 or a_dx_1m0_l[i]==0:
                    t_dx_1m0_l[i] = exp_max
                else:
                    if r_next == 5:
                        t_dx_1m0_l[i] = exponential0(a_dx_1m0_l_new)
                    else:
                        t_dx_1m0_l[i] = (a_dx_1m0_l[i]/a_dx_1m0_l_new)*(t_dx_1m0_l[i]-t_next)+t_next

                a_dx_0c1_l[i] = a_dx_0c1_l_new
                a_dx_0m1_l[i] = a_dx_0m1_l_new
                a_dx_1m0_l[i] = a_dx_1m0_l_new

                a_dx_0c1_r_new = min(w*x0[i*J],w*x1[i])
                a_dx_0m1_r_new = rho(w*x0[i*J],w*x1[i])
                a_dx_1m0_r_new = rho(w*x1[i],w*x0[i*J])

                if a_dx_0c1_r_new==0:
                    t_dx_0c1_r[i] = exp_max
                else:
                    if r_next == 7 or a_dx_0c1_r[i]==0:
                        t_dx_0c1_r[i] = exponential0(a_dx_0c1_r_new)
                    else:
                        t_dx_0c1_r[i] = (a_dx_0c1_r[i]/a_dx_0c1_r_new)*(t_dx_0c1_r[i]-t_next)+t_next

                if a_dx_0m1_r_new ==0:
                    t_dx_0m1_r[i] = exp_max
                else:
                    if r_next == 8 or a_dx_0m1_r[i]==0:
                        t_dx_0m1_r[i] = exponential0(a_dx_0m1_r_new)
                    else:
                        t_dx_0m1_r[i] = (a_dx_0m1_r[i]/a_dx_0m1_r_new)*(t_dx_0m1_r[i]-t_next)+t_next

                if a_dx_1m0_r_new ==0:
                    t_dx_1m0_r[i] = exp_max
                else:
                    if r_next == 9 or a_dx_1m0_r[i]==0:
                        t_dx_1m0_r[i] = exponential0(a_dx_1m0_r_new)
                    else:
                        t_dx_1m0_r[i] = (a_dx_1m0_r[i]/a_dx_1m0_r_new)*(t_dx_1m0_r[i]-t_next)+t_next

                a_dx_0c1_r[i] = a_dx_0c1_r_new
                a_dx_0m1_r[i] = a_dx_0m1_r_new
                a_dx_1m0_r[i] = a_dx_1m0_r_new

                for m in range(J-1):

                    a_dx_0_l_new = w*x0[i*J+m+1]
                    a_dx_0_r_new = w*x0[i*J+m]

                    if a_dx_0_l_new ==0:
                        t_dx_0_l[i*(J-1)+m] = exp_max
                    else:
                        if r_next == 6 or a_dx_0_l[i*(J-1)+m]==0:
                            t_dx_0_l[i*(J-1)+m] = exponential0(a_dx_0_l_new)
                        else:
                            t_dx_0_l[i*(J-1)+m] =(a_dx_0_l[i*(J-1)+m]/a_dx_0_l_new)*(t_dx_0_l[i*(J-1)+m]-t_next)+t_next


                    if a_dx_0_r_new ==0:
                        t_dx_0_r[i*(J-1)+m] = exp_max
                    else:
                        if r_next == 10 or a_dx_0_r[i*(J-1)+m]==0:
                            t_dx_0_r[i*(J-1)+m] = exponential0(a_dx_0_r_new)
                        else:
                            t_dx_0_r[i*(J-1)+m] =(a_dx_0_r[i*(J-1)+m]/a_dx_0_r_new)*(t_dx_0_r[i*(J-1)+m]-t_next)+t_next



                    a_dx_0_l[i*(J-1)+m] = a_dx_0_l_new
                    a_dx_0_r[i*(J-1)+m] = a_dx_0_r_new


        # find min time and it's index

        # generate 11 reactions
        T = array([min(t_r_0c1),min(t_r_0m1),min(t_r_1m0),
         min(t_dx_0c1_l),min(t_dx_0m1_l),min(t_dx_1m0_l),min(t_dx_0_l),
         min(t_dx_0c1_r),min(t_dx_0m1_r),min(t_dx_1m0_r),min(t_dx_0_r)])
        R = array([argmin(t_r_0c1),argmin(t_r_0m1),argmin(t_r_1m0),\
         argmin(t_dx_0c1_l),argmin(t_dx_0m1_l),argmin(t_dx_1m0_l),argmin(t_dx_0_l),\
         argmin(t_dx_0c1_r),argmin(t_dx_0m1_r),argmin(t_dx_1m0_r),argmin(t_dx_0_r)])
        #print("T = "+str(T))
        #print("R = "+str(R))
        t_next = min(T)
        r_next = argmin(T) # which type of reaction
        #print("T = "+str(T))
        #print("time     = "+str(t_next))
        #print("reaction = "+str(r_next))
        l = R[argmin(T)]  # spatial location of reaction
        v0,v1,u0,u1 = get_reaction(N0,N1,l,J,r_next)

        X0[k][:] = x0 + v0
        X1[k][:] = x1 + v1
        Y0[k][:] = y0 + u0
        Y1[k][:] = y1 + u1
        t_grid[k] = t_grid[k-1]+t_next
        #print("t_grid[k] = "+str(t_grid[k]))
        k = k+1

    return X0,X1,Y0,Y1,t_grid

def path_uncoupled(N,w,t_max):
    Np = 100.
    Nt = 1000.
    t_grid = zeros(Nt)

    X0 = zeros((Nt,N))
    Y0 = zeros((Nt,N))

    x0 = Np*ones(N)
    y0 = Np*ones(N)

    events = []

    # reactions
    for i in range(N):
        rate = x0[i]*y0[i]
        stoichiometric_coeffs =(-identity(N)[i],-identity(N)[i])
        system_state = (x0,y0)
        reaction = Reaction(N,i,2,stoichiometric_coeffs,system_state)
        events.append(reaction)

    # diffusions
    for i in range(N-1):

        # left diffusion
        left_diffusion = Diffusion(N,i,i+1,0,x0)
        events.append(left_diffusion)

        # right diffusion
        right_diffusion = Diffusion(N,i+1,i,0,x0)
        events.append(right_diffusion)

    t_grid[0] = 0.
    k = 1

    while k < Nt:

        x0[:] = X0[k-1]
        y0[:] = Y0[k-1]
        system_state=(x0,y0)

        # find event with minimum absolute time
        firing_event = min(events, key=lambda e: e.wait_absolute)
        m = events.index(firing_event)
        # get information needed from event
        delta = firing_event.wait_absolute
        stoichiometric_coeffs = firing_event.stoichiometric_coeffs

        # update events
        firing_event.fire(system_state,delta)
        # update all other events
        for e in events:
            e.no_fire(system_state,delta)


        # update system
        t_grid[i] = t_grid[i-1]+delta
        x0[k][:] = x0 + stoichiometric_coeffs[0]
        y0[k][:] = y0 + stoichiometric_coeffs[1]

        k = k+1
    return X0,Y0,t_grid









class Event:

    def __init__(self,system_state):
        self.time_internal = 0.
        self.wait_internal = exponential0(1.)
        self.update_rate(system_state)
        self.wait_absolute = (self.wait_internal-self.time_internal)/self.rate



    def fire(self,system_state,delta):
        self.update_rate(system_state)
        self.time_internal = self.time_internal + self.rate*delta
        self.wait_internal = exponential0(1.)
        self.wait_absolute = (self.wait_internal-self.time_internal)/self.rate
        return None

    def no_fire(self,system_state,delta):
        self.update_rate(system_state)
        self.time_internal = self.time_internal + self.rate*delta
        # wait_internal remains unchanged
        self.wait_absolute = (self.wait_internal-self.time_internal)/self.rate
        return None


    def update_rate(self,system_state):
        self.rate = 0.
        return None


class Diffusion(Event):
    def __init__(self,mesh,voxel_out,voxel_in,species,system_state):
        self.voxel_out = voxel_out
        self.voxel_in = voxel_in
        self.stoichiometric_coeffs = -identity(mesh)[voxel_out]+identity(mesh)[voxel_in]
        super().__init__(system_state)

    def update_rate(self,system_state):
        self.rate = system_state[self.voxel_out]
        return None

class Reaction(Event):
    def __init__(self,mesh,voxel,order,stoichiometric_coeffs,system_state):
        self.voxel = voxel
        self.order = order
        self.stoichiometric_coeffs = stoichiometric_coeffs
        super().__init__(system_state)

    def update_rate(self,system_state):
        # works for order =1,2
        rate = 1.
        for i in range(self.order):
            rate = rate*system_state[i][self.voxel]
        self.rate = rate
        return None

class Diffusion_SplitComman(Event):
    pass

class Diffusion_SplitFine(Event):
    pass

class Diffusion_SplitCoarse(Event):
    pass


class Reaction_SplitComman(Event):
    pass

class Reaction_SplitFine(Event):
    pass

class Reaction_SplitCoarse(Event):
    pass
