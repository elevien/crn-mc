from pylab import *

def rho(u,v):
    return u - min(u,v)

def exponential0(l):
    exp_max = 100000000000000.
    if (l <= 0):
        return exp_max
    else:
        return exponential(1./l)

def get_reaction(N0,N1,k,J,r_next):
    if r_next == 0:
        # reaction: common
        v0 = -identity(N0)[k]
        v1 = -identity(N1)[k/J]
        u0 = -identity(N0)[k]
        u1 = -identity(N1)[k/J]
    elif r_next == 1:
        # reaction: only fine grid
        v0 = -identity(N0)[k]
        v1 = 0
        u0 = -identity(N0)[k]
        u1 = 0
    elif r_next == 2:
        # reaction: only coarse grid
        v0 = 0.
        v1 = -identity(N1)[k]
        u0 = 0.
        u1 = -identity(N1)[k]
    elif r_next == 3:
        # left diffusion: common
        v0 = identity(N0)[k]-identity(N0)[k+1]
        v1 = identity(N1)[k/J]-identity(N1)[k/J+1]
        u0 = 0
        u1 = 0
    elif r_next == 4:
        # left diffusion: only fine grid
        v0 = identity(N0)[k]-identity(N0)[k+1]
        v1 = 0
        u0 = 0
        u1 = 0
    elif r_next == 5:
        # left diffusion: only coarse grid
        v0 = 0
        v1 = identity(N1)[k]-identity(N1)[k+1]
        u0 = 0
        u1 = 0
    elif r_next == 6:
        # left diffusion: only fine grid
        v0 = identity(N0)[k]-identity(N0)[k+1]
        v1 = 0
        u0 = 0
        u1 = 0
    elif r_next == 7:
        # right diffusion: common
        v0 = -identity(N0)[k]+identity(N0)[k+1]
        v1 = -identity(N1)[k/J]+identity(N1)[k/J+1]
        u0 = 0
        u1 = 0
    elif r_next == 8:
        # right diffusion: only fine grid
        v0 = -identity(N0)[k]+identity(N0)[k+1]
        v1 = 0
        u0 = 0
        u1 = 0
    elif r_next == 9:
        # right diffusion: only coarse grid
        v0 = 0
        v1 = -identity(N1)[k]+identity(N1)[k+1]
        u0 = 0
        u1 = 0
    elif r_next == 10:
        # right diffusion: only fine grid
        v0 = -identity(N0)[k]+identity(N0)[k+1]
        v1 = 0
        u0 = 0
        u1 = 0
    return v0,v1,u0,u1

def path_coupled(N,J,level,w,t_max):
    Np = 20.
    Nt = 10000000.
    t_grid = zeros(Nt)

    N0 = int(N/(pow(J,level)))     # fine grid
    N1 = int(N/(pow(J,level+1))) # coarse grid
    print('N0 = ' + str(N0))
    print('N1 = ' + str(N1))

    X0 = zeros((Nt,N0))
    X1 = zeros((Nt,N1))

    Y0 = zeros((Nt,N0))
    Y1 = zeros((Nt,N1))

    x0 = ones(N0)/(pow(J,level))
    x1 = ones(N1)/(pow(J,level+1))
    y0 = ones(N0)/(pow(J,level))
    y1 = ones(N1)/(pow(J,level+1))

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
    while (t_grid[k]<t_max) and (k<Nt):
        x0[:] = X0[k-1]
        x1[:] = X1[k-1]

        y0[:] = Y0[k-1]
        y1[:] = Y1[k-1]


        for i in range(N1):
            a_sum = 0.
            for m in range(J):
                a_r_0c1_new = min(x0[i*J+m]*y0[i*J+m],x1[i]*y1[i]*pow(J,-1.))
                a_r_0m1_new = rho(x0[i*J+m]*y0[i*J+m],x1[i]*y1[i]*pow(J,-1.))

                t_r_0c1[i*J+m] = (a_r_0c1[i*J+m]/a_r_0c1_new)*(t_r_0c1[i*J+m]-t_next)
                t_r_0m1[i*J+m] = (a_r_0m1[i*J+m]/a_r_0m1_new)*(t_r_0c1[i*J+m]-t_next)

                a_r_0c1[i*J+m] = a_r_0c1_new
                a_r_0m1[i*J+m] = a_r_0m1_new

                a_sum = a_sum +x0[i*J+m]*y0[i*J+m]
            a_r_1m0_new = rho(a_sum,x1[i]*y1[i])
            t_r_1m0[i] = (a_r_1m0[i]/a_r_1m0_new)*(t_r_1m0[i]-t_next)

            a_r_0m1[i] = a_r_0m1_new
            # diffusion
            for i in range(N1-1):
                a_dx_0c1_r_new = min(w*x0[i*J],w*x1[i])
                a_dx_0m1_r_new = rho(w*x0[i*J],w*x1[i])
                a_dx_1m0_r_new = rho(w*x1[i],w*x0[i*J])
                a_dx_0c1_l_new = min(w*x0[(i+1)*J],w*x1[i+1])
                a_dx_0m1_l_new = rho(w*x0[(i+1)*J],w*x1[i+1])
                a_dx_1m0_l_new = rho(w*x1[i+1],w*x0[(i+1)*J])

                t_dx_0c1_r[i] = (a_dx_0c1_r[i]/a_dx_0c1_r_new)*(t_dx_0c1_r[i]-t_next)
                t_dx_0m1_r[i] = (a_dx_0m1_r[i]/a_dx_0m1_r_new)*(t_dx_0m1_r[i]-t_next)
                t_dx_1m0_r[i] = (a_dx_1m0_r[i]/a_dx_1m0_r_new)*(t_dx_1m0_r[i]-t_next)
                t_dx_0c1_l[i] = (a_dx_0c1_l[i]/a_dx_0c1_l_new)*(t_dx_0c1_l[i]-t_next)
                t_dx_0m1_l[i] = (a_dx_0m1_l[i]/a_dx_0m1_l_new)*(t_dx_0m1_l[i]-t_next)
                t_dx_1m0_l[i] = (a_dx_1m0_l[i]/a_dx_1m0_l_new)*(t_dx_1m0_l[i]-t_next)

                a_dx_0c1_r[i] = a_dx_0c1_r_new
                a_dx_0m1_r[i] = a_dx_0m1_r_new
                a_dx_1m0_r[i] = a_dx_1m0_r_new
                a_dx_0c1_l[i] = a_dx_0c1_l_new
                a_dx_0m1_l[i] = a_dx_0m1_l_new
                a_dx_1m0_l[i] = a_dx_1m0_l_new




                for m in range(J-1):
                    a_dx_0_r_new = w*x0[i*J+m]
                    a_dx_0_l_new = w*x0[i*J+m+1]

                    t_dx_0_r[i*(J-1)+m] = (a_dx_0_r[i*(J-1)+m]/a_dx_0_r_new)*(t_dx_0_r[i*(J-1)+m]-t_next)
                    t_dx_0_l[i*(J-1)+m] = (a_dx_0_l[i*(J-1)+m]/a_dx_0_l_new)*(t_dx_0_l[i*(J-1)+m]-t_next)

                    a_dx_0_r[i*(J-1)+m] = a_dx_0_r_new
                    a_dx_0_l[i*(J-1)+m] = a_dx_0_l_new

            # find min time and it's index

            # generate 11 reactions
            T = array([min(t_r_0c1),min(t_r_0m1),min(t_r_1m0),
             min(t_dx_0c1_l),min(t_dx_0m1_l),min(t_dx_1m0_l),min(t_dx_0_l),
             min(t_dx_0c1_r),min(t_dx_0m1_r),min(t_dx_1m0_r),min(t_dx_0_r)])

            R = array([argmin(t_r_0c1),argmin(t_r_0m1),argmin(t_r_1m0),\
             argmin(t_dx_0c1_l),argmin(t_dx_0m1_l),argmin(t_dx_1m0_l),argmin(t_dx_0_l),\
             argmin(t_dx_0c1_r),argmin(t_dx_0m1_r),argmin(t_dx_1m0_r),argmin(t_dx_0_r)])

            t_next = min(T)
            r_next = argmin(T) # which type of reaction
            l = R[argmin(T)] # location of reaction
            v0,v1,u0,u1 = get_reaction(N0,N1,l,J,r_next)

            X0[k][:] = x0 + v0
            X1[k][:] = x1 + v1
            Y0[k][:] = y0 + u0
            Y1[k][:] = y1 + u1
            t_grid[k] = t_next
            k = k+1
    return X0,X1,Y0,Y1
