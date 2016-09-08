from pylab import *

def rho(u,v):
    return u - min(u,v)

def path_coupled(N,J,level,w,T):
    Np = 20.
    Nt = 100000000.
    t_grid = zeros(Nt)

    N0 = N     # fine grid
    N1 = N/(pow(J,level)) # coarse grid

    X0 = zeros((Nt,N0))
    X1 = zeros((Nt,N1))

    Y0 = zeros((Nt,N0))
    Y1 = zeros((Nt,N1))

    x0 = zeros(N0)
    x1 = zeros(N1)
    y0 = zeros(N0)
    y1 = zeros(N1)


    # store absolute times for each reaction channel
    # reaction kinetics
    # there is ony one reaction: X+Y->0

    T_rx_0c1 = zeros(N0)
    T_rx_0m1 = zeros(N0)
    T_rx_1m0 = zeros(N1)

    # no diffusion for Y moleculers
    T_dx_0c1_l = zeros() # left diffusion
    T_dx_0c1_r = zeros() # right diffusions
    T_dx_0m1_l = zeros() #
    T_dx_1m0_l = zeros() #

    # compute reaction rates
    for i in range(N0):
        T_rx_0c1[i] = min(x0[i]*)
        T_rx_0c1[i] = min(x0[i]*)
        T_rx_0c1[i] = min(x0[i]*)
    for i in range(N1):

    # find min time and it's index
    m1 = min(T_rx_0c1)

    #update state of system
    x0[:]
    x0[:]
    x0[:]
    x0[:]

    while (t_grid[i] < T) and (i<Nt):
        i = i+1
        x0[:] = X0[i-1]
        x1[:] = X1[i-1]

        y0[:] = Y0[i-1]
        y1[:] = Y1[i-1]



        for i in range()
        r = rand()
        # find next reaction index
        # using r,x0,x1,y0,y1

        # implement reaction
