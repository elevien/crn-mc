from pylab import *

def rho(u,v):
    return u - min(u,v)

def path_coupled(N,J,level,w,T):
    Np = 20.
    Nt = 100000000.
    t_grid = zeros(Nt)

    N0 = N
    N1 = N/(pow(J,level))

    X0 = zeros((Nt,N0))
    X1 = zeros((Nt,N1))

    Y0 = zeros((Nt,N0))
    Y1 = zeros((Nt,N1))

    x0 = zeros(N0)
    x1 = zeros(N1)
    y0 = zeros(N0)
    y1 = zeros(N1)

    while (t_grid[i] < T) and (i<Nt):
        i = i+1
        x0[:] = X0[i-1]
        x1[:] = X1[i-1]

        y0[:] = Y0[i-1]
        y1[:] = Y1[i-1]

        # compute rate sum
        # using x0,x1,y0,y1
        r = rand()
        # find next reaction index
        # using r,x0,x1,y0,y1

        # implement reaction
