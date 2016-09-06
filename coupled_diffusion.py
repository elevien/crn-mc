from pylab import *

def rho(u,v):
    return u - min(u,v)


def diffusion(x1,x2,r,m_max,N,w):
    m = 0

    # fine grid
    m = m + w*x2[0]
    if r < m/m_max:
        x2[0] = x2[0]-1
        x2[1] = x2[1]+1
        return x1,x2
    m = m + w*x2[1]
    if r < m/m_max:
        x2[0] = x2[0]+1
        x2[1] = x2[1]-1
        return x1,x2

    # coarse grid coupling, right jump
    m = m + w*min(x1[0],x2[1])
    if r < m/m_max:
        x2[1] = x2[1]-1
        x2[2] = x2[2]+1
        x1[0] = x1[0]-1
        x1[1] = x1[1]+1
        return x1,x2

    m = m + w*rho(x1[0],x2[1])
    if r < m/m_max:
        x1[0] = x1[0]-1
        x1[1] = x1[1]+1
        return x1,x2


    m = m + w*rho(x2[0],x1[1])
    if r < m/m_max:
        x2[1] = x2[1]-1
        x2[2] = x2[2]+1
        return x1,x2

    # right boundary

    # fine grid
    m = m + w*x2[2*N-1]
    if r < m/m_max:
        x2[2*N-1] = x2[2*N-1]-1
        x2[2*N-2] = x2[2*N-2]+1
        return x1,x2

    m = m + w*x2[2*N-2]
    if r < m/m_max:
        x2[2*N-1] = x2[2*N-1]+1
        x2[2*N-2] = x2[2*N-2]-1
        return x1,x2


    # coarse grid coupling, left jump
    m = m + w*min(x1[N-1], x2[2*N-2])
    if r < m/m_max:
        x2[2*N-2] = x2[2*N-2]-1
        x2[2*N-3] = x2[2*N-3]+1
        x1[N-1] = x1[N-1]-1
        x1[N-2] = x1[N-2]+1
        return x1,x2
    m = m + w*rho(x1[N-1], x2[2*N-2])
    if r < m/m_max:
        x1[N-1] = x1[N-1]-1
        x1[N-2] = x1[N-2]+1
        return x1,x2
    m = m + w*rho(x2[2*N-2], x1[N-1])
    if r < m/m_max:
        x2[2*N-2] = x2[2*N-2]-1
        x2[2*N-3] = x2[2*N-3]+1
        return x1,x2


    for j in range(1,N-1):


        # fine grid
        m = m + w*x2[2*j]
        if r < m/m_max:
            x2[2*j] = x2[2*j]-1
            x2[2*j+1] = x2[2*j+1]+1
            return x1,x2

        m = m + w*x2[2*j+1]
        if r < m/m_max:
            x2[2*j] = x2[2*j]+1
            x2[2*j+1] = x2[2*j+1]-1
            return x1,x2


        # coarse grid coupling, left jump
        m = m + w*min(x1[j], x2[2*j+1])
        if r < m/m_max:
            x2[2*j+1] = x2[2*j+1]-1
            x2[2*j+2] = x2[2*j+2]+1
            x1[j] = x1[j]-1
            x1[j+1] = x1[j+1]+1
            return x1,x2
        m = m + w*rho(x1[j], x2[2*j+1])
        if r < m/m_max:
            x1[j] = x1[j]-1
            x1[j+1] = x1[j+1]+1
            return x1,x2
        m = m + w*rho(x2[2*j+1], x1[j])
        if r < m/m_max:
            x2[2*j+1] = x2[2*j+1]-1
            x2[2*j+2] = x2[2*j+2]+1
            return x1,x2

        # coarse grid coupling, right jump
        m = m + w*min(x1[j], x2[2*j])
        if r < m/m_max:
            x2[2*j] = x2[2*j]-1
            x2[2*j-1] = x2[2*j-1]+1
            x1[j] = x1[j]-1
            x1[j-1] = x1[j-1]+1
            return x1,x2
        m = m + w*rho(x1[j], x2[2*j])
        if r < m/m_max:
            x1[j] = x1[j]-1
            x1[j-1] = x1[j-1]+1
            return x1,x2
        m = m + w*rho(x2[2*j], x1[j])
        if r < m/m_max:
            x2[2*j] = x2[2*j]-1
            x2[2*j-1] = x2[2*j-1]+1
            return x1,x2



def path_coupled_diffusion(N, M, T, w, Np):


    Np = 20.
    Nt = 100000000.

    N1 = N
    N2 = N*pow(2,M)

    t_grid = zeros(Nt)
    X1 = zeros((Nt,N1))
    X2 = zeros((Nt,N2))

    x1 = zeros(N1)
    x2 = zeros(N2)

    X1[0] = pow(2,M)*ones(N1)
    X2[0] = ones(N2)

    i = 0
    while (t_grid[i] < T) and (i<Nt):
        i = i+1
        x1[:] = X1[i-1]
        x2[:] = X2[i-1]

        m= 0.;

        # compute total reaction rate

        # left boundary

        # fine grid
        m = m + w*x2[0]
        m = m + w*x2[1]

        # coarse grid coupling, right jump
        m = m + w*min(x1[0],x2[1])
        m = m + w*rho(x1[0],x2[1])
        m = m + w*rho(x2[0],x1[1])

        # right boundary

        # fine grid
        m = m + w*x2[2*N1-1]
        m = m + w*x2[2*N1-2]

        # coarse grid coupling, left jump
        m = m + w*min(x1[N1-1], x2[2*N1-2])
        m = m + w*rho(x1[N1-1], x2[2*N1-2])
        m = m + w*rho(x2[2*N1-2], x1[N1-1])


        for j in range(1,N-1):


            # fine grid
            m = m + w*x2[2*j]
            m = m + w*x2[2*j+1]

            # coarse grid coupling, right jump
            m = m + w*min(x1[j],x2[2*j+1])
            m = m + w*rho(x1[j],x2[2*j+1])
            m = m + w*rho(x2[2*j+1],x1[j])

            # coarse grid coupling, left jump
            m = m + w*min(x1[j],x2[2*j])
            m = m + w*rho(x1[j],x2[2*j])
            m = m + w*rho(x2[2*j],x1[j])


        # find next reaction
        m_max = m;
        if m_max >0:
            dt = exponential(1./m_max)
            r = rand()
            x1[:],x2[:] = diffusion(x1,x2,r,m_max,N1,w)

        else:
            dt =0;
        X1[i][:] = x1;
        X2[i][:] = x2;

        t_grid[i] = t_grid[i-1]+dt
    return X1[0:i],X2[0:i],t_grid[0:i]
