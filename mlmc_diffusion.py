import matplotlib.pyplot as plt
from coupled_diffusion import *

N = 10
M = 1
w = 1.
Np = 100
T = 50.


S = 100
K = 5

v_DIFF = zeros(K)
v_CRUDE = zeros(K)
N_range = arange(10,10+K)

for j in range(K):
    DIFF = zeros(S)
    CRUDE = zeros(S)
    for i in range(S):
        N = N_range[j]
        X1,X2,t_grid = path_coupled_diffusion(N, M, T, w, Np)
        DIFF[i] = abs(X1[-1,4]-(X2[-1,8]+X2[-1,9]))
        CRUDE[i] = X1[-1,4]
    v_DIFF[j] = var(DIFF)
    v_CRUDE[j] = var(CRUDE)

plt.plot(N_range,v_DIFF,'r-')
plt.plot(N_range,v_CRUDE,'k-')
#plt.plot(t_grid,X1[:,4],'r-')
#plt.plot(t_grid,X2[:,8]+X2[:,9],'k-')
plt.show()
