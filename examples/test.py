import sys
from crn_mc.mesh import *
from crn_mc.model import *
from crn_mc.simulation.paths import *
from crn_mc.simulation.montecarlo import *
from pylab import *


Nx = 1;L = 1;T = 10.;
mesh = make_lattice1d(Nx,L);systemSize = 100.;m = Model(mesh,systemSize);

m.addspecies("A",exponent=1.);m.addspecies("B",exponent=0.);
m.addspecies("C",exponent=0.);m.addspecies("O",exponent=0.);
m.addreaction([["B",1]],[["B",1],["A",1]],1.3,exponent=1.);
m.addreaction([["A",1]],[["O",1]],1.,exponent=1.);
m.addreaction([["C",1]],[["B",1]],0.5,exponent=0.);
m.addreaction([["B",1]],[["C",1]],0.3,exponent=0.);
# set initial data
ic = [1.,0.,1.,0.];
delta = 1.2;

filename_hybrid = "testpaths-hybrid.txt";filename_exact= "testpaths-exact.txt";
file_hybrid= open(filename_hybrid,"w");file_exact= open(filename_exact,"w");
for i in range(m.dimension):m.systemState[i].value[0]= ic[i];
makepath(m,T,sample_rate = 0.,path_type = 'hybrid',output_file = file_hybrid);
for i in range(m.dimension):m.systemState[i].value[0]= ic[i];
makepath(m,T,sample_rate = 0.,path_type = 'exact',output_file = file_exact);
file_hybrid.close();file_exact.close();

file_hybrid= open(filename_hybrid,"r");file_exact= open(filename_exact,"r");
results_hybrid = json.load(file_hybrid);results_exact = json.load(file_exact);
file_hybrid.close();file_exact.close();

%matplotlib inline
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

clock_hybrid = results_hybrid['results']['clock'];
clock_exact = results_exact['results']['clock'];
plt.plot(clock_hybrid,results_hybrid['results']['path']['A'],'r--',label='hybrid')
plt.plot(clock_exact,results_exact['results']['path']['A'],'k-',label='exact',alpha=0.6)
plt.legend(bbox_to_anchor=(0.95, 0.95), borderaxespad=0.,fontsize=20)
plt.legend(frameon=False)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=20)
plt.xlabel(r'$t$', fontsize=20)
plt.ylabel(r'$Z(t)$', fontsize=20)
