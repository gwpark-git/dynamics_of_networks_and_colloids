
from numpy import *
import matplotlib.pyplot as plt


sf_NP1000 = loadtxt('RF_NP1000_100tau0.dat')[:,1]/(2.*1000.)
t_sf_NP1000 = arange(size(sf_NP1000))/100.

sf_NP1728 = loadtxt('RF_NP1728_30tau0.dat')[:,1]/(2.*1728.)
t_sf_NP1728 = arange(size(sf_NP1728))/100.


plt.close()
plt.figure(figsize=(11,6))
plt.ion()

plt.plot(t_sf_NP1000, sf_NP1000, 'b-', label = r'Np=1000, $\tilde{V}= 10^3$')
plt.plot(t_sf_NP1728, sf_NP1728, 'r-', alpha=0.5, label = r'Np=1728, $\tilde{V}= 12^3$')
plt.legend(loc = 'upper right')
plt.axis([0, 40, -0.2, 0.2])
plt.xlabel(r'dimensionless time / $\tau_0$')
plt.ylabel('virial shear stress function')
plt.grid()
plt.show()
