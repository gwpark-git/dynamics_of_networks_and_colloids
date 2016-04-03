
# from numpy import *
# import matplotlib.pyplot as plt
# import sys
# sys.path.append('../../../post_processing')
# from acf_fcn import *

# sf_NP1728 = loadtxt('RF_NP1728_30tau0.dat')[:,1]/(2.*1728.)
# t_sf_NP1728 = arange(size(sf_NP1728))/100.
# acf_NP1728 = acf_gro(sf_NP1728)
# t_acf_NP1728 = arange(size(acf_NP1728))/100.



sf_NP1000 = loadtxt('RF_NP1000_100tau0.dat')[:size(sf_NP1728),1]/(2.*1000.)
# sf_NP1000 = loadtxt('RF_NP1000_100tau0.dat')[:,1]/(2.*1000.)
t_sf_NP1000 = arange(size(sf_NP1000))/100.
acf_NP1000 = acf_gro(sf_NP1000)
t_acf_NP1000 = arange(size(acf_NP1000))/100.




plt.close()
plt.figure(figsize=(11,6))
plt.ion()

plt.plot(t_acf_NP1000, acf_NP1000/acf_NP1000[0], 'b-', label = r'Np=1000, $\tilde{V}= 10^3$')
plt.plot(t_acf_NP1728, acf_NP1728/acf_NP1728[0], 'r-', label = r'Np=1728, $\tilde{V}= 12^3$')
plt.legend(loc = 'upper right')
plt.axis([0, 1, -0.2, 1.0])
plt.xlabel(r'dimensionless time / $\tau_0$')
plt.ylabel('virial shear stress function')
plt.grid()
plt.show()
