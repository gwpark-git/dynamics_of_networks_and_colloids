
from numpy import *
import matplotlib.pyplot as plt
# rdf_NP0400_EQ = loadtxt('rdf_NP0400_EQ.dat')
# rdf_NP0800_EQ = loadtxt('rdf_NP0800_EQ.dat')
# rdf_NP1000_EQ = loadtxt('rdf_NP1000_EQ.dat')

rdf_NP0400 = loadtxt('rdf_NP0400.dat')
rdf_NP0800 = loadtxt('rdf_NP0800.dat')
rdf_NP1000 = loadtxt('rdf_NP1000.dat')

plt.close()
plt.figure(figsize=(11, 6))
plt.ion()
plt.plot(rdf_NP0400_EQ[:,0], rdf_NP0400_EQ[:,2], 'b-', linewidth=2, alpha=0.2, label = r'$\nu_m$=0.4, repulsive Brownian')
plt.plot(rdf_NP0800_EQ[:,0], rdf_NP0800_EQ[:,2], '-', color='brown', linewidth=2, alpha=0.2, label = r'$\nu_m$=0.4, repulsive Brownian')
plt.plot(rdf_NP1000_EQ[:,0], rdf_NP1000_EQ[:,2], 'r-', linewidth=2, alpha=0.2, label = r'$\nu_m$=0.4, repulsive Brownian')

plt.plot(rdf_NP0400[:,0], rdf_NP0400[:,2], 'b-', label = r'$\nu_m$=0.4, $\tau_0/\tau_B$=10')
plt.plot(rdf_NP0800[:,0], rdf_NP0800[:,2], '-', color='brown', label = r'$\nu_m$=0.8, $\tau_0/\tau_B$=10')
plt.plot(rdf_NP1000[:,0], rdf_NP1000[:,2], 'r-', label = r'$\nu_m$=1.0, $\tau_0/\tau_B$=10')
plt.plot([0, 5], [1, 1], 'k-')
plt.legend(loc = 'upper right')
plt.xlabel('dimensionless distance')
plt.ylabel('isotropic pair-correlation function')
plt.axis([0, 3, -0.1, 2.5])

plt.grid()
plt.show()
