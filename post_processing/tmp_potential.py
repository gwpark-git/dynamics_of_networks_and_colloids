from numpy import *
import matplotlib.pyplot as plt

def U_ij(r):
    return (1./3.)*(2. + r)*(1. - r)**2. 
r = linspace(0.5, 1.5,100)

plt.plot(r, exp(-U_ij(r)), 'b-')
plt.show()
