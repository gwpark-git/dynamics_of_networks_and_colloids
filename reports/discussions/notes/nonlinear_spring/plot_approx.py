from numpy import *
import matplotlib.pyplot as plt

def H_W(H, x):
    return H*1./(1.-x*x)

def H_C(H, x):
    return H*(1.-(1./3.)*x*x)/(1.-x*x)

H = 1.
x_lim = 0.99
dx = x_lim/20.
x = arange(0, x_lim, dx)
xW = arange(0, x_lim, dx)
xC = arange(0.5*dx, x_lim, dx)
# x = linspace(0, x_lim, 10)
# xW = linspace(0, x_lim, 20)
# xC = linspace(dx, x_lim, 20)
plt.clf()
plt.plot(x, H*x, 'k-', label ='Gaussian (set as unity)')
plt.plot(xW, xW*H_W(H, xW), 'bo-', label='Warner approx.')
plt.plot(xC, xC*H_C(H, xC), 'ro-', label='Pade approx.')
plt.legend(loc = 'upper left', numpoints=1)
plt.xlabel('strain')
plt.ylabel('force exerted on a chain')
plt.grid('on')
plt.axis([0, 0.8, 0, 1.6])
plt.savefig('inverse_Langevin.pdf')
# plt.show()

