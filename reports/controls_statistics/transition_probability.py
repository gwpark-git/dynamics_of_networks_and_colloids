
from numpy import *
import matplotlib.pyplot as plt

def UG(r):
    return r**2.0

def Boltzmann_map(r):
    return exp(-UG(r))

def Boltzmann_factor(r, Eb):
    return exp(-(UG(r) - Eb))

def Boltzmann_factor_cut(r, Eb):
    BF = Boltzmann_factor(r, Eb)
    for i,p in enumerate(BF):
        if p > 1.0:
            BF[i] = 1.0
    return BF
    

r = linspace(0, 5, 80)
Eb = 0.5
BM = Boltzmann_map(r)
BF = Boltzmann_factor(r, Eb)
BFC = Boltzmann_factor_cut(r, Eb)
ref_unity = asarray([[0, 1.0],
                     [5, 1.0]])
plt.clf()
plt.ion()
plt.plot(r, BM, 'b-', label = 'Normalized probability')
plt.plot(r, BF, 'r-', label = 'without Normalization')
plt.plot(r, BFC, 'g-', label = 'cut-off probability')
plt.plot(ref_unity[:,0], ref_unity[:,1], 'k--')
plt.grid()

plt.xlabel('test distance')
plt.ylabel('probability')
plt.legend(loc = 'upper right')
plt.show()
