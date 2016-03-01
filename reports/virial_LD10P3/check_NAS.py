
from numpy import *
import matplotlib.pyplot as plt

dat = loadtxt('../NP0400_LD10P3_C100.ener')
plt.ion()
plt.plot(dat[:,4], 'b-')
plt.grid()
plt.show()
