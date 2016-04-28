

from numpy import *
import matplotlib.pyplot as plt

dat = (loadtxt('NP1000_LD10P3_C100.ener')[:,4]/(1000.*10.))[2000:]

# plt.close()
# plt.ion()
# plt.plot(dat, 'b-')
# plt.grid()
# plt.show()
