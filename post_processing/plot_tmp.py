
from numpy import *
import matplotlib.pyplot as plt

dat = loadtxt('debug.traj_test')
dat0 = loadtxt('debug.traj')
# Nt = shape(dat)[0]
# ND = 2
# i=1
# LB = 10.0
# test = []
# for k in range(ND):
#     index_Rik = 2*ND*i + 1 + k
#     for t in range(1, Nt):
#         dR = dat[t, index_Rik] - dat[t-1, index_Rik]
#         if fabs(dR) > 0.9*LB:
#         # if dat[t,0] > 760.5 and dat[t,0] < 760.8:
#             test.append([t, dR])
#             print t, dat[t,0], dR

plt.plot(dat[:,0], dat[:,1], 'b-')
plt.plot(dat0[:,0], dat0[:,1], 'r-')
plt.grid('on')
plt.show()
