
from numpy import *
import matplotlib.pyplot as plt

dat = loadtxt('test_2d_trajectory.dat')

plt.plot(dat[:,0], dat[:,1], 'b-')
plt.grid('on')
plt.show()
