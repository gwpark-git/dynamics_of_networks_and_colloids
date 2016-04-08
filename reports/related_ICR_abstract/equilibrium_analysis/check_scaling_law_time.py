
from numpy import *
import matplotlib.pyplot as plt


dat = asarray([[400, 0.051, 0.162, 0.776],
               [600, 0.077, 0.182, 0.327],
               [700, 0.061, 0.263, 0.199],
               [800, 0.080, 0.183, 0.087]])

plt.clf()
plt.ion()
plt.figure(figsize=(6, 8))
plt.loglog(dat[:,0], dat[:,1], 'bo-', label = 'tc between 0 and 0.01')
plt.loglog(dat[:,0], dat[:,2], 'go-', label = 'tc between 0.1 and 0.2')
plt.loglog(dat[:,0], dat[:,3], 'ro-', label = 'tc between 0.2 and 0.3')
plt.grid()
plt.xticks(dat[:,0], ['%d'% float(val) for val in dat[:,0]])
plt.yticks(dat[:,1:].flatten(), ['%3.2e'% val for val in dat[:,1:].flatten()])
plt.axis([380, 820, 0.05, 0.8])
plt.xlabel('number of particles')
plt.ylabel('initial characteristic time')
plt.legend(loc = 'upper right', numpoints=1)
plt.show()
