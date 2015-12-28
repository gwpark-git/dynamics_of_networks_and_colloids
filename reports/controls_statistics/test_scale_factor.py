


from numpy import *
import matplotlib.pyplot as plt

dat = asarray([[1.0, 6089],
               [1.25, 3384],
               [1.5, 1721],
               [1.75, 741],
               [2.0, 289]])


plt.clf()
plt.ion()
ax = plt.subplot(111)
ax.plot(dat[:,0], dat[:,1], 'bo', markerfacecolor='white', markeredgecolor='blue')
ax.axis([1, 2, 0, 7000])
ax.set_xlabel('scale factor: r=scale_factor*given_r')
ax.set_ylabel('Number of association with the first step')
ax.grid()
ax2 = ax.twinx()
ax2.plot(dat[:,0], dat[:,1]/(25*640), 'b-')
ax2.axis([1, 2, 0./(25.*640.), 7000./(25.*640.)])
ax2.set_ylabel('fraction of NAS to the all associable chain')

plt.show()
