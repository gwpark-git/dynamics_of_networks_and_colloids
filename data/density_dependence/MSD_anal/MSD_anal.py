
from numpy import *
import matplotlib.pyplot as plt

dat = asarray([[0.4, 0.402, 2.487, 4.6],
               [0.6, 0.268, 3.726, 6.6],
               [0.8, 0.157, 6.385, 8.9],
               [1.0, 0.098, 10.176, 11.6],
               [1.2, 0.049, 20.264, 14.8],
               [1.4, 0.028, 35.742, 18.3]])

plt.close()
plt.figure(figsize=(6,8))
plt.ion()
ax = plt.subplot(111)
# plt.semilogy(dat[:,3], dat[:,1], 'bo-', markerfacecolor='white', markeredgecolor='blue', markersize=8)
# plt.yticks(dat[:,1])
# plt.xticks(dat[:,3])
leg_all = []
leg,= ax.semilogy(dat[:,0], dat[:,1], 'bo-', markerfacecolor='white', markeredgecolor='blue', markersize=8, label = 'tc vs. number density for micelle')
leg_all.append(leg)
ax.set_yticks(dat[:,1])
# ax.grid(axis='x')
# ax.set_xticks(dat[:,0])
ax.set_xlabel(r'number density for micelles, $\nu_m$')
# ax.set_ylabel('dimensionless diffusion coefficient (based on topological time)')
ax.set_ylabel(r'characteristic time / $\tau_0$')
# ax.grid()
for i in range(shape(dat)[0]):
    ax.semilogy([dat[i,0], dat[i,0]], [10**-2, dat[i,1]], 'b:', alpha=0.5)
ax2 = ax.twiny()
leg,=ax2.semilogy(dat[:,3], dat[:,1], 'rs-', markerfacecolor='white', markeredgecolor='red', markersize=8, label = 'tc vs. fraction of NAS')
leg_all.append(leg)
ax2.set_xticks(dat[:,3])
for i in range(shape(dat)[0]):
    ax2.semilogy([dat[i,3], dat[i,3]], [10**2, dat[i,1]], 'r:', alpha=0.5)

ax2.set_xlabel('fraction of NAS (%)')

# ax.semilogy(1.0, 0.11, 'bo')
# ax2.semilogy(10.5, 0.11, 'ro')
plt.legend(loc = 'upper right', handles=leg_all, numpoints=1)
# plt.loglog(dat[:,0], dat[:,1], 'o', markerfacecolor='white', markeredgecolor='blue', markersize=8)
# plt.loglog(dat[:2,0], dat[:2,1], 'b-', label = 'slope = %4.1f'%((log10(dat[1,1]) - log10(dat[0,1]))/(log10(dat[1,0]) - log10(dat[0,0]))))
# plt.loglog(dat[1:4,0], dat[1:4,1], 'g-', label = 'slope = %4.1f'%((log10(dat[3,1]) - log10(dat[1,1]))/(log10(dat[3,0]) - log10(dat[1,0]))))
# plt.loglog(dat[3:6,0], dat[3:6,1], 'r-', label = 'slope = %4.1f'%((log10(dat[5,1]) - log10(dat[3,1]))/(log10(dat[5,0]) - log10(dat[3,0]))))

# plt.yticks(dat[:,1])
# plt.axis([0.3, 1.5, 2, 40])
# plt.ylabel(r'dimensionless chracteristic time / $\tau_0$')
# plt.legend(loc = 'upper left')

# plt.xlabel(r'number density for micelles, $\nu_m$')
# plt.xticks(dat[:,0])
# plt.grid()
plt.show()
