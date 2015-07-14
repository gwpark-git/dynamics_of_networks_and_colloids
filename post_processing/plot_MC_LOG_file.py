
from numpy import *
import matplotlib.pyplot as plt
import sys
dat = loadtxt(sys.argv[1], skiprows=1)
Nt = shape(dat)[0]
cut_frac = float(sys.argv[3])
dat = dat[:long(Nt*cut_frac), :]


# plt.clf()
ax = plt.subplot(111)
leg = []
mean_Nt = shape(dat)[0]
mean_ASSO = mean(dat[mean_Nt/2:, 11])/2
leg.append(ax.plot(dat[:,0], dat[:, 12], 'b-', label='#ADD')[0])
leg.append(ax.plot(dat[:,0], dat[:, 13], 'g-', label='#MOV')[0])
leg.append(ax.plot(dat[:,0], dat[:, 14], 'y-', label='#DEL')[0])
ax2 = ax.twinx()

ref_ASSO = asarray([[min(dat[:,0]), mean_ASSO],
                    [max(dat[:,0]), mean_ASSO]])
leg.append(ax2.plot(dat[:,0], dat[:, 11]/2, 'k-', label='#CONNECTED'%(mean_ASSO))[0])
leg.append(ax2.plot(ref_ASSO[:,0], ref_ASSO[:,1], 'k--', linewidth=3, label = 'mean_#CON = %d'%(mean_ASSO))[0])
# leg.append(ax2.plot(dat[:,0], dat[:, 11]/2, 'k-', linewidth=3, label='#CONNECTED (%d)'%(mean_ASSO))[0])
leg.append(ax2.plot(dat[:,0], (mean_ASSO)*(dat[:, 15]/dat[:,0]), 'r--', linewidth= 3, label='%d*(#CANCEL/TIME)'%(mean_ASSO))[0])
# max_const = max(max(dat[:,13]), max(dat[:,12]), max(dat[:,13]), max(dat[:,14]))
# ax.axis([min(dat[:,0]), max(dat[:,0]), 0, 1.1*max_const])
ax.grid('on')
ax.set_xlabel('N_steps')
ax.set_ylabel('Countings for Add, Mov, Del')
ax2.set_ylabel('Countings for Association and Fraction of Cancel')
# ax2 = plt.twinx()
# ax2.plot(dat[:, 0]
plt.legend(handles = leg, loc = 'lower right')
if(size(sys.argv) == 5):
    plt.title(sys.argv[4])

plt.savefig(sys.argv[2])


# plt.show()
