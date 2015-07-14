
from numpy import *
import matplotlib.pyplot as plt
import sys
from matplotlib.font_manager import FontProperties
dat = loadtxt(sys.argv[1], skiprows=1)
Nt = shape(dat)[0]
cut_frac = float(sys.argv[3])
dat = dat[:long(Nt*cut_frac), :]

fontP = FontProperties()
fontP.set_size('x-small')
N_blocks = int(sys.argv[4])
block_dat = range(N_blocks, long(Nt*cut_frac), N_blocks)
N_bd = size(block_dat)
measure_block_MSE = zeros(N_bd)
for i in range(N_bd):
    measure_block_MSE[i] = sum((dat[N_blocks*i:N_blocks*(i+1), 11]/2 - mean(dat[N_blocks*i:N_blocks*(i+1), 11])/2)**2.)/N_bd

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
leg.append(ax2.plot(dat[:,0], (dat[:, 11]/2)/mean_ASSO, 'k-', label='#CONNECTED/mean_eq_#CONNECTED')[0])
leg.append(ax2.plot(ref_ASSO[:,0], ref_ASSO[:,1]/mean_ASSO, 'k--', linewidth=3, label = 'mean_#CON = %d (scaled by 1)'%(mean_ASSO))[0])
leg.append(ax2.plot(block_dat, measure_block_MSE/max(measure_block_MSE), 'b--', linewidth=3, label='bMSE/max(bMSE)')[0])
# leg.append(ax2.plot(dat[:,0], dat[:, 11]/2, 'k-', linewidth=3, label='#CONNECTED (%d)'%(mean_ASSO))[0])
leg.append(ax2.plot(dat[:,0], (dat[:, 15]/dat[:,0]), 'r--', linewidth= 3, label='#CANCEL/TIME')[0])
ax2.axis([0, long(Nt*cut_frac), -0.3, 1.3])
# max_const = max(max(dat[:,13]), max(dat[:,12]), max(dat[:,13]), max(dat[:,14]))
# ax.axis([min(dat[:,0]), max(dat[:,0]), 0, 1.1*max_const])
ax.grid('on')
ax.set_xlabel('N_steps')
ax.set_ylabel('Countings for Add, Mov, Del')
ax2.set_ylabel('Scaled Measure for ASSO, CANCEL, bMSE')
# ax2 = plt.twinx()
# ax2.plot(dat[:, 0]
plt.legend(handles = leg, loc = 'center right', prop=fontP, ncol=2)
if(size(sys.argv) == 6):
    plt.title(sys.argv[5] + '_#mean_eq(ASSO)=922')

plt.savefig(sys.argv[2])


# plt.show()
