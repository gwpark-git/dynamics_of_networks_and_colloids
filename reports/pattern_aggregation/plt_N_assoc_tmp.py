
from numpy import *
import matplotlib.pyplot as plt

def get_N_assoc(dat, Np):
    Nt = int(shape(dat)[0]/Np)
    re = zeros(Nt)
    for ts in range(Nt):
        re[ts] = sum(dat[ts*Np:(ts+1)*Np, 1:])/2.0
    return re
Np = 640

# NAS_G_NC05 = get_N_assoc(loadtxt('../../data_G_NC05/NP640_C100_T3.weight'), Np)
# NAS_G_NC10 = get_N_assoc(loadtxt('../../data_G_NC10/NP640_C100_T3.weight'), Np)
# NAS_G_NC15 = get_N_assoc(loadtxt('../../data_G_NC15/NP640_C100_T3.weight'), Np)
# NAS_G_NC20 = get_N_assoc(loadtxt('../../data_G_NC20/NP640_C100_T3.weight'), Np)
# NAS_G_NC25 = get_N_assoc(loadtxt('../../data_GAUSSIAN/NP640_C100_T3.weight'), Np)

# NAS_GAUSSIAN = get_N_assoc(loadtxt('data_GAUSSIAN/NP640_C100_T3.weight'), Np)
# NAS_FENE_r11 = get_N_assoc(loadtxt('data_FENE_r11/NP640_C100_T3.weight'), Np)
# NAS_FENE_r06 = get_N_assoc(loadtxt('data_FENE_r06/NP640_C100_T3.weight'), Np)

plt.clf()
plt.ion()

ax = plt.subplot(111)
Nt = size(NAS_G_NC05)

mean_NC05 = mean(NAS_G_NC05[size(NAS_G_NC05)/2:])
ref_NC05 = asarray([[0, mean_NC05],
                    [Nt, mean_NC05]])
mean_NC10 = mean(NAS_G_NC10[size(NAS_G_NC10)/2:])
ref_NC10 = asarray([[0, mean_NC10],
                    [Nt, mean_NC10]])
mean_NC15 = mean(NAS_G_NC15[size(NAS_G_NC15)/2:])
ref_NC15 = asarray([[0, mean_NC15],
                    [Nt, mean_NC15]])
mean_NC20 = mean(NAS_G_NC20[size(NAS_G_NC20)/2:])
ref_NC20 = asarray([[0, mean_NC20],
                    [Nt, mean_NC20]])
mean_NC25 = mean(NAS_G_NC25[size(NAS_G_NC25)/2:])
ref_NC25 = asarray([[0, mean_NC25],
                    [Nt, mean_NC25]])



ax.plot(NAS_G_NC05, 'b-', label = 'G, NC05, NAS=%d, f=%3.2f'%(mean_NC05, 2*mean_NC05/Np))
ax.plot(ref_NC05[:,0], ref_NC05[:,1], 'k-')
ax.plot(NAS_G_NC10, 'r-', label = 'G, NC10, NAS=%d, f=%3.2f'%(mean_NC10, 2*mean_NC10/Np))
ax.plot(ref_NC10[:,0], ref_NC10[:,1], 'k-')
ax.plot(NAS_G_NC15, 'g-', label = 'G, NC15, NAS=%d, f=%3.2f'%(mean_NC15, 2*mean_NC15/Np))
ax.plot(ref_NC15[:,0], ref_NC15[:,1], 'k-')
ax.plot(NAS_G_NC20, 'c-', label = 'G, NC20, NAS=%d, f=%3.2f'%(mean_NC20, 2*mean_NC20/Np))
ax.plot(ref_NC20[:,0], ref_NC20[:,1], 'k-')
ax.plot(NAS_G_NC25, 'k-', label = 'G, NC25, NAS=%d, f=%3.2f'%(mean_NC25, 2*mean_NC25/Np))
ax.plot(ref_NC25[:,0], ref_NC25[:,1], 'k-')

# ax.plot(NAS_GAUSSIAN, 'r.-', markerfacecolor='white', markeredgecolor='red', label = 'GAUSSIAN')
# ax.plot(NAS_FENE_r11, 'b.-', markerfacecolor='white', markeredgecolor='blue', label = 'FENE, r11')
# ax.plot(NAS_FENE_r06, 'g.-', markerfacecolor='white', markeredgecolor='green', label = 'FENE, r06')

# Nt = size(NAS_GAUSSIAN)
ax.set_xlabel('time step')
ax.set_ylabel('number of associations (NAS)')
ax.axis([0, Nt, 0, 7000])
ax.legend(loc = 'upper right', numpoints=1)
ax.grid()

ax2 = ax.twinx()
ax2.axis([0, Nt, 0, 2*7000/Np])
# ax2.plot(2*NAS_GAUSSIAN/Np, '.', markerfacecolor='none', markeredgecolor='none')
ax2.set_ylabel('f = 2*NAS/Np')
plt.show()
