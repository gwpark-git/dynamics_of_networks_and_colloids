
from numpy import *
import matplotlib.pyplot as plt

def get_N_assoc(dat, Np):
    Nt = int(shape(dat)[0]/Np)
    re = zeros(Nt)
    for ts in range(Nt):
        re[ts] = sum(dat[ts*Np:(ts+1)*Np, 1:])/2.0
    return re
Np = 640

dat_NR = loadtxt('data_NR/NP640_C100_T3.weight')
dat_ME_EB0_000 = loadtxt('data_ME_EB0_000/NP640_C100_T3.weight')
dat_ME_EB0_100 = loadtxt('data_ME_EB0_100/NP640_C100_T3.weight')
dat_ME_EB0_200 = loadtxt('data_ME_EB0_200/NP640_C100_T3.weight')
dat_ME_EB0_300 = loadtxt('data_ME_EB0_300/NP640_C100_T3.weight')
dat_ME_EB0_400 = loadtxt('data_ME_EB0_400/NP640_C100_T3.weight')
dat_ME_EB0_500 = loadtxt('data_ME_EB0_500/NP640_C100_T3.weight')
dat_ME_EB0_600 = loadtxt('data_ME_EB0_600/NP640_C100_T3.weight')
dat_ME_EB0_700 = loadtxt('data_ME_EB0_700/NP640_C100_T3.weight')
dat_ME_EB0_800 = loadtxt('data_ME_EB0_800/NP640_C100_T3.weight')
dat_ME_EB0_900 = loadtxt('data_ME_EB0_900/NP640_C100_T3.weight')
dat_ME_EB1_100 = loadtxt('data_ME_EB1_100/NP640_C100_T3.weight')
dat_ME_EB1_200 = loadtxt('data_ME_EB1_200/NP640_C100_T3.weight')
dat_ME_EB1_300 = loadtxt('data_ME_EB1_300/NP640_C100_T3.weight')
dat_ME_EB1_400 = loadtxt('data_ME_EB1_400/NP640_C100_T3.weight')
dat_ME_EB1_500 = loadtxt('data_ME_EB1_500/NP640_C100_T3.weight')
dat_ME_EB1_600 = loadtxt('data_ME_EB1_600/NP640_C100_T3.weight')
dat_ME_EB1_700 = loadtxt('data_ME_EB1_700/NP640_C100_T3.weight')
dat_ME_EB1_800 = loadtxt('data_ME_EB1_800/NP640_C100_T3.weight')
dat_ME_EB1_900 = loadtxt('data_ME_EB1_900/NP640_C100_T3.weight')
dat_ME_EB1_000 = loadtxt('data_ME_EB1_000/NP640_C100_T3.weight')
dat_ME_EB2_000 = loadtxt('data_ME_EB2_000/NP640_C100_T3.weight')
dat_ME_EB2_303 = loadtxt('data_ME_EB2_303/NP640_C100_T3.weight')
dat_ME_EB3_000 = loadtxt('data_ME_EB3_000/NP640_C100_T3.weight')
dat_ME_EB4_000 = loadtxt('data_ME_EB4_000/NP640_C100_T3.weight')
dat_ME_EB5_000 = loadtxt('data_ME_EB5_000/NP640_C100_T3.weight')
dat_ME_EB6_000 = loadtxt('data_ME_EB6_000/NP640_C100_T3.weight')
dat_ME_EB7_000 = loadtxt('data_ME_EB7_000/NP640_C100_T3.weight')
dat_ME_EB8_000 = loadtxt('data_ME_EB8_000/NP640_C100_T3.weight')
dat_ME_EB9_000 = loadtxt('data_ME_EB9_000/NP640_C100_T3.weight')
dat_ME_EB10_000 = loadtxt('data_ME_EB10_000/NP640_C100_T3.weight')

dat = asarray([[0.000, get_N_assoc(dat_ME_EB0_000, Np)[0]],
               [0.100, get_N_assoc(dat_ME_EB0_100, Np)[0]],
               [0.200, get_N_assoc(dat_ME_EB0_200, Np)[0]],
               [0.300, get_N_assoc(dat_ME_EB0_300, Np)[0]],
               [0.400, get_N_assoc(dat_ME_EB0_400, Np)[0]],
               [0.500, get_N_assoc(dat_ME_EB0_500, Np)[0]],
               # [0.600, get_N_assoc(dat_ME_EB0_600, Np)[0]],
               # [0.700, get_N_assoc(dat_ME_EB0_700, Np)[0]],
               # [0.800, get_N_assoc(dat_ME_EB0_800, Np)[0]],
               # [0.900, get_N_assoc(dat_ME_EB0_900, Np)[0]],
               [1.000, get_N_assoc(dat_ME_EB1_000, Np)[0]],
               # [1.100, get_N_assoc(dat_ME_EB1_100, Np)[0]],
               # [1.200, get_N_assoc(dat_ME_EB1_200, Np)[0]],
               # [1.300, get_N_assoc(dat_ME_EB1_300, Np)[0]],
               # [1.400, get_N_assoc(dat_ME_EB1_400, Np)[0]],
               [1.500, get_N_assoc(dat_ME_EB1_500, Np)[0]],
               # [1.600, get_N_assoc(dat_ME_EB1_600, Np)[0]],
               # [1.700, get_N_assoc(dat_ME_EB1_700, Np)[0]],
               # [1.800, get_N_assoc(dat_ME_EB1_800, Np)[0]],
               # [1.900, get_N_assoc(dat_ME_EB1_900, Np)[0]],
               [2.000, get_N_assoc(dat_ME_EB2_000, Np)[0]],
               [2.303, get_N_assoc(dat_ME_EB2_303, Np)[0]],
               [3.000, get_N_assoc(dat_ME_EB3_000, Np)[0]],
               [4.000, get_N_assoc(dat_ME_EB4_000, Np)[0]],
               [5.000, get_N_assoc(dat_ME_EB5_000, Np)[0]],
               [6.000, get_N_assoc(dat_ME_EB6_000, Np)[0]],
               [7.000, get_N_assoc(dat_ME_EB7_000, Np)[0]],
               [8.000, get_N_assoc(dat_ME_EB8_000, Np)[0]],
               [9.000, get_N_assoc(dat_ME_EB9_000, Np)[0]],
               [10.000, get_N_assoc(dat_ME_EB10_000, Np)[0]]])

plt.clf()
plt.ion()


plt.plot(dat[:,0], dat[:,1], 'bo-', markerfacecolor='white', markeredgecolor='blue', label='NAS, energy barrier')
i=0
set_x, set_y = dat[i,0]+2, dat[i,1]
plt.annotate('(P, r) = (%4.3f, %4.3f)'%(exp(-dat[i, 0]), dat[i,0]/0.24), xy=(dat[i,0], dat[i,1]), xytext=(set_x, set_y), size='x-small', arrowprops=dict(arrowstyle='-'))

i=4
set_x, set_y = dat[i,0]+2, dat[i,1] + 500
plt.annotate('(P, r) = (%4.3f, %4.3f)'%(exp(-dat[i, 0]), dat[i,0]/0.24), xy=(dat[i,0], dat[i,1]), xytext=(set_x, set_y), size='x-small', arrowprops=dict(arrowstyle='-'))

i=8
set_x, set_y = dat[i,0] - 1, dat[i,1] - 400
plt.annotate('(P, r) = (%4.3f, %4.3f)'%(exp(-dat[i, 0]), dat[i,0]/0.24), xy=(dat[i,0], dat[i,1]), xytext=(set_x, set_y), size='x-small', arrowprops=dict(arrowstyle='-'))

i=11
set_x, set_y = dat[i,0] - 1, dat[i,1] - 1000
plt.annotate('(P, r) = (%4.3f, %4.3f)'%(exp(-dat[i, 0]), dat[i,0]/0.24), xy=(dat[i,0], dat[i,1]), xytext=(set_x, set_y), size='x-small', arrowprops=dict(arrowstyle='-'))

i=12
set_x, set_y = dat[i,0] + 1, dat[i,1] + 200
plt.annotate('(P, r) = (%4.3f, %4.3f)'%(exp(-dat[i, 0]), dat[i,0]/0.24), xy=(dat[i,0], dat[i,1]), xytext=(set_x, set_y), size='x-small', arrowprops=dict(arrowstyle='-'))

i=13
set_x, set_y = dat[i,0] + 1, dat[i,1] + 200
plt.annotate('(P, r) = (%4.3f, %4.3f)'%(exp(-dat[i, 0]), dat[i,0]/0.24), xy=(dat[i,0], dat[i,1]), xytext=(set_x, set_y), size='x-small', arrowprops=dict(arrowstyle='-'))

i=14
set_x, set_y = dat[i,0] + 1, dat[i,1] + 200
plt.annotate('(P, r) = (%4.3f, %4.3f)'%(exp(-dat[i, 0]), dat[i,0]/0.24), xy=(dat[i,0], dat[i,1]), xytext=(set_x, set_y), size='x-small', arrowprops=dict(arrowstyle='-'))

i=15
set_x, set_y = dat[i,0] , dat[i,1] + 500
plt.annotate('(P, r) = (%4.3f, %4.3f)'%(exp(-dat[i, 0]), dat[i,0]/0.24), xy=(dat[i,0], dat[i,1]), xytext=(set_x, set_y), size='x-small', arrowprops=dict(arrowstyle='-'))

    
NAS_NR = get_N_assoc(dat_NR, Np)[0]
ref_NR = asarray([[0., NAS_NR],
                  [max(dat[:,0]), NAS_NR]])
plt.plot(ref_NR[:,0], ref_NR[:,1], 'r--', label='original (w/o barrier)')
# plt.plot(0.0, get_N_assoc(dat_NR, Np)[0], 'ro')
plt.show()

# plt.plot(get_N_assoc(dat_NR, Np), 'k-', label = 'kinetics:normalized')
# plt.plot(get_N_assoc(dat_ME_PB01, Np), 'k.', label = 'kinetics:metropolis, with Eb=0')
# plt.plot(get_N_assoc(dat_ME, Np), 'b-', label = 'kinetics:metropolis, with Eb=1')
# plt.plot(get_N_assoc(dat_ME_PB10, Np), 'r-', label = 'kinetics:metropolis, with Eb=2.30259, PB=10')


# plt.plot(get_N_assoc(dat, Np), 'b-', label = 'data')
# plt.plot(get_N_assoc(dat_long, Np), 'r.', label ='MKL_LONG')
# plt.plot(get_N_assoc(dat_INDEX_FIX, Np), 'r-', label = 'INDEX_FIX')
# plt.plot(get_N_assoc(dat_INDEX_FIX_BOOST, Np), 'b.', label = 'INDEX_FIX_BOOST')
# plt.plot(arange(N1), get_N_assoc(dat_INDEX_FIX_FLAG, Np), 'b-', label = 'FIX_FLAG') # same with previous
# plt.plot(get_N_assoc(dat_INDEX_FIX_FLAG_longer, Np), 'r-', label = 'INDEX_FIX_FLAG_longer')

# plt.plot(arange(N2)*10, get_N_assoc(dat_TEST_longer, Np), 'r.-', label = 'TEST_longer')

# plt.axis([0, 10, -50, 5000])
plt.xlabel('dimensionless energy barrier')
plt.ylabel('number of associations')
plt.legend(loc = 'upper right', numpoints=1)
# plt.grid()
plt.show()
