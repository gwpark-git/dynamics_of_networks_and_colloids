
from numpy import *
import matplotlib.pyplot as plt

def get_N_assoc(dat, Np):
    Nt = int(shape(dat)[0]/Np)
    re = zeros(Nt)
    for ts in range(Nt):
        re[ts] = sum(dat[ts*Np:(ts+1)*Np, 1:])/2.0
    return re
Np = 640

# dat_NR = loadtxt('data_NR/NP640_C100_T3.weight')
# dat_ME_EB0_000 = loadtxt('data_ME_EB0_000/NP640_C100_T3.weight')
# dat_ME_EB0_500 = loadtxt('data_ME_EB0_500/NP640_C100_T3.weight')
# dat_ME_EB1_000 = loadtxt('data_ME_EB1_000/NP640_C100_T3.weight')
# dat_ME_EB2_000 = loadtxt('data_ME_EB2_000/NP640_C100_T3.weight')
# dat_ME_EB2_303 = loadtxt('data_ME_EB2_303/NP640_C100_T3.weight')
# dat_ME_EB3_000 = loadtxt('data_ME_EB3_000/NP640_C100_T3.weight')
# dat_ME_EB4_000 = loadtxt('data_ME_EB4_000/NP640_C100_T3.weight')
# dat_ME_EB5_000 = loadtxt('data_ME_EB5_000/NP640_C100_T3.weight')
# dat_ME_EB6_000 = loadtxt('data_ME_EB6_000/NP640_C100_T3.weight')
# dat_ME_EB7_000 = loadtxt('data_ME_EB7_000/NP640_C100_T3.weight')
# dat_ME_EB8_000 = loadtxt('data_ME_EB8_000/NP640_C100_T3.weight')
# dat_ME_EB9_000 = loadtxt('data_ME_EB9_000/NP640_C100_T3.weight')
# dat_ME_EB10_000 = loadtxt('data_ME_EB10_000/NP640_C100_T3.weight')

# dat = asarray([[0.000, get_N_assoc(dat_ME_EB0_000, Np)[0]],
#                [0.500, get_N_assoc(dat_ME_EB0_500, Np)[0]],
#                [1.000, get_N_assoc(dat_ME_EB1_000, Np)[0]],
#                [2.000, get_N_assoc(dat_ME_EB2_000, Np)[0]],
#                [2.303, get_N_assoc(dat_ME_EB2_303, Np)[0]],
#                [3.000, get_N_assoc(dat_ME_EB3_000, Np)[0]],
#                [4.000, get_N_assoc(dat_ME_EB4_000, Np)[0]],
#                [5.000, get_N_assoc(dat_ME_EB5_000, Np)[0]],
#                [6.000, get_N_assoc(dat_ME_EB6_000, Np)[0]],
#                [7.000, get_N_assoc(dat_ME_EB7_000, Np)[0]],
#                [8.000, get_N_assoc(dat_ME_EB8_000, Np)[0]],
#                [9.000, get_N_assoc(dat_ME_EB9_000, Np)[0]],
#                [10.000, get_N_assoc(dat_ME_EB10_000, Np)[0]]])

dat_e5 = loadtxt('test_1_e5.dat')
dat_e6 = loadtxt('test_2_e6.dat')
dat_e7 = loadtxt('test_3_e7.dat')


plt.clf()
plt.ion()


plt.plot(dat_e5[:,0], dat_e5[:,1], 'b.-', markerfacecolor='white', markeredgecolor='blue', label='NAS, MC_EQ blocks=1e5')
plt.plot(dat_e6[:,0], dat_e6[:,1], 'g.-', markerfacecolor='white', markeredgecolor='green', label='NAS, MC_EQ blocks=1e6')
plt.plot(dat_e7[:,0], dat_e7[:,1], 'r.-', markerfacecolor='white', markeredgecolor='red', label='NAS, MC_EQ blocks=1e7')

ref_5 = asarray([[0, 5000],
                 [10, 5000]])
plt.plot(ref_5[:,0], ref_5[:,1], 'k--', linewidth=2, label = 'ref. approx. previous scheme')

# i=0
# set_x, set_y = dat[i,0]+2, dat[i,1]
# plt.annotate('(P, r) = (%4.3f, %4.3f)'%(exp(-dat[i, 0]), dat[i,0]/0.24), xy=(dat[i,0], dat[i,1]), xytext=(set_x, set_y), size='x-small', arrowprops=dict(arrowstyle='-'))

# i=1
# set_x, set_y = dat[i,0]+2, dat[i,1] + 500
# plt.annotate('(P, r) = (%4.3f, %4.3f)'%(exp(-dat[i, 0]), dat[i,0]/0.24), xy=(dat[i,0], dat[i,1]), xytext=(set_x, set_y), size='x-small', arrowprops=dict(arrowstyle='-'))

# i=4
# set_x, set_y = dat[i,0] - 1, dat[i,1] - 400
# plt.annotate('(P, r) = (%4.3f, %4.3f)'%(exp(-dat[i, 0]), dat[i,0]/0.24), xy=(dat[i,0], dat[i,1]), xytext=(set_x, set_y), size='x-small', arrowprops=dict(arrowstyle='-'))

# i=7
# set_x, set_y = dat[i,0] - 1, dat[i,1] - 1000
# plt.annotate('(P, r) = (%4.3f, %4.3f)'%(exp(-dat[i, 0]), dat[i,0]/0.24), xy=(dat[i,0], dat[i,1]), xytext=(set_x, set_y), size='x-small', arrowprops=dict(arrowstyle='-'))

# i=8
# set_x, set_y = dat[i,0] + 1, dat[i,1] + 200
# plt.annotate('(P, r) = (%4.3f, %4.3f)'%(exp(-dat[i, 0]), dat[i,0]/0.24), xy=(dat[i,0], dat[i,1]), xytext=(set_x, set_y), size='x-small', arrowprops=dict(arrowstyle='-'))

# i=9
# set_x, set_y = dat[i,0] + 1, dat[i,1] + 200
# plt.annotate('(P, r) = (%4.3f, %4.3f)'%(exp(-dat[i, 0]), dat[i,0]/0.24), xy=(dat[i,0], dat[i,1]), xytext=(set_x, set_y), size='x-small', arrowprops=dict(arrowstyle='-'))

# i=10
# set_x, set_y = dat[i,0] + 1, dat[i,1] + 200
# plt.annotate('(P, r) = (%4.3f, %4.3f)'%(exp(-dat[i, 0]), dat[i,0]/0.24), xy=(dat[i,0], dat[i,1]), xytext=(set_x, set_y), size='x-small', arrowprops=dict(arrowstyle='-'))

# i=11
# set_x, set_y = dat[i,0] , dat[i,1] + 500
# plt.annotate('(P, r) = (%4.3f, %4.3f)'%(exp(-dat[i, 0]), dat[i,0]/0.24), xy=(dat[i,0], dat[i,1]), xytext=(set_x, set_y), size='x-small', arrowprops=dict(arrowstyle='-'))

    
# NAS_NR = get_N_assoc(dat_NR, Np)[0]
# ref_NR = asarray([[0., NAS_NR],
#                   [max(dat[:,0]), NAS_NR]])
# plt.plot(ref_NR[:,0], ref_NR[:,1], 'r--', label='original (w/o barrier)')
# # plt.plot(0.0, get_N_assoc(dat_NR, Np)[0], 'ro')
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

plt.axis([0, 10, -50, 7000])
plt.xlabel('dimensionless energy barrier')
plt.ylabel('number of associations')
plt.legend(loc = 'lower left', numpoints=1)
plt.grid()
plt.show()
