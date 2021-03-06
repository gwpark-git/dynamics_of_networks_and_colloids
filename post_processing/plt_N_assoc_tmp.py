
# from numpy import *
# import matplotlib.pyplot as plt

# def get_N_assoc(dat, Np):
#     Nt = int(shape(dat)[0]/Np)
#     re = zeros(Nt)
#     for ts in range(Nt):
#         re[ts] = sum(dat[ts*Np:(ts+1)*Np, 1:])/2.0
#     return re
# # Np = 1350

# # NAS_100 = get_N_assoc(loadtxt('TEST_3D_NP3200_LD20P3_SF15/NP3200_C100_T3.weight'), Np)
# # NAS_050 = get_N_assoc(loadtxt('TEST_3D_NP3200_LD20P3_SF15_half/NP3200_C100_T3.weight'), Np)

# NAS_NP1350 = get_N_assoc(loadtxt('NP1350_LD15P3_C100_SF15/NP1350_LD15P3_C100.weight'), 1350)
# NAS_NP3200 = get_N_assoc(loadtxt('TEST_3D_NP3200_LD20P3_SF15/NP3200_C100_T3.weight'), 3200)

plt.clf()
plt.ion()
plt.figure(figsize=(11,6))
ax = plt.subplot(111)
Nt = size(NAS_NP1350)
dt = 0.001
N_skip = 100

t = arange(Nt)*dt*N_skip
Rt = 100
t = t/Rt

ax.plot(t[:size(NAS_NP1350)], NAS_NP1350, 'r.', label =r'NP1350/15^3 = 0.4')
ax.plot(t, NAS_NP3200, 'b-', label =r'NP3200/20^3 = 0.4')

ax.set_xlabel('Tological dimensionless time')
# ax.set_xlabel('Brownian dimensionless time')
ax.set_ylabel('number of associations')
ax.legend(loc = 'lower right', numpoints=1)
ax.grid()

# ax2 = ax.twinx()
# ax2.axis([0, t[-1], 0, 4000./(25.*3200.)])
# ax2.set_ylabel('fraction = NAS/80,000')
# plt.plot(get_N_assoc(dat_INDEX_FIX, Np), 'r-', label = 'INDEX_FIX')
# plt.plot(get_N_assoc(dat_INDEX_FIX_BOOST, Np), 'b.', label = 'INDEX_FIX_BOOST')
# plt.plot(arange(N1), get_N_assoc(dat_INDEX_FIX_FLAG, Np), 'b-', label = 'FIX_FLAG') # same with previous
# plt.plot(get_N_assoc(dat_INDEX_FIX_FLAG_longer, Np), 'r-', label = 'INDEX_FIX_FLAG_longer')

# plt.plot(arange(N2)*10, get_N_assoc(dat_TEST_longer, Np), 'r.-', label = 'TEST_longer')


plt.show()
