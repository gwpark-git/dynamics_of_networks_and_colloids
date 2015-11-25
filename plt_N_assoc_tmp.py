
from numpy import *
import matplotlib.pyplot as plt

def get_N_assoc(dat, Np):
    Nt = int(shape(dat)[0]/Np)
    re = zeros(Nt)
    for ts in range(Nt):
        re[ts] = sum(dat[ts*Np:(ts+1)*Np, 1:])/2.0
    return re
Np = 640

dat = loadtxt('data/NP640_C100_T3_r11.weight')
dat_long = loadtxt('data_MKL_LONG/NP640_C100_T3_r11.weight')
# dat_INDEX_FIX = loadtxt('data_INDEX_FIX/NP640_C100_T3_r11.weight')
# # dat_INDEX_FIX_BOOST = loadtxt('data_INDEX_FIX_BOOST/NP640_C100_T3_r11.weight')
# dat_INDEX_FIX_FLAG = loadtxt('data_INDEX_FIX_FLAG/NP640_C100_T3_r11.weight')
# N1 = shape(dat_INDEX_FIX_FLAG)[0]/Np
# # dat_INDEX_FIX_FLAG_longer = loadtxt('data_INDEX_FIX_FLAG_longer/NP640_C100_T3_r11.weight')
# dat_TEST_longer = loadtxt('data_TEST_longer/NP640_C100_T3_r11.weight')
# N2 = shape(dat_TEST_longer)[0]/Np


plt.clf()
plt.ion()

plt.plot(get_N_assoc(dat, Np), 'b-', label = 'data')
plt.plot(get_N_assoc(dat_long, Np), 'r.', label ='MKL_LONG')
# plt.plot(get_N_assoc(dat_INDEX_FIX, Np), 'r-', label = 'INDEX_FIX')
# plt.plot(get_N_assoc(dat_INDEX_FIX_BOOST, Np), 'b.', label = 'INDEX_FIX_BOOST')
# plt.plot(arange(N1), get_N_assoc(dat_INDEX_FIX_FLAG, Np), 'b-', label = 'FIX_FLAG') # same with previous
# plt.plot(get_N_assoc(dat_INDEX_FIX_FLAG_longer, Np), 'r-', label = 'INDEX_FIX_FLAG_longer')

# plt.plot(arange(N2)*10, get_N_assoc(dat_TEST_longer, Np), 'r.-', label = 'TEST_longer')


plt.xlabel('time step')
plt.ylabel('number of associations')
plt.legend(loc = 'lower right', numpoints=1)
plt.grid()
plt.show()
