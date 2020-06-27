
from numpy import *

def Wi_from_str(Wi_str):
    identifier = ''
    sign = 1.
    if 'm' in Wi_str:
        identifier = 'm'
        sign = -1.
        # prefactor, exponent = Wi_str.split('m')
    elif 'p' in Wi_str:
        identifier = 'p'
        # prefactor, exponent = Wi_str.split('p')
    else:
        print 'Error, unidentified identifier: ', Wi_str
        return -1
    prefactor, exponent = Wi_str.split(identifier)
    if prefactor == '':
        prefactor = 1.
    return float(prefactor)*(10**(sign*float(exponent)))

def get_Wi_str(fn_str):
    return fn_str[38:41].replace('_','')


# fn_list = []
# with open('list_RT1div25.dat') as f:
#     for line in f:
#         fn_list.append(line.replace('\n', ''))

# N_dat = size(fn_list)
# ener = []
# for fn in fn_list:
#     ener.append(loadtxt(fn))

N_samples = 5
N_cut = int(shape(ener[0])[0]*0.8)

tau_xy_con = []
tau_xy_rep = []
tau_xy = []
Wi = []
for i in range(N_dat/N_samples):
    index = i*N_samples
    Wi.append(Wi_from_str(get_Wi_str(fn_list[index])))
    tau_xy_con.append(mean(ener[index][N_cut:, 27]))
    tau_xy_rep.append(mean(ener[index][N_cut:, 21]))
    tau_xy.append(tau_xy_con[i] + tau_xy_rep[i])

# the following maps in increasing order of Wi
Wi, tau_xy_con, tau_xy_rep, tau_xy = map(list, zip(*sorted(zip(Wi, tau_xy_con, tau_xy_rep, tau_xy))))

    
import matplotlib.pyplot as plt
plt.close()
plt.ion()
plt.plot(Wi, tau_xy_con, 'bo-')
plt.plot(Wi, abs(asarray(tau_xy_rep)), 'ro--', markerfacecolor='white', markeredgecolor='red')
plt.plot(Wi, tau_xy_rep, 'ro-')
plt.plot(Wi, tau_xy, 'ko-')
plt.xscale('log')
plt.yscale('log')
plt.show()

# import matplotlib.pyplot as plt

# plt.close()
# plt.ion()
# ind = 20
# plt.plot(ener[ind][:, 0], ener[ind][:, 27], 'b-')
# plt.plot(ener[ind][:, 0], ener[ind][:, 21], 'r-')
# plt.plot(ener[ind][:, 0], ener[ind][:, 21] + ener[ind][:, 27], 'k-')

# plt.xscale('log')

# plt.show()
