
from numpy import *
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import acf
import sys
sys.path.append('../../../post_processing')
from acf_fcn import *

dat_NP0400 = loadtxt('RF_NP0400_RT1k.dat')[:,1]/(2.*1000.)
acf_NP0400 = acf_gro_BAV(dat_NP0400, 2)
t_acf_NP0400 = arange(size(acf_NP0400))/100.

dat_NP1000 = (loadtxt('RF_NP1000_100tau0.dat')[:,1]/(2.*1000.))[1000:]
t_dat_NP1000 = arange(size(dat_NP1000))
acf_NP1000 = acf_gro_BAV(dat_NP1000, 2)
t_acf_NP1000 = arange(size(acf_NP1000))/100.

from scipy.stats import linregress
reg_BR_NP0400 = linregress(t_acf_NP0400[0:5], log(acf_NP0400[0:5]/acf_NP0400[0]))
reg_BR_NP1000 = linregress(t_acf_NP1000[0:5], log(acf_NP1000[0:5]/acf_NP1000[0]))

reg_NP0400 = linregress(t_acf_NP0400[5:15], log(acf_NP0400[5:15]/acf_NP0400[0]))
reg_NP1000 = linregress(t_acf_NP1000[10:30], log(acf_NP1000[10:30]/acf_NP1000[0]))

t_reg = linspace(0, 1, 100)

y_reg_BR_NP0400 = exp(reg_BR_NP0400[0]*t_reg + reg_BR_NP0400[1])
y_reg_BR_NP1000 = exp(reg_BR_NP1000[0]*t_reg + reg_BR_NP1000[1])

y_reg_NP0400 = exp(reg_NP0400[0]*t_reg + reg_NP0400[1])
y_reg_NP1000 = exp(reg_NP1000[0]*t_reg + reg_NP1000[1])



from matplotlib.font_manager import FontProperties
fontP = FontProperties()
# fontP.set_size('small')

plt.close()
plt.ion()
plt.figure(figsize=(11,8))
ax = plt.subplot(311)
ax.plot(t_acf_NP0400, acf_NP0400/acf_NP0400[0], 'b-', label = 'NP=400')
# ax.plot(t_reg, y_reg_BR_NP0400, 'b--')
# ax.plot(t_reg, y_reg_NP0400, 'b:')
ax.plot(t_acf_NP1000, acf_NP1000/acf_NP1000[0], 'r-', label = 'NP=1000')
# ax.plot(t_reg, y_reg_BR_NP1000, 'r--')
# ax.plot(t_reg, y_reg_NP1000, 'r:')
ax.axis([0, 1, -0.2, 1])
ax.grid()
ax.legend(loc = 'upper right', prop = fontP)
# ax.set_xlabel(r'dimensionless time / $\tau_0$')
ax.set_ylabel(r'normalized stress autocorrelation')
ax.set_title('linear-linear', y=0.85)

ax2 = plt.subplot(312)

ax2.semilogy(t_acf_NP0400, acf_NP0400/acf_NP0400[0], 'b-')
ax2.semilogy(t_reg[:10], y_reg_BR_NP0400[:10], 'b--', label = r'$t_c$(%4.2f, %4.2f) = %4.2f'%(t_acf_NP0400[0], t_acf_NP0400[5], -1./reg_BR_NP0400[0]))
ax2.semilogy(t_reg, y_reg_NP0400, 'b:', label = r'$t_c$(%4.2f, %4.2f) = %4.2f'%(t_acf_NP0400[5], t_acf_NP0400[15], -1./reg_NP0400[0]))

ax2.semilogy(t_acf_NP1000, acf_NP1000/acf_NP1000[0], 'r-')
ax2.semilogy(t_reg[:10], y_reg_BR_NP1000[:10], 'r--', label = r'$t_c$(%4.2f, %4.2f) = %4.2f'%(t_acf_NP1000[0], t_acf_NP1000[5], -1./reg_BR_NP1000[0]))
ax2.semilogy(t_reg, y_reg_NP1000, 'r:', label = r'$t_c$(%4.2f, %4.2f) = %4.2f'%(t_acf_NP1000[10], t_acf_NP1000[30], -1./reg_NP1000[0]))
ax2.axis([0, 1, 10**-4, 10**0])
# ax2.axis([0, 2, -0.1, 1])
ax2.grid()
ax2.legend(loc = 'lower right', prop = fontP, ncol=2)
ax2.set_xlabel(r'dimensionless time / $\tau_0$')

# ax2.set_ylabel(r'normalized stress autocorrelation')
ax2.set_title('linear-log', y=0.85)

ax3 = plt.subplot(313)

def find_val(dat, val):
    for i,x in enumerate(dat):
        if x < val:
            return i
    return i

tsf_NP0400 = t_acf_NP0400[find_val(acf_NP0400/acf_NP0400[0], 0.1)]
tsf_NP1000 = t_acf_NP1000[find_val(acf_NP1000/acf_NP1000[0], 0.1)]

ax3.loglog(t_acf_NP0400/tsf_NP0400, acf_NP0400/acf_NP0400[0], 'b-')
ax3.loglog(t_acf_NP1000/tsf_NP1000, acf_NP1000/acf_NP1000[0], 'r-')
ax3.axis([t_acf_NP0400[2], 10**2, 10**-3, 10**0])
ax3.grid()
ax3.set_xlabel(r'dimensionless time / $\tau_0$')
ax3.set_ylabel('streess autocorrelation')
ax3.set_title('log-log')


plt.show()

# plt.close()
# plt.ion()
# plt.plot(t_dat_NP1000, dat_NP1000_tmp, 'b-')
# plt.show()
