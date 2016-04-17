
from numpy import *
import matplotlib.pyplot as plt
import sys
sys.path.append('../../../post_processing')
from acf_fcn import *


def corr_Gt(dat):
    s_xx, s_xy, s_xz, s_yy, s_yz, s_zz = dat[:,0], dat[:,1], dat[:,2], dat[:, 4], dat[:,5], dat[:,8]
    N_xy = s_xx - s_yy
    N_xz = s_xx - s_zz
    N_yz = s_yy - s_zz
    return acf_gro(s_xy) + acf_gro(s_yz) + acf_gro(s_xz) + (1./6.)*(acf_gro(N_xy) + acf_gro(N_xz) + acf_gro(N_yz))



# dat_NP0400 = (loadtxt('RF_NP0400_tmp.dat')/(2.*1000.))[:,:]
# t_dat_NP0400 = arange(shape(dat_NP0400)[0])
# acf_NP0400 = corr_Gt(dat_NP0400)
# acf_NP0400_simple = acf_gro(dat_NP0400[:,1])
# t_acf_NP0400 = arange(shape(acf_NP0400)[0])/1000.


# dat_NP1000 = (loadtxt('RF_NP1000_5tau0.dat')/(2.*1000.))
# t_dat_NP1000 = arange(shape(dat_NP1000)[0])
# acf_NP1000 = corr_Gt(dat_NP1000)
# acf_NP1000_simple = acf_gro(dat_NP1000[:,1])
# t_acf_NP1000 = arange(shape(acf_NP1000)[0])/1000.



from scipy.stats import linregress
reg_BR_NP0400 = linregress(t_acf_NP0400[0:5], log(acf_NP0400[0:5]/acf_NP0400[0]))
reg_BR_NP1000 = linregress(t_acf_NP1000[0:5], log(acf_NP1000[0:5]/acf_NP1000[0]))
# reg_BR_NP1728 = linregress(t_acf_NP1728[0:5], log(acf_NP1728[0:5]/acf_NP1728[0]))

reg_NP0400 = linregress(t_acf_NP0400[5:15], log(acf_NP0400[5:15]/acf_NP0400[0]))
reg_NP1000 = linregress(t_acf_NP1000[10:30], log(acf_NP1000[10:30]/acf_NP1000[0]))
# reg_NP1728 = linregress(t_acf_NP1728[10:30], log(acf_NP1728[10:30]/acf_NP1728[0]))

t_reg = linspace(0, 1, 100)

y_reg_BR_NP0400 = exp(reg_BR_NP0400[0]*t_reg + reg_BR_NP0400[1])
y_reg_BR_NP1000 = exp(reg_BR_NP1000[0]*t_reg + reg_BR_NP1000[1])
# y_reg_BR_NP1728 = exp(reg_BR_NP1728[0]*t_reg + reg_BR_NP1728[1])

y_reg_NP0400 = exp(reg_NP0400[0]*t_reg + reg_NP0400[1])
y_reg_NP1000 = exp(reg_NP1000[0]*t_reg + reg_NP1000[1])
# y_reg_NP1728 = exp(reg_NP1728[0]*t_reg + reg_NP1728[1])



# from matplotlib.font_manager import FontProperties
# fontP = FontProperties()
# # fontP.set_size('small')

# plt.close()
# plt.ion()
# plt.figure(figsize=(11,8))
# ax = plt.subplot(311)
# ax.plot(t_acf_NP0400, acf_NP0400/acf_NP0400[0], 'b-', label = 'NP=400')
# # ax.plot(t_reg, y_reg_BR_NP0400, 'b--')
# # ax.plot(t_reg, y_reg_NP0400, 'b:')
# ax.plot(t_acf_NP1000, acf_NP1000/acf_NP1000[0], 'g-', label = 'NP=1000')
# # ax.plot(t_acf_NP1000, acf_NP1000_simple/acf_NP1000_simple[0], 'r-', label = 'NP=1000, simple')

# # ax.plot(t_acf_NP1728, acf_NP1728/acf_NP1728[0], 'r-', label = 'NP=1728')

# # ax.plot(t_reg, y_reg_BR_NP1000, 'r--')
# # ax.plot(t_reg, y_reg_NP1000, 'r:')
# ax.axis([0, 2, -0.2, 1])
# ax.grid()
# ax.legend(loc = 'upper right', prop = fontP)
# # ax.set_xlabel(r'dimensionless time / $\tau_0$')
# ax.set_ylabel(r'normalized stress autocorrelation')
# ax.set_title('linear-linear', y=0.85)

# ax2 = plt.subplot(312)

# ax2.semilogy(t_acf_NP0400, acf_NP0400/acf_NP0400[0], 'b-')
# ax2.semilogy(t_reg[:10], y_reg_BR_NP0400[:10], 'b--', label = r'$t_c$(%4.2f, %4.2f) = %4.2f'%(t_acf_NP0400[0], t_acf_NP0400[5], -1./reg_BR_NP0400[0]))
# ax2.semilogy(t_reg, y_reg_NP0400, 'b:', label = r'$t_c$(%4.2f, %4.2f) = %4.2f'%(t_acf_NP0400[5], t_acf_NP0400[15], -1./reg_NP0400[0]))

# ax2.semilogy(t_acf_NP1000, acf_NP1000/acf_NP1000[0], 'g-')
# ax2.semilogy(t_reg[:10], y_reg_BR_NP1000[:10], 'g--', label = r'$t_c$(%4.2f, %4.2f) = %4.2f'%(t_acf_NP1000[0], t_acf_NP1000[5], -1./reg_BR_NP1000[0]))
# ax2.semilogy(t_reg, y_reg_NP1000, 'g:', label = r'$t_c$(%4.2f, %4.2f) = %4.2f'%(t_acf_NP1000[10], t_acf_NP1000[30], -1./reg_NP1000[0]))

# ax2.semilogy(t_acf_NP1000, acf_NP1000_simple/acf_NP1000_simple[0], 'r-')

# ax2.axis([0, 1, 10**-4, 10**0])

# # ax2.semilogy(t_acf_NP1728, acf_NP1728/acf_NP1728[0], 'r-')
# # ax2.semilogy(t_reg[:10], y_reg_BR_NP1728[:10], 'r--', label = r'$t_c$(%4.2f, %4.2f) = %4.2f'%(t_acf_NP1728[0], t_acf_NP1728[5], -1./reg_BR_NP1728[0]))
# # ax2.semilogy(t_reg, y_reg_NP1728, 'r:', label = r'$t_c$(%4.2f, %4.2f) = %4.2f'%(t_acf_NP1728[10], t_acf_NP1728[30], -1./reg_NP1728[0]))
# # ax2.axis([0, 1, 10**-4, 10**0])


# # ax2.axis([0, 2, -0.1, 1])
# ax2.grid()
# ax2.legend(loc = 'lower right', prop = fontP, ncol=2)
# ax2.set_xlabel(r'dimensionless time / $\tau_0$')

# # ax2.set_ylabel(r'normalized stress autocorrelation')
# ax2.set_title('linear-log', y=0.85)

# ax3 = plt.subplot(313)

# ax3.loglog(t_acf_NP0400, acf_NP0400, 'b-')
# ax3.loglog(t_acf_NP1000, acf_NP1000/acf_NP1000[0], 'g-')
# ax3.loglog(t_acf_NP1000, acf_NP1000_simple/acf_NP1000_simple[0], 'r-')

# # ax3.loglog(t_acf_NP1728, acf_NP1728, 'r-')

# # ax3.axis([t_acf_NP0400[2], 10**0, 10**-7, 10**-3])
# ax3.grid()
# ax3.set_xlabel(r'dimensionless time / $\tau_0$')
# ax3.set_ylabel('stress autocorrelation')
# # ax3.set_title(r'log-log, $C_{\tau_{xy}}^{(1000)}(0)/C_{\tau_{xy}}^{(1728)}(0)$ = %4.1f'%(acf_NP1000[0]/acf_NP1728[0]), y=0.7)


# plt.show()

# # plt.close()
# # plt.ion()
# # plt.plot(t_dat_NP1000, dat_NP1000_tmp, 'b-')
# # plt.show()

plt.close()
plt.figure(figsize=(11,6))
plt.ion()
plt.loglog(t_acf_NP0400[:], acf_NP0400[:]/acf_NP0400[0], 'b-', linewidth=2, alpha=0.5, label = r'$\nu_m$ = 0.4')
plt.loglog(t_acf_NP1000[:], acf_NP1000[:]/acf_NP1000[0], 'r-', linewidth=2, alpha=0.5, label = r'$\nu_m$ = 1.0')

plt.loglog(t_acf_NP0400[14], acf_NP0400[14]/acf_NP0400[0], 'b|', markersize=10)
plt.loglog(t_acf_NP1000[26], acf_NP1000[26]/acf_NP1000[0], 'r|', markersize=10)

plt.legend(loc = 'lower left', numpoints=1)
plt.axis([10**-2, 1.1*10**0, 10**-3, 2.*10**0])
plt.xlabel(r'dimensionless time, $\beta_0 t$')
plt.ylabel(r'normalized stress autocorrelation, $C_{\tilde{\tau}_{xy}}(t)/C_{\tilde{\tau}_{xy}}(0)$')
plt.grid()
plt.show()
