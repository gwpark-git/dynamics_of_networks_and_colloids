
from numpy import *
import matplotlib.pyplot as plt
import sys
sys.path.append('../../post_processing')
from acf_fcn import *

N_st = 1000 # it is 10 tau_0
normal_factor = 2.*1000 # 1000 is volume
RF_NP0400 = loadtxt('RF_NP0400_NC25_RT100.dat')[N_st:,1]/normal_factor
acf_NP0400 = acf_gro(RF_NP0400)
t_NP0400 = arange(size(acf_NP0400))/100.

RF_NP0600 = loadtxt('RF_NP0600_NC25_RT100.dat')[N_st:,1]/normal_factor
acf_NP0600 = acf_gro(RF_NP0600)
t_NP0600 = arange(size(acf_NP0600))/100.

RF_NP0800 = loadtxt('RF_NP0800_NC25_RT100.dat')[N_st:,1]/normal_factor
acf_NP0800 = acf_gro(RF_NP0800)
t_NP0800 = arange(size(acf_NP0800))/100.

RF_NP1000 = loadtxt('RF_NP1000_RT100_dt0_01.dat')[N_st:,1]/normal_factor
acf_NP1000 = acf_gro(RF_NP1000)
t_NP1000 = arange(size(acf_NP1000))/100.

RF_NP1200 = loadtxt('RF_NP1200_RT100_dt0_01.dat')[N_st:,1]/normal_factor
acf_NP1200 = acf_gro(RF_NP1200)
t_NP1200 = arange(size(acf_NP1200))/100.

RF_NP1400 = loadtxt('RF_NP1400_RT100_dt0_01.dat')[N_st:,1]/normal_factor
acf_NP1400 = acf_gro(RF_NP1400)
t_NP1400 = arange(size(acf_NP1400))/100.

# RF_NP2000 = loadtxt('RF_NP2000_RT100_dt0_01.dat')[N_st:,1]
# acf_NP2000 = acf_gro(RF_NP2000)
# t_NP2000 = arange(size(acf_NP2000))/100.

# RF_NP4000 = loadtxt('RF_NP4000_RT100_dt0_01.dat')[N_st:,1]
# acf_NP4000 = acf_gro(RF_NP4000)
# t_NP4000 = arange(size(acf_NP4000))/100.


from scipy.stats import linregress
N1 = 5
a_400, y_400, r_400, p_400, std_400 = linregress(t_NP0400[:N1], log(acf_NP0400[:N1]))
a_600, y_600, r_600, p_600, std_600 = linregress(t_NP0600[:N1], log(acf_NP0600[:N1]))
a_800, y_800, r_800, p_800, std_800 = linregress(t_NP0800[:N1], log(acf_NP0800[:N1]))
a_1000, y_1000, r_1000, p_1000, std_1000 = linregress(t_NP1000[:N1], log(acf_NP1000[:N1]))
a_1200, y_1200, r_1200, p_1200, std_1200 = linregress(t_NP1200[:N1], log(acf_NP1200[:N1]))
a_1400, y_1400, r_1400, p_1400, std_1400 = linregress(t_NP1400[:N1], log(acf_NP1400[:N1]))

a_2000, y_2000, r_2000, p_2000, std_2000 = linregress(t_NP2000[:N1], log(acf_NP2000[:N1]))
a_4000, y_4000, r_4000, p_4000, std_4000 = linregress(t_NP4000[:N1], log(acf_NP4000[:N1]))

N2_st = 10
N2 = 50
a2_400, y2_400, r2_400, p2_400, std2_400 = linregress(t_NP0400[N2_st:N2], log(acf_NP0400[N2_st:N2]))
a2_600, y2_600, r2_600, p2_600, std2_600 = linregress(t_NP0600[N2_st:N2], log(acf_NP0600[N2_st:N2]))
a2_800, y2_800, r2_800, p2_800, std2_800 = linregress(t_NP0800[N2_st:N2], log(acf_NP0800[N2_st:N2]))
a2_1000, y2_1000, r2_1000, p2_1000, std2_1000 = linregress(t_NP1000[N2_st:N2], log(acf_NP1000[N2_st:N2]))
# a2_1200, y2_1200, r2_1200, p2_1200, std2_1200 = linregress(t_NP1200[N2_st:N2], log(acf_NP1200[N2_st:N2]))
a2_1400, y2_1400, r2_1400, p2_1400, std2_1400 = linregress(t_NP1400[N2_st:N2], log(acf_NP1400[N2_st:N2]))

# a2_2000, y2_2000, r2_2000, p2_2000, std2_2000 = linregress(t_NP2000[N2_st:N2], log(acf_NP2000[N2_st:N2]))
# a2_4000, y2_4000, r2_4000, p2_4000, std2_4000 = linregress(t_NP4000[N2_st:N2], log(acf_NP4000[N2_st:N2]))

# N3_st = 50
# N3 = 100
# a3_400, y3_400, r3_400, p3_400, std3_400 = linregress(t_NP0400[N3_st:N3], log(acf_NP0400[N3_st:N3]))
# a3_600, y3_600, r3_600, p3_600, std3_600 = linregress(t_NP0600[N3_st:N3], log(acf_NP0600[N3_st:N3]))
# a3_800, y3_800, r3_800, p3_800, std3_800 = linregress(t_NP0800[N3_st:N3], log(acf_NP0800[N3_st:N3]))
# a3_1200, y3_1200, r3_1200, p3_1200, std3_1200 = linregress(t_NP1200[N3_st:N3], log(acf_NP1200[N3_st:N3]))

# N4_st = 50
# N4 = 220
# a3_1000, y3_1000, r3_1000, p3_1000, std3_1000 = linregress(t_NP1000[N4_st:N4], log(acf_NP1000[N4_st:N4]))
# a3_1400, y3_1400, r3_1400, p3_1400, std3_1400 = linregress(t_NP1400[N4_st:N4], log(acf_NP1400[N4_st:N4]))

# a3_2000, y3_2000, r3_2000, p3_2000, std3_2000 = linregress(t_NP2000[N4_st:N4], log(acf_NP2000[N4_st:N4]))
# a3_4000, y3_4000, r3_4000, p3_4000, std3_4000 = linregress(t_NP4000[N4_st:N4], log(acf_NP4000[N4_st:N4]))

from matplotlib.font_manager import FontProperties
fontP = FontProperties()
# fontP.set_size('small')

colP = ['blue', 'green', 'brown', 'red', 'cyan', 'purple', 'gray', 'black']

ref_zero = asarray([[0, 0], [25, 0]])
plt.close()
plt.ion()
plt.figure(figsize=(11,6))
plt.semilogy(t_NP0400, acf_NP0400/acf_NP0400[0], '+-', color=colP[0], label = r'$\nu_0=1.0/R_0^3$, $\lambda(t^{(1)})$ = %4.3f $\tau_0$, $\lambda(t^{(2)})$ = %4.3f $\tau_0$'%(-1./a_400, -1./a2_400))
plt.semilogy(t_NP0600, acf_NP0600/acf_NP0600[0], '+-', color=colP[1], label = r'$\nu_0=1.5/R_0^3$, $\lambda(t^{(1)})$ = %4.3f $\tau_0$, $\lambda(t^{(2)})$ = %4.3f $\tau_0$'%(-1./a_600, -1./a2_600))
plt.semilogy(t_NP0800, acf_NP0800/acf_NP0800[0], '+-', color=colP[2], label = r'$\nu_0=2.0/R_0^3$, $\lambda(t^{(1)})$ = %4.3f $\tau_0$, $\lambda(t^{(2)})$ = %4.3f $\tau_0$'%(-1./a_800, -1./a2_800))

plt.semilogy(t_NP1000, acf_NP1000/acf_NP1000[0], '+-', color=colP[3], label = r'$\nu_0=2.5/R_0^3$, $\lambda(t^{(1)})$ = %4.3f $\tau_0$, $\lambda(t^{(2)})$ = %4.3f $\tau_0$'%(-1./a_1000, -1./a2_1000))
# plt.semilogy(t_NP1200, acf_NP1200/acf_NP1200[0], '+-', color=colP[4], label = r'NP=1200, $\lambda(t^{(1)})$ = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$'%(t_NP1200[N1], -1./a_1200, t_NP1200[N2_st], t_NP1200[N2], -1./a2_1200))
plt.semilogy(t_NP1400, acf_NP1400/acf_NP1400[0], '+-', color=colP[5], label = r'$\nu_0=3.5/R_0^3$, $\lambda(t^{(1)})$ = %4.3f $\tau_0$, $\lambda(t^{(2)})$ = %4.3f $\tau_0$'%( -1./a_1400, -1./a2_1400))

plt.annotate(r'$t^{(1)}\in$(0, %3.2f)'%t_NP0400[N1], xy=(0.01, 0.06))
plt.annotate(r'$t^{(2)}\in$(%3.1f, %3.1f)'%(t_NP0400[N2_st], t_NP0400[N2]), xy=(0.15, 0.05))
# plt.annotate(r'$t^2_1$=%3.1f'%t_NP0400[N2], xy=(0.35, 0.05))
# plt.semilogy(t_NP0400, acf_NP0400/acf_NP0400[0], '-', color=colP[0], label = r'NP=400, tc(0,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$'%(t_NP0400[N1], -1./a_400, t_NP0400[N2_st], t_NP0400[N2], -1./a2_400, t_NP0400[N3_st], t_NP0400[N3], -1./a3_400))
# plt.semilogy(t_NP0600, acf_NP0600/acf_NP0600[0], '-', color=colP[1], label = r'NP=600, tc(0,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$'%(t_NP0600[N1], -1./a_600, t_NP0600[N2_st], t_NP0600[N2], -1./a2_600, t_NP0600[N3_st], t_NP0600[N3], -1./a3_600))
# plt.semilogy(t_NP0800, acf_NP0800/acf_NP0800[0], '-', color=colP[2], label = r'NP=800, tc(0,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$'%(t_NP0800[N1], -1./a_800, t_NP0800[N2_st], t_NP0800[N2], -1./a2_800, t_NP0800[N3_st], t_NP0800[N3], -1./a3_800))

# plt.semilogy(t_NP1000, acf_NP1000/acf_NP1000[0], '-', color=colP[3], label = r'NP=1000, tc(0,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$'%(t_NP1000[N1], -1./a_1000, t_NP1000[N2_st], t_NP1000[N2], -1./a2_1000, t_NP1000[N4_st], t_NP1000[N4], -1./a3_1000))
# plt.semilogy(t_NP1200, acf_NP1200/acf_NP1200[0], '-', color=colP[4], label = r'NP=1200, tc(0,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$'%(t_NP1200[N1], -1./a_1200, t_NP1200[N2_st], t_NP1200[N2], -1./a2_1200, t_NP1200[N3_st], t_NP1200[N3], -1./a3_1200))
# plt.semilogy(t_NP1400, acf_NP1400/acf_NP1400[0], '-', color=colP[5], label = r'NP=1400, tc(0,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$'%(t_NP1400[N1], -1./a_1400, t_NP1400[N2_st], t_NP1400[N2], -1./a2_1400, t_NP1400[N4_st], t_NP1400[N4], -1./a3_1400))


# plt.semilogy(t_NP2000, acf_NP2000/acf_NP2000[0], '-', color=colP[6], label = r'NP=2000, tc(0,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$'%(t_NP2000[N1], -1./a_2000, t_NP2000[N2_st], t_NP2000[N2], -1./a2_2000, t_NP2000[N4_st], t_NP2000[N4], -1./a3_2000))
# plt.semilogy(t_NP4000, acf_NP4000/acf_NP4000[0], '-', color=colP[7], label = r'NP=4000, tc(0,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$'%(t_NP4000[N1], -1./a_4000, t_NP4000[N2_st], t_NP4000[N2], -1./a2_4000, t_NP4000[N4_st], t_NP4000[N4], -1./a3_4000))



# plt.semilogy(t_NP0400, acf_NP0400/acf_NP0400[0], 'b-', label = r'NP=400, $\delta t$ = 0.01 $\tau_0$, tc(0,0.01) = %4.3f $\tau_0$, tc(0.1,0.2) = %4.3f $\tau_0$, tc(0.2,0.3) = %4.3f $\tau_0$'%(-1./a_400, -1./a2_400, -1./a3_400))
# plt.semilogy(t_NP0600, acf_NP0600/acf_NP0600[0], 'g-', label = r'NP=600, $\delta t$ = 0.01 $\tau_0$, tc(0,0.01) = %4.3f $\tau_0$, tc(0.1,0.2) = %4.3f $\tau_0$, tc(0.2,0.3) = %4.3f $\tau_0$'%(-1./a_600, -1./a2_600, -1./a3_600))
# plt.semilogy(t_NP0700, acf_NP0700/acf_NP0700[0], '-', color='brown', label = r'NP=700, $\delta t$ = 0.01 $\tau_0$, tc(0,0.01) = %4.3f $\tau_0$, tc(0.1,0.2) = %4.3f $\tau_0$, tc(0.2,0.3) = %4.3f $\tau_0$'%(-1./a_700, -1./a2_700, -1./a3_700))
# plt.semilogy(t_NP0800, acf_NP0800/acf_NP0800[0], 'r-', label = r'NP=800, $\delta t$ = 0.01 $\tau_0$, tc(0,0.01) = %4.3f $\tau_0$, tc(0.1,0.2) = %4.3f $\tau_0$, tc(0.2,0.3) = %4.3f $\tau_0$'%(-1./a_800, -1./a2_800, -1./a3_800))

plt.semilogy(ref_zero[:,0], ref_zero[:,1], 'k-')

ref_vert_1 = asarray([[t_NP0400[N1], 0.001], [t_NP0400[N1], 1.2]])
ref_vert_2 = asarray([[t_NP0400[N2_st], 0.001], [t_NP0400[N2_st], 1.2]])
ref_vert_3 = asarray([[t_NP0400[N2], 0.001], [t_NP0400[N2], 1.2]])
# ref_vert_4 = asarray([[1, 0.001], [1, 1.2]])
# ref_vert_5 = asarray([[0.27, 0.05], [0.27, 1.2]])
plt.semilogy(ref_vert_1[:,0], ref_vert_1[:,1], 'r-', linewidth=2, alpha=0.3)
plt.semilogy(ref_vert_2[:,0], ref_vert_2[:,1], 'b-', linewidth=2, alpha=0.3)
plt.semilogy(ref_vert_3[:,0], ref_vert_3[:,1], 'b-', linewidth=2, alpha=0.3)
# plt.semilogy(ref_vert_4[:,0], ref_vert_4[:,1], 'k-')
# plt.semilogy(ref_vert_5[:,0], ref_vert_5[:,1], 'b-')

plt.xlabel(r'topological time / $\tau_0$')
plt.ylabel('normalized autocorrelation function')
plt.legend(loc = 'upper right', prop=fontP, numpoints=1)
plt.grid()
plt.axis([0, 2.5, 0.05, 1.0])

# plt.axis([0, 10.0, -0.2, 1.0])


# ax2 = plt.axes([.2, .5, .55, .35])
# ax2.semilogy(t_NP0400, acf_NP0400/acf_NP0400[0], '-', color=colP[0])
# ax2.semilogy(t_NP0600, acf_NP0600/acf_NP0600[0], '-', color=colP[1])
# ax2.semilogy(t_NP0800, acf_NP0800/acf_NP0800[0], '-', color=colP[2])
# ax2.semilogy(t_NP1000, acf_NP1000/acf_NP1000[0], '-', color=colP[3])
# ax2.annotate('zoom in initial part', xy=(1.4, 0.9), size='large')
# ax2.semilogy(ref_zero[:,0], ref_zero[:,1], 'k-')
# ax2.axis([0, 2, -0.05, 1.0])
# ax2.grid()
plt.show()

