
from numpy import *
import matplotlib.pyplot as plt
import sys
sys.path.append('../../../post_processing')
from acf_fcn import *

N_st = 1000 # it is 10 tau_0

# RF_NP0400 = loadtxt('RF_NP0400_NC25_RT100.dat')[N_st:,1]
# acf_NP0400 = acf_gro(RF_NP0400)
# t_NP0400 = arange(size(acf_NP0400))/100.

# RF_NP0600 = loadtxt('RF_NP0600_NC25_RT100.dat')[N_st:,1]
# acf_NP0600 = acf_gro(RF_NP0600)
# t_NP0600 = arange(size(acf_NP0600))/100.

# RF_NP0800 = loadtxt('RF_NP0800_NC25_RT100.dat')[N_st:,1]
# acf_NP0800 = acf_gro(RF_NP0800)
# t_NP0800 = arange(size(acf_NP0800))/100.

# RF_NP1000 = loadtxt('RF_NP1000_RT100_dt0_01.dat')[N_st:,1]
# acf_NP1000 = acf_gro(RF_NP1000)
# t_NP1000 = arange(size(acf_NP1000))/100.
# acf_NP1000_BAV = acf_gro_BAV(RF_NP1000, 100)
from scipy.stats import linregress
N1 = 5
# a_400, y_400, r_400, p_400, std_400 = linregress(t_NP0400[:N1], log(acf_NP0400[:N1]))
# a_600, y_600, r_600, p_600, std_600 = linregress(t_NP0600[:N1], log(acf_NP0600[:N1]))
# a_800, y_800, r_800, p_800, std_800 = linregress(t_NP0800[:N1], log(acf_NP0800[:N1]))
a_1000, y_1000, r_1000, p_1000, std_1000 = linregress(t_NP1000[:N1], log(acf_NP1000[:N1]))

N2_st = 30
N2 = 80
# a2_400, y2_400, r2_400, p2_400, std2_400 = linregress(t_NP0400[N2_st:N2], log(acf_NP0400[N2_st:N2]))
# a2_600, y2_600, r2_600, p2_600, std2_600 = linregress(t_NP0600[N2_st:N2], log(acf_NP0600[N2_st:N2]))
# a2_800, y2_800, r2_800, p2_800, std2_800 = linregress(t_NP0800[N2_st:N2], log(acf_NP0800[N2_st:N2]))
a2_1000, y2_1000, r2_1000, p2_1000, std2_1000 = linregress(t_NP1000[N2_st:N2], log(acf_NP1000[N2_st:N2]))


N3_st = 100
N3 = 250
# a3_400, y3_400, r3_400, p3_400, std3_400 = linregress(t_NP0400[N3_st:N3], log(acf_NP0400[N3_st:N3]))
# a3_600, y3_600, r3_600, p3_600, std3_600 = linregress(t_NP0600[N3_st:N3], log(acf_NP0600[N3_st:N3]))
# a3_800, y3_800, r3_800, p3_800, std3_800 = linregress(t_NP0800[N3_st:N3], log(acf_NP0800[N3_st:N3]))
a3_1000, y3_1000, r3_1000, p3_1000, std3_1000 = linregress(t_NP1000[N3_st:N3], log(acf_NP1000[N3_st:N3]))


from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('small')

colP = ['blue', 'green', 'brown', 'red']

ref_zero = asarray([[0, 0], [100, 0]])
plt.close()
plt.ion()
plt.figure(figsize=(11,6))
# plt.plot(t_NP0400, acf_NP0400/acf_NP0400[0], '-', color=colP[0], label = r'NP=400, tc(0,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$'%(t_NP0400[N1], -1./a_400, t_NP0400[N2_st], t_NP0400[N2], -1./a2_400, t_NP0400[N3_st], t_NP0400[N3], -1./a3_400))
# plt.plot(t_NP0600, acf_NP0600/acf_NP0600[0], '-', color=colP[1], label = r'NP=600, tc(0,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$'%(t_NP0600[N1], -1./a_600, t_NP0600[N2_st], t_NP0600[N2], -1./a2_600, t_NP0600[N3_st], t_NP0600[N3], -1./a3_600))
# plt.plot(t_NP0800, acf_NP0800/acf_NP0800[0], '-', color=colP[2], label = r'NP=800, tc(0,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$'%(t_NP0800[N1], -1./a_800, t_NP0800[N2_st], t_NP0800[N2], -1./a2_800, t_NP0800[N3_st], t_NP0800[N3], -1./a3_800))
# plt.plot(t_NP1000, acf_NP1000_BAV/acf_NP1000_BAV[0], '-', color=colP[3], label = r'NP=1000, tc(0,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$'%(t_NP1000[N1], -1./a_1000, t_NP1000[N2_st], t_NP1000[N2], -1./a2_1000, t_NP1000[N3_st], t_NP1000[N3], -1./a3_1000))
# plt.plot(t_NP1000, acf_NP1000/acf_NP1000[0], '-', color=colP[0], label = r'NP=1000, tc(0,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$, tc(%3.2f,%3.2f) = %4.3f $\tau_0$'%(t_NP1000[N1], -1./a_1000, t_NP1000[N2_st], t_NP1000[N2], -1./a2_1000, t_NP1000[N3_st], t_NP1000[N3], -1./a3_1000))

plt.semilogy(t_NP1000, acf_NP1000_BAV/acf_NP1000_BAV[0], '-', color=colP[3], label = r'NP=1000, $K\delta t = 1\tau_0$')
plt.semilogy(t_NP1000, acf_NP1000/acf_NP1000[0], '-', color=colP[0], label = r'NP=1000, $\delta t = 0.01\tau_0$')
plt.axis([0, 10, 0.001, 1.0])
# plt.semilogy(t_NP1000, acf_NP1000/acf_NP1000[0], 'b-')


# plt.semilogy(t_NP0400, acf_NP0400/acf_NP0400[0], 'b-', label = r'NP=400, $\delta t$ = 0.01 $\tau_0$, tc(0,0.01) = %4.3f $\tau_0$, tc(0.1,0.2) = %4.3f $\tau_0$, tc(0.2,0.3) = %4.3f $\tau_0$'%(-1./a_400, -1./a2_400, -1./a3_400))
# plt.semilogy(t_NP0600, acf_NP0600/acf_NP0600[0], 'g-', label = r'NP=600, $\delta t$ = 0.01 $\tau_0$, tc(0,0.01) = %4.3f $\tau_0$, tc(0.1,0.2) = %4.3f $\tau_0$, tc(0.2,0.3) = %4.3f $\tau_0$'%(-1./a_600, -1./a2_600, -1./a3_600))
# plt.semilogy(t_NP0700, acf_NP0700/acf_NP0700[0], '-', color='brown', label = r'NP=700, $\delta t$ = 0.01 $\tau_0$, tc(0,0.01) = %4.3f $\tau_0$, tc(0.1,0.2) = %4.3f $\tau_0$, tc(0.2,0.3) = %4.3f $\tau_0$'%(-1./a_700, -1./a2_700, -1./a3_700))
# plt.semilogy(t_NP0800, acf_NP0800/acf_NP0800[0], 'r-', label = r'NP=800, $\delta t$ = 0.01 $\tau_0$, tc(0,0.01) = %4.3f $\tau_0$, tc(0.1,0.2) = %4.3f $\tau_0$, tc(0.2,0.3) = %4.3f $\tau_0$'%(-1./a_800, -1./a2_800, -1./a3_800))

plt.semilogy(ref_zero[:,0], ref_zero[:,1], 'k-')

# ref_vert_1 = asarray([[0.1, 0.001], [0.1, 1.2]])
# ref_vert_2 = asarray([[0.2, 0.001], [0.2, 1.2]])
# ref_vert_3 = asarray([[0.5, 0.001], [0.5, 1.2]])
# ref_vert_4 = asarray([[1, 0.001], [1, 1.2]])
# # ref_vert_5 = asarray([[0.27, 0.05], [0.27, 1.2]])
# plt.semilogy(ref_vert_1[:,0], ref_vert_1[:,1], 'k-')
# plt.semilogy(ref_vert_2[:,0], ref_vert_2[:,1], 'k-')
# plt.semilogy(ref_vert_3[:,0], ref_vert_3[:,1], 'k-')
# plt.semilogy(ref_vert_4[:,0], ref_vert_4[:,1], 'k-')
# plt.semilogy(ref_vert_5[:,0], ref_vert_5[:,1], 'b-')

plt.xlabel('topological dimensionless time')
plt.ylabel('normalized autocorrelation function')
plt.legend(loc = 'upper right', prop=fontP)
plt.grid()
# plt.axis([0, 1.5, 0.001, 1.0])

# plt.axis([0, 25.0, -0.4, 1.2])


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

