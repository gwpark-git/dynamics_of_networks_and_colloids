
def find_val(dat, val):
    for i,x in enumerate(dat/dat[0]):
        if x < val:
            return i
    return i

# from numpy import *
# import matplotlib.pyplot as plt

# def corr_Gt(dat):
#     s_xx, s_xy, s_xz, s_yy, s_yz, s_zz = dat[:,0], dat[:,1], dat[:,2], dat[:, 4], dat[:,5], dat[:,8]
#     N_xy = s_xx - s_yy
#     N_xz = s_xx - s_zz
#     N_yz = s_yy - s_zz
#     return acf_gro(s_xy) + acf_gro(s_yz) + acf_gro(s_xz) + (1./6.)*(acf_gro(N_xy) + acf_gro(N_xz) + acf_gro(N_yz))

# import sys
# sys.path.append('../../../../post_processing')
# from acf_fcn import *
# normal_factor = 2.*1000.

# sf_NP0400 = loadtxt('RF_NP0400_NC25_RT100.dat')[1000:,:]
# t_sf_NP0400 = arange(size(sf_NP0400))/100.
# acf_NP0400 = corr_Gt(sf_NP0400/normal_factor)
# t_acf_NP0400 = arange(size(acf_NP0400))/100.

# sf_NP0600 = loadtxt('RF_NP0600_NC25_RT100.dat')[1000:,:]
# t_sf_NP0600 = arange(size(sf_NP0600))/100.
# acf_NP0600 = corr_Gt(sf_NP0600/normal_factor)
# t_acf_NP0600 = arange(size(acf_NP0600))/100.

# sf_NP0800 = loadtxt('RF_NP0800_NC25_RT100.dat')[1000:,:]
# t_sf_NP0800 = arange(size(sf_NP0800))/100.
# acf_NP0800 = corr_Gt(sf_NP0800/normal_factor)
# t_acf_NP0800 = arange(size(acf_NP0800))/100.


# sf_NP1000 = loadtxt('../RF_NP1000_first_run.dat')[1000:,:]
# t_sf_NP1000 = arange(size(sf_NP1000))/100.
# normal_factor = 2.*1000.
# acf_NP1000 = corr_Gt(sf_NP1000/normal_factor)
# t_acf_NP1000 = arange(size(acf_NP1000))/100.

# sf_NP1200 = loadtxt('../../longer_computation/virial_stress_longer/RF_NP1200_combined_4.dat')[1000:,:]
# t_sf_NP1200 = arange(size(sf_NP1200))/100.
# acf_NP1200 = corr_Gt(sf_NP1200/normal_factor)
# t_acf_NP1200 = arange(size(acf_NP1200))/100.

# sf_NP1400 = loadtxt('RF/RF_NP1400_RT100_dt0_01_longer_combined.dat')[1000:,:]
# t_sf_NP1400 = arange(size(sf_NP1400))/100.
# acf_NP1400 = corr_Gt(sf_NP1400/normal_factor)
# t_acf_NP1400 = arange(size(acf_NP1400))/100.

colP = ['blue', 'green', 'brown', 'red', 'cyan', 'purple']
symP = ['o', 's', 'D', '^', 'v', 'h']
cond = [0.4, 0.6, 0.8, 1.0, 1.2, 1.4]
peak_info = []
t_peak_NP0400, acf_peak_NP0400 = t_acf_NP0400[90:100][argmax(acf_NP0400[90:100])], max(acf_NP0400[90:100])
peak_info.append([10, t_peak_NP0400, acf_peak_NP0400, acf_NP0400[1]])
t_peak_NP0600, acf_peak_NP0600 = t_acf_NP0600[100:120][argmax(acf_NP0600[100:120])], max(acf_NP0600[100:120])
peak_info.append([15, t_peak_NP0600, acf_peak_NP0600, acf_NP0600[0]])
t_peak_NP0800, acf_peak_NP0800 = t_acf_NP0800[100:160][argmax(acf_NP0800[100:160])], max(acf_NP0800[100:160])
peak_info.append([20, t_peak_NP0800, acf_peak_NP0800, acf_NP0800[0]])
t_peak_NP1000, acf_peak_NP1000 = t_acf_NP1000[200:430][argmax(acf_NP1000[200:430])], max(acf_NP1000[200:430])
peak_info.append([25, t_peak_NP1000, acf_peak_NP1000, acf_NP1000[0]])
t_peak_NP1200, acf_peak_NP1200 = t_acf_NP1200[600:700][argmax(acf_NP1200[600:700])], max(acf_NP1200[600:700])
peak_info.append([30, t_peak_NP1200, acf_peak_NP1200, acf_NP1200[0]])
t_peak_NP1400, acf_peak_NP1400 = t_acf_NP1400[1200:1400][argmax(acf_NP1400[1200:1400])], max(acf_NP1400[1200:1400])
peak_info.append([35, t_peak_NP1400, acf_peak_NP1400, acf_NP1400[0]])
peak_info = asarray(peak_info)

plt.close()
plt.ion()
plt.figure(figsize=(11,6))
# plt.loglog(t_acf_NP0400[:100], acf_NP0400[:100]/acf_NP0400[0], '-', linewidth=2, alpha=0.7, color=colP[0], label = r'$\nu_m$=%4.1f, $\tilde{G}(0)$=%4.1e'%(cond[0], peak_info[0,3]))
plt.loglog(t_acf_NP0400[:100], acf_NP0400[:100]/acf_NP0400[0], '-', linewidth=2, alpha=0.7, color=colP[0], label = r'$\nu_m$=%4.1f'%(cond[0]))

plt.loglog(t_acf_NP0600[:120], acf_NP0600[:120]/acf_NP0600[0], '-', linewidth=2, alpha=0.7, color=colP[1], label = r'$\nu_m$=%4.1f'%(cond[1]))

plt.loglog(t_acf_NP0800[:160], acf_NP0800[:160]/acf_NP0800[0], '-', linewidth=2, alpha=0.7, color=colP[2], label = r'$\nu_m$=%4.1f'%(cond[2]))

plt.loglog(t_acf_NP1000[:430], acf_NP1000[:430]/acf_NP1000[0], '-', linewidth=2, alpha=0.7, color=colP[3], label = r'$\nu_m$=%4.1f'%(cond[3]))

plt.loglog(t_acf_NP1200[:700], acf_NP1200[:700]/acf_NP1200[0], '-', linewidth=2, alpha=0.7, color=colP[4], label = r'$\nu_m$=%4.1f'%(cond[4]))

plt.loglog(t_acf_NP1400[:1400], acf_NP1400[:1400]/acf_NP1400[0], '-', linewidth=2, alpha=0.7, color=colP[5], label = r'$\nu_m$=%4.1f'%(cond[5]))


# def get_text_scientifi\nu_motation(val):
#     return '%4.1f



plt.xlabel(r'dimensionless time, $\beta_0 t$')
plt.ylabel(r'normalized stress autocorrelation, $C_{\tilde{\tau}_{xy}}(t)/C_{\tilde{\tau}_{xy}}(0)$')
# plt.axis([0, 20, 10**-8, 10**-2])
plt.axis([0, 5, 10**-3, 2])
from matplotlib.font_manager import FontProperties
fontP = FontProperties()
# fontP.set_size('small')
plt.grid()

plt.legend(loc = 'lower right', numpoints=1, ncol=1, prop=fontP)
# # plt.loglog(t_sf_NP1000, sf_NP1000, 'b-')

fig = plt.axes([.195, .15, .32, .4])
# peak_info[:,0] = peak_info[:,0]**(2./3.)
peak_info[:, 3] /= peak_info[0, 3]
plt.loglog(cond[:], peak_info[:,3], 'k-')
cnt=0
plt.loglog(cond[cnt], peak_info[cnt,3], marker=symP[cnt], color=colP[cnt])
cnt+=1
plt.loglog(cond[cnt], peak_info[cnt,3], marker=symP[cnt], color=colP[cnt])
cnt+=1 
plt.loglog(cond[cnt], peak_info[cnt,3], marker=symP[cnt], color=colP[cnt])
cnt+=1 
plt.loglog(cond[cnt], peak_info[cnt,3], marker=symP[cnt], color=colP[cnt])
cnt+=1 
plt.loglog(cond[cnt], peak_info[cnt,3], marker=symP[cnt], color=colP[cnt])
cnt+=1 
plt.loglog(cond[cnt], peak_info[cnt,3], marker=symP[cnt], color=colP[cnt])
plt.axis([0.38, 1.5, 0.9*10**0, 1.1*10**1])
plt.xticks(cond[:], ['%4.1f'% val for val in cond])
# plt.xlabel(r'micelle number density, $\nu_m$', y=-0.5)
plt.text(0.55, 0.9, r'micelle number density, $\nu_m$')
# plt.text(0.4, 2.*10**-3, r'$C_{\tilde{\tau}_{xy}}(0)$', rotation='vertical')
plt.ylabel(r'$C_{\tilde{\tau}_{xy}}(0)/C^{(\nu_m=0.4)}_{\tilde{\tau}_{xy}}(0)$')
# plt.xlabel(r'micelle number density, $\nu_m$')
# plt.xlabel(r'$C_{\tilde{\tau}_{xy}}(0)$')
# plt.ylabel(r'$\tau_{eff}/\tau_0$')
plt.grid(axis='x')
# plt.loglog(
plt.show()


# dat_time = asarray([t_acf_NP0400[find_val(acf_NP0400, 0.1)], t_acf_NP0600[find_val(acf_NP0600, 0.1)], t_acf_NP0800[find_val(acf_NP0800, 0.1)], t_acf_NP1000[find_val(acf_NP1000, 0.1)], t_acf_NP1200[find_val(acf_NP1200, 0.1)], t_acf_NP1400[find_val(acf_NP1400, 0.1)]])
# plt.close()
# plt.loglog(cond, dat_time, 'bo-')
# plt.grid()
# plt.show()
