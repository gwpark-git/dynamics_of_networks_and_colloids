
# from numpy import *
# import matplotlib.pyplot as plt

# import sys
# sys.path.append('../../../post_processing')
# from acf_fcn import *
# normal_factor = 2.*1000.

# sf_NP0400 = loadtxt('RF_NP0400_NC25_RT100.dat')[1000:,1]
# t_sf_NP0400 = arange(size(sf_NP0400))/100.
# acf_NP0400 = acf_gro(sf_NP0400/normal_factor)
# t_acf_NP0400 = arange(size(acf_NP0400))/100.

# sf_NP0600 = loadtxt('RF_NP0600_NC25_RT100.dat')[1000:,1]
# t_sf_NP0600 = arange(size(sf_NP0600))/100.
# acf_NP0600 = acf_gro(sf_NP0600/normal_factor)
# t_acf_NP0600 = arange(size(acf_NP0600))/100.

# sf_NP0800 = loadtxt('RF_NP0800_NC25_RT100.dat')[1000:,1]
# t_sf_NP0800 = arange(size(sf_NP0800))/100.
# acf_NP0800 = acf_gro(sf_NP0800/normal_factor)
# t_acf_NP0800 = arange(size(acf_NP0800))/100.


# sf_NP1000 = loadtxt('../RF_NP1000_first_run.dat')[1000:,1]
# t_sf_NP1000 = arange(size(sf_NP1000))/100.
# normal_factor = 2.*1000.
# acf_NP1000 = acf_gro(sf_NP1000/normal_factor)
# t_acf_NP1000 = arange(size(acf_NP1000))/100.

# sf_NP1200 = loadtxt('../../longer_computation/virial_stress_longer/RF_NP1200_combined_4.dat')[1000:,1]
# t_sf_NP1200 = arange(size(sf_NP1200))/100.
# acf_NP1200 = acf_gro(sf_NP1200/normal_factor)
# t_acf_NP1200 = arange(size(acf_NP1200))/100.

# sf_NP1400 = loadtxt('RF/RF_NP1400_RT100_dt0_01_longer_combined.dat')[1000:,1]
# t_sf_NP1400 = arange(size(sf_NP1400))/100.
# acf_NP1400 = acf_gro(sf_NP1400/normal_factor)
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
plt.loglog(t_acf_NP0400[:100], acf_NP0400[:100], '-', linewidth=2, alpha=0.7, color=colP[0], label = r'$\nu_m$=%4.1f, $C_{\tilde{\tau}_{xy}}(0)=%4.1f\times 10^{%d}$'%(cond[0], 10**(log10(peak_info[0,3]) - int(log10(peak_info[0,3]) - 1)), int(log10(peak_info[0,3])) - 1))

plt.loglog(t_acf_NP0600[:120], acf_NP0600[:120], '-', linewidth=2, alpha=0.7, color=colP[1], label = r'$\nu_m$=%4.1f, $C_{\tilde{\tau}_{xy}}(0)=%4.1f\times 10^{%d}$'%(cond[1], 10**(log10(peak_info[1,3]) - int(log10(peak_info[1,3]) - 1)), int(log10(peak_info[1,3])) - 1))

plt.loglog(t_acf_NP0800[:160], acf_NP0800[:160], '-', linewidth=2, alpha=0.7, color=colP[2], label = r'$\nu_m$=%4.1f, $C_{\tilde{\tau}_{xy}}(0)=%4.1f\times 10^{%d}$'%(cond[2], 10**(log10(peak_info[2,3]) - int(log10(peak_info[2,3]) - 1)), int(log10(peak_info[2,3])) - 1))

plt.loglog(t_acf_NP1000[:430], acf_NP1000[:430], '-', linewidth=2, alpha=0.7, color=colP[3], label = r'$\nu_m$=%4.1f, $C_{\tilde{\tau}_{xy}}(0)=%4.1f\times 10^{%d}$'%(cond[3], 10**(log10(peak_info[3,3]) - int(log10(peak_info[3,3]) - 1)), int(log10(peak_info[3,3])) - 1))

plt.loglog(t_acf_NP1200[:700], acf_NP1200[:700], '-', linewidth=2, alpha=0.7, color=colP[4], label = r'$\nu_m$=%4.1f, $C_{\tilde{\tau}_{xy}}(0)=%4.1f\times 10^{%d}$'%(cond[4], 10**(log10(peak_info[4,3]) - int(log10(peak_info[4,3]) - 1)), int(log10(peak_info[4,3])) - 1))

plt.loglog(t_acf_NP1400[:1400], acf_NP1400[:1400], '-', linewidth=2, alpha=0.7, color=colP[5], label = r'$\nu_m$=%4.1f, $C_{\tilde{\tau}_{xy}}(0)=%4.1f\times 10^{%d}$'%(cond[5], 10**(log10(peak_info[5,3]) - int(log10(peak_info[5,3]) - 1)), int(log10(peak_info[5,3])) - 1))


# def get_text_scientifi\nu_motation(val):
#     return '%4.1f



# plt.xlabel(r'dimensionless time, $\beta_0 t$')
# plt.ylabel(r'normalized stress autocorrelation, $C_{\tilde{\tau}_{xy}}(t)/C_{\tilde{\tau}_{xy}}(0)$')
# # plt.axis([0, 20, 10**-8, 10**-2])
# plt.axis([0, 20, 10**-4, 2])
# from matplotlib.font_manager import FontProperties
# fontP = FontProperties()
# # fontP.set_size('small')
# plt.grid()

# plt.legend(loc = 'lower right', numpoints=1, ncol=1, prop=fontP)
# # plt.plot(t_sf_NP1000, sf_NP1000, 'b-')

# fig = plt.axes([.195, .15, .32, .4])
# # peak_info[:,0] = peak_info[:,0]**(2./3.)
# plt.loglog(cond[:], peak_info[:,3], 'k-')
# cnt=0
# plt.loglog(cond[cnt], peak_info[cnt,3], marker=symP[cnt], color=colP[cnt])
# cnt+=1
# plt.loglog(cond[cnt], peak_info[cnt,3], marker=symP[cnt], color=colP[cnt])
# cnt+=1 
# plt.loglog(cond[cnt], peak_info[cnt,3], marker=symP[cnt], color=colP[cnt])
# cnt+=1 
# plt.loglog(cond[cnt], peak_info[cnt,3], marker=symP[cnt], color=colP[cnt])
# cnt+=1 
# plt.loglog(cond[cnt], peak_info[cnt,3], marker=symP[cnt], color=colP[cnt])
# cnt+=1 
# plt.loglog(cond[cnt], peak_info[cnt,3], marker=symP[cnt], color=colP[cnt])
# plt.axis([0.38, 1.5, 3.*10**-4, 10**-2])
# plt.xticks(cond[:], ['%4.1f'% val for val in cond])
# plt.text(0.55, 3.5*10**-4, r'micelle number density, $\nu_m$')
# # plt.text(0.4, 2.*10**-3, r'$C_{\tilde{\tau}_{xy}}(0)$', rotation='vertical')
# plt.ylabel(r'$C_{\tilde{\tau}_{xy}}(0)$')
# # plt.xlabel(r'micelle number density, $\nu_m$')
# # plt.xlabel(r'$C_{\tilde{\tau}_{xy}}(0)$')
# # plt.ylabel(r'$\tau_{eff}/\tau_0$')
# plt.grid(axis='x')
# # plt.loglog(
# plt.show()
