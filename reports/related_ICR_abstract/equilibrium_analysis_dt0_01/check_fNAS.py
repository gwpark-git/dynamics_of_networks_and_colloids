

# from numpy import *
# import matplotlib.pyplot as plt

# # data before modify time step
# # mean_fNAS_dt01 = []
# # fNAS0400 = 100.*loadtxt('move_path/NP0400_LD10P3_RT100_block100/NP0400_LD10P3_C100.ener')[:,4]/(400.*25.)
# # t_fNAS0400 = arange(size(fNAS0400))/1000.
# # mean_fNAS0400 = mean(fNAS0400[size(fNAS0400)/2:])
# # mean_fNAS.append(mean_fNAS0400)

# # fNAS0600 = 100.*loadtxt('move_path/NP0600_LD10P3_RT100_NB100/NP0600_LD10P3_C100.ener')[:,4]/(600.*25.)
# # t_fNAS0600 = arange(size(fNAS0600))/1000.
# # mean_fNAS0600 = mean(fNAS0600[size(fNAS0600)/2:])
# # mean_fNAS.append(mean_fNAS0600)

# # fNAS0800 = 100.*loadtxt('move_path/NP0800_LD10P3_RT100_block100/NP0800_LD10P3_C100.ener')[:,4]/(800.*25.)
# # t_fNAS0800 = arange(size(fNAS0800))/1000.
# # mean_fNAS0800 = mean(fNAS0800[size(fNAS0800)/2:])
# # mean_fNAS.append(mean_fNAS0800)

# # fNAS1000 = 100.*loadtxt('../longer_computation/NP1000_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP1000_LD10P3_C100.ener')[:,4]/(1000.*25.)
# # t_fNAS1000 = arange(size(fNAS1000))/1000.
# # mean_fNAS1000 = mean(fNAS1000[size(fNAS1000)/2:])
# # mean_fNAS.append(mean_fNAS1000)

# # fNAS1200 = 100.*loadtxt('../longer_computation/NP1200_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP1200_LD10P3_C100.ener')[:,4]/(1200.*25.)
# # t_fNAS1200 = arange(size(fNAS1200))/1000.
# # mean_fNAS1200 = mean(fNAS1200[size(fNAS1200)/2:])
# # mean_fNAS.append(mean_fNAS1200)

# # fNAS1000 = 100.*loadtxt('../longer_computation/NP1000_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP1000_LD10P3_C100.ener')[:,4]/(1000.*25.)
# # t_fNAS1000 = arange(size(fNAS1000))/1000.
# # mean_fNAS1000 = mean(fNAS1000[size(fNAS1000)/2:])
# # mean_fNAS.append(mean_fNAS1000)

# data after modify time step
mean_fNAS_dt01 = []

move_path='/Volumes/Task_REPO_MAC_PRO/works/Brownian_simulation/data/stochastic_simulation/density_dependence'

fNAS0400_dt01 = 100.*loadtxt('%s/NP0400_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP0400_LD10P3_C100.ener'%move_path)[:,4]/(400.*25.)
t_fNAS0400_dt01 = arange(size(fNAS0400_dt01))/100.
mean_fNAS0400_dt01 = mean(fNAS0400_dt01[size(fNAS0400_dt01)/2:])
mean_fNAS_dt01.append(mean_fNAS0400_dt01)

fNAS0600_dt01 = 100.*loadtxt('%s/NP0600_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP0600_LD10P3_C100.ener'%move_path)[:,4]/(600.*25.)
t_fNAS0600_dt01 = arange(size(fNAS0600_dt01))/100.
mean_fNAS0600_dt01 = mean(fNAS0600_dt01[size(fNAS0600_dt01)/2:])
mean_fNAS_dt01.append(mean_fNAS0600_dt01)

fNAS0800_dt01 = 100.*loadtxt('%s/NP0800_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP0800_LD10P3_C100.ener'%move_path)[:,4]/(800.*25.)
t_fNAS0800_dt01 = arange(size(fNAS0800_dt01))/100.
mean_fNAS0800_dt01 = mean(fNAS0800_dt01[size(fNAS0800_dt01)/2:])
mean_fNAS_dt01.append(mean_fNAS0800_dt01)

fNAS1000_dt01 = 100.*loadtxt('../longer_computation/NP1000_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP1000_LD10P3_C100.ener')[:,4]/(1000.*25.)
t_fNAS1000_dt01 = arange(size(fNAS1000_dt01))/100.
mean_fNAS1000_dt01 = mean(fNAS1000_dt01[size(fNAS1000_dt01)/2:])
mean_fNAS_dt01.append(mean_fNAS1000_dt01)

fNAS1200_dt01 = 100.*loadtxt('../longer_computation/NP1200_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP1200_LD10P3_C100.ener')[:,4]/(1200.*25.)
t_fNAS1200_dt01 = arange(size(fNAS1200_dt01))/100.
mean_fNAS1200_dt01 = mean(fNAS1200_dt01[size(fNAS1200_dt01)/2:])
mean_fNAS_dt01.append(mean_fNAS1200_dt01)

fNAS1400_dt01 = 100.*loadtxt('../longer_computation/NP1400_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP1400_LD10P3_C100.ener')[:,4]/(1400.*25.)
t_fNAS1400_dt01 = arange(size(fNAS1400_dt01))/100.
mean_fNAS1400_dt01 = mean(fNAS1400_dt01[size(fNAS1400_dt01)/2:])
mean_fNAS_dt01.append(mean_fNAS1400_dt01)

fNAS2000_dt01 = 100.*loadtxt('../longer_computation/NP2000_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP2000_LD10P3_C100.ener')[:,4]/(2000.*25.)
t_fNAS2000_dt01 = arange(size(fNAS2000_dt01))/100.
mean_fNAS2000_dt01 = mean(fNAS2000_dt01[size(fNAS2000_dt01)/2:])
mean_fNAS_dt01.append(mean_fNAS2000_dt01)

fNAS4000_dt01 = 100.*loadtxt('../longer_computation/NP4000_LD10P3_RT100_NB10_dt0_01_NSKIP100/NP4000_LD10P3_C100.ener')[:,4]/(4000.*25.)
t_fNAS4000_dt01 = arange(size(fNAS4000_dt01))/100.
mean_fNAS4000_dt01 = mean(fNAS4000_dt01[size(fNAS4000_dt01)/2:])
mean_fNAS_dt01.append(mean_fNAS4000_dt01)



# colP = ['blue', 'green', 'brown', 'red', 'cyan', 'purple', 'gray', 'black']
# Np = [400, 600, 800, 1000, 1200, 1400]
# plt.close()
# plt.figure(figsize=(6,8))
# plt.ion()
# plt.plot(t_fNAS0400_dt01, fNAS0400_dt01, '-', color=colP[0], label = 'NP=400, mean=%3.2f(%%)'%(mean_fNAS0400_dt01))
# plt.plot(t_fNAS0600_dt01, fNAS0600_dt01, '-', color=colP[1], label = 'NP=600, mean=%3.2f(%%)'%(mean_fNAS0600_dt01))
# plt.plot(t_fNAS0800_dt01, fNAS0800_dt01, '-', color=colP[2], label = 'NP=800, mean=%3.2f(%%)'%(mean_fNAS0800_dt01))
# plt.plot(t_fNAS1000_dt01, fNAS1000_dt01, '-', color=colP[3], label = 'NP=1000, mean=%3.2f(%%)'%(mean_fNAS1000_dt01))
# plt.plot(t_fNAS1200_dt01, fNAS1200_dt01, '-', color=colP[4], label = 'NP=1200, mean=%3.2f(%%)'%(mean_fNAS1200_dt01))
# plt.plot(t_fNAS1400_dt01, fNAS1400_dt01, '-', color=colP[5], label = 'NP=1400, mean=%3.2f(%%)'%(mean_fNAS1400_dt01))
# plt.plot(t_fNAS2000_dt01, fNAS2000_dt01, '-', color=colP[6], label = 'NP=2000, mean=%3.2f(%%)'%(mean_fNAS2000_dt01))
# plt.plot(t_fNAS4000_dt01, fNAS4000_dt01, '-', color=colP[7], label = 'NP=4000, mean=%3.2f(%%)'%(mean_fNAS4000_dt01))

# plt.grid()
# plt.yticks(mean_fNAS_dt01, ['%3.2f' % val for val in mean_fNAS_dt01])
# # plt.xticks(Np)
# plt.axis([0, 50, 0, 35])

# from matplotlib.font_manager import FontProperties
# fontP = FontProperties()
# fontP.set_size('small')
# plt.legend(loc = 'upper right', ncol=1, prop=fontP)
# plt.xlabel('topological dimensionless time')
# plt.ylabel('fraction of NAS (%)')
# plt.show()


NP = [400, 600, 800, 1000, 1200, 1400, 2000, 4000]

NP_dat=NP[:-2]
mean_fNAS_dt01_dat = mean_fNAS_dt01[:-2]
# mean_fNDAS_dt01_dat = asarange([])
plt.close()
plt.figure(figsize=(6,8))
plt.ion()
plt.plot(NP_dat, mean_fNAS_dt01_dat, 'b-')
for i in range(size(NP_dat)):
    plt.plot(NP_dat[i], mean_fNAS_dt01_dat[i], 'o', markersize=6, markerfacecolor=colP[i], label = 'NP_dat=%d'%(NP_dat[i]))
plt.legend(loc = 'upper left', numpoints=1)
plt.xlabel('Number of particles')
plt.ylabel('fraction of NAS (%)')
plt.xticks(NP_dat)
plt.yticks(mean_fNAS_dt01_dat)
plt.grid()
plt.axis([300, 1500, 3, 20])
plt.show()
