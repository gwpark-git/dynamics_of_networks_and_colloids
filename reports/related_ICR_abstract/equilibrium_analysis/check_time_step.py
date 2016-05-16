

# from numpy import *
# import matplotlib.pyplot as plt

# NAS_NP0400 = loadtxt('../NP0400_LD10P3_RT100_block100/NP0400_LD10P3_C100.ener')[:,4]/(400.*25.)
# t_NP0400 = arange(size(NAS_NP0400))/1000.
mean_fNAS_NP0400 = mean(NAS_NP0400[size(NAS_NP0400)/2:])
# NAS_NP0800 = loadtxt('../NP0800_LD10P3_RT100_block100/NP0800_LD10P3_C100.ener')[:,4]/(800.*25.)
# t_NP0800 = arange(size(NAS_NP0800))/1000.
mean_fNAS_NP0800 = mean(NAS_NP0800[size(NAS_NP0800)/2:])

NAS_NP0400_dt = loadtxt('../NP0400_LD10P3_RT100_NB10_dt0_01/NP0400_LD10P3_C100.ener')[:,4]/(400.*25.)
t_NP0400_dt = arange(size(NAS_NP0400_dt))/100.
mean_fNAS_NP0400_dt = mean(NAS_NP0400_dt[size(NAS_NP0400_dt)/2:])
NAS_NP0800_dt = loadtxt('../NP0800_LD10P3_RT100_NB10_dt0_01/NP0800_LD10P3_C100.ener')[:,4]/(800.*25.)
t_NP0800_dt = arange(size(NAS_NP0800_dt))/100.
mean_fNAS_NP0800_dt = mean(NAS_NP0800_dt[size(NAS_NP0800_dt)/2:])

plt.clf()
plt.ion()
plt.figure(figsize=(11,6))
plt.plot(t_NP0400, NAS_NP0400, 'b-', linewidth=3, alpha=0.15, label = r'Np=400, dt=0.001, mean=%3.2f (%%)'%(mean_fNAS_NP0400*100.))
plt.plot(t_NP0800, NAS_NP0800, 'r-', linewidth=3, alpha=0.15, label = r'Np=800, dt=0.001, mean=%3.2f (%%)'%(mean_fNAS_NP0800*100.))
plt.plot(t_NP0400_dt, NAS_NP0400_dt, 'b-', label = r'Np=400, dt=0.01, mean=%3.2f (%%)'%(mean_fNAS_NP0400_dt*100.))
plt.plot(t_NP0800_dt, NAS_NP0800_dt, 'r-', label = r'Np=800, dt=0.01, mean=%3.2f (%%)'%(mean_fNAS_NP0800_dt*100.))


from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('small')
plt.legend(loc = 'lower right', ncol=2)
mean_fNAS = [0, mean_fNAS_NP0400, mean_fNAS_NP0800, 0.1]
plt.yticks(mean_fNAS, ['%3.2f'% val for val in asarray(mean_fNAS)*100.])
plt.xlabel('topological dimensionless time')
plt.ylabel('fraction of NAS (%)')
plt.grid()
plt.axis([0, 15, 0, 0.1])
plt.show()


# Np = [400, 600, 700, 800]

# plt.clf()
# plt.ion()
# plt.figure(figsize=(6,8))
# plt.plot(Np, mean_fNAS, 'bo-')
# plt.grid()
# plt.xticks(Np)
# plt.yticks(mean_fNAS, ['%4.3f'%val for val in asarray(mean_fNAS)*100.])
# plt.xlabel('number of particles')
# plt.ylabel('fraction of NAS (%)')
# plt.show()
