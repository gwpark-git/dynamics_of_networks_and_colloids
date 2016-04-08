

# from numpy import *
# import matplotlib.pyplot as plt

# NAS_NP0400 = loadtxt('../NP0400_LD10P3_RT100_block100/NP0400_LD10P3_C100.ener')[:,4]
# NAS_NP0600 = loadtxt('../NP0600_LD10P3_RT100_NB100/NP0600_LD10P3_C100.ener')[:,4]
# NAS_NP0700 = loadtxt('../NP0700_LD10P3_RT100_block100/NP0700_LD10P3_C100.ener')[:,4]
# NAS_NP0800 = loadtxt('../NP0800_LD10P3_RT100_block100/NP0800_LD10P3_C100.ener')[:,4]

# t_NP0400 = arange(size(NAS_NP0400))/1000.
# t_NP0600 = arange(size(NAS_NP0600))/1000.
# t_NP0700 = arange(size(NAS_NP0700))/1000.
# t_NP0800 = arange(size(NAS_NP0800))/1000.

# mean_fNAS_NP0400 = mean(NAS_NP0400[size(NAS_NP0400)/2:])/(400.*25.)
# mean_fNAS_NP0600 = mean(NAS_NP0600[size(NAS_NP0600)/2:])/(600.*25.)
# mean_fNAS_NP0700 = mean(NAS_NP0700[size(NAS_NP0700)/2:])/(700.*25.)

# mean_fNAS_NP0800 = mean(NAS_NP0800[size(NAS_NP0800)/2:])/(800.*25.)

# mean_fNAS = [mean_fNAS_NP0400, mean_fNAS_NP0600, mean_fNAS_NP0700, mean_fNAS_NP0800]

plt.clf()
plt.ion()
plt.figure(figsize=(6,8))
plt.plot(t_NP0400, NAS_NP0400/(400.*25.), 'b-', label = 'Np=400')
plt.plot(t_NP0600, NAS_NP0600/(600.*25.), 'g-', label = 'Np=600')
plt.plot(t_NP0700, NAS_NP0700/(700.*25.), '-', color='brown', label = 'Np=700')
plt.plot(t_NP0800, NAS_NP0800/(800.*25.), 'r-', label = 'Np=800')
plt.legend(loc = 'lower right')
plt.yticks(mean_fNAS, ['%3.2f'% val for val in asarray(mean_fNAS)*100.])
plt.xlabel('topological dimensionless time')
plt.ylabel('fraction of NAS (%)')
plt.grid()
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
