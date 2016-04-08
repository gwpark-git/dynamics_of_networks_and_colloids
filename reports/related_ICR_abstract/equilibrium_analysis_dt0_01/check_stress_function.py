
from numpy import *
import matplotlib.pyplot as plt

mean_sf = []

sf_NP0400 = loadtxt('RF_NP0400_NC25_RT100.dat')[:,1]
t_NP0400 = arange(size(sf_NP0400))/100.
mean_sf_NP0400 = mean(sf_NP0400[size(sf_NP0400)/2:])
mean_sf.append(mean_sf_NP0400)

sf_NP0600 = loadtxt('RF_NP0600_NC25_RT100.dat')[:,1]
t_NP0600 = arange(size(sf_NP0600))/100.
mean_sf_NP0600 = mean(sf_NP0600[size(sf_NP0600)/2:])
mean_sf.append(mean_sf_NP0600)

sf_NP0800 = loadtxt('RF_NP0800_NC25_RT100.dat')[:,1]
t_NP0800 = arange(size(sf_NP0800))/100.
mean_sf_NP0800 = mean(sf_NP0800[size(sf_NP0800)/2:])
mean_sf.append(mean_sf_NP0800)

sf_NP1000 = loadtxt('RF_NP1000_RT100_dt0_01.dat')[:,1]
t_NP1000 = arange(size(sf_NP1000))/100.
mean_sf_NP1000 = mean(sf_NP1000[size(sf_NP1000)/2:])
mean_sf.append(mean_sf_NP1000)

sf_NP1200 = loadtxt('RF_NP1200_RT100_dt0_01.dat')[:,1]
t_NP1200 = arange(size(sf_NP1200))/100.
mean_sf_NP1200 = mean(sf_NP1200[size(sf_NP1200)/2:])
mean_sf.append(mean_sf_NP1200)

sf_NP1400 = loadtxt('RF_NP1400_RT100_dt0_01.dat')[:,1]
t_NP1400 = arange(size(sf_NP1400))/100.
mean_sf_NP1400 = mean(sf_NP1400[size(sf_NP1400)/2:])
mean_sf.append(mean_sf_NP1400)

# sf_NP2000 = loadtxt('RF_NP2000_RT100_dt0_01.dat')[:,1]
# t_NP2000 = arange(size(sf_NP2000))/100.
# mean_sf_NP2000 = mean(sf_NP2000[size(sf_NP2000)/2:])
# mean_sf.append(mean_sf_NP2000)

# sf_NP4000 = loadtxt('RF_NP4000_RT100_dt0_01.dat')[:,1]
# t_NP4000 = arange(size(sf_NP4000))/100.
# mean_sf_NP4000 = mean(sf_NP4000[size(sf_NP4000)/2:])
# mean_sf.append(mean_sf_NP4000)


colP = ['blue', 'green', 'brown', 'red', 'cyan', 'purple', 'gray', 'black']

plt.close()
plt.ion()
plt.figure(figsize=(11,6))

plt.plot(t_NP0400, sf_NP0400, '-', color=colP[0], label = 'Np=400, mean = %3.2f'%(mean_sf_NP0400))
plt.plot(t_NP0600, sf_NP0600 + 300, '-', color=colP[1], label = 'Np=600, shifted by 300, mean = %3.2f'%(mean_sf_NP0600))
plt.plot(t_NP0800, sf_NP0800 + 600, '-', color=colP[2], label = 'Np=800, shifted by 600, mean = %3.2f'%(mean_sf_NP0800))
plt.plot(t_NP1000, sf_NP1000 + 900, '-', color=colP[3], label = 'Np=1000, shifted by 900, mean = %3.2f'%(mean_sf_NP1000))
plt.plot(t_NP1200, sf_NP1200 + 1200, '-', color=colP[4], label = 'Np=1200, shifted by 1200, mean = %3.2f'%(mean_sf_NP1200))
plt.plot(t_NP1400, sf_NP1400 + 1500, '-', color=colP[5], label = 'Np=1400, shifted by 1500, mean = %3.2f'%(mean_sf_NP1400))


# plt.plot(t_NP1000, sf_NP1000, '-', color=colP[0], label = 'Np=400, mean = %3.2f'%(mean_sf_NP1000))
# plt.plot(t_NP1000, sf_NP1000 + 900, '-', color=colP[3], label = 'Np=1000, shifted by 900, mean = %3.2f'%(mean_sf_NP1000))

mean_sf_shifted = [val + 300*i for i,val in enumerate(mean_sf)]
text_y = []
for i in range(size(mean_sf)):
    text_y.append('%4.0f\n(%3.2f)'%(mean_sf_shifted[i], mean_sf[i]))
plt.yticks(mean_sf_shifted, text_y)
plt.grid()
from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('small')
plt.legend(loc = 'lower right', prop=fontP)
plt.xlabel('topological dimensionless time')
plt.ylabel('shear stress function')
plt.axis([0, 50, -300, 2000])
plt.title('shear stress function')
plt.show()
