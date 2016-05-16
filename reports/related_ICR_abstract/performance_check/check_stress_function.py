
from numpy import *
import matplotlib.pyplot as plt

mean_sf = []

sf_NP1000_dt001_NB10 = loadtxt('../../equilibrium_analysis_dt0_01/RF_NP1000_RT100_dt0_01.dat')[:,1]
t_NP1000_dt001_NB10 = arange(size(sf_NP1000_dt001_NB10))/100.
mean_sf_NP1000_dt001_NB10 = mean(sf_NP1000_dt001_NB10[size(sf_NP1000_dt001_NB10)/2:])
mean_sf.append(mean_sf_NP1000_dt001_NB10)

sf_NP1000_dt01_NB1 = loadtxt('../RF_NP1000_dt0_1_NB1.dat')[:,1]
t_NP1000_dt01_NB1 = arange(size(sf_NP1000_dt01_NB1))/100.
mean_sf_NP1000_dt01_NB1 = mean(sf_NP1000_dt01_NB1[size(sf_NP1000_dt01_NB1)/2:])
mean_sf.append(mean_sf_NP1000_dt01_NB1)

sf_NP1000_dt01_NB10 = loadtxt('../RF_NP1000_dt0_1_NB10.dat')[:,1]
t_NP1000_dt01_NB10 = arange(size(sf_NP1000_dt01_NB10))/100.
mean_sf_NP1000_dt01_NB10 = mean(sf_NP1000_dt01_NB10[size(sf_NP1000_dt01_NB10)/2:])
mean_sf.append(mean_sf_NP1000_dt01_NB10)



colP = ['blue', 'green', 'red']

plt.close()
plt.ion()
plt.figure(figsize=(11,6))

plt.plot(t_NP1000_dt01_NB10, sf_NP1000_dt01_NB10 + 600, '-', color=colP[2], label = r'shifted by 600, $dt=10^{-1} \tau_B = 10^{-3} \tau_0$, $\delta t = \tau_B = 10^{-2} \tau_0$, mean=%3.2f'%(mean_sf_NP1000_dt01_NB10))
plt.plot(t_NP1000_dt01_NB1, sf_NP1000_dt01_NB1 + 300, '-', color=colP[1], label = r'shifted by 300, $dt=10^{-1} \tau_B = 10^{-3} \tau_0$, $\delta t = 10^{-1} \tau_B = 10^{-3} \tau_0$, mean=%3.2f'%(mean_sf_NP1000_dt01_NB1))

plt.plot(t_NP1000_dt001_NB10, sf_NP1000_dt001_NB10, '-', color=colP[0], label = r'shifted by 0, $dt=10^{-2} \tau_B = 10^{-4} \tau_0$, $\delta t = 10^{-1} \tau_B = 10^{-3} \tau_0$, mean=%3.2f'%(mean_sf_NP1000_dt001_NB10))


mean_sf_shifted = [val + 300*i for i,val in enumerate(mean_sf)]
text_y = []
for i in range(size(mean_sf)):
    text_y.append('%4.0f\n(%3.2f)'%(mean_sf_shifted[i], mean_sf[i]))
plt.yticks(mean_sf_shifted, text_y)
plt.grid()
from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('small')
plt.legend(loc = 'upper right', prop=fontP)
plt.xlabel(r'topological time / $\tau_0$')
plt.ylabel('shear stress function')
plt.axis([0, 100, -300, 1200])
plt.title('shear stress function')
plt.show()
