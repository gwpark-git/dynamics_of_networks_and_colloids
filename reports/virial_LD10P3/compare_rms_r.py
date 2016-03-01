
from numpy import *
import matplotlib.pyplot as plt

mean_Rrms_arr = []

Rrms_0400 = sqrt(loadtxt('dist_NP0400_NC25.dat')[:,0])
mean_Rrms_0400 = mean(Rrms_0400[size(Rrms_0400)/2:])
mean_Rrms_arr.append(mean_Rrms_0400)

Rrms_0500 = sqrt(loadtxt('dist_NP0500_NC25.dat')[:,0])
mean_Rrms_0500 = mean(Rrms_0500[size(Rrms_0500)/2:])
mean_Rrms_arr.append(mean_Rrms_0500)

Rrms_0600 = sqrt(loadtxt('dist_NP0600_NC25.dat')[:,0])
mean_Rrms_0600 = mean(Rrms_0600[size(Rrms_0600)/2:])
mean_Rrms_arr.append(mean_Rrms_0600)

Rrms_0700 = sqrt(loadtxt('dist_NP0700_NC25.dat')[:,0])
mean_Rrms_0700 = mean(Rrms_0700[size(Rrms_0700)/2:])
mean_Rrms_arr.append(mean_Rrms_0700)

plt.clf()
plt.ion()
plt.plot(Rrms_0400, 'b-', label = 'Np=400')
plt.plot(Rrms_0500 + 0.1, 'r-', label = 'Np=500, shifted')
plt.plot(Rrms_0600 + 0.2, 'g-', label = 'Np=600, shifted')
plt.plot(Rrms_0700 + 0.3, 'k-', label = 'Np=700, shifted')
plt.grid()
# plt.yticks(
plt.yticks([mean_Rrms_0400, mean_Rrms_0500 + 0.1, mean_Rrms_0600 + 0.2, mean_Rrms_0700 + 0.3], ['%.3f'%val for val in asarray(mean_Rrms_arr)])
plt.legend(loc = 'upper right')
plt.xlabel('time index (1 = topological update)')
plt.ylabel('average RMS for distance of bridges')
plt.axis([0, 10000, 0.96, 1.1 + 0.3])
plt.show()
