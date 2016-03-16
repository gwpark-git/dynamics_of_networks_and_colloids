
from numpy import *
import matplotlib.pyplot as plt

mean_beta_arr = []

beta_0400 = loadtxt('beta_NP0400_NC25.dat')
mean_beta_0400 = mean(beta_0400[size(beta_0400)/2:])
mean_beta_arr.append(mean_beta_0400)
beta_0500 = loadtxt('beta_NP0500_NC25.dat')
mean_beta_0500 = mean(beta_0500[size(beta_0500)/2:])
mean_beta_arr.append(mean_beta_0500)

beta_0600 = loadtxt('beta_NP0600_NC25.dat')
mean_beta_0600 = mean(beta_0600[size(beta_0600)/2:])
mean_beta_arr.append(mean_beta_0600)

beta_0700 = loadtxt('beta_NP0700_NC25.dat')
mean_beta_0700 = mean(beta_0700[size(beta_0700)/2:])
mean_beta_arr.append(mean_beta_0700)
plt.clf()
plt.ion()
plt.plot(beta_0400, 'b-', label = 'Np=400')
plt.plot(beta_0500, 'r-', label = 'Np=500')
plt.plot(beta_0600, 'g-', label = 'Np=600')
plt.plot(beta_0700, 'k-', label = 'Np=700')
plt.grid()
# plt.yticks(
plt.yticks(mean_beta_arr, ['%.2f'%val for val in asarray(mean_beta_arr)])
plt.legend(loc = 'lower right')
plt.xlabel('time index (1 = topological update)')
plt.ylabel('average detachment frequency')
plt.show()
