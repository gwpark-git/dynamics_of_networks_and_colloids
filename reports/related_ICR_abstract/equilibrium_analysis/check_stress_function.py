
# from numpy import *
# import matplotlib.pyplot as plt

# sf_NP0400 = loadtxt('RF_NP0400_NC25_longer_TMP.dat')[:,1]
# t_NP0400 = arange(size(sf_NP0400))/1000.
# mean_sf_NP0400 = mean(sf_NP0400[size(sf_NP0400)/2:])
# sf_NP0600 = loadtxt('RF_NP0600_NC25_longer_TMP.dat')[:,1]
# t_NP0600 = arange(size(sf_NP0600))/1000.
# mean_sf_NP0600 = mean(sf_NP0600[size(sf_NP0600)/2:])

# sf_NP0800 = loadtxt('RF_NP0800_NC25_longer_TMP2.dat')[:,1]
# t_NP0800 = arange(size(sf_NP0800))/1000.
# mean_sf_NP0800 = mean(sf_NP0800[size(sf_NP0800)/2:])

# mean_arr_sf = asarray([mean_sf_NP0400, mean_sf_NP0600, mean_sf_NP0800])
# mean_arr_sf_lab = asarray([-200, mean_sf_NP0400, mean_sf_NP0600 + 300, mean_sf_NP0800 + 600, 1000])

plt.clf()
plt.ion()
plt.figure(figsize=(11,6))
plt.plot(t_NP0400, sf_NP0400, 'b-', label = 'Np=400')
plt.plot(t_NP0600, sf_NP0600 + 300, 'g-', label = 'Np=600, shifted by 300')
plt.plot(t_NP0800, sf_NP0800 + 600, 'r-', label = 'Np=800, shifted by 600')
plt.yticks(mean_arr_sf_lab)
plt.legend(loc = 'upper right')
plt.xlabel('topological dimensionless time')
plt.ylabel('shear stress function (shifted)')
plt.grid()
plt.show()
