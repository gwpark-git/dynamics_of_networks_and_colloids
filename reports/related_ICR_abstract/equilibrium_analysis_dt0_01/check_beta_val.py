


from numpy import *
import matplotlib.pyplot as plt

dat = []

dist_NP1000 = loadtxt('dist_NP1000_RT100_new.dat')
beta_NP1000 = mean(dist_NP1000[:,1][-100:])
dat.append([1000, beta_NP1000, 1./beta_NP1000, 0.409, 2.687])
t_NP1000 = arange(shape(dist_NP1000)[0])/100.

dist_NP1200 = loadtxt('dist_NP1200_RT100_new.dat')
beta_NP1200 = mean(dist_NP1200[:,1][-100:])
dat.append([1200, beta_NP1200, 1./beta_NP1200, 0.422, 0.207])
t_NP1200 = arange(shape(dist_NP1200)[0])/100.

dist_NP1400 = loadtxt('dist_NP1400_RT100_new.dat')
beta_NP1400 = mean(dist_NP1400[:,1][-100:])
dat.append([1400, beta_NP1400, 1./beta_NP1400, 0.526, 1.114])
t_NP1400 = arange(shape(dist_NP1400)[0])/100.


dat = asarray(dat)

plt.close()
plt.figure(figsize=(6,8))
plt.ion()

colP = ['red', 'cyan', 'purple']
NP = [1000, 1200, 1400]
av_maximum = [mean(dist_NP1000[:,2]), mean(dist_NP1200[:,2]) , mean(dist_NP1400[:,2]) ]

plt.plot(NP, av_maximum, 'b-')
for i in range(size(NP)):
    plt.plot(NP[i], av_maximum[i], 'o', markersize=6, markerfacecolor=colP[i], label = 'Np=%d'%NP[i])
plt.xticks(NP)
plt.yticks(av_maximum)

# plt.plot(t_NP1000, dist_NP1000[:,2], '-', color=colP[0], label = 'NP=1000')
# plt.plot(t_NP1200, dist_NP1200[:,2] + 1, '-', color=colP[1], label = 'NP=1200, shifted by 1')
# plt.plot(t_NP1400, dist_NP1400[:,2] + 2, '-', color=colP[2], label = 'NP=1400, shifted by 2')
# plt.grid()
# av_maximum = [mean(dist_NP1000[:,2]), mean(dist_NP1200[:,2]) + 1, mean(dist_NP1400[:,2]) + 2]
# av_maximum_text = []
# for i in range(size(av_maximum)):
#     av_maximum_text.append('%3.2f\n(%3.2f)'%(av_maximum[i], av_maximum[i] - 1*i))

# from matplotlib.font_manager import FontProperties
# fontP = FontProperties()
# fontP.set_size('small')
# plt.yticks(av_maximum, av_maximum_text)
# plt.legend(loc = 'upper right', prop=fontP)
# plt.xlabel('topological dimensionless time')
# plt.ylabel('shifted maximum distance for bridges')
# plt.show()

