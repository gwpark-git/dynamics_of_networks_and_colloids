
from numpy import *
import matplotlib.pyplot as plt


def plot_acf_Rsca(acf_Rsca, colP, label):
    # Nt = shape(acf_Rsca)[0]
    # t = arange(Nt)/100.
    plt.plot(acf_Rsca[:,0], acf_Rsca[:,1], '-', color=colP, alpha=0.8, label = label+' C[0]=%3.2e'%(acf_Rsca[0,1]))

# dat = loadtxt('tracking_individual_chain/acf_Rsca_000000.dat')

# plt.close()
# plt.ion()
# plt.plot(dat[:,0], dat[:,1]/dat[0,1], 'b-')
# plt.grid()
# plt.show()
colP = ['blue', 'cyan', 'green', 'brown', 'red']

plt.close()
plt.ion()
plt.figure(figsize=(6,8))
cnt = 0
plot_acf_Rsca(loadtxt('tracking_individual_chain/acf_Rsca_%06d.dat'%cnt), colP[cnt/1000], label = 'chain index = %04d'%cnt); cnt+=1000;
plot_acf_Rsca(loadtxt('tracking_individual_chain/acf_Rsca_%06d.dat'%cnt), colP[cnt/1000], label = 'chain index = %04d'%cnt); cnt+=1000;
plot_acf_Rsca(loadtxt('tracking_individual_chain/acf_Rsca_%06d.dat'%cnt), colP[cnt/1000], label = 'chain index = %04d'%cnt); cnt+=1000;
plot_acf_Rsca(loadtxt('tracking_individual_chain/acf_Rsca_%06d.dat'%cnt), colP[cnt/1000], label = 'chain index = %04d'%cnt); cnt+=1000;
plot_acf_Rsca(loadtxt('tracking_individual_chain/acf_Rsca_%06d.dat'%cnt), colP[cnt/1000], label = 'chain index = %04d'%cnt); cnt+=1000;

plt.legend(loc = 'upper right')
plt.xlabel(r'dimensionless time, $\beta_0t$')
plt.ylabel(r'autocovariance for individual chain')
plt.axis([-1, 21, -0.05, 0.4])
plt.show()
