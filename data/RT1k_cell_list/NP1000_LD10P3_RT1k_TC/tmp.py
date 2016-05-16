

from numpy import *
import matplotlib.pyplot as plt
from scipy.linalg import norm

# Rvec = loadtxt('tracking_individual_chain/Rvec_000000.dat')

def plot_Rsca(Rvec, colP, label):
    Nt = shape(Rvec)[0]
    t = arange(Nt)/100.
    Rsca = zeros(Nt)
    for i in range(Nt):
        Rsca[i] = norm(Rvec[i, :])
    plt.plot(t, Rsca, '-', color=colP, alpha=0.8, label = label)

    
plt.close()
plt.ion()
plt.figure(figsize=(6,6))
colP = ['blue', 'cyan', 'green', 'brown', 'red']
cnt = 0
plot_Rsca(loadtxt('tracking_individual_chain/Rvec_%06d.dat'%cnt), colP[cnt/1000], label = 'chain index = %d'%cnt); cnt+=1000;
plot_Rsca(loadtxt('tracking_individual_chain/Rvec_%06d.dat'%cnt), colP[cnt/1000], label = 'chain index = %d'%cnt); cnt+=1000;
plot_Rsca(loadtxt('tracking_individual_chain/Rvec_%06d.dat'%cnt), colP[cnt/1000], label = 'chain index = %d'%cnt); cnt+=1000;
plot_Rsca(loadtxt('tracking_individual_chain/Rvec_%06d.dat'%cnt), colP[cnt/1000], label = 'chain index = %d'%cnt); cnt+=1000;
plot_Rsca(loadtxt('tracking_individual_chain/Rvec_%06d.dat'%cnt), colP[cnt/1000], label = 'chain index = %d'%cnt); cnt+=1000;

# plot_Rsca(loadtxt('tracking_individual_chain/Rvec_00000%d.dat'%cnt), colP[cnt]); cnt+=1;
# plot_Rsca(loadtxt('tracking_individual_chain/Rvec_00000%d.dat'%cnt), colP[cnt]); cnt+=1;
# plot_Rsca(loadtxt('tracking_individual_chain/Rvec_00000%d.dat'%cnt), colP[cnt]); cnt+=1;

cut_vert = asarray([[20, 0],
                    [20, 1.8]])
plt.plot(cut_vert[:,0], cut_vert[:,1], 'k--', linewidth=3)

plt.legend(loc = 'lower right')
plt.xlabel(r'dimensionless time, $\beta_0t$')
plt.ylabel(r'norm of relative vector, $|\mathbf{r}_{k}|$')
# plt.grid()
plt.show()
