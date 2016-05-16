

from numpy import *
import matplotlib.pyplot as plt

Np             = 1000
Nc             = 10
N_tot          = Nc*Np

# t = arange(100*100)/100.
# RR_flat = zeros([size(t), 9])
# for k in range(N_tot):
#     RR_k = loadtxt('virial_expression/RR_%06d.dat'%k)
#     RR_flat[:shape(RR_k)[0], :] += RR_k
# RxRy_chain = RR_flat[:,1]/max(RR_flat[:,1])

# RxRy_virial = loadtxt('RF_NP1000_LD10P3_C100.dat')[2000:,1]
# RxRy_virial /= max(RxRy_virial)
# t_virial = arange(size(RxRy_virial))/100.
# plt.close()
# plt.ion()
# plt.figure(figsize=(11,6))
# plt.plot(t, RxRy_chain, 'b.', label = 'shear stress (using the previous post-processing)')
# plt.plot(t, RxRy_chain, 'b-', linewidth=1, alpha=0.1)
# plt.plot(t_virial, RxRy_virial, 'r.', markersize=2, label = 'shear stress sum over all the chain')

# ref_unity = asarray([[0, 0],
#                      [40, 0]])
# plt.plot(ref_unity[:,0], ref_unity[:,1], 'k-')
# plt.axis([0, 40, -1, 1])
# plt.legend(loc = 'upper right', numpoints=1)
# plt.xlabel(r'dimensionless time, $\beta_0 t$')
# plt.ylabel('normalized shear stress')

# plt.grid()
# plt.show()

plt.close()
plt.ion()
plt.figure(figsize=(6,6))
plt.plot(RxRy_virial[:4000], RxRy_chain[:4000], 'b-')
plt.grid()
plt.axis([-1, 1, -1, 1])
plt.xlabel('normalized shear stress (original)')
plt.ylabel('normalized shear stress (sum over chains)')
plt.show()
