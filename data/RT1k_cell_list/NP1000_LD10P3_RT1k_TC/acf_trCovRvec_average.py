
# from numpy import *
# import matplotlib.pyplot as plt

# Np      = 1000
# Nc      = 10

# N_tot   = Nc*Np

# t = arange(50*100)/100.
# acf_av_TrCovRvec = zeros(size(t))
# for k in range(N_tot):
#     acf_TrCovRvec_k = loadtxt('tracking_individual_chain/cross_corr_%06d.dat'%k)[:,1]
#     acf_av_TrCovRvec[:size(acf_TrCovRvec_k)] += acf_TrCovRvec_k
# acf_av_TrCovRvec /= float(N_tot)

# ref_fraction = asarray([[-1, 0.1041],
#                         [21, 0.1041]])

# from scipy.stats import linregress
# val_reg = linregress(t[:100], log(acf_av_TrCovRvec[:100]))
# y_reg = exp(val_reg[0]*t)

# plt.close()
# plt.ion()
# plt.figure(figsize=(6, 6))
# plt.plot(t, acf_av_TrCovRvec/acf_av_TrCovRvec[0], 'r-', linewidth=2, label = 'Np=1000')
# plt.plot(t, y_reg, 'r--', label = r'reg up to 1$\tau_0$, tc=%3.2f$\tau_0$'%(-1./val_reg[0]))
# # plt.plot(ref_fraction[:,0], ref_fraction[:,1], 'k:', linewidth=3, label = 'fNAS=10.41%')
# plt.fill_between(t, acf_av_TrCovRvec/acf_av_TrCovRvec[0], 0, color='red', alpha=0.2)
# plt.legend(loc = 'upper right')
# plt.grid()
# plt.xlabel(r'dimensionless time, $\beta_0 t$')
# plt.ylabel('averaged autocorrelation')
# plt.axis([-0.1, 3, -0.05, 1.2])
# plt.show()

# dat = zeros([1000, 2])
# dat[:,0] = t[:1000]
# dat[:,1] = acf_av_TrCovRvec[:1000]
# savetxt('acf_av_TrCovRvec_NP1000.dat', dat)
