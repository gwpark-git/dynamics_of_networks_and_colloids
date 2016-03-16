

# from numpy import *
# import matplotlib.pyplot as plt
# from acf_fcn import *

# dat = loadtxt('RF_NP0400_NC25_longer.dat')[:,1]
# Nt = shape(dat)[0]

def smooth(dat):
    re = zeros(size(dat))
    re[0] = dat[0]
    re[-1] = dat[-1]
    for i in range(1, size(dat) - 1):
        re[i] = mean(dat[i-1:i+2])
    return re

t = arange(Nt)


# dat_s1 = smooth(dat)
# dat_s2 = smooth(dat_s1)

# acf_dat = acf_gro(dat)
# acf_dat_s1 = acf_gro(dat_s1)
# acf_dat_s2 = acf_gro(dat_s2)


# # dat_s = zeros(Nt)
# # dat_s[0] = dat[0]
# # dat_s[-1] = dat[-1]
# # for i in range(1, Nt-1):
# #     dat_s[i] = mean(dat[i-1:i+2])
    

# plt.clf()
# plt.ion()
# # plt.plot(t, dat, 'b-', alpha = 0.3)
# # plt.plot(t, dat_s1, 'r-', alpha = 0.3)
# # plt.plot(t, dat_s2, 'k-')

# plt.plot(acf_dat, 'b-', alpha = 0.3)
# plt.plot(acf_dat_s1, 'r-', alpha = 0.3)
# plt.plot(acf_dat_s2, 'k-')
# plt.grid()
# plt.show()


from scipy.fftpack import *
# sig = zeros([Nt, 2])
# sig[:,0] = t
# sig[:,1] = dat

fft_dat = rfft(dat)
w = rfftfreq(size(dat), d=1.)
N_cut = 22000
fft_dat[-N_cut:] = zeros(N_cut)
# fft_dat[:N_cut] = zeros(N_cut)
inv_dat = irfft(fft_dat)
acf_invdat_22000 = acf_gro(inv_dat)


fft_dat = rfft(dat)
w = rfftfreq(size(dat), d=1.)
N_cut = 23000
fft_dat[-N_cut:] = zeros(N_cut)
# fft_dat[:N_cut] = zeros(N_cut)
inv_dat = irfft(fft_dat)
acf_invdat_23000 = acf_gro(inv_dat)

fft_dat = rfft(dat)
w = rfftfreq(size(dat), d=1.)
N_cut = 24000
fft_dat[-N_cut:] = zeros(N_cut)
# fft_dat[:N_cut] = zeros(N_cut)
inv_dat = irfft(fft_dat)
acf_invdat_24000 = acf_gro(inv_dat)

fft_dat = rfft(dat)
w = rfftfreq(size(dat), d=1.)
N_cut = 21000
fft_dat[-N_cut:] = zeros(N_cut)
# fft_dat[:N_cut] = zeros(N_cut)
inv_dat = irfft(fft_dat)
acf_invdat_21000 = acf_gro(inv_dat)


# plt.clf()
# plt.ion()
# # plt.plot(w, fft_dat, 'b-')
# plt.plot(t, dat, 'b-', label = 'given data')
# plt.plot(t, inv_dat, 'r-', linewidth=2, label = 'filtered in Fourier domain')
# plt.axis([5000, 7000, -100, 150])
# plt.legend(loc = 'upper right')
# plt.xlabel('time index (1 = topological update)')
# plt.ylabel('shear stress function')
# # plt.plot(log10(w), log10(abs(fft_dat)), 'b-')
# plt.grid()
# plt.show()

# acf_invdat = acf_gro(inv_dat)
