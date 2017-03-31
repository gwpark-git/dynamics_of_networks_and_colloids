
from numpy import *
from acf_fcn import *
import sys

def tauij_acf_con(dat):
    tau_12 = dat[:, 21]
    tau_13 = dat[:, 22]
    tau_23 = dat[:, 23]
    return (corr_FFT(tau_12, tau_12) + corr_FFT(tau_13, tau_13) + corr_FFT(tau_23, tau_23))/5.

def Nij_con(dat, i, j):
    con_st = 17
    tau_ii = dat[:, con_st + i]
    tau_jj = dat[:, con_st + j]
    N_ij = tau_ii - tau_jj
    return N_ij

def Nij_acf_con(dat):
    N_12 = Nij_con(dat, 1, 2)
    N_13 = Nij_con(dat, 1, 3)
    N_23 = Nij_con(dat, 2, 3)
    return (corr_FFT(N_12, N_12) + corr_FFT(N_13, N_13) + corr_FFT(N_23, N_23))/30.


if size(sys.argv) < 3:
    print 'Get stress autocorrelation from given .ener file for associating telechelic polymers'
    print 'USAGE :'
    print 'argv[1] == .ener file name'
    print 'argv[2] == output file name'
    print 'argv[3] (optional) == number of cutting variables (HALF will cut half). Default = 0.'
    print 'argv[4] (optional) == stride'
    print 'argv[?] (disabled) == number of divisor (average over).'
else:
    if (sys.argv[1] == sys.argv[2]):
        print 'ERR: input and output have the same file name'
    else:
        print 'processing stress autocorrelation from %s to %s'%(sys.argv[1], sys.argv[2])
        ener = loadtxt(sys.argv[1])
        Nt = shape(ener)[0]
        N_cut = 0
        if size(sys.argv) >= 4:
            if(sys.argv[3] == 'HALF'):
                N_cut = Nt/2
            else:
                N_cut = int(sys.argv[3])
        stride = 1
        if size(sys.argv) == 5:
            stride = int(sys.argv[4])
        print 'compute stress autocorrelation for number of time steps %ld out of %ld with stride %ld'%(N_cut, Nt, stride)                
        if size(sys.argv) > 5:
            N_div = int(sys.argv[5])
            print 'average over %d blocks'%(N_div)
            Nt_block = (size(ener[:, 0]) - N_cut)/2 # it will the raw data range
            N_st = N_cut
            Nt_inc = Nt_block/N_div # it will increment of data
            print '\taverage over %ld blocks:'%(N_div)
            print '\t### %ld out of %ld is being processing: from %ld to %ld'%(0, N_div, N_st, N_st + Nt_block)
            corr_con = tauij_acf_con(ener[N_st:N_st+Nt_block, :]) + Nij_acf_con(ener[N_st:N_st+Nt_block, :])
            for i in range(1, N_div):
                N_st = N_cut + i*Nt_inc
                print '\t### %ld out of %ld is being processing: from %ld to %ld'%(i, N_div, N_st, N_st + Nt_block)
                corr_con += tauij_acf_con(ener[N_st:N_st+Nt_block, :]) + Nij_acf_con(ener[N_st:N_st+Nt_block, :])
            corr_con /= float(N_div)
            # N_div = int(sys.argv[4])
            # print 'average over %d blocks'%(N_div)            
            # Nt_block = Nt/N_div

            # corr_con = corr(ener[N_cut:N_cut+Nt_block, 27], ener[N_cut:N_cut+Nt_block, 27])
            # corr_rep = corr(ener[N_cut:N_cut+Nt_block, 21], ener[N_cut:N_cut+Nt_block, 21])
            # corr_con_rep = corr(ener[N_cut:N_cut+Nt_block, 27], ener[N_cut:N_cut+Nt_block, 21])
            # corr_rep_con = corr(ener[N_cut:N_cut+Nt_block, 21], ener[N_cut:N_cut+Nt_block, 27])

            # for i in range(1, N_div):
            #     corr_con = corr(ener[N_cut + i*Nt_block:N_cut+(i+1)*Nt_block, 27], ener[N_cut + i*Nt_block:N_cut+(i+1)*Nt_block, 27])
            #     corr_rep = corr(ener[N_cut + i*Nt_block:N_cut+(i+1)*Nt_block, 21], ener[N_cut + i*Nt_block:N_cut+(i+1)*Nt_block, 21])
            #     corr_con_rep = corr(ener[N_cut + i*Nt_block:N_cut+(i+1)*Nt_block, 27], ener[N_cut + i*Nt_block:N_cut+(i+1)*Nt_block, 21])
            #     corr_rep_con = corr(ener[N_cut + i*Nt_block:N_cut+(i+1)*Nt_block, 21], ener[N_cut + i*Nt_block:N_cut+(i+1)*Nt_block, 27])

            # corr_con /= float(N_div)
            # corr_rep /= float(N_div)
            # corr_con_rep /= float(N_div)
            # corr_rep_con /= float(N_div)
            
        else:
            Nt = shape(ener)[0]
            dat_range = arange(N_cut, Nt, stride)
            corr_con = tauij_acf_con(ener[dat_range, :]) + Nij_acf_con(ener[dat_range, :])
            # corr_con = (corr(ener[dat_range,21], ener[dat_range,21]) + corr(ener[dat_range,22], ener[dat_range,22]) + corr(ener[dat_range,23], ener[dat_range,23]))/3.
        dat = zeros([size(corr_con), 2])
        dat[:,0] = ener[:size(corr_con),0] - ener[0, 0]
        dat[:,0] *= float(stride)
        dat[:,1] = corr_con
        # dat[:,1] = corr_con + corr_rep + corr_con_rep + corr_rep_con
        # dat[:,2] = corr_con
        # dat[:,3] = corr_rep
        # dat[:,4] = corr_con_rep
        # dat[:,5] = corr_rep_con
        savetxt(sys.argv[2], dat)
                


