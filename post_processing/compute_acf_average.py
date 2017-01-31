
from numpy import *
from acf_fcn import *
import sys

def Nij_con(dat, i, j):
    con_st = 23
    tau_ii = dat[:, con_st + i]
    tau_jj = dat[:, con_st + j]
    N_ij = tau_ii - tau_jj
    return N_ij

def Nij_acf_con(dat):
    N_12 = Nij_con(dat, 1, 2)
    N_13 = Nij_con(dat, 1, 3)
    N_23 = Nij_con(dat, 2, 3)
    return (corr(N_12, N_12) + corr(N_13, N_13) + corr(N_23, N_23))/30.


def tauij_acf_con(dat):
    tau_12 = dat[:, 27]
    tau_13 = dat[:, 28]
    tau_23 = dat[:, 29]
    return (corr(tau_12, tau_12) + corr(tau_13, tau_13) + corr(tau_23, tau_23))/5.


def Nij_rep(dat, i, j):
    con_st = 17
    tau_ii = dat[:, con_st + i]
    tau_jj = dat[:, con_st + j]
    N_ij = tau_ii - tau_jj
    return N_ij

def Nij_acf_rep(dat):
    N_12 = Nij_rep(dat, 1, 2)
    N_13 = Nij_rep(dat, 1, 3)
    N_23 = Nij_rep(dat, 2, 3)
    return (corr(N_12, N_12) + corr(N_13, N_13) + corr(N_23, N_23))/30.


def tauij_acf_rep(dat):
    tau_12 = dat[:, 21]
    tau_13 = dat[:, 22]
    tau_23 = dat[:, 23]
    return (corr(tau_12, tau_12) + corr(tau_13, tau_13) + corr(tau_23, tau_23))/5.

def tauij_corr_con_rep(dat):
    tau_12_con = dat[:, 27]
    tau_13_con = dat[:, 28]
    tau_23_con = dat[:, 29]
    tau_12_rep = dat[:, 21]
    tau_13_rep = dat[:, 22]
    tau_23_rep = dat[:, 23]
    return (corr(tau_12_con, tau_12_rep) + corr(tau_13_con, tau_13_rep) + corr(tau_23_con, tau_23_rep))/5.

def Nij_corr_con_rep(dat):
    N_12_con = Nij_con(dat, 1, 2)
    N_13_con = Nij_con(dat, 1, 3)
    N_23_con = Nij_con(dat, 2, 3)
    N_12_rep = Nij_rep(dat, 1, 2)
    N_13_rep = Nij_rep(dat, 1, 3)
    N_23_rep = Nij_rep(dat, 2, 3)
    return (corr(N_12_con, N_12_rep) + corr(N_13_con, N_13_rep) + corr(N_23_con, N_23_rep))/30.

def tauij_corr_rep_con(dat):
    tau_12_con = dat[:, 27]
    tau_13_con = dat[:, 28]
    tau_23_con = dat[:, 29]
    tau_12_rep = dat[:, 21]
    tau_13_rep = dat[:, 22]
    tau_23_rep = dat[:, 23]
    return (corr(tau_12_rep, tau_12_con) + corr(tau_13_rep, tau_13_con) + corr(tau_23_rep, tau_23_con))/5.

def Nij_corr_rep_con(dat):
    N_12_con = Nij_con(dat, 1, 2)
    N_13_con = Nij_con(dat, 1, 3)
    N_23_con = Nij_con(dat, 2, 3)
    N_12_rep = Nij_rep(dat, 1, 2)
    N_13_rep = Nij_rep(dat, 1, 3)
    N_23_rep = Nij_rep(dat, 2, 3)
    return (corr(N_12_rep, N_12_con) + corr(N_13_rep, N_13_con) + corr(N_23_rep, N_23_con))/30.

if size(sys.argv) < 3:
    print 'Get stress autocorrelation from given .ener file for associating telechelic polymers'
    print 'USAGE :'
    print 'argv[1] == .ener file name'
    print 'argv[2] == output file name'
    print 'argv[3] (optional) == number of cutting variables (HALF will cut half). Default = 0.'
    print 'argv[4] (optional) == stride'
    print 'argv[5] (optional) == number of divisor (average over).'
    print '       averaging scheme is used the following'
    print '       |-----------|-----------| : data range is sectionized by half (say cut=HALF)'
    print '       |--|--|--|--|-----------| : initial data range when (say N_div=4)'
    print '       |-----------|--|--|--|--| : data coverage with each section'

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
            # |-----------|-----------| : data range is sectionized by half
            # |--|--|--|--|-----------| : initial data range when N_div=4
            # |-----------|--|--|--|--| : data coverage with each section
            # Nt_block = Nt/N_div

            print '\taverage over %ld blocks:'%(N_div)
            print '\t### %ld out of %ld is being processing'%(0, N_div)

            corr_con = tauij_acf_con(ener[N_st:N_st+Nt_block, :]) + Nij_acf_con(ener[N_st:N_st+Nt_block, :])
            corr_rep = tauij_acf_rep(ener[N_st:N_st+Nt_block, :]) + Nij_acf_rep(ener[N_st:N_st+Nt_block, :])
            corr_con_rep = tauij_corr_con_rep(ener[N_st:N_st+Nt_block, :]) + Nij_corr_con_rep(ener[N_st:N_st+Nt_block, :])
            corr_rep_con = tauij_corr_rep_con(ener[N_st:N_st+Nt_block, :]) + Nij_corr_rep_con(ener[N_st:N_st+Nt_block, :])
            
            # corr_con = (corr(ener[N_st:N_st+Nt_block, 27], ener[N_st:N_st+Nt_block, 27])
            #             + corr(ener[N_st:N_st+Nt_block, 28], ener[N_st:N_st+Nt_block, 28])
            #             + corr(ener[N_st:N_st+Nt_block, 29], ener[N_st:N_st+Nt_block, 29]))/3.
            # corr_rep = (corr(ener[N_st:N_st+Nt_block, 21], ener[N_st:N_st+Nt_block, 21])
            #             + corr(ener[N_st:N_st+Nt_block, 22], ener[N_st:N_st+Nt_block, 22])
            #             + corr(ener[N_st:N_st+Nt_block, 23], ener[N_st:N_st+Nt_block, 23]))/3.
            # corr_con_rep = (corr(ener[N_st:N_st+Nt_block, 27], ener[N_st:N_st+Nt_block, 21])
            #                 + corr(ener[N_st:N_st+Nt_block, 28], ener[N_st:N_st+Nt_block, 22])
            #                 + corr(ener[N_st:N_st+Nt_block, 29], ener[N_st:N_st+Nt_block, 23]))/3.
            # corr_rep_con = (corr(ener[N_st:N_st+Nt_block, 21], ener[N_st:N_st+Nt_block, 27])
            #                 + corr(ener[N_st:N_st+Nt_block, 22], ener[N_st:N_st+Nt_block, 28])
            #                 + corr(ener[N_st:N_st+Nt_block, 23], ener[N_st:N_st+Nt_block, 29]))/3.

            for i in range(1, N_div):
                print '\t### %ld out of %ld is being processing'%(i, N_div)                
                N_st = N_cut + i*Nt_inc
                corr_con += tauij_acf_con(ener[N_st:N_st+Nt_block, :]) + Nij_acf_con(ener[N_st:N_st+Nt_block, :])
                corr_rep += tauij_acf_rep(ener[N_st:N_st+Nt_block, :]) + Nij_acf_rep(ener[N_st:N_st+Nt_block, :])
                corr_con_rep += tauij_corr_con_rep(ener[N_st:N_st+Nt_block, :]) + Nij_corr_con_rep(ener[N_st:N_st+Nt_block, :])
                corr_rep_con += tauij_corr_rep_con(ener[N_st:N_st+Nt_block, :]) + Nij_corr_rep_con(ener[N_st:N_st+Nt_block, :])
                
                # corr_con += (corr(ener[N_st:N_st+Nt_block, 27], ener[N_st:N_st+Nt_block, 27])
                #              + corr(ener[N_st:N_st+Nt_block, 28], ener[N_st:N_st+Nt_block, 28])
                #              + corr(ener[N_st:N_st+Nt_block, 29], ener[N_st:N_st+Nt_block, 29]))/3.
                # corr_rep += (corr(ener[N_st:N_st+Nt_block, 21], ener[N_st:N_st+Nt_block, 21])
                #              + corr(ener[N_st:N_st+Nt_block, 22], ener[N_st:N_st+Nt_block, 22])
                #              + corr(ener[N_st:N_st+Nt_block, 23], ener[N_st:N_st+Nt_block, 23]))/3.
                # corr_con_rep += (corr(ener[N_st:N_st+Nt_block, 27], ener[N_st:N_st+Nt_block, 21])
                #                  + corr(ener[N_st:N_st+Nt_block, 28], ener[N_st:N_st+Nt_block, 22])
                #                  + corr(ener[N_st:N_st+Nt_block, 29], ener[N_st:N_st+Nt_block, 23]))/3.
                # corr_rep_con += (corr(ener[N_st:N_st+Nt_block, 21], ener[N_st:N_st+Nt_block, 27])
                #                  + corr(ener[N_st:N_st+Nt_block, 22], ener[N_st:N_st+Nt_block, 28])
                #                  + corr(ener[N_st:N_st+Nt_block, 23], ener[N_st:N_st+Nt_block, 29]))/3.

            corr_con /= float(N_div)
            corr_rep /= float(N_div)
            corr_con_rep /= float(N_div)
            corr_rep_con /= float(N_div)
            
        else:
            Nt = shape(ener)[0]
            dat_range = arange(N_cut, Nt, stride)

            corr_con = tauij_acf_con(ener[dat_range, :]) + Nij_acf_con(ener[dat_range, :])
            corr_rep = tauij_acf_rep(ener[dat_range, :]) + Nij_acf_rep(ener[dat_range, :])
            corr_con_rep = tauij_corr_con_rep(ener[dat_range, :]) + Nij_corr_con_rep(ener[dat_range, :])
            corr_rep_con = tauij_corr_rep_con(ener[dat_range, :]) + Nij_corr_rep_con(ener[dat_range, :])
            
            # corr_con = (corr(ener[dat_range,27], ener[dat_range,27]) + corr(ener[dat_range,28], ener[dat_range,28]) + corr(ener[dat_range,29], ener[dat_range,29]))/3.
            # corr_rep = (corr(ener[dat_range,21], ener[dat_range,21]) + corr(ener[dat_range,22], ener[dat_range,22]) + corr(ener[dat_range,23], ener[dat_range,23]))/3.
            # corr_con_rep = (corr(ener[dat_range,27], ener[dat_range,21]) + corr(ener[dat_range,28], ener[dat_range,22]) + corr(ener[dat_range,29], ener[dat_range,23]))/3.
            # corr_rep_con = (corr(ener[dat_range,21], ener[dat_range,27]) + corr(ener[dat_range,22], ener[dat_range,28]) + corr(ener[dat_range,23], ener[dat_range,29]))/3.
        dat = zeros([size(corr_con), 6])
        dat[:,0] = ener[:size(corr_con),0] - ener[0, 0]
        dat[:,0] *= float(stride)
        dat[:,1] = corr_con + corr_rep + corr_con_rep + corr_rep_con
        dat[:,2] = corr_con
        dat[:,3] = corr_rep
        dat[:,4] = corr_con_rep
        dat[:,5] = corr_rep_con
        savetxt(sys.argv[2], dat)
                


