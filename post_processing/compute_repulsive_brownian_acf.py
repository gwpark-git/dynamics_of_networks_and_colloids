
from numpy import *
from acf_fcn import *
import sys

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
            # corr_con = (corr(ener[N_st:N_st+Nt_block, 27], ener[N_st:N_st+Nt_block, 27])
            #             + corr(ener[N_st:N_st+Nt_block, 28], ener[N_st:N_st+Nt_block, 28])
            #             + corr(ener[N_st:N_st+Nt_block, 29], ener[N_st:N_st+Nt_block, 29]))/3.
            corr_rep = (corr(ener[N_st:N_st+Nt_block, 21], ener[N_st:N_st+Nt_block, 21])
                        + corr(ener[N_st:N_st+Nt_block, 22], ener[N_st:N_st+Nt_block, 22])
                        + corr(ener[N_st:N_st+Nt_block, 23], ener[N_st:N_st+Nt_block, 23]))/3.
            # corr_con_rep = (corr(ener[N_st:N_st+Nt_block, 27], ener[N_st:N_st+Nt_block, 21])
            #                 + corr(ener[N_st:N_st+Nt_block, 28], ener[N_st:N_st+Nt_block, 22])
            #                 + corr(ener[N_st:N_st+Nt_block, 29], ener[N_st:N_st+Nt_block, 23]))/3.
            # corr_rep_con = (corr(ener[N_st:N_st+Nt_block, 21], ener[N_st:N_st+Nt_block, 27])
            #                 + corr(ener[N_st:N_st+Nt_block, 22], ener[N_st:N_st+Nt_block, 28])
            #                 + corr(ener[N_st:N_st+Nt_block, 23], ener[N_st:N_st+Nt_block, 29]))/3.

            for i in range(1, N_div):
                print '\t### %ld out of %ld is being processing'%(i, N_div)                
                N_st = N_cut + i*Nt_inc
                # corr_con += (corr(ener[N_st:N_st+Nt_block, 27], ener[N_st:N_st+Nt_block, 27])
                #              + corr(ener[N_st:N_st+Nt_block, 28], ener[N_st:N_st+Nt_block, 28])
                #              + corr(ener[N_st:N_st+Nt_block, 29], ener[N_st:N_st+Nt_block, 29]))/3.
                corr_rep += (corr(ener[N_st:N_st+Nt_block, 21], ener[N_st:N_st+Nt_block, 21])
                             + corr(ener[N_st:N_st+Nt_block, 22], ener[N_st:N_st+Nt_block, 22])
                             + corr(ener[N_st:N_st+Nt_block, 23], ener[N_st:N_st+Nt_block, 23]))/3.
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
            # corr_con = (corr(ener[dat_range,27], ener[dat_range,27]) + corr(ener[dat_range,28], ener[dat_range,28]) + corr(ener[dat_range,29], ener[dat_range,29]))/3.
            corr_rep = (corr(ener[dat_range,21], ener[dat_range,21]) + corr(ener[dat_range,22], ener[dat_range,22]) + corr(ener[dat_range,23], ener[dat_range,23]))/3.
            # corr_con_rep = (corr(ener[dat_range,27], ener[dat_range,21]) + corr(ener[dat_range,28], ener[dat_range,22]) + corr(ener[dat_range,29], ener[dat_range,23]))/3.
            # corr_rep_con = (corr(ener[dat_range,21], ener[dat_range,27]) + corr(ener[dat_range,22], ener[dat_range,28]) + corr(ener[dat_range,23], ener[dat_range,29]))/3.
        dat = zeros([size(corr_rep), 6])
        dat[:,0] = ener[:size(corr_rep),0] - ener[0, 0]
        dat[:,0] *= float(stride)
        dat[:,1] = corr_rep
        # dat[:,2] = corr_con
        dat[:,3] = corr_rep
        # dat[:,4] = corr_con_rep
        # dat[:,5] = corr_rep_con
        savetxt(sys.argv[2], dat)
                


