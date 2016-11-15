
from numpy import *
from acf_fcn import *
import sys

if size(sys.argv) < 3:
    print 'Get stress autocorrelation from given .ener file for associating telechelic polymers'
    print 'USAGE :'
    print 'argv[1] == .ener file name'
    print 'argv[2] == output file name'
    print 'argv[3] == number of cutting variables (HALF will cut half). Default = 0'
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
                N_cut = int(sys.argv[2])
        print 'compute stress autocorrelation for number of time steps %ld out of %ld'%(N_cut, Nt)
        corr_con = corr(ener[:,27], ener[:,27])
        corr_rep = corr(ener[:,21], ener[:,21])
        corr_con_rep = corr(ener[:,27], ener[:,21])
        corr_rep_con = corr(ener[:,21], ener[:,27])
        dat = zeros([size(corr_con), 6])
        dat[:,0] = ener[:size(corr_con),0] - ener[0, 0]
        dat[:,1] = corr_con + corr_rep + corr_con_rep + corr_rep_con
        dat[:,2] = corr_con
        dat[:,3] = corr_rep
        dat[:,4] = corr_con_rep
        dat[:,5] = corr_rep_con
        savetxt(sys.argv[2], dat)


