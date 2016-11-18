
from numpy import *
from acf_fcn import *
import sys

if size(sys.argv) < 3:
    print 'Get stress autocorrelation from given .ener file for associating telechelic polymers'
    print 'USAGE :'
    print 'argv[1] == .ener file name'
    print 'argv[2] == output file name'
    print 'argv[3] (optional) == number of cutting variables (HALF will cut half). Default = 0.'
    print 'argv[4] (optional) == number of divisor (average over).'
else:
    if (sys.argv[1] == sys.argv[2]):
        print 'ERR: input and output have the same file name'
    else:
        print 'processing stress autocorrelation from %s to %s'%(sys.argv[1], sys.argv[2])
        ener = loadtxt(sys.argv[1])
        Nt = shape(ener)[0]
        N_cut = 0
        if size(sys.argv) == 4:
            if(sys.argv[3] == 'HALF'):
                N_cut = Nt/2
            else:
                N_cut = int(sys.argv[2])
        print 'compute stress autocorrelation for number of time steps %ld out of %ld'%(N_cut, Nt)                
        if size(sys.argv) == 5:

            N_div = int(sys.argv[4])
            print 'average over %d blocks'%(N_div)            
            Nt_block = Nt/N_div

            corr_con = corr(ener[N_cut:N_cut+Nt_block, 27], ener[N_cut:N_cut+Nt_block, 27])
            corr_rep = corr(ener[N_cut:N_cut+Nt_block, 21], ener[N_cut:N_cut+Nt_block, 21])
            corr_con_rep = corr(ener[N_cut:N_cut+Nt_block, 27], ener[N_cut:N_cut+Nt_block, 21])
            corr_rep_con = corr(ener[N_cut:N_cut+Nt_block, 21], ener[N_cut:N_cut+Nt_block, 27])

            for i in range(1, N_div):
                corr_con = corr(ener[N_cut + i*Nt_block:N_cut+(i+1)*Nt_block, 27], ener[N_cut + i*Nt_block:N_cut+(i+1)*Nt_block, 27])
                corr_rep = corr(ener[N_cut + i*Nt_block:N_cut+(i+1)*Nt_block, 21], ener[N_cut + i*Nt_block:N_cut+(i+1)*Nt_block, 21])
                corr_con_rep = corr(ener[N_cut + i*Nt_block:N_cut+(i+1)*Nt_block, 27], ener[N_cut + i*Nt_block:N_cut+(i+1)*Nt_block, 21])
                corr_rep_con = corr(ener[N_cut + i*Nt_block:N_cut+(i+1)*Nt_block, 21], ener[N_cut + i*Nt_block:N_cut+(i+1)*Nt_block, 27])

            corr_con /= float(N_div)
            corr_rep /= float(N_div)
            corr_con_rep /= float(N_div)
            corr_rep_con /= float(N_div)
            
        else:
            corr_con = (corr(ener[:,27], ener[:,27]) + corr(ener[:,28], ener[:,28]) + corr(ener[:,29], ener[:,29]))/3.
            corr_rep = (corr(ener[:,21], ener[:,21]) + corr(ener[:,22], ener[:,22]) + corr(ener[:,23], ener[:,23]))/3.
            corr_con_rep = (corr(ener[:,27], ener[:,21]) + corr(ener[:,28], ener[:,22]) + corr(ener[:,29], ener[:,23]))/3.
            corr_rep_con = (corr(ener[:,21], ener[:,27]) + corr(ener[:,22], ener[:,28]) + corr(ener[:,23], ener[:,29]))/3.
        dat = zeros([size(corr_con), 6])
        dat[:,0] = ener[:size(corr_con),0] - ener[0, 0]
        dat[:,1] = corr_con + corr_rep + corr_con_rep + corr_rep_con
        dat[:,2] = corr_con
        dat[:,3] = corr_rep
        dat[:,4] = corr_con_rep
        dat[:,5] = corr_rep_con
        savetxt(sys.argv[2], dat)
                


