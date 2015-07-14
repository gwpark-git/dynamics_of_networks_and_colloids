
from numpy import *
import matplotlib.pyplot as plt
import sys
from scipy.stats import linregress

if size(sys.argv) == 5:
    dat = loadtxt(sys.argv[1])
    cut_frac_st = float(sys.argv[3])
    cut_frac_end = float(sys.argv[4])
    Nt = shape(dat)[0]
    val = linregress(dat[Nt*cut_frac_st:Nt*cut_frac_end, 0], dat[Nt*cut_frac_st:Nt*cut_frac_end, 1])
    t_reg = linspace(dat[0, 0], dat[-1, 0], 100)
    y_reg = val[0]*t_reg + val[1]

    plt.clf()
    plt.plot(dat[:,0], dat[:,1], 'b-', label = 'MSD')
    plt.plot(t_reg, y_reg, 'r-', label = 'lin. reg. slope=%6.3e'%(val[0]))
    plt.grid('on')
    plt.xlabel(r'reduced time, $\tilde{t}$')
    plt.ylabel('mean-square displacement, MSD')
    plt.legend(loc = 'upper left')
    plt.savefig(sys.argv[2])
    # plt.show()
else:
    print 'USAGE for plotting MSD:'
    print 'argv[1] == input MSD file'
    print 'argv[2] == output pdf file'
    print 'argv[3] == cut fraction'
