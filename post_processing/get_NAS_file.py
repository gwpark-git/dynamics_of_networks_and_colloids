

from numpy import *
import matplotlib.pyplot as plt
import sys


def get_N_assoc(dat, Np):
    Nt = int(shape(dat)[0]/Np)
    re = zeros(Nt)
    for ts in range(Nt):
        re[ts] = sum(dat[ts*Np:(ts+1)*Np, 1:])/2.0
    return re

if size(sys.argv) < 4:
    print 'USAGE of calculation NAS:'
    print 'argv[1] == given weight function'
    print 'argv[2] == given traj file for time'
    print 'argv[3] == output file name'
    print 'argv[4] == number of particles'
else:
    dat = loadtxt(sys.argv[1])
    t = loadtxt(sys.argv[2], usecols=(0,))
    Np = int(sys.argv[4])
    NAS = get_N_assoc(dat, Np)
    re = zeros([size(NAS), 2])
    re[:,0] = t
    re[:,1] = NAS
    savetxt(sys.argv[3], re)
