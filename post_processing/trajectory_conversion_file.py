
from numpy import *
import sys

def sign(x):
    if x < 0.:
        return -1.
    return 1.

if size(sys.argv) < 6:
    print 'USAGE: Converting Trajectory without PBC'
    print 'Note: this is purpose to measure dynamic qunatities'
    print 'argv[1] == input trajectory'
    print 'argv[2] == output data'
    print 'argv[3] == spatial dimension'
    print 'argv[4] == number of particles'
    print 'argv[5] == box dimension'
else:
    dat = loadtxt(sys.argv[1])
    ND = int(sys.argv[3])
    NP = long(sys.argv[4])
    LB = float(sys.argv[5])
    # ND = 2
    # NP = 100
    Nt = shape(dat)[0]
    # LB = 10.0

    # dat_conv = zeros([Nt, 2*ND*NP + 1])


    def inv_PBC(x_now, x_next, LB):
        dX = x_next - x_now
        if abs(dX) > 0.5*LB:
            return inv_PBC(x_now, x_next - sign(dX)*LB, LB)
        return x_next

    for i in range(NP):
        for k in range(ND):
            index_Rik = 2*ND*i + 1 + k
            for t in range(1, Nt):
                dat[t, index_Rik] = inv_PBC(dat[t-1, index_Rik], dat[t, index_Rik], LB)
                # dR = dat[t, index_Rik] - dat[t-1, index_Rik]
                # if abs(dR) > 0.5*LB:
                #     print i, t, abs(dat[t-1, index_Rik] - dat[t-2, index_Rik]), abs(dR)

                # if abs(dat[t, index_Rik]) > LB:
                #     print t, i, k, dat[t, index_Rik]
                # if fabs(dR) > 0.5*LB :
                #     dat[t, index_Rik] -= sign(dR)*LB
                #     if i==0:
                #         dR_2 = dat[t, index_Rik] - dat[t-1, index_Rik]
                #         if abs(dR_2) > 0.5*LB:
                #             dat[t, index_Rik] -= sign(dR_2)*LB
                #             dR_3 = dat[t, index_Rik] - dat[t-1, index_Rik]
                #             if abs(dR_3) > 0.5*LB:
                #                 print i, t, abs(dat[t-1, index_Rik] - dat[t-2, index_Rik]), abs(dR), abs(dR_2)

    savetxt(sys.argv[2], dat)
