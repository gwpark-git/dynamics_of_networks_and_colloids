
from numpy import *
import sys
from read_file import *

if size(sys.argv) < 6:
    print 'USAGE: GETTING MSD'
    print 'argv[1] == given trajectory'
    print 'argv[2] == output file'
    print 'argv[3] == number of spatial dimension'
    print 'argv[4] == number of particles'
    print 'argv[5] == half option (TRUE/FALSE, default FALSE)'
    # print 'argv[5] == skip frequency from given trajectory'
else:
    ND = int(sys.argv[3])
    NP = int(sys.argv[4])
    # skip_w = int(sys.argv[5])
    half_option = sys.argv[5]
    # traj = loadtxt(sys.argv[1])
    # traj = get_seq_data_float(sys.argv[1])
    traj = loadtxt(sys.argv[1])
    Nt = shape(traj)[0]

    Nst = 0
    if (half_option == 'TRUE'):
        Nst = int(Nt/2.)
        Nt = int(Nt/2.)
    
    tr_RR = zeros([Nt, 2])
    tr_RR[:,0] = traj[Nst:Nst+Nt,0]
    for i in range(NP):
        print 'processing particle %d'%(i)
        for k in range(ND):
            index_Rik = 2*ND*i + 1 + k
            tr_RR[:,1] += (traj[Nst:Nst+Nt, index_Rik] - traj[Nst, index_Rik])**2.
    tr_RR[:,1] /= float(NP)
    savetxt(sys.argv[2], tr_RR)

    # traj = loadtxt(sys.argv[1])
    # ND = 2
    # NP = 100
    # Nt = shape(traj)[0]

    # tr_RR = zeros([Nt, 2])
    # tr_RR[:,0] = traj[:,0]
    # for i in range(NP):
    #     for k in range(ND):
    #         index_Rik = 2*ND*i + 1 + k
    #         tr_RR[:,1] += (traj[:, index_Rik] - traj[0, index_Rik])**2.
    # tr_RR[:,1] /= float(NP)
    # savetxt(sys.argv[2], tr_RR)
