
import sys
from numpy import *
import matplotlib.pyplot as plt

def get_rpdist(fn_rpdist):
    rpdist = []
    with open(fn_rpdist, 'r') as f:
        for line in f:
            tmp_str = line.split('\t')[:-1]
            # if tmp_str[-1][-1] <> 'e' or tmp_str[-1][-2] <> 'e':
            # print line
            for x in tmp_str:
                rpdist.append(float(x))
    return rpdist

if size(sys.argv) > 1:

    fn_base = sys.argv[1]
    Np = int(fn_base.split('NP')[1][:4])
    Ld = float(fn_base.split('_LD')[1][:2])
    Nc = int(fn_base.split('_NC')[1][:2])
    C_rep = float(fn_base.split('_C')[1][:3])
    Rt = float(fn_base.split('_RT')[1][:3].split('_')[0])

    # fn_base = 'NP0400_LD10P3_C100_RT1_alpha150_lc006'
    ener = loadtxt('%s.ener'%(fn_base))
    # Nc = 10.
    # Np = 400

    rbdist = loadtxt('%s.rbdist'%(fn_base))[:-1]
    rbdist = rbdist[size(rbdist)/2:]

    rpdist = get_rpdist('%s.rpdist'%(fn_base))
    rpdist = rpdist[size(rpdist)/2:]

    # C_rep = 100.
    # Rt = 1.

    plt.close()
    plt.ion()
    plt.subplot(411)
    plt.plot(ener[:,0]/(C_rep*Rt), 100.*ener[:,4]/(Nc*Np), 'b-')
    plt.plot(ener[:,0]/(C_rep*Rt), 100.*ener[:,30]/(Nc*Np), 'r-')
    plt.xlabel(r'time / $\tau_0$')
    plt.ylabel('fNAS (%)')
    plt.title(fn_base)
    plt.subplot(412)
    plt.plot(ener[:,0]/(C_rep*Rt), ener[:,27], 'b-', label = 'elastic shear stress')
    plt.plot(ener[:,0]/(C_rep*Rt), ener[:,21], 'r-', label = 'repulsive shear stress', alpha=0.2)
    plt.xlabel(r'time / $\tau_0$')
    plt.ylabel('shear stress')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.subplot(413)
    plt.hist(rbdist, bins=50, color='blue', alpha=0.5, normed=True)
    plt.xlabel(r'distacne / $R_0$')
    plt.ylabel('bridge dist.')

    plt.subplot(414)
    plt.hist(rpdist, bins=50, color='red', alpha=0.5, normed=True)
    plt.xlabel(r'distacne / $R_0$')
    plt.ylabel('micelle dist.')
    get_axis = asarray(plt.axis())
    get_axis[1] = Ld/2.
    plt.axis(get_axis)
    plt.savefig('%s_stats.png'%(fn_base), dpi=300, bbox_inches='tight')
    # plt.show()
else:
    print 'USAGE: plot history of system with fraction of active bridge, stresses, bridge distribution, and micelle distance distributions'
    print 'argv[1] == base file name'
    print 'note that the output will be automatically with _stats.png'
    
