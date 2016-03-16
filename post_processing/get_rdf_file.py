import os, sys
lib_path = os.path.abspath(os.path.join('..','post_processing'))
sys.path.append(lib_path)
from lib_rdf import *

from numpy import *
import matplotlib.pyplot as plt

if size(sys.argv) < 10:
    print 'USAGE:'
    print 'argv[1] == trajectory file name'
    print 'argv[2] == output rdf file name'
    print 'argv[3] == starting time step'
    print 'argv[4] == number of blocks for time average'
    print 'argv[5] == factor to multiplicate blocks'
    print 'argv[6] == N_dimension'
    print 'argv[7] == box_dimension'
    print 'argv[8] == number of particles'
    print 'argv[9] == dr'
else:

    fn_traj = sys.argv[1]
    fn_rdf_out = sys.argv[2]
    # fn_pdf = sys.argv[3]
    t_st = int(sys.argv[3])
    N_t_block = int(sys.argv[4])
    fac_t = int(sys.argv[5])
    # dt = float(sys.argv[6])
    traj = loadtxt(fn_traj)
    Nt = shape(traj)[0]
    N_dimension = int(sys.argv[6])
    box_dimension = float(sys.argv[7])
    half_bd = box_dimension/2.0
    Np = int(sys.argv[8])
    dr = float(sys.argv[9])
    rho = 0.4
    cut_ratio = 0.5
    fac_t = 1
    t = t_st + arange(N_t_block)*fac_t
    rdf, rho_local = get_rdf_ref(traj, t, dr, Np, N_dimension, box_dimension, cut_ratio)
    print rho_local
    savetxt(fn_rdf_out, rdf)
    # ref_unity = asarray([[0, 1], [cut_ratio*box_dimension, 1]])
    # kt = [0,1,2,3,4,5,6,7,8,9,10,15,20]
    # plt.clf()
    # plt.ion()
    # plt.plot(rdf[:,0], rdf[:,2], 'b-', label = 'RDF, ts=[%d, %d]*%d'%(t[0], t[-1], dt))
    # plt.plot(ref_unity[:,0], ref_unity[:,1], 'k--', linewidth=2, label = r'unity, ref. $\rho_D\approx$ %4.3f, $\rho_D/\rho\approx$ %4.3f'%(rho_local, rho_local/0.4))
    # plt.legend(loc = 'upper right')
    # plt.xlabel('dimensionless distance')
    # plt.ylabel('radial distribution function')
    # plt.axis([0, cut_ratio*box_dimension, 0.8, 1.2])
    # plt.xticks(kt)
    # plt.grid()
    # plt.savefig('rdf_FENE_r11.pdf')
    # plt.show()

