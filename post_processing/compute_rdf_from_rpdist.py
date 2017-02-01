from numpy import *
import matplotlib.pyplot as plt
import sys
# sys.path.append('/Users/parkgunwoo/current_processing/stochastic_simulation_current/post_processing/')
import lib_rdf as lr

if size(sys.argv) < 5:
    print 'USAGE:'
    print 'argv[1] == filename for rpdist'
    print 'argv[2] == number of effective time. Note that HALF option will be used'
    print 'argv[3] == dr'
    print 'argv[4] == cut ratio. In this case, cutoff_length/box_dimension should be add'
else:
    fn_rpdist = sys.argv[1]
    rpdist = loadtxt(fn_rpdist)
    rpdist = (rpdist.flatten())[size(rpdist)/2:]

    Nt= int(sys.argv[2])/2
    Np = 512
    Np = int(fn_rpdist.split('NP')[1][0:4])
    dr = float(sys.argv[3])
    Ld = float(fn_rpdist.split('_LD')[1][0:2])
    cut_ratio = float(sys.argv[4])
    rdf, rho_local = lr.get_rdf_from_rpdist(rpdist, dr, Np, Nt, 3, Ld, cut_ratio)
    savetxt('%s.rdf'%(fn_rpdist), rdf)
    # rdf[:,2] *= 2.0 # compensate duplication counting on ddf
    # rho_local *= 2.0

    # plt.close()
    # plt.ion()
    # plt.plot(rdf[:,0], rdf[:,2], 'b-')
    # plt.grid()
    # plt.show()
