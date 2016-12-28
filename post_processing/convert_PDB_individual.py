
from numpy import *
import sys

if size(sys.argv) > 1:

    def get_px_index(index_bead):
        return index_bead*6 + 1


    # fn_base = 'NP0400_LD10P3_C025'
    fn_base = sys.argv[1]
    Np = 400
    Nd = 3
    traj_particles = zeros([Np, Nd])
    time_count = 0
    with open ('%s.traj'%(fn_base), 'r') as f_traj:
        for line in f_traj:
            line_str = line.replace('\t\n', '').split('\t')
            time_count += 1


    Np = (size(line_str) - 1)/6
    with open ('%s_LAST.pdb'%(fn_base), 'w') as f:
        # f.write('%d\n'%Np)
        # f.write('HEUR Micelle from stochastic_simulation_HEUR developed by Gun Woo Park\n')
        f.write('COMPND    HEUR. INP_FN_BASE: %s\n'%(fn_base))
        f.write('AUTHOR    Gun Woo Park\n')
        for i in range(Np):
            index_st = get_px_index(i)
            f.write('HETATM %4d PS1  PSD P%4d       %4.3f   %4.3f   %4.3f  1.00  0.00      PSDOPS\n'%(i, i, float(line_str[index_st + 0]), float(line_str[index_st + 1]), float(line_str[index_st + 2])))
            # f.write('MIC\t%s\t%s\t%s\n'%(line_str[index_st + 0], line_str[index_st + 1], line_str[index_st + 2]))
        with open ('%s.hash'%(fn_base), 'r') as f_hash:
            cnt = 0
            for line in f_hash:
                if cnt >= (time_count-1)*Np:
                    adj_list = line.replace('\t\n', '').split('\t')
                    if size(adj_list) > 1:
                        f.write('CONNECT ')
                        for i in range(size(adj_list)):
                            f.write('%4d '%(float(adj_list[i])))
                        f.write('\n')
                cnt += 1
        f.write('END')

else:
    print 'USAGE: to convert single time step of trajectory and hash into PDB file format'
    print 'argv[1] == base filename'
    print 'Note that the results will write with the same base file anme with _LAST.pdb '
    
