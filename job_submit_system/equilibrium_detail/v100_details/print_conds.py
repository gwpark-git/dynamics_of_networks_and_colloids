
NP_arr = ['0512']
LD_arr = ['LD08P3']
C_arr = ['C025', 'C020', 'C015', 'C010', 'C005', 'C001']
RT_arr = ['RT1_', 'RT10_']
alpha_arr = ['alpha200']
lc_arr = ['lc000', 'lc006', 'lc012']
NC_arr = ['NC05', 'NC06', 'NC07', 'NC08', 'NC09', 'NC10']

for RT in RT_arr:
    for i_NP, NP in enumerate(NP_arr):
        # for LD in LD_arr:
        LD = LD_arr[i_NP]
        for C in C_arr:
            for alpha in alpha_arr:
                for lc in lc_arr:
                    for NC in NC_arr:
                        print 'bash automation_submit_recasna.sh equilibrium_detail/v100_details/NP%s_%s_%s_%s%s_%s_%s.inp $N_runs $init_sample_number ../%sVER/EQ_NP%s_%s/NP%s_%s.traj'%(NP, LD, C, RT, alpha, lc, NC, LD, NP, LD, NP, LD)
