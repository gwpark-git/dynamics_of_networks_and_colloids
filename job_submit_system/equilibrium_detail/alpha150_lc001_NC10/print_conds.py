
subpath = 'alpha150_lc001_NC10'
Crep = '025'
NP_arr = ['0512', '0324', '0432', '0375']
LD_arr = ['LD08P3', 'LD06P3', 'LD06P3', 'LD05P3']
RT_arr = ['RT1_', 'RT10_']
alpha_arr = ['alpha150']
lc_arr = ['lc001']
NC_arr = ['NC10']

for i_NP, NP in enumerate(NP_arr):
    # for LD in LD_arr:
    LD = LD_arr[i_NP]
    for RT in RT_arr:
        for alpha in alpha_arr:
            for lc in lc_arr:
                for NC in NC_arr:
                    print 'bash automation_submit_recasna.sh equilibrium_detail/%s/NP%s_%s_C%s_%s%s_%s_%s.inp $N_runs $init_sample_number ../%sVER/EQ_NP%s_%s/NP%s_%s.traj'%(subpath, NP, LD, Crep, RT, alpha, lc, NC, LD, NP, LD, NP, LD)
