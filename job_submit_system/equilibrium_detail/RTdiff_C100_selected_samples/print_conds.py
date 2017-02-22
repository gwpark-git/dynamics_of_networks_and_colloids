
subpath = 'RTdiff_C100_selected_samples'
Crep = '100'
NP_arr = ['0407', '0400', '0600', '0512', '0375']
LD_arr = ['LD16P3', 'LD10P3', 'LD10P3', 'LD08P3', 'LD05P3']
RT_arr = ['RT1div25_', 'RT10_', 'RT100_']
alpha_arr = ['alpha150']
lc_arr = ['lc001']
NC_arr = ['NC05']

for i_NP, NP in enumerate(NP_arr):
    # for LD in LD_arr:
    LD = LD_arr[i_NP]
    for RT in RT_arr:
        for alpha in alpha_arr:
            for lc in lc_arr:
                for NC in NC_arr:
                    print 'bash automation_submit_recasna.sh equilibrium_detail/%s/NP%s_%s_C%s_%s%s_%s_%s.inp $N_runs $init_sample_number ../%sVER/EQ_NP%s_%s/NP%s_%s.traj'%(subpath, NP, LD, Crep, RT, alpha, lc, NC, LD, NP, LD, NP, LD)
