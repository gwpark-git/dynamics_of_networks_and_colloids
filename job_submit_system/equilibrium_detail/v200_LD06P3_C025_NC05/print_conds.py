
# NP_arr = ['0407', '0400', '0600', '0512', '0375']
NP_arr = ['0432']
LD_arr = ['LD06P3']
# LD_arr = ['LD16P3', 'LD10P3', 'LD10P3', 'LD08P3', 'LD05P3']
RT_arr = ['RT1_', 'RT10_']
alpha_arr = ['alpha100', 'alpha110', 'alpha120', 'alpha130', 'alpha140', 'alpha150', 'alpha160', 'alpha170', 'alpha180', 'alpha190', 'alpha200']
# lc_arr = ['lc000', 'lc001', 'lc003', 'lc006', 'lc012']
lc_arr = ['lc000', 'lc001', 'lc003']
NC_arr = ['NC05']

for i_NP, NP in enumerate(NP_arr):
    # for LD in LD_arr:
    LD = LD_arr[i_NP]
    for RT in RT_arr:
        for alpha in alpha_arr:
            for lc in lc_arr:
                for NC in NC_arr:
                    print 'bash automation_submit.sh equilibrium_detail/v200_LD06P3_C025_NC05/NP%s_%s_C025_%s%s_%s_%s.inp $N_runs $init_sample_number ../%sVER/EQ_NP%s_%s/NP%s_%s.traj'%(NP, LD, RT, alpha, lc, NC, LD, NP, LD, NP, LD)
