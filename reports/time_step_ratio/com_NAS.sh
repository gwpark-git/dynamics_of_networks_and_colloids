module_py=../../post_processing/get_NAS_file.py
fn=NP1350_LD15P3_C100
given_path=(NP1350_LD15P3_C100_SF15_RT010 NP1350_LD15P3_C100_SF15_RT050 NP1350_LD15P3_C100_SF15_RT100 NP1350_LD15P3_C100_SF15_RT200 NP1350_LD15P3_C100_SF15_RT300 NP1350_LD15P3_C100_SF15_RT500 NP1350_LD15P3_C100_SF15_RT1k)

for gp in ${given_path[@]}
do
    echo $gp 'is processing'
    python $module_py ../../$gp/$fn.weight ../../$gp/$fn.traj $gp.NAS 1350
    echo $gp 'is done'
done

