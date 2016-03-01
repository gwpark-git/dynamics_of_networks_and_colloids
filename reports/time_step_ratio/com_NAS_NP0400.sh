module_py=../../post_processing/get_NAS_file.py
fn=NP0400_LD15P3_C100
given_path=(NP0400_LD15P3_C100_SF15_RT010 NP0400_LD15P3_C100_SF15_RT050 NP0400_LD15P3_C100_SF15_RT100 NP0400_LD15P3_C100_SF15_RT200 NP0400_LD15P3_C100_SF15_RT300 NP0400_LD15P3_C100_SF15_RT500 NP0400_LD15P3_C100_SF15_RT1k)

for gp in ${given_path[@]}
do
    echo $gp 'is processing'
    python $module_py ../../$gp/$fn.weight ../../$gp/$fn.traj NP0400/$gp.NAS 1350
    echo $gp 'is done'
done

