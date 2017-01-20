rename -vs RT1 RT1div25 *
for fn in $(ls *.inp)
do
    sed -i '' 's/RT1_NC03_FENE_05_inherit/RT1div25_NC03_FENE_05_inherit/g' $fn
    sed -i '' 's/Rt=1/Rt=0.04/g' $fn
done

# rename -vs lc012 lc006 *
# for fn in $(ls *.inp)
# do
#     sed -i '' 's/lc012/lc006/g' $fn
#     sed -i '' 's/l_cap=0.12/l_cap=0.06/g' $fn
# done
# for fn in $(ls *.inp)
# do
#     sed -i '' 's/FENE_06/FENE_05/g' $fn
#     sed -i '' 's/cutoff_connection=6/cutoff_connection=5/g' $fn
#     sed -i '' 's/ratio_RM_R0=6/ratio_RM_R0=5/g' $fn
# done

# rename -vs C100 C025 *
# rename -vs NA03 NC05 *
# rename -vs alpha150 alpha100 *
# rename -vs lc000 lc012 *

# for fn in $(ls *.inp)
# do
#     sed -i '' 's/C100/C025/g' $fn
#     sed -i '' 's/repulsion_coefficient=100/repulsion_coefficient=25/g' $fn
#     sed -i '' 's/NA03/NC05/g' $fn
#     sed -i '' 's/N_chains_per_particle=10/N_chains_per_particle=5/g' $fn
#     sed -i '' 's/alpha150/alpha100/g' $fn
#     sed -i '' 's/scale_factor_chain=1.5/scale_factor_chain=1.0/g' $fn
#     sed -i '' 's/lc000/lc012/g' $fn
#     sed -i '' 's/l_cap=0.00/l_cap=0.12/g' $fn
    
# done
