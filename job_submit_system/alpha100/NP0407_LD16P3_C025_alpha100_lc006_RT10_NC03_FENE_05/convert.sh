for fn in $(ls *.inp)
do
    sed -i '' 's/FN=NP0407_LD16P3_C025_alpha100_lc006_RT10_NC03/FN=NP0407_LD16P3_C025_alpha100_lc006_RT1_NC03/g' $fn
done
# rename -vs RT1 RT10 *
# for fn in $(ls *.inp)
# do

#     sed -i '' 's/RT1/RT10/g' $fn
#     sed -i '' 's/Rt=1/Rt=10/g' $fn

#     sed -i '' 's/Nt=2000000/Nt=10000000/g' $fn
#     sed -i '' 's/N_skip_ener=10/N_skip_ener=50/g' $fn
#     sed -i '' 's/N_skip_file=1000/N_skip_file=5000/g' $fn
#     sed -i '' 's/N_skip_rdist=10000/N_skip_rdist=50000/g' $fn
# done

# for fn in $(ls *.inp)
# do
# done

# rename -vs lc012 lc006 *
# for fn in $(ls *.inp)
# do
#     sed -i '' 's/lc012/lc006/g' $fn
#     sed -i '' 's/l_cap=0.12/l_cap=0.06/g' $fn
# done

# rename -vs NC05 NC03 *
# for fn in $(ls *.inp)
# do
#     sed -i '' 's/NC05/NC03/g' $fn
#     sed -i '' 's/N_chains_per_particle=5/N_chains_per_particle=3/g' $fn
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
