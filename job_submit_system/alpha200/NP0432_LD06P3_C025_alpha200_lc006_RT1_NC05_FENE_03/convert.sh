rename -vs 512 432 *
rename -vs LD08P3 LD06P3 *

for fn in $(ls *.inp)
do
    sed -i '' 's/512/432/g' $fn
    sed -i '' 's/LD08P3/LD06P3/g' $fn
    sed -i '' 's/box_dimension=8.0/box_dimension=6.0/g' $fn

    sed -i '' 's/FENE_04/FENE_03/g' $fn
    sed -i '' 's/cutoff_connection=4/cutoff_connection=3/g' $fn
    sed -i '' 's/ratio_RM_R0=4/ratio_RM_R0=3/g' $fn
done

# rename -vs lc012 lc003 *
# for fn in $(ls *.inp)
# do
#     sed -i '' 's/lc012/lc003/g' $fn
#     sed -i '' 's/l_cap=0.12/l_cap=0.03/g' $fn
# done

# for fn in $(ls *.inp)
# do
#     sed -i '' 's/FENE_05/FENE_04/g' $fn
#     sed -i '' 's/cutoff_connection=5/cutoff_connection=4/g' $fn
#     sed -i '' 's/ratio_RM_R0=5/ratio_RM_R0=4/g' $fn
# done

# rename -vs NP0600 NP0512 *
# rename -vs LD10P3 LD08P3 *
# rename -vs alpha100 alpha200 *


# for fn in $(ls *.inp)
# do
#     sed -i '' 's/600/512/g' $fn
#     sed -i '' 's/LD10P3/LD08P3/g' $fn
#     sed -i '' 's/box_dimension=10.0/box_dimension=8.0/g' $fn

#     sed -i '' 's/alpha100/alpha200/g' $fn
#     sed -i '' 's/scale_factor_chain=1.0/scale_factor_chain=2.0/g' $fn
# done

# rename -vs NP0407 NP0600 *
# rename -vs LD16P3 LD10P3 *

# for fn in $(ls *.inp)
# do
#     sed -i '' 's/407/600/g' $fn
#     sed -i '' 's/LD16P3/LD10P3/g' $fn
#     sed -i '' 's/box_dimension=16.0/box_dimension=10.0/g' $fn
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
