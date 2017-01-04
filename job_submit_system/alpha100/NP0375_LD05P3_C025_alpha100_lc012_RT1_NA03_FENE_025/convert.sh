rename -vs NP0400 NP0375 *
rename -vs LD10P3 LD05P3 *
rename -vs alpha200 alpha100 *
rename -vs lc006 lc012 *
for fn in $(ls *.inp)
do
    sed -i '' 's/400/375/g' $fn
    sed -i '' 's/LD10P3/LD05P3/g' $fn
    sed -i '' 's/box_dimension=10/box_dimension=5/g' $fn
    sed -i '' 's/alpha200/alpha100/g' $fn
    sed -i '' 's/scale_factor_chain=2/scale_factor_chain=1/g' $fn
    sed -i '' 's/FENE_03/FENE_025/g' $fn
    sed -i '' 's/cutoff_connection=3/cutoff_connection=2.5/g' $fn
    sed -i '' 's/ratio_RM_R0=3/ratio_RM_R0=2.5/g' $fn
    sed -i '' 's/lc006/lc012/g' $fn
    sed -i '' 's/l_cap=0.06/l_cap=0.12/g' $fn
done
