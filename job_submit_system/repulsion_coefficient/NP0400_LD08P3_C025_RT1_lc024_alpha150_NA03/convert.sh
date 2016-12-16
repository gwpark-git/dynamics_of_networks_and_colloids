rename -vs LD10P3 LD08P3 *
rename -vs RT1k RT1 *
rename -vs alpha150_ alpha150_lc024_ *

for fn in $(ls *.inp)
do

    sed -i '' 's/NP0400_LD10P3_C025_alpha150_RT100_NA03/NP0400_LD08P3_C025_RT1_lc024_alpha150_NA03/g' $fn
    sed -i '' 's/LD10P3/LD08P3/g' $fn
    sed -i '' 's/box_dimension=10.0/box_dimension=8.0/g' $fn

    sed -i '' 's/RT1k/RT1/g' $fn
    sed -i '' 's/Rt=1000/Rt=1/g' $fn

    sed -i '' 's/l_cap=0.12/l_cap=0.24/g' $fn
    sed -i '' 's/alpha150_RT1/alpha150_lc024_RT1/g' $fn

done
