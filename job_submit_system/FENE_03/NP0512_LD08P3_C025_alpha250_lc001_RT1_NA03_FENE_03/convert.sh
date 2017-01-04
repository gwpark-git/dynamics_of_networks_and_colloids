rename -vs NP0400_LD10P3 NP0512_LD08P3 *
rename -vs alpha200 alpha300 *

for fn in $(ls *.inp)
do
    sed -i '' 's/Nt=2000000/Nt=200000/g' $fn
    sed -i '' 's/_inherit/_short_inherit/g' $fn

    sed -i '' 's/400/512/g' $fn
    sed -i '' 's/LD10P3/LD08P3/g' $fn
    sed -i '' 's/box_dimension=10/box_dimension=8/g' $fn

    sed -i '' 's/alpha200/alpha300/g' $fn
    sed -i '' 's/scale_factor_chain=2.0/scale_factor_chain=3.0/g' $fn
done
