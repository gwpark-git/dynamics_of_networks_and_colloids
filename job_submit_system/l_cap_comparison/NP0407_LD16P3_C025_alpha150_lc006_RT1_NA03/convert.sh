rename -vs NP0400 NP0407 *
rename -vs LD10P3 LD16P3 *

for fn in $(ls *.inp)
do
    sed -i '' 's/NP0400/NP0407/g' $fn
    sed -i '' 's/Np=400/Np=407/g' $fn
    sed -i '' 's/LD10P3/LD16P3/g' $fn
    sed -i '' 's/box_dimension=10.0/box_dimension=16.0/g' $fn

done
