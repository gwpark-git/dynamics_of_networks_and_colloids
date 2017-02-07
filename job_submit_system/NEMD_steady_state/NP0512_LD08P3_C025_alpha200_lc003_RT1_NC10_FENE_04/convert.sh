
rename -vs NP0600 NP0512 *
rename -vs LD10 LD08 *

for fn in $(ls *.inp)
do
    sed -i '' 's/600/512/g' $fn
    sed -i '' 's/LD10/LD08/g' $fn
    sed -i '' 's/box_dimension=10/box_dimension=8/g' $fn

done
