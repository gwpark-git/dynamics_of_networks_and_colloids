rename -vs NP0432_LD06P3 NP1728_LD06LX24 *
for fn in $(ls *.inp)
do
    sed -i '' 's/NP0432/NP1728/g' $fn
    sed -i '' 's/Np=432/Np=1728/g' $fn
    sed -i '' 's/LBx=6/LBx=24/g' $fn
    sed -i '' 's/LD06P3/LD06LX24/g' $fn
done
