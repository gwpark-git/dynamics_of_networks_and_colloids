rename -vs NP0432_LD06P3 NP3456_LD06LX48 *
for fn in $(ls *.inp)
do
    sed -i '' 's/NP0432/NP3456/g' $fn
    sed -i '' 's/Np=432/Np=3456/g' $fn
    sed -i '' 's/LBx=6/LBx=48/g' $fn
    sed -i '' 's/LD06P3/LD06LX48/g' $fn
done
