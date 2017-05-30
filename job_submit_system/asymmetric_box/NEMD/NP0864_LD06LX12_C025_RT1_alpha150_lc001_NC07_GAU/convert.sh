for fn in $(ls *.inp)
do
    sed -i '' 's/NP0432/NP0864/g' $fn
    sed -i '' 's/Np=432/Np=864/g' $fn
    sed -i '' 's/LBx=6/LBx=12/g' $fn
    sed -i '' 's/LD06P3/LD06LX12/g' $fn
done
