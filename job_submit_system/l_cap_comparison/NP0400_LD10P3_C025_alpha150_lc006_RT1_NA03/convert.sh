# rename -vs RT100 RT1 *
# rename -vs alpha150_RT1 alpha150_lc006_RT1 *

for fn in $(ls *.inp)
do
    sed -i '' 's/Nt=2000000/Nt=200000/g' $fn
    sed -i '' 's/RT100/RT1/g' $fn
    sed -i '' 's/Rt=100/Rt=1/g' $fn

    sed -i '' 's/alpha150_RT1/alpha150_lc006_RT1/g' $fn
    sed -i '' 's/l_cap=0.12/l_cap=0.06/g' $fn

done
