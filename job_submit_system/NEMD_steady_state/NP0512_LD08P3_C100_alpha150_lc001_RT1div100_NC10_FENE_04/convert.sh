rename -vs C025 C100 *
rename -vs alpha200 alpha150 *
rename -vs lc000 lc001 *

for fn in $(ls *.inp)
do
    sed -i '' 's/C025/C100/g' $fn
    sed -i '' 's/repulsion_coefficient=25/repulsion_coefficient=100/g' $fn
    sed -i '' 's/alpha200/alpha150/g' $fn
    sed -i '' 's/scale_factor_chain=2.0/scale_factor_chain=1.5/g' $fn
    sed -i '' 's/lc000/lc001/g' $fn
    sed -i '' 's/l_cap=0.00/l_cap=0.01/g' $fn

    sed -i '' 's/Nt=10000000/Nt=2000000/g' $fn
    sed -i '' 's/N_skip_ener=100/N_skip_ener=20/g' $fn
    sed -i '' 's/N_skip_file=10000/N_skip_file=2000/g' $fn
    sed -i '' 's/N_skip_rdist=100000/N_skip_rdist=20000/g' $fn
done
