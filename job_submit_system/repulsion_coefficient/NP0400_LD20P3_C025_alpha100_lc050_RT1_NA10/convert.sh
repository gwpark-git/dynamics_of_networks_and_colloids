rename -vs C010 C025 *
rename -vs alpha200 alpha100 *
rename -vs alpha100_RT1 alpha100_lc050_RT1 *
rename -vs _inherit _short_inherit *
rename -vs RT1 RT1_NA10 *

for fn in $(ls *.inp)
do
    sed -i '' 's/C010/C025/g' $fn
    sed -i '' 's/repulsion_coefficient=10/repulsion_coefficient=25/g' $fn

    sed -i '' 's/alpha200/alpha100/g' $fn
    sed -i '' 's/scale_factor_chain=2/scale_factor_chain=1/g' $fn

    sed -i '' 's/alpha100_RT1/alpha100_lc050_RT1/g' $fn
    sed -i '' 's/l_cap=0.12/l_cap=0.50/g' $fn
    
    sed -i '' 's/_inherit/_short_inherit/g' $fn
    sed -i '' 's/Nt=2000000/Nt=200000/g' $fn

    sed -i '' 's/RT1/RT1_NA10/g' $fn
done
