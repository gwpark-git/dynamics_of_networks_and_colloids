rename -vs alpha200 alpha150 *
for fn in $(ls *.inp)
do
    sed -i '' 's/alpha200/alpha150/g' $fn
    sed -i '' 's/scale_factor_chain=2.0/scale_factor_chain=1.5/g' $fn
done
