# rename -vs alpha200 alpha150 *
# for fn in $(ls *.inp)
# do
#     sed -i '' 's/alpha200/alpha150/g' $fn
#     sed -i '' 's/scale_factor_chain=2.0/scale_factor_chain=1.5/g' $fn
# done
rename -vs NP0512 NP0324 *
rename -vs LD08P3 LD06P3 *
for fn in $(ls *.inp)
do
    sed -i '' 's/512/324/g' $fn
    sed -i '' 's/LD08P3/LD06P3/g' $fn

    sed -i '' 's/cutoff_connection=4.0/cutoff_connection=3.0/g' $fn
    sed -i '' 's/box_dimension=8.0/box_dimension=6.0/g' $fn
done