rename -vs 324 432 *
for fn in $(ls *.inp)
do
    sed -i '' 's/324/432/g' $fn
done
