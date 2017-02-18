rename -vs 216 324 *
for fn in $(ls *.inp)
do
    sed -i '' 's/216/324/g' $fn
done
