for fn in $(ls *.inp)
do
    sed -i '' 's/NA03/NC05/g' $fn
done
