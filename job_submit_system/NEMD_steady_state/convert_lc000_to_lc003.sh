rename -vs lc000 lc003 *
for fn in $(ls *.inp)
do
    sed -i '' 's/lc000/lc003/g' $fn
    sed -i '' 's/l_cap=0.00/l_cap=0.03/g' $fn
   
done
