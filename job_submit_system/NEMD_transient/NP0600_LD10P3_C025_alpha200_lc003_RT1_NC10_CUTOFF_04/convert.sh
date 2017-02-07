
rename -vs NP0400 NP0600 *
for fn in $(ls *.inp)
do
    sed -i '' 's/400/600/g' $fn
    

done
