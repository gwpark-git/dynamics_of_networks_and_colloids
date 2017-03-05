rename -vs C100 C025 *
rename -vs NC10 NC05 *
for fn in $(ls *.inp)
do
    sed -i '' 's/C100/C025/g' $fn
    sed -i '' 's/repulsion_coefficient=100/repulsion_coefficient=25/g' $fn
    sed -i '' 's/short_inherit/inherit/g' $fn
done

for fn in $(ls *.inp)
do
    sed -i '' 's/NC10/NC05/g' $fn
    sed -i '' 's/N_chains_per_particle=10/N_chains_per_particle=5/g' $fn
done
