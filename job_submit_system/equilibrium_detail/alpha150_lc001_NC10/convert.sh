for fn in $(ls *.inp)
do
    sed -i 's/NC07/NC10/g' $fn
    sed -i 's/N_chains_per_particle=7/N_chains_per_particle=10/g' $fn
    fn_new=$(echo $fn | sed 's/NC07/NC10/g')
    mv $fn $fn_new
done