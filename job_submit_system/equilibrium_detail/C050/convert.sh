for fn in $(ls *.inp)
do
    sed -i '' 's/N_THREADS_BD=6/N_THREADS_BD=8/g' $fn
done
   

# rename -vs C025 C050 *
# for fn in $(ls *.inp)
# do
#     sed -i '' 's/C025/C050/g' $fn
#     sed -i '' 's/repulsion_coefficient=25/repulsion_coefficient=50/g' $fn
# done
