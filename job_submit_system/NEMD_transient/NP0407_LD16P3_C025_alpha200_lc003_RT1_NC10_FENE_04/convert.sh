# for fn in $(ls *.inp)
# do
#     sed -i '' 's/inherit/transient/g' $fn
#     sed -i '' 's/Nt=10000000/Nt=100000/g' $fn
#     sed -i '' 's/N_skip_ener=100/N_skip_ener=1/g' $fn
#     sed -i '' 's/tauR=0.0/tauR=0.00/g' $fn
# done
         
for fn in $(ls *.inp)
do
    sed -i '' 's/N_skip_file=10000/N_skip_file=1000/g' $fn
    sed -i '' 's/N_skip_rdist=100000/N_skip_rdist=10000/g' $fn
done
          
