# for fn in $(ls *.inp)
# do
#     fn_new=$(echo $fn | sed 's/RT1_/RT10_/g')
#     cp $fn $fn_new
#     sed -i '' 's/RT1_/RT10_/g' $fn_new
#     sed -i '' 's/Rt=1/Rt=10/g' $fn_new
# done

# for fn in $(ls *.inp)
# do
#     sed -i '' 's/N_THREADS_BD=6/N_THREADS_BD=8/g' $fn
#     sed -i '' 's#/Users/parkgunwoo/current_processing/stochastic_simulation_current/LD08P3VER/EQ_NP0512_LD08P3/NP0512_LD08P3#/NP0512_LD08P3#g' $fn
# done
