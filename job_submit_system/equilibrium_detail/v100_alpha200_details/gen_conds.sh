# for fn in $(ls *.inp)
# do
#     rep_str1=
# done

# for fn in $(ls *.inp)
# do
#     #RT 10
#     fn_new=$(echo $fn | sed 's/RT1/RT10/g')
#     cp $fn $fn_new
#     sed -i '' 's/RT1/RT10/g' $fn_new
#     sed -i '' 's/Rt=1/Rt=10/g' $fn_new
# done

# for fn in $(ls *.inp)
# do
#     # alpha 110
#     fn_new=$(echo $fn | sed 's/alpha100/alpha110/g')
#     cp $fn $fn_new
#     sed -i '' 's/alpha100/alpha110/g' $fn_new
#     sed -i '' 's/scale_factor_chain=1.0/scale_factor_chain=1.1/g' $fn_new

#         # alpha 120
#     fn_new=$(echo $fn | sed 's/alpha100/alpha120/g')
#     cp $fn $fn_new
#     sed -i '' 's/alpha100/alpha120/g' $fn_new
#     sed -i '' 's/scale_factor_chain=1.0/scale_factor_chain=1.2/g' $fn_new

#     # alpha 130
#     fn_new=$(echo $fn | sed 's/alpha100/alpha130/g')
#     cp $fn $fn_new
#     sed -i '' 's/alpha100/alpha130/g' $fn_new
#     sed -i '' 's/scale_factor_chain=1.0/scale_factor_chain=1.3/g' $fn_new

#     # alpha 140
#     fn_new=$(echo $fn | sed 's/alpha100/alpha140/g')
#     cp $fn $fn_new
#     sed -i '' 's/alpha100/alpha140/g' $fn_new
#     sed -i '' 's/scale_factor_chain=1.0/scale_factor_chain=1.4/g' $fn_new

#     # alpha 150
#     fn_new=$(echo $fn | sed 's/alpha100/alpha150/g')
#     cp $fn $fn_new
#     sed -i '' 's/alpha100/alpha150/g' $fn_new
#     sed -i '' 's/scale_factor_chain=1.0/scale_factor_chain=1.5/g' $fn_new


#     # alpha 160
#     fn_new=$(echo $fn | sed 's/alpha100/alpha160/g')
#     cp $fn $fn_new
#     sed -i '' 's/alpha100/alpha160/g' $fn_new
#     sed -i '' 's/scale_factor_chain=1.0/scale_factor_chain=1.6/g' $fn_new

#     # alpha 170
#     fn_new=$(echo $fn | sed 's/alpha100/alpha170/g')
#     cp $fn $fn_new
#     sed -i '' 's/alpha100/alpha170/g' $fn_new
#     sed -i '' 's/scale_factor_chain=1.0/scale_factor_chain=1.7/g' $fn_new

#     # alpha 180
#     fn_new=$(echo $fn | sed 's/alpha100/alpha180/g')
#     cp $fn $fn_new
#     sed -i '' 's/alpha100/alpha180/g' $fn_new
#     sed -i '' 's/scale_factor_chain=1.0/scale_factor_chain=1.8/g' $fn_new

#     # alpha 190
#     fn_new=$(echo $fn | sed 's/alpha100/alpha190/g')
#     cp $fn $fn_new
#     sed -i '' 's/alpha100/alpha190/g' $fn_new
#     sed -i '' 's/scale_factor_chain=1.0/scale_factor_chain=1.9/g' $fn_new

#     # alpha 200
#     fn_new=$(echo $fn | sed 's/alpha100/alpha200/g')
#     cp $fn $fn_new
#     sed -i '' 's/alpha100/alpha200/g' $fn_new
#     sed -i '' 's/scale_factor_chain=1.0/scale_factor_chain=2.0/g' $fn_new
    
    
# done

# # for fn in $(ls *.inp)
# # do
# #     for alpha in $(cat alpha_cond.txt)
# #     do
# #         rep_str='s/alpha100/alpha'$alpha'/g'
# #         fn_new=$(echo $fn | sed $rep_str)
# #         # echo $fn_new
# #         cp $fn $fn_new
# #         sed -i '' $rep_str $fn_new

# #     done
# # done

# # for alpha in $(cat alpha_cond.txt)
# # do

    
# # done

# # for fn in $(ls *.inp)
# # do
# #     # for NC20
# #     fn_NC20=$(echo $fn | sed 's/NC10/NC20/g')
# #     cp $fn $fn_NC20
# #     sed -i '' 's/NC10/NC20/g' $fn_NC20
# #     sed -i '' 's/N_chains_per_particle=10/N_chains_per_particle=20/g' $fn_NC20
    
# #     # for NC05
# #     fn_NC05=$(echo $fn | sed 's/NC10/NC05/g')
# #     cp $fn $fn_NC05
# #     sed -i '' 's/NC10/NC05/g' $fn_NC05
# #     sed -i '' 's/N_chains_per_particle=10/N_chains_per_particle=5/g' $fn_NC05
    
# # done


# # for fn in $(ls *.inp)
# # do
# #     # for lc006
# #     fn_lc006=$(echo $fn | sed 's/lc012/lc006/g')
# #     cp $fn $fn_lc006
# #     sed -i '' 's/lc012/lc006/g' $fn_lc006
# #     sed -i '' 's/l_cap=0.12/l_cap=0.06/g' $fn_lc006
    
# #     # for lc003
# #     fn_lc003=$(echo $fn | sed 's/lc012/lc003/g')
# #     cp $fn $fn_lc003
# #     sed -i '' 's/lc012/lc003/g' $fn_lc003
# #     sed -i '' 's/l_cap=0.12/l_cap=0.03/g' $fn_lc003

# #     # for lc001
# #     fn_lc001=$(echo $fn | sed 's/lc012/lc001/g')
# #     cp $fn $fn_lc001
# #     sed -i '' 's/lc012/lc001/g' $fn_lc001
# #     sed -i '' 's/l_cap=0.12/l_cap=0.01/g' $fn_lc001

# #     # for lc000
# #     fn_lc000=$(echo $fn | sed 's/lc012/lc000/g')
# #     cp $fn $fn_lc000
# #     sed -i '' 's/lc012/lc000/g' $fn_lc000
# #     sed -i '' 's/l_cap=0.12/l_cap=0.00/g' $fn_lc000

# # done
