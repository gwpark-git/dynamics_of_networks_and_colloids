for fn in $(ls *.inp)
do
    sed -i '' 's/Nt=2000000/Nt=200000/g' $fn
    sed -i '' 's/_inherit/_short_inherit/g' $fn
    # sed -i '' 's/alpha100/alpha200/g' $fn
    # sed -i '' 's/scale_factor_chain=1.0/scale_factor_chain=2.0/g' $fn
    # sed -i '' 's/NC05/NC05_ZNRC/g' $fn
    # sed -i '' 's/FENE_05/FENE_04/g' $fn
    # sed -i '' 's/cutoff_connection=5/cutoff_connection=4/g' $fn
    # sed -i '' 's/ratio_RM_R0=5/ratio_RM_R0=4/g' $fn
    # cat ZNRC.txt >> $fn
done
