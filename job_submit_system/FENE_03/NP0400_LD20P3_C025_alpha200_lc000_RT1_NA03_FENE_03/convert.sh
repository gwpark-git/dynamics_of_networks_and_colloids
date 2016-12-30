for fn in $(ls *.inp)
do
    sed -i '' 's/FN=NP0400_LD20P3_C025_alpha200_lc000_RT1_NA03/FN=NP0400_LD20P3_C025_alpha200_lc000_RT1_NA03_FENE_03/g' $fn
    done
