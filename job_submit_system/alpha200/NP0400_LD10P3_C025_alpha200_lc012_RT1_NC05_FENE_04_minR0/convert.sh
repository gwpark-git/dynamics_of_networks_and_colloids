# rename -vs 512 400 *
# rename -vs LD08 LD10 *
# rename -vs NC05_short NC05_minR0_short *
for fn in $(ls *.inp)
do
    # sed -i '' 's/512/400/g' $fn
    # sed -i '' 's/LD08/LD10/g' $fn
    # sed -i '' 's/box_dimension=8.0/box_dimension=10.0/g' $fn
    # sed -i '' 's/_inherit/_minR0_inherit/g' $fn
    # sed -i '' 's/transition_probability=DISSOCIATION/transition_probability=MINIMUM_R0_DISSOCIATION/g' $fn
    # sed -i '' 's/FN=NP0512_LD08P3_C025_alpha200_lc012_RT1_NC05/FN=NP0512_LD08P3_C025_alpha200_lc012_RT1_NC05_minR0/g' $fn
    sed -i '' 's/FN=NP0400_LD10P3_C025_alpha200_lc012_RT1_NC05/FN=NP0400_LD10P3_C025_alpha200_lc012_RT1_NC05_minR0/g' $fn
done
