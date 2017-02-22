rename -vs RT1 RT1div100 *
for fn in $(ls *.inp)
do
    sed -i '' 's/RT1_NC10_FENE_04_inherit/RT1div100_NC10_FENE_04_inherit/g' $fn
    sed -i '' 's/Rt=1/Rt=0.01/g' $fn
done
