# for fn in $(ls *.inp)
# do
#     sed -i '' 's/RT1_NC10_FENE_04_inherit/RT10_NC10_FENE_04_inherit/g' $fn    
#     sed -i '' 's/inherit/short_inherit/g' $fn
#     sed -i '' 's/Rt=1/Rt=10/g' $fn
# done

# rename -vs RT1 RT10 *


fn_m3='NP0600_LD10P3_C100_alpha150_lc001_RT10_NC10_short_inherit_Wi_m3_shear.inp'
fn_m4=$(echo $fn_m3 | sed 's/m3/m4/g')
cp $fn_m3 $fn_m4
sed -i '' 's/m3/m4/g' $fn_m4
sed -i '' 's/Wi_tau_C=0.001/Wi_tau_C=0.0001/g' $fn_m4

fn_2m4=$(echo $fn_m3 | sed 's/m3/2m4/g')
cp $fn_m3 $fn_2m4
sed -i '' 's/m3/2m4/g' $fn_2m4
sed -i '' 's/Wi_tau_C=0.001/Wi_tau_C=0.0002/g' $fn_2m4

fn_4m4=$(echo $fn_m3 | sed 's/m3/4m4/g')
cp $fn_m3 $fn_4m4
sed -i '' 's/m3/4m4/g' $fn_4m4
sed -i '' 's/Wi_tau_C=0.001/Wi_tau_C=0.0004/g' $fn_4m4

fn_6m4=$(echo $fn_m3 | sed 's/m3/6m4/g')
cp $fn_m3 $fn_6m4
sed -i '' 's/m3/6m4/g' $fn_6m4
sed -i '' 's/Wi_tau_C=0.001/Wi_tau_C=0.0006/g' $fn_6m4



