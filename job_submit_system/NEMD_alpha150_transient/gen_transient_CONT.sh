# fn='NP0512_LD08P3_C025_alpha150_lc001_RT1_NC05_FENE_03'
fn=$1
fn_submit='submit_'$fn'_transient'
rm ${fn_submit}_gen.sh
# cat control_random_stream.txt > ${fn_submit}_gen.sh
for cnt in $(seq 0 1 99)
do
    fn_new=$(echo ${fn}_CONT$(printf "%03d" $cnt))
    mkdir $fn_new
    cp *.inp $fn_new/.
    # cp ${fn}.inp $fn_new/$fn_new.inp
    str_sed='s/CONTINUATION_STEP=-1/CONTINUATION_STEP='$cnt'/g'
    str_sed2='s/transient/transient_CONT'$(printf "%03d" $cnt)'/g'
    for fn_inp in $fn_new/*.inp
    do
        # echo $fn_inp
        sed -i $str_sed $fn_inp
        sed -i $str_sed2 $fn_inp
    done
    cp $fn_submit.sh ${fn_submit}_tmp.sh
    str_sed_submit='s/CONT=-1/CONT='$(printf "%03d" $cnt)'/g'
    sed -i $str_sed_submit ${fn_submit}_tmp.sh
    echo $fn_new is being processing
    echo 'echo '$fn_new >> ${fn_submit}_gen.sh
    bash ${fn_submit}_tmp.sh >> ${fn_submit}_gen.sh
done
rm ${fn_submit}_tmp.sh


