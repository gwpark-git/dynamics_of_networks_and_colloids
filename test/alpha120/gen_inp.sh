fn='EQ_TOP_NP0800_LD20P3_C025_alpha120_NANC.inp'
C_init='C025'
for i in $(seq 50 25 200)
do
    C_index=C$(printf %03d $i)
    fn_str_rep='s/'$C_init'/'$C_index'/g'
    # echo $fn_str_rep
    fn_new=$(echo $fn | sed $fn_str_rep)
    cp $fn $fn_new
    echo $fn_new
    replace_string_1='s/repulsion_coefficient=25/repulsion_coefficient='$i'/g'
    replace_string_2='s/filename_base=EQ_TOP_NP0800_LD20P3_C025_alpha120_NANC/filename_base=EQ_TOP_NP0800_LD20P3_'$C_index'_alpha120_NANC/g'
    # echo $replace_string_1
    # echo $replace_string_2
    sed -i '' $replace_string_1 $fn_new
    sed -i '' $replace_string_2 $fn_new
done
