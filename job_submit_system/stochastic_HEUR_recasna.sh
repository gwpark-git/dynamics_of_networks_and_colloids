#!/bin/sh

inp_file=$1
out_file=$2

export OMP_NUM_THREADS=$3
source /opt/exp_soft/unina.it/intel/composer_xe_2013_sp1.3.174/mkl/bin/mklvars.sh intel64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/exp_soft/unina.it/gsl-2.1/lib/
export execute_path=$(pwd)
chmod u+x $execute_path/stochastic_simulation
# if [ $# -lt 4 ]
# then
$execute_path/stochastic_simulation $inp_file
# else
#     for (( i=4; i<=$#; i++))
#     do
# 	eval arg_inp=\$$i
# 	inp_file=$inp_file' '$arg_inp
#     done
#     $execute_path/stochastic_simulation $inp_file
# fi    

export LFC_HOST=lfc02.scope.unina.it
export VO_UNINA_IT_DEFAULT_SE=se01.scope.unina.it
tar czvf $out_file *
lcg-cr $out_file --vo unina.it -l lfn://grid/unina.it/gpark_dir/$out_file


