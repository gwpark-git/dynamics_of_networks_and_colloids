#!/bin/sh
inp_file=$1
out_file=$2

export OMP_NUM_THREADS=2
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/exp_soft/unina.it/gsl-2.1/lib/
export execute_path=$(pwd)
chmod u+x $execute_path/stochastic_simulation
$execute_path/stochastic_simulation $inp_file

export LFC_HOST=lfc02.scope.unina.it
tar czvf $out_file *
lcg-cr $out_file --vo unina.it -l lfn://grid/unina.it/gpark_dir/$out_file


