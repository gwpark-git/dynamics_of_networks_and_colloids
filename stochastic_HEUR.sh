#!/bin/sh
out_path='result_HEUR'
out_file='simulation_result.tgz'
inp_file='run.inp'

export OMP_NUM_THREADS=2

export HEUR_PATH=/home/gpark/stochastic_HEUR
./stochastic_simulation run.inp

tar czvf /ustre/home/$USER/$out_file $out_path/*
mv /lustre/home/$USER/$out_file .
lcg-cr $out_file --vo unina.it -l lfn://grid/unina.it/gpark_dir/$out_file
