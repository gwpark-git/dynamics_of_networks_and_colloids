
if [ $# -lt 2 ]
then
    echo "Script for job submitting procedure"
    echo "1: given input file"
    echo "2: number of different runs (>1 means the seed will be re-generated)"
    echo "3: init seed number"
    echo "4 > given data"
else
    server=emi2-ce0[1-2].scope.unina.it:8443/cream-pbs-unina_hpc
    # server=recasna-ce01.unina.it:8443/cream-pbs-recas-normal
    file_inp=$1 # input file 
    N_runs=$2
    init_seed_number=$3
    if [ $N_runs -ne 1 ];
    then
	echo Number of runs is set with number $N_runs. Seed numbers \in .inp file will be ignored.
    fi

    rm $file_inp.job/*
    mkdir $file_inp.job
    cp ../stochastic_simulation $file_inp.job/stochastic_simulation
    cp stochastic_HEUR.sh $file_inp.job/stochastic_HEUR.sh
    arg_str=''
    if [ $# -gt 2 ]
    then
	for (( i=4; i<=$#; i++))
    # it will count and copy the given data files into the job submittion folder
	do
	    eval arg_fn=\$$i
	    base_fn=$(basename $arg_fn)
	    cp $arg_fn $file_inp.job/$base_fn
	    arg_str=$arg_str' '$base_fn
	    echo $arg_str
	done
    fi
    echo $arg_string
    echo '* JOB DESCRIPTION:' >> $file_inp.job/jobs.org
    cat $file_inp >> $file_inp.job/jobs.org
    python gen_job_files.py $file_inp $N_runs $init_seed_number stochastic_HEUR.sh $server $arg_str >> $file_inp.job/jobs.org
    
    cp execute_submit.sh $file_inp.job/execute_submit.sh
    cd $file_inp.job
    source execute_submit.sh run
    cd ..
fi