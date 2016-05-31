if [[ $1 == "run" ]]
then 
    for fn in *.jdl
    do
	glite-wms-job-submit -o jobID.txt -a $fn
    done
else
    echo 'Script for job submitting'
    echo 'the argument is used for preventing accidently called'
fi
