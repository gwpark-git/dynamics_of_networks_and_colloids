for fn in $(ls *CONT*/*job/jobID.txt);
do
    idval=$(glite-wms-job-status -i $fn | grep Current | grep Aborted | wc -l)
    if [ $idval -eq 1 ]
    then
	CONT=$(echo $fn | awk -F/ '{ print $1}' | awk -FCONT '{ print $2}')
	path=$(echo $fn | awk -F/ '{ print $2}' | sed 's/.job//g')
	# identification of if phrase
	echo 'CONT='$CONT', path='$path
	cat *CONT${CONT}*.sh | grep $path >> resubmit_NEMD_alpha150_transient.sh
    fi
done