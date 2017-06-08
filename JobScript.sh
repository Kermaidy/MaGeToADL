#!/bin/bash

SCRIPT_PATH=/path/to/testadl
ERROR_PATH=${SCRIPT_PATH}/JobErrors/
OUTPUT_PATH=${SCRIPT_PATH}/JobOutput/

LAUNCH_DATE=`echo $(date +"%Y%m%d%H%M")`
LAUNCH_DATE=${LAUNCH_DATE}
JOBNAME=submit_at_$LAUNCH_DATE

cd ${SCRIPT_PATH}

TIER2PATH=${SCRIPT_PATH}/Tier2

declare -i iter

if [ $1 = "1" ];
then
	    JOBNAME_ID=${JOBNAME}_$iter
	    echo "Processing simulation : ${iter}"
	    qsub -P short -e $ERROR_PATH -o $OUTPUT_PATH -N $JOBNAME_ID ${SCRIPT_PATH}/LaunchJobSimulation.sh ${iter} ${SCRIPT_PATH}
elif [ $1 = "2" ];
then
	iter=0
	cd ${TIER2PATH}
	files=$(find -name "*.root")

	for file in ${files}
	do
    		iter=iter+1
    		JOBNAME_ID=${JOBNAME}_$iter
    		echo "Processing file : ${file}"
    		qsub -P short -e $ERROR_PATH -o $OUTPUT_PATH -N $JOBNAME_ID ${SCRIPT_PATH}/LaunchExecModuleIni.sh ${file} ${SCRIPT_PATH}
  	done
fi
