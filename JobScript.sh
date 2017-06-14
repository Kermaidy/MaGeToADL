#!/bin/bash

SCRIPT_PATH=/lfs/l3/gerda/kermaidy/Analysis/software/src/MaGeToADL
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
	cd ${SCRIPT_PATH}/../mage/macros/
	files=$(find -maxdepth 1 -name "Gerda_Gerda$2_17.root")
    	iter=0
	for file in ${files}
        do
	    iter=iter+1
	    JOBNAME_ID=${JOBNAME}_$iter
	    echo "Processing simulation : ${iter} ${file}"
	    qsub -P short -e $ERROR_PATH -o $OUTPUT_PATH -N $JOBNAME_ID ${SCRIPT_PATH}/LaunchJobSimulation.sh ${file} ${SCRIPT_PATH}
	done
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
