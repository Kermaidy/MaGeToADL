#!/bin/bash

SCRIPT_PATH=/lfs/l1/gerda/vibothe/software/src/MaGeToADL
ERROR_PATH=${SCRIPT_PATH}/JobErrors/
OUTPUT_PATH=${SCRIPT_PATH}/JobOutput/

LAUNCH_DATE=`echo $(date +"%Y%m%d%H%M")`
LAUNCH_DATE=${LAUNCH_DATE}
JOBNAME=submit_at_$LAUNCH_DATE

cd ${SCRIPT_PATH}
TIER=$1
CONFIGFILE=${SCRIPT_PATH}/config/BEGE.txt

RAWDATAPATH=${SCRIPT_PATH}/RawData
RAWPULSESPATH=${SCRIPT_PATH}/RawPulses
TIER1PATH=${SCRIPT_PATH}/Tier1
TIER2PATH=${SCRIPT_PATH}/Tier2

declare -i iter

if [ ${TIER} = "0" ];
then
	cd ${RAWDATAPATH}
	files=$(find -maxdepth 1 -name "*.root")

	iter=0
	for file in ${files}
        do
	    iter=iter+1
	    JOBNAME_ID=${JOBNAME}_$iter
	    echo "Processing simulation : ${iter} ${file}"
	    qsub -P short -e $ERROR_PATH -o $OUTPUT_PATH -N $JOBNAME_ID ${SCRIPT_PATH}/LaunchGenericJobSimulation.sh -i ${file} -config ${CONFIGFILE}
#	    sh ${SCRIPT_PATH}/LaunchGenericJobSimulation.sh -i ${file} -config ${CONFIGFILE}
	done
elif [ ${TIER} = "1" ];
then
	cd ${RAWPULSESPATH}
      	files=$(find -maxdepth 1 -name "*.root")

	iter=0
	for file in ${files}
        do
            iter=iter+1
            JOBNAME_ID=${JOBNAME}_$iter
            echo "Processing simulation : ${iter} ${file}"
            qsub -e $ERROR_PATH -o $OUTPUT_PATH -N $JOBNAME_ID ${SCRIPT_PATH}/LaunchGenericJobConvolution.sh -i ${file} 
        done
elif [ ${TIER} = "2" ];
then
	cd ${TIER1PATH}
	files=$(find -maxdepth 1 -name "*.root")

	iter=0

	for file in ${files}
	do
    		iter=iter+1
    		JOBNAME_ID=${JOBNAME}_$iter
    		echo "Processing file : ${file}"
    		qsub -e $ERROR_PATH -o $OUTPUT_PATH -N $JOBNAME_ID ${SCRIPT_PATH}/LaunchGenericExecModuleIni.sh ${file} ${SCRIPT_PATH} ${TIER2PATH} ${TIER1PATH}
  	done
fi
