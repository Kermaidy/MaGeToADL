#!/bin/bash

SCRIPT_PATH=/lfs/l3/gerda/kermaidy/Analysis/software/src/MaGeToADL
ERROR_PATH=${SCRIPT_PATH}/JobErrors/
OUTPUT_PATH=${SCRIPT_PATH}/JobOutput/

LAUNCH_DATE=`echo $(date +"%Y%m%d%H%M")`
LAUNCH_DATE=${LAUNCH_DATE}
JOBNAME=submit_at_$LAUNCH_DATE

cd ${SCRIPT_PATH}

ISCAL=$2
DET=$3
SOURCE=$4
SUBDIR=$5

#
# If DET = 0, simulate pulses for the whole GERDA array, if not, only a single detector
#

#TIER2PATH=${SCRIPT_PATH}/Tier2_CCD_0
#TIER1PATH=${SCRIPT_PATH}/Tier1_CCD_0
if [ ${DET} = "ORTEC" ];
then
	TIER2PATH=${SCRIPT_PATH}/Tier2_ORTEC
	TIER1PATH=${SCRIPT_PATH}/Tier1_ORTEC
else
	TIER2PATH=${SCRIPT_PATH}/Tier2
	TIER1PATH=${SCRIPT_PATH}/Tier1
fi


declare -i iter

if [ $1 = "1" ];
then
	if [ ${DET} = "ORTEC" ];
	then
		cd RawData_ORTEC/
		files=$(find -maxdepth 1 -name "*.root")
	elif [ ${ISCAL} = "0" ];
	then
		cd RawData/
		files=$(find -maxdepth 1 -name "Gerda*nbb_*.root")
	elif [ ${ISCAL} = "2" ];
	then
		cd /lfs/l3/gerda/kermaidy/Analysis/GERDA_ARRAY_BKG_SIMULATION_Analysis/${SOURCE}/${SUBDIR}
		files=$(find -maxdepth 1 -name "*${SOURCE}*.root")
	else
		cd RawData/
		files=$(find -maxdepth 1 -name "Gerda_*.root")
    	fi
	iter=0
	for file in ${files}
        do
	    iter=iter+1
	    JOBNAME_ID=${JOBNAME}_$iter
	    echo "Processing simulation : ${iter} ${file}"
	    qsub -P short -e $ERROR_PATH -o $OUTPUT_PATH -N $JOBNAME_ID ${SCRIPT_PATH}/LaunchJobSimulation.sh ${file} ${SCRIPT_PATH} ${ISCAL} ${DET} ${SOURCE} ${SUBDIR}
	done
elif [ $1 = "2" ];
then
	if [ ${DET} = "ORTEC" ];
        then
		cd RawPulses_ORTEC/
		files=$(find -maxdepth 1 -name "*.root")
        elif [ ${ISCAL} = "0" ];
        then
	        cd RawPulses/
            	files=$(find -maxdepth 1 -name "Gerda*nbb_*.root")
        elif [ ${ISCAL} = "2" ];
        then        
		mkdir /lfs/l3/gerda/kermaidy/Analysis/software/src/MaGeToADL/Tier1_BKG_${SOURCE}
		cd RawPulses_BKG_${SOURCE}/
            	files=$(find -maxdepth 1 -name "*.root")
	else
                cd RawPulses/
	    	files=$(find -maxdepth 1 -name "Gerda_*.root")
        fi
	iter=0
	for file in ${files}
        do
            iter=iter+1
            JOBNAME_ID=${JOBNAME}_$iter
            echo "Processing simulation : ${iter} ${file}"
            qsub -e $ERROR_PATH -o $OUTPUT_PATH -N $JOBNAME_ID ${SCRIPT_PATH}/LaunchJobConvolution.sh ${file} ${SCRIPT_PATH} ${iter} ${SOURCE}
        done
elif [ $1 = "3" ];
then
	iter=0
	if [ -n "${SOURCE}" ];
	then
		mkdir /lfs/l3/gerda/kermaidy/Analysis/software/src/MaGeToADL/Tier2_BKG_${SOURCE}
		cd ${TIER1PATH}_BKG_${SOURCE}
	else
		cd ${TIER1PATH}
	fi
	files=$(find -maxdepth 1 -name "*.root")

	for file in ${files}
	do
    		iter=iter+1
    		JOBNAME_ID=${JOBNAME}_$iter
    		echo "Processing file : ${file}"
    		qsub -e $ERROR_PATH -o $OUTPUT_PATH -N $JOBNAME_ID ${SCRIPT_PATH}/LaunchExecModuleIni.sh ${file} ${SCRIPT_PATH} ${SOURCE}
  	done
fi
