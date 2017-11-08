#!/bin/bash

FILE=$1
SOURCE=$3

SWMOD_PATH=/lfs/l3/gerda/kermaidy/Analysis/software/src/swmod

. ${SWMOD_PATH}/bin/swmod.sh init

swmod load root@6.05.02
swmod load clhep@2.1.3.1
swmod load gerda@master

cd $2

echo "execModuleIni -o ${FILE%-*}-Tier2.root summary_gerda-like.ini ${FILE}"
if [ -n "${SOURCE}" ];
then
	execModuleIni -o Tier2_BKG_${SOURCE}/${FILE%-*}-Tier2.root summary_gerda-like.ini Tier1_BKG_${SOURCE}/${FILE}
else
#	execModuleIni -o Tier2/${FILE%-*}-Tier2.root summary_gerda-like.ini Tier1/${FILE}
#	execModuleIni -o Tier2_CCD_00025/${FILE%-*}-Tier2.root summary_gerda-like.ini Tier1_CCD_00025/${FILE}
	execModuleIni -o Tier2_ORTEC/${FILE%-*}-Tier2.root summary_gerda-like.ini Tier1_ORTEC/${FILE}
fi
