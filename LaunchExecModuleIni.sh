#!/bin/bash

FILE=$1

SWMOD_PATH=/lfs/l3/gerda/kermaidy/Analysis/software/src/swmod

. ${SWMOD_PATH}/bin/swmod.sh init

swmod load root@6.05.02
swmod load clhep@2.1.3.1
swmod load gerda@master

echo "execModuleIni -o ${FILE%-*}-Tier2.root summary_gerda-like.ini ${FILE}"
execModuleIni -o Tier2/${FILE%-*}-Tier2.root summary_gerda-like.ini Tier1/${FILE}
