#!/bin/bash

FILE=$1

SWMOD_PATH=/Path/To/swmod

. ${SWMOD_PATH}/bin/swmod.sh init

swmod load root@6.05.02
swmod load clhep@2.1.3.1
swmod load gerda@master

cd $2

echo "execModuleIni -o $3/${FILE%-*}-Tier2.root summary_gerda-like.ini $4/${FILE}"
execModuleIni -o $3/${FILE%-*}-Tier2.root summary_gerda-like.ini $4/${FILE}
