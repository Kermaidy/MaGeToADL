#!/bin/bash

SWMOD_PATH=/Path/To/swmod

. ${SWMOD_PATH}/bin/swmod.sh init

swmod load root@6.05.02
swmod load clhep@2.1.3.1
swmod load gerda@master

cd $2

###
#
# SIMULATION OF THE FULL GERDA ARRAY OR SINGLE DETECTOR
#

./SimulatePulse $4/ $5/ $1 0 $3 0


