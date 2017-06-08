#!/bin/bash

SWMOD_PATH=/path/to/swmod

. ${SWMOD_PATH}/bin/swmod.sh init

swmod load root@6.05.02
swmod load clhep@2.1.3.1
swmod load gerda@master

./SimulatePulse $1 $2 $3
