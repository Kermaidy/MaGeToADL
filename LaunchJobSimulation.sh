#!/bin/bash

SWMOD_PATH=/lfs/l3/gerda/kermaidy/Analysis/software/src/swmod

. ${SWMOD_PATH}/bin/swmod.sh init

swmod load root@6.05.02
swmod load clhep@2.1.3.1
swmod load gerda@master

cd $2

./SimulatePulse ../mage/macros/ $1
