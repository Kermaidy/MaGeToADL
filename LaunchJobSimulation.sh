#!/bin/bash

SWMOD_PATH=/lfs/l3/gerda/kermaidy/Analysis/software/src/swmod

. ${SWMOD_PATH}/bin/swmod.sh init

swmod load root@6.05.02
swmod load clhep@2.1.3.1
swmod load gerda@master

cd $2

#./SimulatePulse RawData_DL/ RawPulses_DL/ $1 $3 0 0
./SimulatePulse RawData_ORTEC/ RawPulses_ORTEC/ $1 $3 1001 0
#./SimulatePulse /lfs/l3/gerda/kermaidy/Analysis/GERDA_ARRAY_BKG_SIMULATION_Analysis/$5/$6/ RawPulses_BKG_$5/ $1 $3 $4 1


#./SimulatePulse RawData/ RawPulses_CCD_0/ $1 $3 $4 0

###
#
# SIMULATION OF THE FULL GERDA ARRAY OR SINGLE DETECTOR
#

#./SimulatePulse RawData/ RawPulses/ $1 $3 $4 0


