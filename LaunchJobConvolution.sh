#!/bin/bash

SWMOD_PATH=/lfs/l3/gerda/kermaidy/Analysis/software/src/swmod

. ${SWMOD_PATH}/bin/swmod.sh init

swmod load root@6.05.02
swmod load clhep@2.1.3.1
swmod load gerda@master

cd $2

#2-poles circuit like convolution wo Noise
./ConvolutePulses RawPulses_ORTEC/$1 Tier1_ORTEC/$1 $3 2 4 10 100 3

#RC circuit like convolution
#./ConvolutePulses RawPulses/$1 Tier1/$1 $3 1 0 -1 -1 0

#RC circuit like convolution
#./ConvolutePulses RawPulses/$1 Tier1/$1 $3 1 0

#RC circuit like convolution for CCD study wo Noise
#./ConvolutePulses RawPulses_CCD_0/$1 Tier1_CCD_0/$1 $3 0 0

#Test pulser from data convolution 
#./ConvolutePulses RawPulses/$1 Tier1/$1 $3 1 1

#RC circuit like convolution for BKG simulation
#./ConvolutePulses RawPulses_BKG_$4/$1 Tier1_BKG_$4/$1 $3 1 2 0.065
