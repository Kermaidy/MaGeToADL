#!/bin/bash

SWMOD_PATH=/Path/To/swmod

. ${SWMOD_PATH}/bin/swmod.sh init

swmod load root@6.05.02
swmod load clhep@2.1.3.1
swmod load gerda@master

cd $2

#RC circuit like convolution
# Execute "./ConvolutePules" to get informations about options
./ConvolutePulses $4/$1 $5/$1 $3 1 0 -1 -1 0
