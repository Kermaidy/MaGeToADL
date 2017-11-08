MaGeToADL manual

1.	Pre-requisite:
- MGDO / GELATIO / ROOT / SWMOD
- Git clone ADL4 @ mppmu GitHub (https://github.com/mppmu/ADL4.git) 
and compile it


2.	Git clone MaGeToADL @ Kermaidy GitHub account
https://github.com/Kermaidy/MaGeToADL.git


3.	Add to your ${HOME}/.bashrc
- export ADLDIR="/path/to/adl-4.2/examples/GERDA"
- export MAGETOADLDIR="/path/to/MaGeToADL"
- export LD_LIBRARY_PATH="${MAGETOADLDIR}/lib:$LD_LIBRARY_PATH" 



4.	- Set SWMOD_PATH in LaunchExecModuleIni.sh, LaunchJobConvolution.sh and LaunchJobSimulation.sh
- Set SCRIPT_PATH in JobScript.sh
- Set GERDA, GELATIO and MGDO paths in the Makefile


5.	Create following repositories in /Path/To/MaGeToADL:
- bin/ build/ lib/ JobErrors/ JobOutput/  -> Output of compilation and jobs
- RawData/ 	 -> output ROOT files from MaGe
- RawPulses/ -> output of ADL raw pulses processing files
- Tier1/		 -> output of pulses convolution – Noise – Signal decay – E.R.
- Tier2/		 -> output of GELATIO 


6.	Compile MaGeToADL: “make”


7.	Check for individual pieces of code (should return informative message):
- ./SimulatePulses
- ./ConvolutePulses
- execModuleIni
