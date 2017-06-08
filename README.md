# ADL_exercise
Show how to use ADL with 2 detectors in parallel and how to launch jobs.

1.	Pre-requisite:
- MGDO / GELATIO / ROOT / SWMOD

LET’S SET THE ENVIRONMENT
2.	Create a /path/to/MunichWorkshop repository

3.	Get (https://github.com/mppmu/ADL4.git) and compile ADL4 (make all) in /path/to/MunichWorkshop/ADL4/ repository

4.	Add export LD_LIBRARY_PATH="/path/to/MunichWorkshop/testadl/lib:$LD_LIBRARY_PATH" 
to your ${HOME}/.bashrc

5.	Get testadl folder and do “make” in 
/path/to/MunichWorkshop/testadl/ repository
-> this should produce the executable SimulatePulse

6.	Set SWMOD_PATH in LaunchExecModuleIni.sh and LaunchJobSimulation.sh
And SCRIPT_PATH in JobScript.sh

LET’S COMPUTE 1 SIMULATION
7.	Run ./SimulatePulse
-> this should start a loop of 1000 events and produce Tier1 “ADLTest-Tier1.root” output file

8.	Run sh LaunchExecModuleIni.sh ADLTest-Tier1.root .
-> this should start GELATIO and produce Tier2 “ADLTest-Tier2.root” output file

NOW WITH JOBS
9.	Run sh JobScript 1
-> this should produce 10 Tier1 “ADLTest_#-Tier1.root” output files

10.	Run sh JobScript 2
  -> this should produce 10 Tier2 “ADLTest_#-Tier2.root” output files
