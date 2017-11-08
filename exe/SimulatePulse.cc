
#include "ADLOutput.hh"
#include "ADLDetectorTrace.hh"
#include "ADLCluster.hh"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

using namespace std;

int main(int argc, const char* argv[])
{
 string inputrootpath;
 string inputrootfilename;
 string rootfilename;
 string outputrootpath;
 string outputrootfilename;
 string isCal;
 int whichDet;
 int resetPos;

 if(argc>6){
   inputrootpath = argv[1];
   outputrootpath = argv[2];
   rootfilename = argv[3];
   isCal = argv[4];
   whichDet = atoi(argv[5]);
   resetPos = atoi(argv[6]);
   inputrootfilename = inputrootpath + rootfilename;
   outputrootfilename = outputrootpath + rootfilename;
 }
 else{
     cout << "You did not provide enough arguments " << endl;
     cout << "Try: ./SimulatePulse #INPUTPATHNAME #OUTPUTPATHNAME #INPUTFILENAME $ISCAL[0 or 1] #WHICHDET[0=GerdaArray or 1=SAGE_HADES or 2=ORTEC_MPIK] #RESETPOS[0 or 1 if MaGe det. pos. not centered around 0]" << endl;
     return 0;
 }

 ADLOutput ADLoutput;

 cout << "Initialize the simulation " << endl;
 ADLoutput.DefineSchema(outputrootfilename,whichDet,resetPos);
    
 cout << "\r Run the simulation " << endl;
 TFile* fOutputFile = ADLoutput.RunSimulation(inputrootfilename,isCal,whichDet);
 
 std::cout << " " << std::endl;
 std::cout << "Open outputs ROOT file" << endl;

 if(ADLoutput.WriteOutputs(fOutputFile)){
   cout << "Simulation ended properly " << endl;
   cout << " " << std::endl;
 }
 else{  
   cout << "Data not written. Exit. " << endl;
   cout << " " << std::endl;
 }
}
