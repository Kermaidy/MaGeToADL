#include "ADLOutput.hh"
#include "ADLDetectorTrace.hh"
#include "ADLCluster.hh"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

using namespace std;

int main(int argc, const char* argv[])
{
 string intputrootpath;
 string intputrootfilename;
 string rootfilename;
 string outputrootfilename;

 if(argc>2){
   intputrootpath = argv[1];
   rootfilename = argv[2];
   intputrootfilename = intputrootpath + rootfilename;
   outputrootfilename = "./Tier1/" + rootfilename;
 }
 else{
     cout << "You did not provide enough arguments " << endl;
     cout << "Try: ./SimulatePulse #INPUTFILENAME #TYPE #SIMID " << endl;
     return 0;
 }

 ADLOutput ADLoutput;

 cout << "Initialize the simulation " << endl;
 ADLoutput.DefineSchema(outputrootfilename);
    
 cout << "\r Run the simulation " << endl;
 TFile* fOutputFile = ADLoutput.RunSimulation(intputrootfilename);
 
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
