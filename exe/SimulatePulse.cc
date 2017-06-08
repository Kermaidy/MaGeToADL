#include "ADLOutput.hh"
#include "ADLDetectorTrace.hh"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

using namespace std;

int main(int argc, const char* argv[])
{
 string intputrootfilename;
 string outputrootfilename;
 string type;
 string simID;

 if(argc>3){
     intputrootfilename = argv[1];
     type = argv[2];
     simID = argv[3];
	outputrootfilename = "./Tier1/" + type + "_" + simID + "-Tier1.root";
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
