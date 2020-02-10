
#include "ADLOutput.hh"
#include "ADLDetectorTrace.hh"
#include "ADLCluster.hh"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

using namespace std;

void PrintHelp()
{
  fprintf(stderr,"  \n");
  fprintf(stderr,"    HELP : \n");
  fprintf(stderr,"    To execute:./SimulatePulse -option1 xxx -option2 yyy .... \n");
  fprintf(stderr,"    -i *	 ->input root filename \n");
  fprintf(stderr,"    -ip *      ->input root path (default:.../MaGeToADL/RawData \n");
  fprintf(stderr,"    -op *      ->output root path (default:.../MaGeToADL/RawPulses \n");
  fprintf(stderr,"    -d  *      ->Debug 0=Off, 1=On (Default: 0) \n");
  fprintf(stderr,"    -config *  -> provide a config file (default: default.txt) \n");  
  fprintf(stderr,"    -pos *     -> centering detector position (1=GERDA MaGe(default),2:HADES) \n");
  fprintf(stderr,"    -n *       -> No. of events simulated (0(default):All, else =*) \n");
  fprintf(stderr,"    -h         -> this help \n");

  fprintf(stderr,"  \n");

  exit(1);
}

int main(int argc, const char* argv[])
{
 string inputrootpath="./RawData/";

 string rootfilename;
 string outputrootpath="./RawPulses/";

 string configfile = "config/default.txt";
 int debug =0;
 string Nentries="0";
 int resetPos=1;   //Mage detectors not centred around zero, correct for it.
 int i=0;
 while(++i<argc){
   const char* opt=argv[i];
   
   if(strcmp(opt,"-d")==0)       sscanf(argv[++i],"%d",&debug);      //debug mode 0=off, 1=On 
   else if(strcmp(opt,"-ip")==0) inputrootpath=argv[++i];                   // specify input path
   else if(strcmp(opt,"-op")==0) outputrootpath=argv[++i];                  // specify output path
   else if(strcmp(opt,"-i")==0) rootfilename=argv[++i];                    // specify input filename
   else if(strcmp(opt,"-h")==0)  PrintHelp();                             
   else if(strcmp(opt,"-config")==0) configfile=argv[++i];            //input a config file
   else if(strcmp(opt,"-pos")==0) sscanf(argv[i++],"%d",&resetPos);   //centre the detector position
   else if(strcmp(opt,"-n")==0) Nentries=argv[++i];   //No. of events simulated ("0"(default):All, else:(1/500)th of all)
   else {
   cout << "You did not provide filename and/or configfile " << endl;
//   return 0;
   }
 }
 string inputrootfilename=inputrootpath+rootfilename;
 string outputrootfilename=outputrootpath+rootfilename;


 ADLOutput ADLoutput(debug);

 cout << "Initialize the simulation " << endl;
 ADLoutput.DefineSchema(configfile,outputrootfilename,resetPos);
    
 cout << "\r Run the simulation " << endl;
 TFile* fOutputFile = ADLoutput.RunSimulation(configfile,inputrootfilename,Nentries);
 
 std::cout << " " << std::endl;
 std::cout << "Open outputs ROOT file" << endl;

 if(ADLoutput.WriteOutputs(fOutputFile)){
   cout << "Simulation ended properly " << endl;
   cout << " " << endl;
 }
 else{  
   cout << "Data not written. Exit. " << endl;
   cout << " " << endl;
 }
 return 0;
}
