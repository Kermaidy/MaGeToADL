#include "ADLDetectorTrace.hh"
#include "TRandom3.h"
#include <fstream>
#include <stdlib.h>     /* getenv */
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

using namespace std;

//---------------------------------------------------------------------------//
ADLDetectorTrace::ADLDetectorTrace(int debug, int channel){

  debugADL = debug;
  SetADLDebug(debugADL);
  detector_channel = channel;
}

//---------------------------------------------------------------------------//

ADLDetectorTrace::~ADLDetectorTrace() {
}

//---------------------------------------------------------------------------//

void GetCustomDetPosition(int NDET, std::vector<double> &x0, std::vector<double> &y0, std::vector<double> &z0)
{
  // Only 1-detector case

  std::string header, StringTmp, xtmp, ytmp, ztmp;
  std::ifstream File("CustomDetectorPosition.txt");

  if(File)  // Check if the file exists
    {
      getline(File,header);
      File >> StringTmp  >> xtmp >> ytmp >> ztmp;
    }
  else{  // if the file doesn't exist
    cerr << "CustomDetectorPosition.txt not found. Try again" << std::endl;
    exit(1);
  }

  File.close();  // Close file

  for(int i = 0;i<NDET;i++){
    x0.push_back(atof(xtmp.c_str())); 
    y0.push_back(atof(ytmp.c_str())); 
    z0.push_back(atof(ztmp.c_str())); 
  }
}

void GetMaGeDetPosition(std::vector<double> &x0, std::vector<double> &y0, std::vector<double> &z0)
{
  // Full GERDA detector array case

  std::string header, StringTmp, xtmp, ytmp, ztmp;
  std::ifstream File("GERDADetectorPosition.txt");

  if(File)  // Check if the file exists
    {
      getline(File,header);
      while(File.good())  // Loop until reach the eof
        {
          File >> StringTmp  >> xtmp >> ytmp >> ztmp;
          x0.push_back(atof(xtmp.c_str()));  
          y0.push_back(atof(ytmp.c_str()));  
          z0.push_back(atof(ztmp.c_str()));  
        }
    }
  else{  // if the file doesn't exist
    cerr << "GERDADetectorPosition.txt not found. Try again" << std::endl;
    exit(1);
  }

  File.close();  // Close file

  x0.pop_back();
  y0.pop_back();
  z0.pop_back();
}

void ADLDetectorTrace::SetSetupFile(int channel){

  detector_channel = channel;

  std::ostringstream oss;
  if(channel <10) oss << 0;
  else oss.str("");
  oss << channel;

  char* envPath = getenv("MAGETOADLDIR");

  //detector_setupfile = envPath;
  detector_setupfile = "";

  if(channel == 1000){       detector_setupfile += "configfiles/Det_HADES/ICOAX.txt"; detector_channel = 0;}
  else if(channel == 1001) { detector_setupfile += "configfiles/Det_ORTEC/ICOAX.txt"; detector_channel = 0;}
  else detector_setupfile += "configfiles/Det_" + oss.str() + "/BEGe_" + oss.str() + ".txt";
}

std::string ADLDetectorTrace::GetSetupFile() {return detector_setupfile;}

void ADLDetectorTrace::ConfigureADL(std::string setupfile, int resetPos)
{
  if(debugADL) cout<<"\r Setup file : " << setupfile<<endl;

  if(debugADL) cout<<"Setup ADL"<<endl;

  ADLSetup(const_cast<char*>(setupfile.c_str()));

  if(debugADL) cout<<"Setup fields"<<endl;

  SetupFields(const_cast<char*>(setupfile.c_str()));

  if(debugADL) cout<<"Store fields"<<endl;

  ADL_Epot[detector_channel] = GetADLEpot();
  ADL_Wpot[detector_channel] = GetADLWpot();
  ADL_Stru[detector_channel] = GetADLStru();

  GridSize[detector_channel] = GetSimionGridSize();
  Center[detector_channel]   = GetSimionCenter();
  if(detector_channel >= 0)  Height[detector_channel] = 0; // No need to correct for z-offset for single detectors
  else                       Height[detector_channel] = GetSimionHeight()/2.;

  if(resetPos == 1)
    GetMaGeDetPosition(x0,y0,z0); // If MaGe detectors are not centered around zero, correct for this.    
  else if(resetPos == 2)
    GetCustomDetPosition(NDET,x0,y0,z0);
  else{
    x0.assign(NDET,0.);
    y0.assign(NDET,0.);
    z0.assign(NDET,0.);
  }

  if(debugADL) 
    std::cout << detector_channel << "       " << 
    ADL_Epot[detector_channel]->h.nx  << " " << 
    ADL_Epot[detector_channel]->h.ny  << " " << 
    ADL_Epot[detector_channel]->h.nz  << " " <<  
    ADL_Wpot[detector_channel]->h.nx  << " " << 
    ADL_Wpot[detector_channel]->h.ny  << " " << 
    ADL_Wpot[detector_channel]->h.nz  << " " <<
    ADL_Stru[detector_channel]->h.nx  << " " << 
    ADL_Stru[detector_channel]->h.ny  << " " << 
    ADL_Stru[detector_channel]->h.nz  << std::endl;
  
  if(debugADL) StatusFields();
}

void ADLDetectorTrace::SetPotentials(std::string setupfile)
{
  //ADLSetup(const_cast<char*>(setupfile.c_str()));

    SetADLEpot(ADL_Epot[detector_channel]);
    SetADLWpot(ADL_Wpot[detector_channel]);
    SetADLStru(ADL_Stru[detector_channel]);
    
    SetSimionGridSize(GridSize[detector_channel]);
    SetSimionCenter(Center[detector_channel]);
    SetSimionHeight(Height[detector_channel]);
}

void ADLDetectorTrace::SetPositionOffset(double inv){

  inverted = inv;

  gridsize = GetSimionGridSize();
  xcenter  = Center[detector_channel] - x0[detector_channel];              //center of detector in cm
  ycenter  = Center[detector_channel] - y0[detector_channel];              //center of detector in cm
  height   = Height[detector_channel] - z0[detector_channel];              //bottom of detector in cm

  if(debugADL){
    std::cout << "DEBUG: ADL detector gridsize : " << gridsize          << std::endl;
    std::cout << "DEBUG: ADL detector center   : " << GetSimionCenter() << std::endl;
    std::cout << "DEBUG: ADL detector height   : " << GetSimionHeight() << std::endl;
    std::cout << "DEBUG: Position offset set " << std::endl;
    std::cout << "DEBUG: ADL to MaGe xcenter   : " << xcenter << std::endl;
    std::cout << "DEBUG: ADL to MaGe ycenter   : " << ycenter << std::endl;
    std::cout << "DEBUG: ADL to MaGe height    : " << height  << std::endl;
  }
}

void ADLDetectorTrace::CreateADLevent(){
  ADL_evt = new_event();
  ADL_evt->HP.T0= 1.000;               //Time the interaction occurs in the trace
}

void ADLDetectorTrace::DeleteADLevent(){
  delete ADL_evt;
}

void ADLDetectorTrace::SetWaveformTimeOffset(double wfpretrigger,double auxwfpretrigger){
  wfPreTrigger = wfpretrigger;
  AuxwfPreTrigger = auxwfpretrigger;
}

double ADLDetectorTrace::SetADLhits(int hits_totnum, std::vector<double> &hits_edep, std::vector<double> &hits_xpos, std::vector<double> &hits_zpos, std::vector<double> &hits_iddet){
  
  double ETotDet = 0;
  int j = 0;
   for(Int_t i = 0;i<hits_totnum;i++){
    ADL_evt->HP.Eint[j]   = 0;             //Energy of interaction
    ADL_evt->HP.Pos[j][0] = 0;	  //Position where this interaction occures in the ADL referential
    ADL_evt->HP.Pos[j][1] = 0;
    ADL_evt->HP.Pos[j][2] = 0;
  } 
  for(Int_t i = 0;i<hits_totnum;i++){
    
    if(hits_iddet[i] == detector_channel){ // Consider only hits in the given detector.
      //Fill in the Hit Pattern (HP):
      double P0[4]={0,hits_xpos[i] + xcenter,ycenter,inverted*(hits_zpos[i] + inverted*height)};
      
      if(IsInDetector(P0)){
	if(debugADL) std::cout << "Hits in detector " << detector_channel << std::endl;
	ETotDet += hits_edep[i];
	ADL_evt->HP.Eint[j]  =hits_edep[i];             //Energy of interaction
	ADL_evt->HP.Pos[j][0]=hits_xpos[i] + xcenter;	  //Position where this interaction occures in the ADL referential
	ADL_evt->HP.Pos[j][1]=ycenter;
	ADL_evt->HP.Pos[j][2]=inverted*(hits_zpos[i] + inverted*height);
	
	j++; // Only iterate on points in the detector
      }
      else if(debugADL) std::cout << "Hits not in detector " << detector_channel << std::endl;
      
      if(debugADL) std::cout << "    MaGe cluster hits position : " << hits_xpos[i] << " " << hits_zpos[i] << std::endl;
      if(debugADL) std::cout << "    ADL  cluster hits position : " << P0[1] << " " << P0[3] << std::endl;
    }
  }
  return ETotDet;
}

double ADLDetectorTrace::SetADLhits(int hits_totnum, Float_t*hits_edep, Float_t*hits_xpos, Float_t*hits_ypos, Float_t*hits_zpos, Int_t*hits_iddet){
  
  double ETotDet = 0;
  int j = 0;
   for(Int_t i = 0;i<hits_totnum;i++){
    ADL_evt->HP.Eint[j]   = 0;             //Energy of interaction
    ADL_evt->HP.Pos[j][0] = 0;	  //Position where this interaction occures in the ADL referential
    ADL_evt->HP.Pos[j][1] = 0;
    ADL_evt->HP.Pos[j][2] = 0;
  } 
  for(Int_t i = 0;i<hits_totnum;i++){
    
    if(hits_iddet[i] == detector_channel){ // Consider only hits in the given detector.
      //Fill in the Hit Pattern (HP):
      double P0[4]={0,hits_xpos[i] + xcenter,hits_ypos[i] + ycenter,inverted*(hits_zpos[i] + inverted*height)};
      
      if(IsInDetector(P0)){
	if(debugADL) std::cout << "Hits in detector " << detector_channel << std::endl;
	ETotDet += hits_edep[i];
	ADL_evt->HP.Eint[j]  =hits_edep[i];             //Energy of interaction
	ADL_evt->HP.Pos[j][0]=hits_xpos[i] + xcenter;	  //Position where this interaction occures in the ADL referential
	ADL_evt->HP.Pos[j][1]=hits_ypos[i] + ycenter;
	ADL_evt->HP.Pos[j][2]=inverted*(hits_zpos[i] + inverted*height);
	
	j++; // Only iterate on points in the detector
      }
      else if(debugADL) std::cout << "Hits not in detector " << detector_channel << std::endl;
      
      if(debugADL) std::cout << "    MaGe ref. hits position : " << hits_xpos[i] << " " << hits_ypos[i] << " " << hits_zpos[i] << std::endl;
      if(debugADL) std::cout << "    ADL  ref. hits position : " << P0[1] << " " << P0[2] << " " << P0[3] << std::endl;
    }
  }
  return ETotDet;
}

double ADLDetectorTrace::SetADLhits(int hits_totnum, Float_t*hits_edep, Float_t*hits_xpos, Float_t*hits_ypos, Float_t*hits_zpos, Int_t*hits_iddet, Float_t*hits_ADLpos, Int_t*hits_isOut){

  double ETotDet = 0;
  int j = 0;
  for(Int_t i = 0;i<hits_totnum;i++){
    ADL_evt->HP.Eint[j]   = 0;             //Energy of interaction
    ADL_evt->HP.Pos[j][0] = 0;	  //Position where this interaction occures in the ADL referential
    ADL_evt->HP.Pos[j][1] = 0;
    ADL_evt->HP.Pos[j][2] = 0;
  }	

  for(Int_t i = 0;i<hits_totnum;i++){
    if(hits_iddet[i] == detector_channel){
      //Fill in the Hit Pattern (HP):
      double P0[4]={0,hits_xpos[i] + xcenter,hits_ypos[i] + ycenter,inverted*(hits_zpos[i] + inverted*height)};
      hits_ADLpos[i] = GetDetectorPos(P0);
      if(IsInDetector(P0)){
	if(debugADL) std::cout << "Hits in detector " << detector_channel << std::endl;
	hits_isOut[i] = 0;
	ETotDet += hits_edep[i];
	ADL_evt->HP.Eint[j]  =hits_edep[i];             //Energy of interaction
	ADL_evt->HP.Pos[j][0]=hits_xpos[i] + xcenter;	  //Position where this interaction occures in the ADL referential
	ADL_evt->HP.Pos[j][1]=hits_ypos[i] + ycenter;
	ADL_evt->HP.Pos[j][2]=inverted*(hits_zpos[i] + inverted*height);
	
	j++; // Only iterate on points in the detector
      }
      else{
	// Hit not in detector. Set outer position to hit position 
	hits_isOut[i] = 1;
	if(debugADL) std::cout << "Hits not in detector " << detector_channel << std::endl;
	if(debugADL) std::cout << "    Hits position " << hits_ADLpos[i] << std::endl;
	if(debugADL) std::cout << "    MaGe ref. hits position : " << hits_xpos[i] << " " << hits_ypos[i] << " " << hits_zpos[i] << std::endl;
	if(debugADL) std::cout << "    ADL ref. hits position  : " << P0[1] << " " << P0[2] << " " << P0[3] << std::endl;
      }
    }
    else if(hits_isOut[i] < 0) hits_isOut[i] = 2; // If hits_isOut has not been set yet and is not related to detector channel
  }
  return ETotDet;
}

int ADLDetectorTrace::CalculateTrace(){
  if(debugADL) StatusTraces(ADL_evt);
  CalculateTraces(ADL_evt);

  return 0;
}

int ADLDetectorTrace::SetADLWaveform(MGTWaveform* waveform) {

  int TraceLength = GetDIMT();
  int adl_iter = 0;

  if(TraceLength>0){
    for (size_t i=0;i<waveform->GetLength();i++){
      if(i<wfPreTrigger) (*waveform)[i] = 0.;
      else if(i<wfPreTrigger + TraceLength/4 - 5){ 
	(*waveform)[i] = (ADL_evt->TD).Tr[0][adl_iter]
	               + (ADL_evt->TD).Tr[0][adl_iter+1] 
	               + (ADL_evt->TD).Tr[0][adl_iter+2] 
	               + (ADL_evt->TD).Tr[0][adl_iter+3];
	adl_iter += 4;
      }
      else (*waveform)[i] = 4. * (ADL_evt->TD).Tr[0][TraceLength-5];
    }
  }
  else return 1;
  return 0;
}

int ADLDetectorTrace::SetADLauxWaveform(MGTWaveform* waveform) {
  
  int TraceLength = GetDIMT();
  int adl_iter = 0;

  if(TraceLength>0){
    for (size_t i=0;i<waveform->GetLength();i++){
      if(i<AuxwfPreTrigger) (*waveform)[i] = 0.;
      else if(i<AuxwfPreTrigger + TraceLength){ 
	if(debugADL) std::cout << adl_iter 
				    << " " << ADL_evt->HP.Eint[0]
				    << " " << (ADL_evt->TD).Tr[0][adl_iter] 
				    << " " << GetNUMRES_XYZh()[adl_iter][1] 
				    << " " << GetNUMRES_XYZh()[adl_iter][2] 
				    << " " << GetNUMRES_XYZh()[adl_iter][3] 
				    << " " << GetWeight(0,GetNUMRES_XYZh()[adl_iter])
				    << " " << GetWeight(0,GetNUMRES_XYZe()[adl_iter])
				    << std::endl;
	(*waveform)[i] = (ADL_evt->TD).Tr[0][adl_iter];
	adl_iter ++;
	//	std::cout<< "DEBUG :  Trace amplitudes : " <<  waveform->At(i) << " " << AuxAmplitude << " " << (ADL_evt->TD).Tr[0][adl_iter] << std::endl;

      }
      else (*waveform)[i] = (ADL_evt->TD).Tr[0][TraceLength-1];
    }
  }
  else return 1;
  return 0;
}

int ADLDetectorTrace::GetTraceDim(){return GetDIMT();}
int ADLDetectorTrace::GetCenter(){return GetSimionCenter();}
double** ADLDetectorTrace::GetElectronPath(){return GetNUMRES_XYZe();}
double** ADLDetectorTrace::GetHolePath(){return GetNUMRES_XYZh();}
