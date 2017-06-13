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
  
  Noise = GetNoise(0);
  AuxNoise = GetNoise(1);
}

//---------------------------------------------------------------------------//

ADLDetectorTrace::~ADLDetectorTrace() {
}

//---------------------------------------------------------------------------//

void ADLDetectorTrace::SetSetupFile(int channel){

  detector_channel = channel;

  std::ostringstream oss;
  if(channel <10) oss << 0;
  else oss.str("");
  oss << channel;

  char* envPath = getenv("MAGEDIR");

  detector_setupfile = envPath;

  if(channel == 666) detector_setupfile += "/adl/configfiles/Det_HADES/ICOAX.txt";
  else detector_setupfile += "/adl/configfiles/Det_" + oss.str() + "/Det.txt";
}

std::string ADLDetectorTrace::GetSetupFile() {return detector_setupfile;}

void ADLDetectorTrace::ConfigureADL(std::string setupfile)
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
  Height[detector_channel]   = GetSimionHeight();

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
  xcenter  = Center[detector_channel];              //center of detector in cm
  ycenter  = Center[detector_channel];              //center of detector in cm
  height   = Height[detector_channel];              //bottom of detector in cm
  
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

void ADLDetectorTrace::SetWaveformAttribute(double wfpretrigger, double baseline, double amplitude, double rms_noise){
  wfPreTrigger = wfpretrigger;
  Baseline = baseline;
  Amplitude = amplitude;
  RMS_noise = rms_noise;
}

void ADLDetectorTrace::SetAuxWaveformAttribute(double wfpretrigger, double baseline, double amplitude, double rms_noise){
  AuxwfPreTrigger = wfpretrigger;
  AuxBaseline = baseline;
  AuxAmplitude = amplitude;
  AuxRMS_noise = rms_noise;
}

double ADLDetectorTrace::SetADLhits(int hits_totnum, std::vector<double> &hits_edep, std::vector<double> &hits_xpos, std::vector<double> &hits_zpos, std::vector<double> &hits_iddet){
  
  double ETotDet = 0;
  int j = 0;
  
  for(Int_t i = 0;i<hits_totnum;i++){
    
    if(hits_iddet[i] == detector_channel){ // Consider only hits in the given detector.
      //Fill in the Hit Pattern (HP):
      double P0[4]={0,hits_xpos[i] + xcenter,ycenter,inverted*(hits_zpos[i] + height)};
      
      if(IsInDetector(P0)){
	if(debugADL) std::cout << "Hits in detector " << detector_channel << std::endl;
	ETotDet += hits_edep[i];
	ADL_evt->HP.Eint[j]  =hits_edep[i];             //Energy of interaction
	ADL_evt->HP.Pos[j][0]=hits_xpos[i] + xcenter;	  //Position where this interaction occures in the ADL referential
	ADL_evt->HP.Pos[j][1]=ycenter;
	ADL_evt->HP.Pos[j][2]=inverted*(hits_zpos[i] + height);
	
	j++; // Only iterate on points in the detector
      }
      else if(debugADL) std::cout << "Hits not in detector " << detector_channel << std::endl;
      
      if(debugADL) std::cout << "    MaGe ref. hits position : " << hits_xpos[i] << " " << hits_zpos[i] << std::endl;
      if(debugADL) std::cout << "    ADL  ref. hits position : " << P0[1] << " " << P0[3] << std::endl;
    }
  }
  return ETotDet;
}

double ADLDetectorTrace::SetADLhits(int hits_totnum, Float_t*hits_edep, Float_t*hits_xpos, Float_t*hits_ypos, Float_t*hits_zpos, Int_t*hits_iddet){
  
  double ETotDet = 0;
  int j = 0;
  
  for(Int_t i = 0;i<hits_totnum;i++){
    
    if(hits_iddet[i] == detector_channel){ // Consider only hits in the given detector.
      //Fill in the Hit Pattern (HP):
      double P0[4]={0,hits_xpos[i] + xcenter,hits_ypos[i] + ycenter,inverted*(hits_zpos[i] + height)};
      
      if(IsInDetector(P0)){
	if(debugADL) std::cout << "Hits in detector " << detector_channel << std::endl;
	ETotDet += hits_edep[i];
	ADL_evt->HP.Eint[j]  =hits_edep[i];             //Energy of interaction
	ADL_evt->HP.Pos[j][0]=hits_xpos[i] + xcenter;	  //Position where this interaction occures in the ADL referential
	ADL_evt->HP.Pos[j][1]=hits_ypos[i] + ycenter;
	ADL_evt->HP.Pos[j][2]=inverted*(hits_zpos[i] + height);
	
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
    if(hits_iddet[i] == detector_channel){
      //Fill in the Hit Pattern (HP):
      double P0[4]={0,hits_xpos[i] + xcenter,hits_ypos[i] + ycenter,inverted*(hits_zpos[i] + height)};
      hits_ADLpos[i] = GetDetectorPos(P0);
      if(IsInDetector(P0)){
	if(debugADL) std::cout << "Hits in detector " << detector_channel << std::endl;
	hits_isOut[i] = 0;
	ETotDet += hits_edep[i];
	ADL_evt->HP.Eint[j]  =hits_edep[i];             //Energy of interaction
	ADL_evt->HP.Pos[j][0]=hits_xpos[i] + xcenter;	  //Position where this interaction occures in the ADL referential
	ADL_evt->HP.Pos[j][1]=hits_ypos[i] + ycenter;
	ADL_evt->HP.Pos[j][2]=inverted*(hits_zpos[i] + height);
	
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

  gRandom->SetSeed(time(NULL));
  int randomline = int(gRandom->Uniform(0,Noise.size()));

  if(debugADL){
    if(Noise[randomline][detector_channel].size() != waveform->GetLength()){ 
      std::cout<< "DEBUG :  No proper noise found in library for channel " << detector_channel << std::endl; 
      std::cout<< "DEBUG :  Noise size : " <<  Noise.size() << "x" << Noise[0].size() << "x" << Noise[randomline][detector_channel].size() << "/" << waveform->GetLength() << std::endl;
      return 1; }
  }
  if(TraceLength>0){
    for (size_t i=0;i<waveform->GetLength();i++){
      if(i<wfPreTrigger) (*waveform)[i] = 4. * Baseline + Noise[randomline][detector_channel][i];
      else if(i<wfPreTrigger + TraceLength/4 - 5){ 
	(*waveform)[i] = 4. * Baseline + Noise[randomline][detector_channel][i] + Amplitude * ((ADL_evt->TD).Tr[0][adl_iter]
											+ (ADL_evt->TD).Tr[0][adl_iter+1] 
											+ (ADL_evt->TD).Tr[0][adl_iter+2] 
											+ (ADL_evt->TD).Tr[0][adl_iter+3]);
	adl_iter += 4;
      }
     else (*waveform)[i] = 4. * Baseline + Noise[randomline][detector_channel][i] + 4. * Amplitude * (ADL_evt->TD).Tr[0][TraceLength];
    }
  }
  else return 1;
  return 0;
}

int ADLDetectorTrace::SetADLauxWaveform(MGTWaveform* waveform) {
  
  int TraceLength = GetDIMT();
  int adl_iter = 0;

  gRandom->SetSeed(time(NULL));
  int randomline = int(gRandom->Uniform(0,AuxNoise.size()));

  if(debugADL){
    if(AuxNoise[randomline][detector_channel].size() != waveform->GetLength()){ 
      std::cout<< "DEBUG :  No proper aux noise found in library for channel " << detector_channel << std::endl; 
      std::cout<< "DEBUG :  Noise size : " <<  AuxNoise.size() << "x" << AuxNoise[0].size() << "x" << AuxNoise[randomline][detector_channel].size() << "/" << waveform->GetLength() << std::endl;
      return 1; 
    }
  }
  if(AuxNoise[randomline][detector_channel].size() != waveform->GetLength()) return 1;

  if(TraceLength>0){
    for (size_t i=0;i<waveform->GetLength();i++){
      if(i<AuxwfPreTrigger) (*waveform)[i] = AuxBaseline + AuxNoise[randomline][detector_channel][i];
      else if(i<AuxwfPreTrigger + TraceLength){ 
	if(debugADL) std::cout << adl_iter 
				    << " " << (ADL_evt->TD).Tr[0][adl_iter] 
				    << " " << GetNUMRES_XYZh()[adl_iter][1] 
				    << " " << GetNUMRES_XYZh()[adl_iter][2] 
				    << " " << GetNUMRES_XYZh()[adl_iter][3] 
				    << std::endl;
	(*waveform)[i] = AuxBaseline + AuxNoise[randomline][detector_channel][i] + AuxAmplitude * (ADL_evt->TD).Tr[0][adl_iter];
	adl_iter ++;
      }
      else (*waveform)[i] = AuxBaseline + AuxNoise[randomline][detector_channel][i] + AuxAmplitude * (ADL_evt->TD).Tr[0][TraceLength];
    }
  }
  else return 1;
  return 0;
}

/*
std::vector<std::vector<std::vector<double> > > ADLDetectorTrace::GetNoise(int isAux)
{
  std::vector<std::vector<std::vector<double> > > Data;
  std::vector<std::vector<double> > VectorTmp2;
  std::vector<double>* VectorTmp = 0;
  std::vector<double> VectorTmp1;

  string treename = "noiseTree";

  for(int channel = 0;channel<40;channel++){
    ostringstream oss;
    if(channel < 10) oss << 0;
    oss << channel;

    VectorTmp2.reserve(2000);

    string branchname = "ch" + oss.str();
    if(isAux) branchname += "_aux";
    
    for(int filenumber = 0;filenumber<200;filenumber++){
      oss.str("");
      oss << filenumber;
      
      string noisefilename = "/lfs/l3/gerda/akirsch/gerda-BEGesimulation/run-58-59-60-61-62-63-64_Results/output/noise_library/noiseLibrary_all_";
      noisefilename += oss.str() + ".root";
      
      TFile* fNoiseFile = new TFile(noisefilename.c_str(),"READ");
      TTree* NoiseTree = 0;
      TBranch* BranchTmp = 0;
            
      if(fNoiseFile->GetListOfKeys()->Contains(treename.c_str()))
      	{
	  std::cout << "\r Get noise from channel " << channel << " and file " <<  filenumber << std::flush;
	  
	  if(debugADL)  std::cout << "Noise library " << noisefilename << " opened" << std::endl;
	  fNoiseFile->GetObject(treename.c_str(),NoiseTree);
	  if(NoiseTree->GetListOfBranches()->Contains(branchname.c_str())){
	   	    
	    NoiseTree->SetBranchAddress(branchname.c_str(),&VectorTmp,&BranchTmp);

	    NoiseTree->GetEntry(0);
	    VectorTmp2.push_back((*VectorTmp));
	    VectorTmp->clear();
	    
	    //for(int line = 0;line<NoiseTree->GetEntries();line++){
	    //NoiseTree->GetEntry(line);
	    //VectorTmp2.push_back((*VectorTmp));	  	      
	    //VectorTmp->clear();
	    //}
	  }
	  else{
	    VectorTmp->push_back(0);
	    VectorTmp2.push_back((*VectorTmp));	  
	    VectorTmp->clear();
	    k++;  
	    if(debugADL) std::cout << "\r Branch " << branchname << " not found in  " << treename << " tree" << std::endl;
	  }
	}
      else std::cout << "\r Tree " << treename << " not found in  " << noisefilename << std::endl;
      fNoiseFile->Close();
    }
    Data.push_back(VectorTmp2);
    VectorTmp2.clear();
  }
  
  return Data;
}
*/

std::vector<std::vector<std::vector<double> > > ADLDetectorTrace::GetNoise(int isAux)
{
  std::vector<std::vector<std::vector<double> > > Data;
  std::vector<std::vector<double> > VectorTmp2;
  std::vector<double>* VectorTmp = 0;
  std::vector<double> VectorTmp1;

  string treename = "noiseTree";

  for(int filenumber = 0;filenumber<1;filenumber++){
    ostringstream oss;
    oss << filenumber;
    
    string noisefilename = "/lfs/l3/gerda/akirsch/gerda-BEGesimulation/run-58-59-60-61-62-63-64_Results/output/noise_library/noiseLibrary_all_";
    noisefilename += oss.str() + ".root";
    
    TFile* fNoiseFile = new TFile(noisefilename.c_str(),"READ");
    TTree* NoiseTree = 0;
    TBranch* BranchTmp = 0;
    
    VectorTmp2.reserve(200);

    if(fNoiseFile->GetListOfKeys()->Contains(treename.c_str()))
      {
	for(int channel = 0;channel<40;channel++){
	  oss.str("");
	  if(channel < 10) oss << 0;
	  oss << channel;
	  	  
	  int k = 0;
	  string branchname = "ch" + oss.str();
	  if(isAux) branchname += "_aux";
	  
	  std::cout << "\r Get noise from channel " << channel << " and file " <<  filenumber << std::flush;
	  
	  if(debugADL)  std::cout << "Noise library " << noisefilename << " opened" << std::endl;
	  fNoiseFile->GetObject(treename.c_str(),NoiseTree);
	  if(NoiseTree->GetListOfBranches()->Contains(branchname.c_str())){
	   	    
	    NoiseTree->SetBranchAddress(branchname.c_str(),&VectorTmp,&BranchTmp);

	    NoiseTree->GetEntry(0);
	    VectorTmp2.push_back((*VectorTmp));	  	      
	    VectorTmp->clear();
	    /*
	    for(int line = 0;line<NoiseTree->GetEntries();line++){
	      NoiseTree->GetEntry(line);
	      VectorTmp2.push_back((*VectorTmp));	  	      
	      VectorTmp->clear();
	    }
	    */
	  }
	  else{
	    VectorTmp->push_back(0);
	    VectorTmp2.push_back((*VectorTmp));	  
	    VectorTmp->clear();
	    k++;  
	    if(debugADL) std::cout << "\r Branch " << branchname << " not found in  " << treename << " tree" << std::endl;
	  }
	}
	Data.push_back(VectorTmp2);
	VectorTmp2.clear();
      }
    else std::cout << "\r Tree " << treename << " not found in  " << noisefilename << std::endl;
  fNoiseFile->Close();
  }

  
  return Data;
}

int ADLDetectorTrace::GetTraceDim(){return GetDIMT();}
int ADLDetectorTrace::GetCenter(){return GetSimionCenter();}
double** ADLDetectorTrace::GetElectronPath(){return GetNUMRES_XYZe();}
double** ADLDetectorTrace::GetHolePath(){return GetNUMRES_XYZh();}
