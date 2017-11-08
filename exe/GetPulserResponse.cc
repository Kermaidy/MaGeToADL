#include <fstream>
#include <sstream>
#include <string>
#include <typeinfo>
#include <time.h>

// ROOT includes
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TObject.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TApplication.h"

// MGDO includes
#include "MGTWaveform.hh"
#include "MGTEvent.hh"
#include "MGTRun.hh"

// GELATIO includes
#include "GETGERDADigitizerData.hh"

using namespace std;

double AuxSamplingFrequency = 0.100;  //GHz
double SamplingFrequency = 0.025;     //GHz

MGTWaveform *Differentiate(MGTWaveform* kAnInput)
{  
  int fEndSample = kAnInput->GetLength();

  MGTWaveform *fAnOutput = new MGTWaveform(NULL,0,0.1,0.0,MGWaveform::kADC,0);
  fAnOutput->SetLength(fEndSample);
  (*fAnOutput)[0] = 0;

  for(int i = 1;i<fEndSample;i++) (*fAnOutput)[i] = kAnInput->At(i)-kAnInput->At(i-1);
  
  return fAnOutput;
}

std::vector<double> GetTestPulser(int channel,int isAux)
{
  int debugPULSER = 0;
  std::vector<double> testpulser;

  TTree* tree = new TTree();
  MGTEvent* event = new MGTEvent();
  MGTWaveform* wf;
  if(!isAux) wf = new MGTWaveform(NULL,0,SamplingFrequency,0.0,MGWaveform::kADC,0);
  else       wf = new MGTWaveform(NULL,0,AuxSamplingFrequency,0.0,MGWaveform::kADC,0);  

  std::string PulseFilename = "/lfs/l3/gerda/Daq/data-phaseII/blind/v03.00/gen/tier1/ged/cal/run0060/gerda-run0060-20160325T162341Z-cal-ged-tier1.root";

  TFile* RootFile = new TFile(PulseFilename.c_str(),"READ");
  uint64_t time=0;
  
  if(debugPULSER) std::cout << "\r Open pulser data file " << PulseFilename << std::endl;

  tree=(TTree*)RootFile->Get("MGTree");
  tree->SetBranchAddress("event", &event);

  if(debugPULSER) std::cout << "Get number of entries " << std::endl;

  tree->GetEntry();
  
  int Nentries = tree->GetEntries();
  int Npulser = 100;
  int iter = 0;
  double Norm = 0.;

  if(debugPULSER) std::cout << "Number of entries " << Nentries << std::endl;
  if(debugPULSER) std::cout << " " << std::endl;
  
  for(int i=0;i<Nentries;i++){
    tree->GetEntry(i);
    if(debugPULSER) std::cout << "\r Search for pulser data : " << i << " " << event->GetEventType() << std::flush;
    if(event->GetEventType() == 1){ //If test pulser
      std::cout << "\r   #" << iter+1 << "/" << Npulser << std::flush;
      if(!isAux && iter == 0) wf->operator=(*event->GetWaveformID(channel));
      else if(!isAux)         wf->operator+=(*event->GetWaveformID(channel));
      else if(iter == 0)      wf->operator=(*event->GetAuxWaveformID(channel));
      else                    wf->operator+=(*event->GetAuxWaveformID(channel));
      iter++;
      if(iter == Npulser) break;
    }
  }

  RootFile->Close();

  if(debugPULSER) std::cout << "Fill pulser data size : " << wf->GetLength() << std::endl;
  if(debugPULSER) std::cout << "Fill pulser data in vector " << std::endl;
  
  (*wf) *= -1.;       // inverse the pulser wf
  (*wf) -= (*wf)[0];  // remove the wf baseline

  
  wf->operator=(*Differentiate(wf));
  testpulser = wf->GetVectorData();
  
  if(debugPULSER) std::cout << "Normalize filtering function" << std::endl;
  // Normalize filtering function
  for( size_t k = 0; k < testpulser.size(); k ++ ) if(Norm<testpulser[k]) Norm = testpulser[k];//Norm += testpulser[k];
  for( size_t k = 0; k < testpulser.size(); k ++ ) testpulser[k] /= Norm;

  if(debugPULSER){
  //if(channel == 4 && isAux){    
    TGraph* gr = new TGraph();
    for(int i = 0;i<testpulser.size();i++) gr->SetPoint(i,i*10,testpulser[i]);
    TApplication* myApp = new TApplication("mAapp",0,0);
    TCanvas* c = new TCanvas("c","Pulser sample",700,500);
    c->cd();
    gr->Draw();

    myApp->Run();
  }

  if(debugPULSER) std::cout << "Return pulser " << std::endl;
  return testpulser;
  
}

int main(int argc, const char* argv[])
{
  std::cout << " " << std::endl;
  std::cout << "Enter in GetPulserResponse " << std::endl;
  std::cout << " " << std::endl;

  std::vector<std::vector<double> > Pulser, auxPulser;

  std::cout << "Fill data in vectors " << std::endl;
  for(int channel = 0;channel<40;channel++){
    std::cout << "\r Get pulser response from channel " << channel << std::endl;
    Pulser.push_back(GetTestPulser(channel,0));
    auxPulser.push_back(GetTestPulser(channel,1));
  }

  std::cout << "\r Store vectors in ROOT file" << std::endl;
  std::cout << "Open ROOT file" << std::endl;
  std::string OutFilename = "/lfs/l3/gerda/kermaidy/Analysis/software/src/MaGeToADL/PulserLibrary.root";
  TFile* OutRootFile = new TFile(OutFilename.c_str(),"RECREATE");

  TTree* Tree = new TTree();
  Tree->Branch("LF", &Pulser);
  Tree->Branch("HF", &auxPulser);

  std::cout << "Fill TTree" << std::endl;
  Tree->Fill();
  std::cout << "Write TTree" << std::endl;
  Tree->Write("Pulser",TObject::kWriteDelete);
  std::cout << "Normal termination" << std::endl;
  OutRootFile->Close();
}
