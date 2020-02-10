#ifndef _ADLOUTPUT_HH
#define _ADLOUTPUT_HH

// MGDO includes
#include "MGTWaveform.hh"
#include "MGTEvent.hh"
#include "MGTRun.hh"

// GELATIO includes
#include "GETGERDADigitizerData.hh"

// ADL includes
#include "ADLDetectorTrace.hh"
#include "ADLCluster.hh"

//---------------------------------------------------------------------------//
// ROOT declarations
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TBranch.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"

//---------------------------------------------------------------------------//


class ADLOutput
{
public:
  //default constructor
  ADLOutput(int);

  //destructor
  ~ADLOutput();

  void DefineSchema(std::string,std::string,int);
  TFile* RunSimulation(std::string,std::string,std::string);
  int WriteOutputs(TFile*);

  std::vector<double> x0,y0,z0;
  std::vector<int> channels;
  int nchannel; 

 //private  members
private:

  double SimulatePulse(int);   // ADL-4.2
  void SimulateEmptyWaveform(int);   // simulates empty wf

  static const int MAX_NTRACE=2000;
  static const int MAX_NHITS=1000;
  static const int NDET=40;

  // hits
  Int_t    hits_totnum;
  Float_t  hits_tote;
  Float_t  hits_edep[MAX_NHITS];
  Float_t  hits_xpos[MAX_NHITS];
  Float_t  hits_ypos[MAX_NHITS];
  Float_t  hits_zpos[MAX_NHITS];
  Int_t    hits_iddet[MAX_NHITS];  // which Ge detector this hit is in

  int traceCalculated[NDET]; // Flag to calculate trace in Ge det only once per event

  // hits : ADL informations
  Float_t  hits_ADLpos[MAX_NHITS];
  Int_t  hits_isOut[MAX_NHITS];

  Int_t    trace_totnum;
  Int_t    trace_iddet;
  Int_t    eventnumber;
  Float_t  trace_xposE[MAX_NTRACE];
  Float_t  trace_yposE[MAX_NTRACE];
  Float_t  trace_zposE[MAX_NTRACE];
  Float_t  trace_xposH[MAX_NTRACE];
  Float_t  trace_yposH[MAX_NTRACE];
  Float_t  trace_zposH[MAX_NTRACE];

  bool fSimulateTraces;
  bool fClusterization;
  bool fRecordADLOutPos;
  bool fRecordADLTraces;

  /***********************/
  /*** ADL simulation  ***/
  /***********************/

  ADLDetectorTrace* ADLDetector;

  int debugADL;

  int setupADLdetPos;    // Used to define the detectors positon offset in ADL only for the 1st event
  int testPulser;

  TFile* fOutputFile;
  TTree* MGTree;
  MGTRun* theRun;
  MGTEvent* event;
  MGTWaveform* waveform;
  MGTWaveform* auxwaveform;
  GETGERDADigitizerData* digiData;

  size_t fPreTrigger;
  unsigned long long fTimestamp;
  unsigned long long fDecimalTimestamp;
  bool IsPulseInverted;
  unsigned int TriggerNumber;
  int fEventNumber;
  bool fMuVetoed;
  unsigned int fMuVetoSample;
  int wfLength;
  int wfPreTrigger;
  int auxwfLength;
  int auxwfPreTrigger;

  Double_t theTime;
  double SamplingFrequency; // GHz
  double AuxSamplingFrequency; // GHz
  size_t iChannel;
  int validChannelCounter;

  /***********************/

};

#endif

