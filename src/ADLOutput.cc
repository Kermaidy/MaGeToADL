// include files for ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TObject.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TVector2.h"

#include <sstream>
#include <string>
#include <typeinfo>
#include <time.h>

//---------------------------------------------------------------------------//

#include "ADLOutput.hh"
using namespace std;

//---------------------------------------------------------------------------//
ADLOutput::ADLOutput():
      fSimulateTraces(true),
      fRecordADLOutPos(true),
      fRecordADLTraces(true)
{
  fOutputFile = 0;
  MGTree = 0;
  theRun = 0;
  event = 0;
  waveform = 0;
  auxwaveform = 0;
  digiData = 0;
}

ADLOutput::~ADLOutput()
{
 
}

void ADLOutput::DefineSchema(string outputrootfilename)
{
      setupADLdetPos = 0;
      
      // Set waveform amplitudes/noise
      FEP_kev = 2.6145;
      Baseline = 5604.0;
      FEP_ADC = 6046;
      RMS_noise = 5.3;

      // ADL detector initialization
      ADLDetector = new  ADLDetectorTrace(0,0);

      cout << "ADL detector created " << endl;

      for(int channel=3; channel < NDET; channel++){
	cout << "\r Setup channel " << channel << flush;
        ADLDetector->SetSetupFile(channel);
        ADLDetector->ConfigureADL(ADLDetector->GetSetupFile()); // ADL-4.2
      }
    
      fOutputFile = new TFile(outputrootfilename.c_str(),"RECREATE");

        // Tier1 tree definition
      MGTree = new TTree("MGTree","MGTree");
      
      // Waveform branch
      MGTree->Branch("event",&event, 32000,-99);
      
      // Run info
      MGRunType::RunType fRunType;
      fRunType = MGRunType::kData;
      std::string theDAQ = "Struck";
      int fRunNumber = 1;
      time_t fStartTime = 1479394213;
      time_t fStopTime = 1479394213.000160000;
      std::string theRunDescription = "No description";
      
      theRun = new MGTRun();
      theRun->SetRunType(fRunType);
      theRun->SetRunNumber(fRunNumber);
      theRun->SetRunDescription(theRunDescription);
      theRun->SetStartTime(fStartTime);
      theRun->SetStopTime(fStopTime);
      //The Stop time is available at the end of run
      theRun->SetParentDAQLabel(theDAQ);
      
      // Set MGDO digitizer
      fPreTrigger = 100;
      fTimestamp = 1479394213;
      fDecimalTimestamp = 0;
      IsPulseInverted = false;
      TriggerNumber = 1;
      fEventNumber = 1;
      fMuVetoed = false;
      fMuVetoSample = 0;
      wfLength = 4096;
      wfPreTrigger = 1900;
      auxwfLength = 1000;
      auxwfPreTrigger = 300;
      
      theTime = 1479394213;
      SamplingFrequency = 0.025; // GHz
      AuxSamplingFrequency = 0.1; // GHz
      iChannel = 0;
      validChannelCounter = 0;
}

TFile* ADLOutput::RunSimulation(string InputFilename){
    
  std::cout << " " << std::endl;
  std::string treename = "fTree";
  TTree* Tree = 0;
  
  std::cout<< " " <<std::endl;

  TFile* fInputFile = new TFile(InputFilename.c_str(),"READ");
  
  if(fInputFile->GetListOfKeys()->Contains(treename.c_str()))
    {
      fInputFile->GetObject(treename.c_str(),Tree);
      
      if(Tree->GetNbranches() > 0){
	
	//std::cout << "Number of branches in " << treename << " tree is " << Tree->GetNbranches() << std::endl;
	
	Tree->SetBranchAddress("hits_totnum",&hits_totnum);
	Tree->SetBranchAddress("hits_tote",&hits_tote);
	Tree->SetBranchAddress("hits_iddet",hits_iddet);
	Tree->SetBranchAddress("hits_edep",hits_edep);
	Tree->SetBranchAddress("hits_xpos",hits_xpos);
	Tree->SetBranchAddress("hits_ypos",hits_ypos);
	Tree->SetBranchAddress("hits_zpos",hits_zpos);
        
	clock_t start, end;
	double cpu_time_used;

	start = clock();

	//	for(int i =0;i<Tree->GetEntries();i++){
	for(int i =0;i<1001;i++){
	  if(i % 1000 == 0){
	    end = clock();
	    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	    cout << i << " events treated in " << cpu_time_used << " seconds" << endl;
	  }
	  Tree->GetEntry(i);
	  
	  //std::cout << "\r   Event : " << i+1 << "/" << Tree->GetEntries() 
	  //<< " containg " << hits_totnum << " hits" 
	  //<< " with " << hits_tote << " keV" << std::flush;
	  
	  // Set MGDO event structure
	  event = new MGTEvent();
	  MGEventType::EventType fEventType;
	  fEventType = MGEventType::kBaseline;
	  event->SetAuxWaveformArrayStatus(true);
	  event->SetTotalNumberOfIDs(NDET);
	  event->SetAllActiveIDs(true);
	  event->InitializeArrays("GETGERDADigitizerData",1);
	  event->SetBypassStreamerArray(kFALSE); //! Event type info
	  event->SetEventType(fEventType);
	  event->SetTime(theTime);
	  event->SetWFEncScheme(MGTWaveform::kDiffVarInt);
	  event->SetAuxWFEncScheme(MGTWaveform::kDiffVarInt);
	  event->SetETotal(hits_tote);
	  validChannelCounter = 0;
	  
	  for(int j=0;j<NDET;j++) traceCalculated[j] = 0;
	  
	  //Simulate pulses with ADL once per event
	  for (int j=0; j<hits_totnum; j++) {
                    if(hits_iddet[j] >=3 && traceCalculated[hits_iddet[j]] == 0){
		      SimulatePulse(hits_iddet[j]);
		      traceCalculated[hits_iddet[j]] = 1;
                    }
	  }
	}
	std::cout << " " << std::endl;
	std::cout << "\r Simulation finished " << std::endl;
      }
      else std::cout << "Tree " << treename << " is empty" << std::endl;
    }
  else{
    std::cout << "Tree " << treename << " not found in  " << InputFilename << std::endl;
    std::cout << "Exit" << std::endl;
  }
  fInputFile->Close();
  
  return fOutputFile;
}

void ADLOutput::SimulatePulse(int channel){
 
  int debugADL = 0;

  double ETotDet = 0;

  if(debugADL) std::cout << "DEBUG: Simulate pulse in channel " << channel << std::endl;

  ADLDetector->SetSetupFile(channel);
  if(debugADL) std::cout << "DEBUG: Setup potentials" << std::endl;
  ADLDetector->SetPotentials(ADLDetector->GetSetupFile()); // ADL-4.2
  if(debugADL) std::cout << "DEBUG: Potentials set" << std::endl;
  if(channel == 9 || channel == 14 ||channel == 16 ||channel == 20 ||channel == 22 ||channel == 33) // Consider inverted BEGe
    ADLDetector->SetPositionOffset(-1);
  else ADLDetector->SetPositionOffset(1);

  ADLDetector->CreateADLevent();
  if(debugADL) std::cout << "DEBUG: ADL event created" << std::endl;

  if(fRecordADLOutPos)
    ETotDet = ADLDetector->SetADLhits(hits_totnum,hits_edep,hits_xpos,hits_ypos,hits_zpos,hits_iddet,hits_ADLpos,hits_isOut);
  else
    ETotDet = ADLDetector->SetADLhits(hits_totnum,hits_edep,hits_xpos,hits_ypos,hits_zpos,hits_iddet);
  if(debugADL) std::cout << "DEBUG: ADL hits set" << std::endl;
  if(debugADL) std::cout << "DEBUG: Deposited energy in channel " << channel << " is " << ETotDet << "/" << hits_tote << " MeV" << std::endl;

  if(hits_tote > 1.0){

    if(ADLDetector->CalculateTrace(ADLDetector->GetSetupFile())) std::cerr<< "Failed to calculate trace" <<std::endl;
    if(debugADL) std::cout << "DEBUG: ADL calculate trace" << std::endl;
    
    MGWaveformTag::EWaveformTag fWaveformTag = MGWaveformTag::kNormal;
    
    //Fill digiData. 
    digiData = new ((*(event->GetDigitizerData()))[validChannelCounter]) GETGERDADigitizerData(); 
    digiData->SetClockFrequency(SamplingFrequency);
    
    digiData->SetPretrigger(fPreTrigger);
    digiData->SetTimeStamp(fTimestamp);
    digiData->SetDecimalTimeStamp(fDecimalTimestamp);
    digiData->SetIsInverted(IsPulseInverted);
    digiData->SetTriggerNumber(TriggerNumber);
    digiData->SetID(channel);
    digiData->SetEnergy(ETotDet);
    digiData->SetEventNumber(eventnumber);
    digiData->SetIsMuVetoed(fMuVetoed);
    digiData->SetMuVetoSample(fMuVetoSample);
    digiData->SetWaveformTag(fWaveformTag);
    if(debugADL) std::cout << "DEBUG: ADL digitizer set" << std::endl;
    
    // Set MGDO waveform 
    waveform = new ((*(event->GetWaveforms()))[validChannelCounter]) MGTWaveform(NULL,0,SamplingFrequency,0.0,MGWaveform::kADC,0);
    waveform->SetLength(wfLength);
    waveform->SetTOffset(0.);
    waveform->SetID(channel);
    
    ADLDetector->SetWaveformAttribute(wfPreTrigger,Baseline,FEP_ADC/FEP_kev,RMS_noise);
    if(debugADL) std::cout << "DEBUG: Waveform attribute set" << std::endl;
    
    if(ADLDetector->SetADLWaveform(waveform)) std::cerr<< "Failed to set waveform" <<std::endl;
    if(debugADL) std::cout << "DEBUG: Waveform set" << std::endl;
    
    if (event->GetAuxWaveformArrayStatus()){
      auxwaveform = new ((*(event->GetAuxWaveforms()))[validChannelCounter])
	MGTWaveform(NULL,0,AuxSamplingFrequency,0.0,MGWaveform::kADC,0);     
      auxwaveform->SetSamplingFrequency(AuxSamplingFrequency);
      auxwaveform->SetID(channel);
      auxwaveform->SetTOffset(0.);      
      auxwaveform->SetLength(auxwfLength);
      
      ADLDetector->SetAuxWaveformAttribute(auxwfPreTrigger,Baseline,FEP_ADC/FEP_kev,RMS_noise);
      if(ADLDetector->SetADLauxWaveform(auxwaveform)) std::cerr<< "Failed to set aux waveform" <<std::endl;
    }

    validChannelCounter++;

    if(fRecordADLTraces){
      double** ePath = ADLDetector->GetElectronPath();  // Get matrix containing e- path in (x,y,z) coord.
      double** hPath = ADLDetector->GetHolePath();      // Get matrix containing h  path in (x,y,z) coord.
      trace_totnum = ADLDetector->GetTraceDim();
      trace_iddet = channel;
      for(int i = 0;i<trace_totnum;i++){
          trace_xposE[i] = ePath[i][1] - ADLDetector->GetCenter(); // center e/h traces around 0 in the (x,y) plane
          trace_yposE[i] = ePath[i][2] - ADLDetector->GetCenter();
          trace_zposE[i] = ePath[i][3];
          trace_xposH[i] = hPath[i][1] - ADLDetector->GetCenter();
          trace_yposH[i] = hPath[i][2] - ADLDetector->GetCenter();
          trace_zposH[i] = hPath[i][3];
      }
    }
    fOutputFile->cd();
    MGTree->Fill();
    if(debugADL) std::cout << "Print MGTree : " << MGTree->GetEntries() << " events recorded" << std::endl;
  }
  ADLDetector->DeleteADLevent();
}

int ADLOutput::WriteOutputs(TFile* fOutputFile)
{
  fOutputFile->cd();
  if(fOutputFile->IsOpen()){
    MGTree->Write();

    fOutputFile->Close();
    return 1;
  }
  else
    return 0;
}
