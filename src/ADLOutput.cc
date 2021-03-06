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

#include <fstream>
#include <sstream>
#include <string>
#include <typeinfo>
#include <time.h>

//---------------------------------------------------------------------------//

#include "ADLOutput.hh"
#include "ADLCluster.hh"
#include "ADLDetectorTrace.hh"
using namespace std;

//---------------------------------------------------------------------------//
ADLOutput::ADLOutput(int debug):
  fSimulateTraces(true),
  fClusterization(false),
  fRecordADLOutPos(false),
  fRecordADLTraces(false)
{
  fOutputFile = 0;
  MGTree = 0;
  theRun = 0;
  event = 0;
  waveform = 0;
  auxwaveform = 0;
  digiData = 0;
  debugADL = debug;
  for (int i=0; i<NDET; i++)  traceCalculated[i]=2;  //to not to calculate traces for excluded detectors
}

ADLOutput::~ADLOutput()
{
 
}

void ADLOutput::DefineSchema(string configfile, string outputrootfilename, int resetPos)
{
  setupADLdetPos = 0;
  // ADL detector initialization
  ADLDetector = new  ADLDetectorTrace(debugADL,0);
    
  cout << "\r ADL detector created " << endl;

  channels=ADLDetector->ReadConfigfile(configfile);
  nchannel=channels.size();
    
  cout<<"nchannel="<<nchannel<<endl;
      
  //    for (int channel=0; channel<nchannel; ++channel){
  for(int ch=0; ch < nchannel; ch++){

    if(debugADL)cout << "\r Setup channel " << channels[ch]  << endl;
    ADLDetector->SetSetupFile(channels[ch]);
    if(debugADL) std::cout << "Set Setup file:"<<channels[ch]<<endl;	  
    ADLDetector->ConfigureADL(ADLDetector->GetSetupFile(),resetPos); // ADL-4.2
  }
    
  // Open output ROOT file
  fOutputFile = new TFile(outputrootfilename.c_str(),"RECREATE");
  if(debugADL)cout<<"create output file "<<endl;  
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

TFile* ADLOutput::RunSimulation(std::string configfile,std::string InputFilename,std::string Nentries){
  
  std::cout << " " << std::endl;

  std::string treename = "fTree";
  TTree* Tree = 0;
  int traces = 0;
  
  std::cout<< " " <<std::endl;

  TFile* fInputFile = new TFile(InputFilename.c_str(),"READ");

  if(fInputFile->GetListOfKeys()->Contains(treename.c_str()))
    {
      fInputFile->GetObject(treename.c_str(),Tree);
      
      if(Tree->GetNbranches() > 0){
	//	std::cout << "Number of branches in " << treename << " tree is " << Tree->GetNbranches() << std::endl;
	
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

	int Nevents;
	if(Nentries == "0") Nevents = Tree->GetEntries();
	else  Nevents =atoi(Nentries.c_str());// (Tree->GetEntries())/500;

	for(int i =0;i<Nevents;i++){
	  if(i % 1000 == 0){
	    end = clock();
	    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	    cout << i << " events / " << Nevents << " treated in " << cpu_time_used << " seconds" << endl;
	  }
	  Tree->GetEntry(i);
          double energy[nchannel],SimuEnergy[nchannel];

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
	  //Simulate pulses with ADL once per event
     
	  if(!hits_totnum){
	    SimulateEmptyWaveform(0);
	    traces++;
	  }
	  else{
	    for(int j=0;j<nchannel;j++){
            traceCalculated[channels[j]] = 0;
	    energy[j]=0.0;
            SimuEnergy[j]=0.0;  
	    }
	    for(int j=0;j<hits_totnum; j++) {
	      if( traceCalculated[hits_iddet[j]] == 0){
                energy[hits_iddet[j]]+=hits_edep[j];

		if(hits_tote > 0.5)SimuEnergy[hits_iddet[j]]= SimulatePulse(hits_iddet[j]);
		else               SimulateEmptyWaveform(hits_iddet[j]);
		
		traceCalculated[hits_iddet[j]] = 1;
		traces++;
		
	      }
	    }   
	  }
	    fOutputFile->cd();
	  MGTree->Fill();
	  if(debugADL) std::cout << "Print MGTree : " << MGTree->GetEntries() << " events recorded" << std::endl;
	  
      
	}
	std::cout << " " << std::endl;
	std::cout << traces << " traces processed in total" << std::endl;
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

//void ADLOutput::SimulatePulse(int channel){
double ADLOutput::SimulatePulse(int channel){
  
  double edepThrs = 0.9999; // Threshold requiring a certain fraction of energy deposition in clusters as compared to individual hits
  double edepFlag = 0;     // Energy deposition ratio between detectors and clusters
  double ETotDet = 0;
  
  if(debugADL) std::cout << "DEBUG: Simulate pulse in channel " << channel << std::endl;

  ADLDetector->SetSetupFile(channel);
  if(debugADL) std::cout << "DEBUG: Setup potentials" << std::endl;
  ADLDetector->SetPotentials(ADLDetector->GetSetupFile()); // ADL-4.2
  if(debugADL) std::cout << "DEBUG: Potentials set" << std::endl;
  
  ADLDetector->SetPositionOffset(channel);
  
  
  ADLDetector->CreateADLevent();
  if(debugADL) std::cout << "DEBUG: ADL event created" << std::endl;

  if(fClusterization){
    //////////////////////////////////////////
    //
    //  Make cluster out of hits to save CPU time
    //
    //////////////////////////////////////////

    if(debugADL) std::cout << "DEBUG: ADL clusterization of " << hits_totnum << " hits" << std::endl;
    //for(int i = 0;i<hits_totnum;i++)  printf("    Hits position : %d %.03f %.03f %.03f %.03f \n",hits_iddet[i], hits_xpos[i],hits_ypos[i],hits_zpos[i],hits_edep[i]);

    ADLCluster HitsCluster(channel);
   
    edepFlag = HitsCluster.LaunchClustering(hits_totnum,hits_xpos,hits_ypos,hits_zpos,hits_edep,hits_iddet);
   
    //if(channel == 34) std::cout << " " << edepFlag << " " << hits_tote;

    if(edepFlag>=edepThrs && edepFlag <= 1.00001){
      if(debugADL) printf("DEBUG: Good clustering : %.05f  \n", edepFlag);
      std::vector<std::vector<double> > clusters = HitsCluster.GetClusters();
     
      //for(int i = 0;i<clusters[0].size();i++) printf("   Clusters position : %.03f %.03f %.03f %.05f \n", edepFlag, clusters[0][i],clusters[1][i],clusters[3][i]);
      
      ETotDet = ADLDetector->SetADLhits(clusters[0].size(),clusters[3],clusters[0],clusters[1],clusters[2]);
    }
    else 
      if(debugADL) printf("   Bad clustering : %.05f  \n", edepFlag);
    ///////////////////////////////////////////////
  }

  if(edepFlag < edepThrs || edepFlag > 1.00001){
    if(fRecordADLOutPos)
      ETotDet = ADLDetector->SetADLhits(hits_totnum,hits_edep,hits_xpos,hits_ypos,hits_zpos,hits_iddet,hits_ADLpos,hits_isOut);
  
    else
      ETotDet = ADLDetector->SetADLhits(hits_totnum,hits_edep,hits_xpos,hits_ypos,hits_zpos,hits_iddet);
    //if(channel == 34) std::cout << " NO cluster energy " << ETotDet << std::endl ;
  }
  if(debugADL) std::cout << "DEBUG: ADL hits set" << std::endl;
  if(debugADL) std::cout << "DEBUG: Deposited energy in channel " << channel << " is " << ETotDet << "/" << hits_tote << " MeV" << std::endl;
  if(ADLDetector->CalculateTrace()) std::cerr<< "Failed to calculate trace" <<std::endl;
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
  if(channel >= 1000) digiData->SetID(0);
  else  digiData->SetID(channel);
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
  
  ADLDetector->SetWaveformTimeOffset(wfPreTrigger,auxwfPreTrigger);
  if(debugADL) std::cout << "DEBUG: Waveform time offset set" << std::endl;
  
  if(ADLDetector->SetADLWaveform(waveform)) if(debugADL) std::cerr<< "Failed to set waveform" <<std::endl;
  if(debugADL) std::cout << "DEBUG: Waveform set" << std::endl;
  
  if (event->GetAuxWaveformArrayStatus()){
    auxwaveform = new ((*(event->GetAuxWaveforms()))[validChannelCounter])
      MGTWaveform(NULL,0,AuxSamplingFrequency,0.0,MGWaveform::kADC,0);     
    auxwaveform->SetSamplingFrequency(AuxSamplingFrequency);
    auxwaveform->SetID(channel);
    auxwaveform->SetTOffset(0.);      
    auxwaveform->SetLength(auxwfLength);
      
    if(ADLDetector->SetADLauxWaveform(auxwaveform)) if(debugADL) std::cerr<< "Failed to set aux waveform" <<std::endl;
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
  /*
    MGTree->Fill();
    if(debugADL) std::cout << "Print MGTree : " << MGTree->GetEntries() << " events recorded" << std::endl;
  */
  ADLDetector->DeleteADLevent();
  return ETotDet;
}

//**************To create empty waveforms (so that Exec_Module_ini does not crash*************************************************

void ADLOutput::SimulateEmptyWaveform(int channel){

  MGWaveformTag::EWaveformTag fWaveformTag = MGWaveformTag::kNormal;
  double ETotDet = 0;

  //Fill digiData. 
 
  digiData = new ((*(event->GetDigitizerData()))[validChannelCounter]) GETGERDADigitizerData(); 
  digiData->SetClockFrequency(SamplingFrequency);
  
  digiData->SetPretrigger(fPreTrigger);
  digiData->SetTimeStamp(fTimestamp);
  digiData->SetDecimalTimeStamp(fDecimalTimestamp);
  digiData->SetIsInverted(IsPulseInverted);
  digiData->SetTriggerNumber(TriggerNumber);
  if(channel >= 1000) digiData->SetID(0);
  else  digiData->SetID(channel);
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
  for (size_t i=0;i<waveform->GetLength();i++) (*waveform)[i] = 0.;
  
  if (event->GetAuxWaveformArrayStatus()){
    auxwaveform = new ((*(event->GetAuxWaveforms()))[validChannelCounter])
      MGTWaveform(NULL,0,AuxSamplingFrequency,0.0,MGWaveform::kADC,0);     
    auxwaveform->SetSamplingFrequency(AuxSamplingFrequency);
    auxwaveform->SetID(channel);
    auxwaveform->SetTOffset(0.);      
    auxwaveform->SetLength(auxwfLength);
    for (size_t i=0;i<auxwaveform->GetLength();i++) (*auxwaveform)[i] = 0.;
  }
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
