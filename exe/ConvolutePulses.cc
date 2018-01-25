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
#include "MGTypes.hh"

// GELATIO includes
#include "GETGERDADigitizerData.hh"

using namespace std;
using namespace MGWaveformTag;

int debugADL = 0;
int display = 0;

int NDET = 40;
double AuxSamplingFrequency = 0.100;  //GHz
double SamplingFrequency = 0.025;     //GHz

TTree* MGTree;

void GetFitParameters(int pole, std::vector<std::vector<double> > &vec1,std::vector<std::vector<double> > &vec2)
{
  std::string channel, a_decay, tau1_decay, tau2_decay, tau1_RC, tau2_RC, header;
  std::string filename;
  if(pole == 1) filename = "./FitElecRespOutput_1pole.txt";
  else filename = "./FitElecRespOutput_2poles.txt";

  std::cout << "\r Get BEGe detectors optimized E.R. parameters " << pole << std::endl;

  std::ifstream File(filename.c_str());

  if(File)  // Check if the file exists
    {
      getline(File,header);
      while(File.good())  // Loop until reach the eof
        {
          File >> channel  >> a_decay >> tau1_decay >> tau2_decay >> tau1_RC >> tau2_RC ;
          vec1[0].push_back(atof(a_decay.c_str()));   
          vec1[1].push_back(atof(tau1_decay.c_str()));   
          vec1[2].push_back(atof(tau2_decay.c_str()));   
	  vec2[0].push_back(atof(tau1_RC.c_str()));   
          vec2[1].push_back(atof(tau2_RC.c_str()));   
         }
    }
  else{  // if the file doesn't exist
    cerr << filename << " not found. Try again" << std::endl;
    exit(1);
  }

  File.close();  // Close file

  vec1[0].pop_back();
  vec1[1].pop_back();
  vec1[2].pop_back();
  vec2[0].pop_back();
  vec2[1].pop_back();
}

std::vector<double> GetDetectorPulseShape(std::string flag)
{
  std::vector<double> vec;
  std::string StringTmp, data, header;
  std::string filename = "./GERDADetectorGain.txt";

  std::cout << "\r Get BEGe detector " << flag << std::endl;

  std::ifstream File(filename.c_str());

  if(File)  // Check if the file exists
    {
      getline(File,header);
      while(File.good())  // Loop until reach the eof
        {
          if(flag == "baseline") File >> StringTmp  >> data >> StringTmp >> StringTmp;
          else if(flag == "gain") File >> StringTmp  >> StringTmp >> StringTmp >> data;
          vec.push_back(atof(data.c_str()));   // Store path & file name in a std::vector of std::strings
        }
    }
  else{  // if the file doesn't exist
    cerr << filename << " not found. Try again" << std::endl;
    exit(1);
  }

  File.close();  // Close file

  vec.pop_back();

  return vec;

}

std::vector<std::vector<std::vector<double> > > GetNoise(int isAux, std::string idsimu)
{
  std::vector<std::vector<std::vector<double> > > Data;
  std::vector<std::vector<double> > VectorTmp2;
  std::vector<double>* VectorTmp = 0;
  std::vector<double> VectorTmp1;

  //Rescale channels of the noise library (include Nat. Ge detectors)
  for(int channel = 0;channel<3;channel++){
    VectorTmp1.push_back(0);
    VectorTmp2.push_back(VectorTmp1);	  
    VectorTmp1.clear();
    Data.push_back(VectorTmp2);
    VectorTmp2.clear();
  }

  string treename = "noiseTree";

  //  for(int filenumber = 0;filenumber<1;filenumber++){
  ostringstream oss;
    //    oss << filenumber;
    
    string noisefilename = "/lfs/l3/gerda/akirsch/gerda-BEGesimulation/run-58-59-60-61-62-63-64_Results/output/noise_library/noiseLibrary_all_";
    noisefilename += idsimu + ".root";
    
    TFile* fNoiseFile = new TFile(noisefilename.c_str(),"READ");
    TTree* NoiseTree = 0;
    TBranch* BranchTmp = 0;
    
    if(fNoiseFile->GetListOfKeys()->Contains(treename.c_str()))
      {
	for(int channel = 0;channel<37;channel++){
	  oss.str("");
	  if(channel < 10) oss << 0;
	  oss << channel;
	  	  
	  int k = 0;
	  string branchname = "ch" + oss.str();
	  if(isAux) branchname += "_aux";
	  
	  std::cout << "\r Get noise from channel " << channel << " and file " << idsimu << " Branch : " << branchname << std::flush;
	  
	  if(debugADL)  std::cout << "Noise library " << noisefilename << " opened" << std::endl;
	  fNoiseFile->GetObject(treename.c_str(),NoiseTree);
	  if(NoiseTree->GetListOfBranches()->Contains(branchname.c_str())){
	   	    
	    NoiseTree->SetBranchAddress(branchname.c_str(),&VectorTmp,&BranchTmp);

	    for(int line = 0;line<NoiseTree->GetEntries();line++){
	    //for(int line = 0;line<1;line++){
	      NoiseTree->GetEntry(line);
	      VectorTmp2.push_back((*VectorTmp));	  	      
	      VectorTmp->clear();
	    }
	    
	  }
	  else{
	    VectorTmp->push_back(0);
	    VectorTmp2.push_back((*VectorTmp));	  
	    VectorTmp->clear();
	    k++;  
	    if(debugADL) std::cout << "\r Branch " << branchname << " not found in  " << treename << " tree" << std::endl;
	  }
	  Data.push_back(VectorTmp2);
	  VectorTmp2.clear();
	}
      }
    else std::cout << "\r Tree " << treename << " not found in  " << noisefilename << std::endl;
  fNoiseFile->Close();
  //  }

  
  return Data;
}


MGTWaveform* SetNoise(MGTWaveform* waveform,std::vector<std::vector<std::vector<double> > > &Noise, double Amplitude, double Baseline, double Scale, int channel) {
    
  int wfSize = waveform->GetLength();
  
  MGTWaveform *wf = new MGTWaveform(NULL,0,SamplingFrequency,0.0,MGWaveform::kADC,0);
  wf->SetLength(wfSize);
  
  gRandom->SetSeed(time(NULL));
  int randomline = int(gRandom->Uniform(0,Noise[channel].size()));
  
  if(Noise[channel][randomline].size() != wfSize){std::cerr << "Channel : " << channel << " Noise entry " << randomline << " doesn't have the proper size " << Noise.size() << "x" << Noise[channel].size() << "x" << Noise[channel][randomline].size() << " " << wfSize << std::endl; return wf;}
  
  for (size_t i=0;i<wfSize;i++) (*wf)[i] = Amplitude * Scale *  waveform->At(i) + 4.*Baseline + Noise[channel][randomline][i];
  
  return wf;
}

MGTWaveform* SetAuxNoise(MGTWaveform* waveform,std::vector<std::vector<std::vector<double> > > &AuxNoise, double AuxAmplitude, double AuxBaseline,double AuxScale,int channel) {
    
  int wfSize = waveform->GetLength();

  MGTWaveform *wf = new MGTWaveform(NULL,0,AuxSamplingFrequency,0.0,MGWaveform::kADC,0);
  wf->SetLength(wfSize);
  
  gRandom->SetSeed(time(NULL));
  int randomline = int(gRandom->Uniform(0,AuxNoise[channel].size()));
  
  if(AuxNoise[channel][randomline].size() != wfSize){std::cerr << "Channel : " << channel << " AuxNoise entry " << randomline << " doesn't have the proper size " << AuxNoise.size() << "x" << AuxNoise[channel].size() << "x" << AuxNoise[channel][randomline].size() << " " << wfSize << std::endl; return wf;}
  
  for (size_t i=0;i<wfSize;i++){
    //    std::cout<< "DEBUG :  Trace amplitudes : " <<  waveform->At(i) << " " << AuxBaseline << " " << AuxNoise[detector_channel][randomline][i] << std::endl;
    (*wf)[i] = AuxAmplitude * AuxScale * waveform->At(i) + AuxBaseline + AuxNoise[channel][randomline][i];
  }  

  return wf;
}

std::vector<std::vector<double> > GetTestPulser(int isAux)
{
  int debugPULSER = 0;
  std::vector<std::vector<double> > VectorTmp2;
  std::vector<std::vector<double> >* VectorTmp = 0;

  TTree* PulserTree = 0;
  TBranch* BranchTmp = 0;
  
  std::string PulseFilename = "/lfs/l3/gerda/kermaidy/Analysis/software/src/MaGeToADL/PulserLibrary.root";
  std::string treename = "Pulser";

  TFile* RootFile = new TFile(PulseFilename.c_str(),"READ");
  
  if(debugPULSER) std::cout << "\r Open pulser data file " << PulseFilename << std::endl;
  
  if(RootFile->GetListOfKeys()->Contains(treename.c_str()))
    {	
      int k = 0;
      string branchname = "LF";
      if(isAux) branchname = "HF";
      
      if(debugPULSER)  std::cout << "Pulse library " << PulseFilename << " opened" << std::endl;
      RootFile->GetObject(treename.c_str(),PulserTree);
      if(PulserTree->GetListOfBranches()->Contains(branchname.c_str())){
	
	PulserTree->SetBranchAddress(branchname.c_str(),&VectorTmp,&BranchTmp);
	
	PulserTree->GetEntry();
	VectorTmp2 = (*VectorTmp);	  	      
	VectorTmp->clear();
      }
    }
  else std::cout << "\r Tree " << treename << " not found in  " << PulseFilename << std::endl;
  RootFile->Close();
  
  if(debugPULSER){
  //if(channel == 4 && isAux){    
    TGraph* gr = new TGraph();
    for(int i = 0;i<VectorTmp2[0].size();i++) gr->SetPoint(i,i,VectorTmp2[0][i]);
    TApplication* myApp = new TApplication("mAapp",0,0);
    TCanvas* c = new TCanvas("c","Pulser sample",700,500);
    c->cd();
    gr->Draw();

    myApp->Run();
  }

  if(debugPULSER) std::cout << "Return pulser " << std::endl;
  return VectorTmp2;
  
 }
 
 MGTWaveform* ConvoluteWF(MGTWaveform*kAnInput,std::vector<double> &fFilter,double &Scale,int isAux)
{
  int debugCONV = 0;

  if(debugCONV) std::cout << "Enter ConvoluteWF" << std::endl;
  
  int shift;  // To get convolved waveform rising edge approximately at the center
  MGTWaveform *kAnOutput;
  if(!isAux) {kAnOutput = new MGTWaveform(NULL,0,SamplingFrequency,0.0,MGWaveform::kADC,0); shift = 800;}
  else       {kAnOutput = new MGTWaveform(NULL,0,AuxSamplingFrequency,0.0,MGWaveform::kADC,0); shift = 199;}

  int fWaveformSize = kAnInput->GetLength();
  int fFilterSize = fFilter.size();
  int fOutputSize = fWaveformSize;

  if(debugCONV) std::cout << "Set output WF size" << std::endl;
  kAnOutput->SetLength(fOutputSize);
  if(debugCONV) std::cout << "Input WF size  : " << fWaveformSize << std::endl;
  if(debugCONV) std::cout << "Output WF size : " << fOutputSize  << std::endl;

  if(debugCONV) std::cout << "Start filtering" << std::endl;
  int k;
  int tMax = 0;
  double AmpMax = -1;
  
  if( fFilterSize > 0){
    for( int i = 0; i < 0.80*fOutputSize; i ++ ){
      (*kAnOutput)[i] = 0.;
      for( int j = 0; j <= fFilterSize/2; j ++ ){
	if(j<=i) (*kAnOutput)[i] += (kAnInput->At(i+shift-j))*fFilter[j+fFilterSize/4];
      }
      if(kAnOutput->At(i) > AmpMax) {AmpMax = kAnOutput->At(i); tMax = i;}
    }
    for( size_t i = tMax; i < fWaveformSize; i++ ) (*kAnOutput)[i] = kAnOutput->At(tMax);
  }

  Scale = (kAnInput->At(fWaveformSize-1)-kAnInput->At(1))/kAnOutput->At(tMax);

  if(debugCONV && isAux){
    TApplication* myApp = new TApplication("mAapp",0,0);
    TGraph* grFil = new TGraph();  for(int i = 0;i<fFilter.size();i++)         grFil->SetPoint(i,i,fFilter[i]);
    TGraph* grIn  = new TGraph();  for(int i = 0;i<kAnInput->GetLength();i++)  grIn->SetPoint(i,i,kAnInput->At(i));
    TGraph* grOut = new TGraph();  for(int i = 0;i<kAnOutput->GetLength();i++) grOut->SetPoint(i,i,kAnOutput->At(i));
    
    std::cout << kAnInput->GetLength() << " " << kAnOutput->GetLength() << std::endl;
    
    TCanvas* cFil = new TCanvas("cFil","Filter response",700,500);
    cFil->cd();
    grFil->Draw("ALP");

    TCanvas* cIn = new TCanvas("cIn","Input pulse",700,500);
    cIn->cd();
    grIn->Draw("ALP");

    TCanvas* cOut = new TCanvas("cOut","Output pulse",700,500);
    cOut->cd();
    grOut->Draw("ALP");
    
    myApp->Run();
  }
  if(debugCONV) std::cout << "End filtering" << std::endl;

  return kAnOutput;
}

MGTWaveform* ApplyRC(MGTWaveform*kAnInput, double RCtau,int isAux)
{
  MGTWaveform *kAnOutput;
  double timeStep = 1./(AuxSamplingFrequency*1e3);
  if(!isAux) kAnOutput = new MGTWaveform(NULL,0,SamplingFrequency,0.0,MGWaveform::kADC,0);    
  else       kAnOutput = new MGTWaveform(NULL,0,AuxSamplingFrequency,0.0,MGWaveform::kADC,0);
  
  int i,j;
  int fWaveformSize = kAnInput->GetLength();
  
  kAnOutput->SetLength(fWaveformSize);
  
  //Apply a simple RC circuit preamp response if time constant is not zero (Y. Kermaidic May 2017)
  if(RCtau>0){
    double integral;    
    for(j = 0; j<fWaveformSize;j++) (*kAnOutput)[j] = kAnInput->At(j);
    
    for(j = 1; j<fWaveformSize;j++){
      integral = kAnOutput->At(j-1) + (kAnInput->At(j-1) - kAnOutput->At(j-1))/RCtau * timeStep;
      (*kAnOutput)[j] = integral;
    }
  }
  else {std::cerr << "RC time constant not defined. Exit." << std::endl; exit(1);}
  
  return kAnOutput;
}

MGTWaveform* Apply2poles(MGTWaveform*kAnInput, double RCtau1, double RCtau2,int isAux)
{
  MGTWaveform *kAnOutput;
  double timeStep = 10;//1./(AuxSamplingFrequency*1e3);
  if(!isAux) kAnOutput = new MGTWaveform(NULL,0,SamplingFrequency,0.0,MGWaveform::kADC,0);    
  else       kAnOutput = new MGTWaveform(NULL,0,AuxSamplingFrequency,0.0,MGWaveform::kADC,0);
  
    kAnOutput->SetLength(kAnInput->GetLength());

    double a1 = timeStep/RCtau2;
    double b1 = 1-timeStep/(RCtau1+RCtau2);
    double b2 = -pow(timeStep,2)/(RCtau2*RCtau1);
    
    double nextpoint;
    
    for(int j = 0; j<kAnInput->GetLength();j++)
        (*kAnOutput)[j] = kAnInput->At(j);
    
    for(int j = 2; j<kAnInput->GetLength();j++){
        nextpoint = (a1*kAnInput->At(j-1) + b1*kAnOutput->At(j-1) + b2*kAnOutput->At(j-2)) / (a1+b1+b2);
	(*kAnOutput)[j] = nextpoint;
//        if(nextpoint < 5) (*kAnOutput)[j] = nextpoint;
//        else (*kAnOutput)[j] = 0;
    }
    
    return kAnOutput;
}

MGTWaveform* ApplyRCdecay(MGTWaveform*kAnInput,double A,double tau1, double tau2,int isAux)
{

    std::vector<double> kAnOutput1(kAnInput->GetLength(),0);
    std::vector<double> kAnOutput2(kAnInput->GetLength(),0);

    MGTWaveform *kAnOutput;
    double timeStep = 1./(AuxSamplingFrequency*1e3);
    if(!isAux) kAnOutput = new MGTWaveform(NULL,0,SamplingFrequency,0.0,MGWaveform::kADC,0);
    else       kAnOutput = new MGTWaveform(NULL,0,AuxSamplingFrequency,0.0,MGWaveform::kADC,0);
    
    kAnOutput->SetLength(kAnInput->GetLength());

    //1st decay constant parameters
    double x1 = exp(-0.01/tau1);
    double a10 = (1+x1)/2.;
    double a11 = -(1+x1)/2.;
    double b11 = x1;
    
    //2nd decay constant parameters
    double x2 = exp(-0.01/tau2);
    double a20 = (1+x2)/2.;
    double a21 = -(1+x2)/2.;
    double b21 = x2;
    
    for(int j = 1; j<kAnInput->GetLength();j++){
        kAnOutput1[j] = a10*kAnInput->At(j) + a11*kAnInput->At(j-1) + b11*kAnOutput1[j-1];
        kAnOutput2[j] = a20*kAnInput->At(j) + a21*kAnInput->At(j-1) + b21*kAnOutput2[j-1];
        (*kAnOutput)[j] = A*kAnOutput1[j] + (1-A)*kAnOutput2[j];
    }


    return kAnOutput;

//    return kAnInput;
}

int WriteOutputs(TFile* fOutputFile)
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

int main(int argc, const char* argv[])
{
  std::cout << " " << std::endl;
  std::cout << "Enter in ConvolutePulses " << std::endl;
  std::cout << " " << std::endl;

  if(argc < 6){
    std::cerr << "Not enough arguments given "         << std::endl;
    std::cerr << "Try ./ConvolutePulses "              << std::endl;
    std::cerr << " $INFILENAME "                       << std::endl;
    std::cerr << " $OUTFILENAME "                      << std::endl;
    std::cerr << " $IDSIMU "                           << std::endl;
    std::cerr << " $ISNOISE[0 or 1=Gerda or 2=Gaus] "  << std::endl;
    std::cerr << " $RMSNOISE[if ISNOISE=2, set AuxWf rms noise and x4 for Wf] " << std::endl;
    std::cerr << " $GAIN[-1 for the GERDA array] "     << std::endl;
    std::cerr << " $BASELINE[-1 for the GERDA array] " << std::endl;
    std::cerr << " $ISFILTER[0 or 1(pulser) or 2(1-pole circuit) or 3(2-poles circuit)] " << std::endl;
    std::cerr << " $RCtimeCnst[mus] (if manual E.R.)"  << std::endl;
    std::cerr << " $RCdecaytimeCnst[mus] (if manual E.R.)" << std::endl;
    exit(0);
  }

  std::string wfInRootFilename = argv[1];
  std::string wfOutRootFilename = argv[2];
  std::string idSimu = argv[3];
  int isNoise = atoi(argv[4]);
  int RMSnoise = atof(argv[5]);
  double gain = atof(argv[6]);
  double baseline = atof(argv[7]);
  int filter = atoi(argv[8]);

  double scale = 1., auxscale = 1.;

  // Set waveform amplitudes/noise
  std::cout << "Get detector pulse informations (baseline/gain) " << std::endl;
  std::vector<double> Baseline;
  if(baseline == -1) Baseline = GetDetectorPulseShape("baseline");
  else Baseline.push_back(baseline);
  std::vector<double> Gain;
  if(gain == -1) Gain = GetDetectorPulseShape("gain");
  else Gain.push_back(gain);

  std::vector<std::vector<std::vector<double> > > Noise, AuxNoise;

  if(isNoise == 1){
    std::cout << " Get noise data " << std::endl;
    Noise = GetNoise(0,idSimu);
    AuxNoise = GetNoise(1,idSimu);
  }
  else if(isNoise == 2){
    std::cout << " Set gaussian noise" << std::endl;
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0,4.*RMSnoise);
    
    std::vector<std::vector<double> > noise2; 
    std::vector<double> noise1;
    
    int NnoiseSamples= 1000;

    for(int i = 0; i < NnoiseSamples;i++){
      for(int j = 0; j < 4096;j++)
	noise1.push_back(distribution(generator));
      noise2.push_back(noise1);
      noise1.clear();
    }

    Noise.assign(NDET,noise2);
    noise2.clear();
    
    std::normal_distribution<double> auxdistribution(0.0,RMSnoise);
    
    for(int i = 0; i < NnoiseSamples;i++){
      for(int j = 0; j < 1000;j++)
	noise1.push_back(distribution(generator));
      noise2.push_back(noise1);
      noise1.clear();
    }
    
    AuxNoise.assign(NDET,noise2);
    noise2.clear();
  }
  else{
    std::cout << " Set noise data to 0 " << std::endl;
    std::vector<double> noise1; noise1.assign(4096,0);
    std::vector<std::vector<double> > noise2; noise2.assign(1,noise1);
    Noise.assign(NDET,noise2);
    std::vector<double> auxnoise1; auxnoise1.assign(1000,0);
    std::vector<std::vector<double> > auxnoise2; auxnoise2.assign(1,auxnoise1);
    AuxNoise.assign(NDET,auxnoise2);
  }
  
  std::vector<std::vector<double> > Pulser, AuxPulser;

  if(filter == 1){
    std::cout << "\r Get pulser data " << std::endl;
    Pulser = GetTestPulser(0);
    AuxPulser = GetTestPulser(1);
  }

  std::vector<std::vector<double> > DecayCnst(3,std::vector<double>(0)); // 1-amplitude / 2-time constants
  std::vector<std::vector<double> > ERCnst(2,std::vector<double>(0));    // 2-time constants

  // Recover optimized BEGe detector pulses parameters from E.R. fit routine output
  if(filter == 2){
	if(argc == 11) for(int i = 0;i<NDET;i++){ ERCnst[0].push_back(atof(argv[9])); DecayCnst[0].push_back(1.); DecayCnst[0].push_back(atof(argv[10])); DecayCnst[0].push_back(0);}
	else if(argc == 10) for(int i = 0;i<NDET;i++) ERCnst[0].push_back(atof(argv[9]));
  	else if(argc == 9) GetFitParameters(1, DecayCnst, ERCnst);
  }
  else GetFitParameters(2, DecayCnst, ERCnst);

  // Parameter of the ORTEC detector
  //    - RC decay: 0.982818,50.0869,0.451906
  //    - ER      : 2.06571,140

  if(display) std::cout << "Start initializing output " << std::endl;
  // Output ROOT file
  TFile* fOutputFile = new TFile(wfOutRootFilename.c_str(),"RECREATE");
     
  // Output Tier1 tree definition
  MGTree = new TTree("MGTree","MGTree");
  
  MGTEvent* eventOut = 0;
  MGTWaveform* waveform;
  MGTWaveform* auxwaveform;
  GETGERDADigitizerData* digiData;

  // Output Waveform branch
  MGTree->Branch("event",&eventOut, 32000,-99);

  // Run info
  MGRunType::RunType fRunType;
  fRunType = MGRunType::kData;
  std::string theDAQ = "Struck";
  int fRunNumber = 1;
  time_t fStartTime = 1479394213;
  time_t fStopTime = 1479394213.000160000;
  std::string theRunDescription = "No description";
  
  MGTRun* theRun = new MGTRun();
  theRun->SetRunType(fRunType);
  theRun->SetRunNumber(fRunNumber);
  theRun->SetRunDescription(theRunDescription);
  theRun->SetStartTime(fStartTime);
  theRun->SetStopTime(fStopTime);
  //The Stop time is available at the end of run
  theRun->SetParentDAQLabel(theDAQ);
  
  // Set MGDO digitizer
  size_t fPreTrigger = 100;
  unsigned long long fTimestamp = 1479394213;
  unsigned long long fDecimalTimestamp = 0;
  bool IsPulseInverted = false;
  unsigned int TriggerNumber = 1;
  int fEventNumber = 1;
  bool fMuVetoed = false;
  unsigned int fMuVetoSample = 0;
  int wfLength = 4096;
  int auxwfLength = 1000;
  
  Double_t theTime = 1479394213;
  size_t iChannel = 0;
  
  if(display) std::cout << "End initializing output " << std::endl;

  std::cout << "Start recovering raw pulses " << std::endl;
  // Input ROOT file
  TFile* wfRootFile = new TFile(wfInRootFilename.c_str(),"READ");
  
  if(!wfRootFile->IsOpen()){std::cerr << wfInRootFilename << " not opened. Exit." << std::endl; exit(1);}

  MGTEvent* eventIn = 0;

  if(!wfRootFile->GetListOfKeys()->Contains("MGTree")){std::cerr << wfInRootFilename << " does not contain MGTree. Exit." << std::endl; exit(1);}

  TTree* tree=(TTree*)wfRootFile->Get("MGTree");
  tree->SetBranchAddress("event", &eventIn);
  tree->GetEntry();
  
  int Nentries = tree->GetEntries();
  clock_t start = clock(), end;
  double cpu_time_used;

  if(display) std::cout << "Tree entries : " << Nentries << std::endl;
  
  for(int i=1;i<Nentries;i++){
    if(i % 1000 == 0){
      end = clock();
      cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
      cout << i << " events / " << Nentries << " treated in " << cpu_time_used << " seconds" << endl;
    }
    
    tree->GetEntry(i);
    
    if(display) std::cout << "Event number : " << eventIn->GetEventNumber() << std::endl;
    if(display) std::cout << "Waveform number : " << eventIn->GetNWaveforms() << std::endl;

    // Set MGDO event structure
    eventOut = new MGTEvent();
    MGEventType::EventType fEventType;
    fEventType = MGEventType::kBaseline;
    eventOut->SetAuxWaveformArrayStatus(true);
    eventOut->SetTotalNumberOfIDs(NDET);
    eventOut->SetAllActiveIDs(true);
    eventOut->InitializeArrays("GETGERDADigitizerData",1);
    eventOut->SetBypassStreamerArray(kFALSE); //! Event type info
    eventOut->SetEventType(fEventType);
    eventOut->SetTime(theTime);
    eventOut->SetWFEncScheme(MGTWaveform::kDiffVarInt);
    eventOut->SetAuxWFEncScheme(MGTWaveform::kDiffVarInt);
    eventOut->SetETotal(eventIn->GetETotal());

    //    MGWaveformTag::EWaveformTag fWaveformTag = MGWaveformTag::kNormal;

    if(debugADL) std::cout << "Event : " << i << std::endl;

    int channel = -1;
    for(int validChannelCounter = 0;validChannelCounter < eventIn->GetNWaveforms();validChannelCounter++){
      //Fill digiData. 
      for(int testchannel = channel+1;testchannel<NDET;testchannel++){
	if(debugADL) std::cout << "Counter : " << validChannelCounter << " Mwf : " << eventIn->GetNWaveforms() << " Find ID index : " << testchannel << " : " << eventIn->FindIDIndex(testchannel)  << std::endl;
	if(eventIn->FindIDIndex(testchannel) >= 0){
	  channel = testchannel;
	  break;
	}
      }
      if(channel == -1){ std::cerr << "\r No valid channel found for event " << i << std::endl; exit(1);}
      
      //if(channel != 20) continue;
      
      if(debugADL) std::cout << i << " Counter : " << validChannelCounter << " Channel : " << channel << " ETot : " << eventIn->GetETotal() << std::endl;

      GETGERDADigitizerData* digiData = new ((*(eventOut->GetDigitizerData()))[validChannelCounter]) GETGERDADigitizerData(); 
      digiData->SetClockFrequency(SamplingFrequency);      
      digiData->SetPretrigger(fPreTrigger);
      digiData->SetTimeStamp(fTimestamp);
      digiData->SetDecimalTimeStamp(fDecimalTimestamp);
      digiData->SetIsInverted(IsPulseInverted);
      digiData->SetTriggerNumber(TriggerNumber);
      digiData->SetID(channel);
      digiData->SetEnergy(eventIn->GetETotal());
      digiData->SetEventNumber(eventIn->GetEventNumber());
      digiData->SetIsMuVetoed(fMuVetoed);
      digiData->SetMuVetoSample(fMuVetoSample);
      digiData->SetWaveformTag(MGWaveformTag::kNormal);
      if(debugADL) std::cout << "DEBUG: ADL digitizer set" << std::endl;
      
      // Set MGDO waveform 
      waveform = new ((*(eventOut->GetWaveforms()))[validChannelCounter]) MGTWaveform(NULL,0,SamplingFrequency,0.0,MGWaveform::kADC,0);
      waveform->SetLength(wfLength);
      waveform->SetTOffset(0.);
      waveform->SetID(channel);
      
      if (eventOut->GetAuxWaveformArrayStatus()){
	auxwaveform = new ((*(eventOut->GetAuxWaveforms()))[validChannelCounter])
	  MGTWaveform(NULL,0,AuxSamplingFrequency,0.0,MGWaveform::kADC,0);     
	auxwaveform->SetID(channel);
	auxwaveform->SetTOffset(0.);      
	auxwaveform->SetLength(auxwfLength);
      }
     
      if(debugADL) std::cout << "Waveform length : " 
			     << eventIn->GetWaveformID(channel)->GetLength() << " & " 
			     << eventIn->GetAuxWaveformID(channel)->GetLength() << std::endl;
      
      waveform->operator=(*ApplyRCdecay(eventIn->GetWaveformID(channel),DecayCnst[0][channel],DecayCnst[1][channel],DecayCnst[2][channel],0));
      auxwaveform->operator=(*ApplyRCdecay(eventIn->GetAuxWaveformID(channel),DecayCnst[0][channel],DecayCnst[1][channel],DecayCnst[2][channel],1));

      if(filter == 1){
	if(debugADL) std::cout << "Convolute raw pulses with test pulser response " << std::endl;
	waveform->operator=(*ConvoluteWF(waveform,Pulser[channel],scale,0));
	auxwaveform->operator=(*ConvoluteWF(auxwaveform,AuxPulser[channel],auxscale,1));
      }
      else if(filter == 2){
	if(debugADL) std::cout << "Convolute raw pulses with RC circuit response " << channel << " " << ERCnst[0][channel]  << std::endl;
	waveform->operator=(*ApplyRC(waveform,ERCnst[0][channel],0));
	auxwaveform->operator=(*ApplyRC(auxwaveform,ERCnst[0][channel],1));
      }
      else if(filter == 3){
	if(debugADL) std::cout << "Convolute raw pulses with 2 poles circuit response " << channel << "  " << ERCnst[0][channel] << "  " << ERCnst[1][channel] << std::endl;
	waveform->operator=(*Apply2poles(waveform,ERCnst[0][channel],ERCnst[1][channel],0));
	auxwaveform->operator=(*Apply2poles(auxwaveform,ERCnst[0][channel],ERCnst[1][channel],1));
      }
      
      if(filter != 0){
	if(debugADL) std::cout << "Add noise/gain/baseline to pulses " << std::endl;
	waveform->operator=(*SetNoise(waveform,Noise,Gain[channel],Baseline[channel],scale,channel));
	auxwaveform->operator=(*SetAuxNoise(auxwaveform,AuxNoise,Gain[channel],Baseline[channel],auxscale,channel));
      }
      else{
	waveform->operator=(*SetNoise(eventIn->GetWaveformID(channel),Noise,Gain[channel],Baseline[channel],scale,channel));
	auxwaveform->operator=(*SetAuxNoise(eventIn->GetAuxWaveformID(channel),AuxNoise,Gain[channel],Baseline[channel],auxscale,channel));
      }

      if(display){
	TApplication* myApp = new TApplication("mAapp",0,0);
	TGraph* gr = new TGraph();     
	for(int i = 0;i<waveform->GetLength();i++)  gr->SetPoint(i,i*40,(*waveform)[i]);
	TGraph* grAux = new TGraph();  
	for(int i = 0;i<auxwaveform->GetLength();i++)grAux->SetPoint(i,i*10,(*auxwaveform)[i]);

	std::cout << waveform->GetLength() << " " << auxwaveform->GetLength() << std::endl;
	
	TCanvas* cAux = new TCanvas("cAux","aux Pulser sample",700,500);
	cAux->cd();
	grAux->Draw();
	TCanvas* c = new TCanvas("c","Pulser sample",700,500);
	c->cd();
	gr->Draw();
	
	myApp->Run();
      }
    }
    fOutputFile->cd();
    MGTree->Fill(); 
  }
  wfRootFile->Close();

  WriteOutputs(fOutputFile);

  return 0;
}
