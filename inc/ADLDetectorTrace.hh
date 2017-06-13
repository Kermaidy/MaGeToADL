#ifndef _ADLDETECTORTRACE_HH
#define _ADLDETECTORTRACE_HH

#include <sstream>      // std::ostringstream

// MGDO includes
#include "MGTWaveform.hh"

// ADL includes
#include "ADL.h"

class ADLDetectorTrace
{
public:
  //default constructor
  ADLDetectorTrace(int,int);

  //copy constructor
  ADLDetectorTrace(const ADLDetectorTrace &);

  //destructor
  ~ADLDetectorTrace();

  // public functions
  std::string GetSetupFile();
  void SetSetupFile(int);
  void ConfigureADL(std::string);
  void SetPotentials(std::string);
  void SetPositionOffset(double);
  void CreateADLevent();
  void DeleteADLevent();
  void SetWaveformAttribute(double,double,double,double);
  void SetAuxWaveformAttribute(double,double,double,double);
  double SetADLhits(int, Float_t*, Float_t*, Float_t*, Float_t*, Int_t*, Float_t*, Int_t*);
  double SetADLhits(int, Float_t*, Float_t*, Float_t*, Float_t*, Int_t*);
  double SetADLhits(int, std::vector<double> &,std::vector<double> &,std::vector<double> &,std::vector<double> &);
//  double SetADLhits(int, Float_t&, Float_t&, Float_t&, Float_t&, Int_t&);
  int CalculateTrace();
  std::vector<std::vector<std::vector<double> > > GetNoise(int);
  int SetADLWaveform(MGTWaveform*);
  int SetADLauxWaveform(MGTWaveform*);
  int GetTraceDim();
  int GetCenter();
  double** GetElectronPath();
  double** GetHolePath();

private:

  struct SIMION_PA *ADL_Epot[40];
  struct SIMION_PA *ADL_Wpot[40];
  struct SIMION_PA *ADL_Stru[40];

  double GridSize[40];
  double Center[40];
  double Height[40];

  std::string detector_setupfile; // ADL configuration file
  int detector_channel; // Detector channel
  int debugADL;

  double inverted; // Correct ADL hits position for inverted detectors
  double xoffset;  // Transpose MaGe coordinate
  double yoffset;  //      into the ADL referential
  double zoffset;

  double xcenter;   // Detector center in cm in ADL coordinate 
  double ycenter;   // Detector center in cm in ADL coordinate 
  double height;   // Detector height in cm 
  double gridsize; // grid size used for det. potentials computation

  struct ADL_EVENT *ADL_evt; // Define ADL event to set hits position and calculate carriers path
  
  std::vector<std::vector<std::vector<double> > > Noise;    // Store WF noise library
  std::vector<std::vector<std::vector<double> > > AuxNoise; // Store aux WF noise library

  double Amplitude;     // Signal amplitude in ADC
  double Baseline;      // Signal baseline
  double RMS_noise;     // Baseline noise amplitude (assumed to be gaussian)
  double wfPreTrigger;  // Time before trigger (usually half of the signal length)

  double AuxAmplitude;
  double AuxBaseline;
  double AuxRMS_noise;
  double AuxwfPreTrigger;
};

#endif
