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
  std::vector<int> ReadConfigfile(std::string);
  void SetSetupFile(int);
  void ConfigureADL(std::string,int);
  void SetPotentials(std::string);
  void SetPositionOffset(int);
  void CreateADLevent();
  void DeleteADLevent();
  void SetWaveformTimeOffset(double,double);
  double SetADLhits(int, Float_t*, Float_t*, Float_t*, Float_t*, Int_t*, Float_t*, Int_t*);
  double SetADLhits(int, Float_t*, Float_t*, Float_t*, Float_t*, Int_t*);
  double SetADLhits(int, std::vector<double> &,std::vector<double> &,std::vector<double> &,std::vector<double> &);
//  double SetADLhits(int, Float_t&, Float_t&, Float_t&, Float_t&, Int_t&);
  int CalculateTrace();
  int SetADLWaveform(MGTWaveform*);
  int SetADLauxWaveform(MGTWaveform*);
  int GetTraceDim();
  int GetCenter();
  double** GetElectronPath();
  double** GetHolePath();
  std::vector<double> x0,y0,z0;
  std::vector<int> channels,inv;
private:


  static const int NDET=40;

  struct SIMION_PA *ADL_Epot[NDET];
  struct SIMION_PA *ADL_Wpot[NDET];
  struct SIMION_PA *ADL_Stru[NDET];

  double GridSize[NDET];
  double Center[NDET];
  double Height[NDET];

//  std::vector<double> x0,y0,z0;

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
  
  double wfPreTrigger;  // Time before trigger (usually half of the signal length)
  double AuxwfPreTrigger;
};

#endif
