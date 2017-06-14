#ifndef _ADLCLUSTER_HH
#define _ADLCLUSTER_HH

#include <vector>

//---------------------------------------------------------------------------//
//   Clustering algorithm from https://rosettacode.org/wiki/Category:C       //
//---------------------------------------------------------------------------//


class ADLCluster
{
public:
  //default constructor
  ADLCluster(int);
  
  //copy constructor
  ADLCluster(const ADLCluster &);
  
  //destructor
  ~ADLCluster();

  void SortHitsId(int,int,int*);
  int CheckClusters(int,int,int,float*,float*, double);
  void SetClusterEnergy(int,int,float*,float*,float*,double);
  double CheckEdep(int, int, float*);
  double LaunchClustering(int, float*, float*,float*,float*,int*);
  int GetDetHits(int,int*);
  void GetRadCoord(int,float*,float*);
  std::vector<std::vector<double> > GetClusters();
  

  //private  members
private:

    static const int MAX_NHITS=1000;    
    
    float hr[MAX_NHITS];
    
    std::vector<std::vector<double> > clustersPos;
    
    int detId;
};

#endif

