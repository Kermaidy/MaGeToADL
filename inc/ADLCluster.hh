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

  std::vector<int> SortHitsId(int,int,int*);
  int CheckClusters(int,int,int,float*,float*,std::vector<int> &, double);
  void SetClusterEnergy(int,int,float*,float*,std::vector<int> &,float*,double);
  double CheckEdep(int, int, float*, std::vector<int> &);
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

