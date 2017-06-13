#ifndef _ADLCLUSTER_HH
#define _ADLCLUSTER_HH

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
  
  double LaunchClustering(int, float*, float*,float*,float*,int*);
  void GetRadCoord(int,float*,float*);
  std::vector<std::vector<double> > GetClusters();
  

  //private  members
private:

    static const int MAX_NHITS=1000;    
    
    float hr[MAX_NHITS];
    
    std::vector<std::vector<double> > clustersPos(4);
    
    int detId;
};

#endif

