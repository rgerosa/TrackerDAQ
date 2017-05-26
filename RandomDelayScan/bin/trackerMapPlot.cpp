using namespace std;
#include "CommonTools/TrackerMap/interface/TrackerMap.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


int main(int argc, char**argv){

  if(argc < 7){
    std::cerr<<"Error in parsing parameters --> 6 are required "<<std::endl;
    return -1;
  }

  string inputFile1 = string(argv[1]);
  string title = string(argv[2]);
  int    isLogScale = atoi(argv[3]);
  float  min = atof(argv[4]);
  float  max = atof(argv[5]);
  string outputname = string(argv[6]);

  //declare the tracker map object
  edm::ParameterSet pset;
  pset.addUntrackedParameter<bool>("logScale",isLogScale);
  TrackerMap themap(pset);
  themap.setTitle(title);
  
  std::ifstream input1(inputFile1.c_str());
  map<uint32_t,float> map1;
  if(input1.is_open()){
    while(!input1.eof()){
      uint32_t detid;
      float    value;
      input1 >> detid >> value;
      map1[detid] = value;
    }
  }
  input1.close();

  cout<<map1.size()<<endl;

  for(auto imap: map1){
    themap.fill_current_val(imap.first,imap.second);              
  }
  
  themap.save(true,min,max,outputname+".png",1400,800);
  themap.save(true,min,max,outputname+".pdf",1400,800);
}
