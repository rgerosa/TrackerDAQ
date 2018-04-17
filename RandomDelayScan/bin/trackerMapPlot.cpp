using namespace std;
#include "CommonTools/TrackerMap/interface/TrackerMap.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


int main(int argc, char**argv){

  if(argc < 6){
    std::cerr<<"Error in parsing parameters --> 5 are required "<<std::endl;
    return -1;
  }

  string inputFile1 = string(argv[1]);
  int    isLogScale = atoi(argv[2]);
  float  min = atof(argv[3]);
  float  max = atof(argv[4]);
  string outputname = string(argv[5]);
  string zaxisname  = "Delay (ns)";

  //declare the tracker map object
  edm::ParameterSet pset;
  pset.addUntrackedParameter<bool>("logScale",isLogScale);
  TrackerMap themap(pset);

  TString title ("Delay per module");
  title.ReplaceAll("-"," ");
  themap.setTitle(string(title.Data()));
  
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
  
  themap.save(true,min,max,outputname+".png",1500,800);
  themap.save(true,min,max,outputname+".pdf",1500,800);
  themap.save(true,min,max,outputname+".root",1500,800);
}
