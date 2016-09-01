using namespace std;
#include "CommonTools/TrackerMap/interface/TrackerMap.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


int main(int argc, char**argv){

  if(argc < 8){
    std::cerr<<"Error in parsing parameters --> 7 are required "<<std::endl;
    return -1;
  }

  string inputFile1 = string(argv[1]);
  string inputFile2 = string(argv[2]);
  string title = string(argv[3]);
  int    isLogScale = atoi(argv[4]);
  float  min = atof(argv[5]);
  float  max = atof(argv[6]);
  string outputname = string(argv[7]);

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
      if(value == 0) continue;
      map1[detid] = value;
    }
  }
  input1.close();

  map<uint32_t,float> map2;
  if(inputFile2 != "NULL" ){
    std::ifstream input2(inputFile2.c_str());
    if(input2.is_open()){
      while(!input2.eof()){
	uint32_t detid;
	float    value;
	input2 >> detid >> value;
	if(value == 0) continue;
	map2[detid] = value;
      }
      input2.close();
    }
  }
  
  if(map2.size() != 0){
    for(auto imap: map1)
      themap.fill_current_val(imap.first,imap.second-map2[imap.first]);          
  }
  else{
    for(auto imap: map1){
      themap.fill_current_val(imap.first,imap.second);              
    }
  }

  themap.save(true,min,max,outputname+".png",1400,800);
  themap.save(true,min,max,outputname+".pdf",1400,800);
}
