#include <algorithm>

class BadChannel {
  
public:

  BadChannel(){};
  ~BadChannel(){};
  BadChannel(const uint32_t & detid, const float & fraction):
    detid_(detid),
    fraction_(fraction){
  }

  bool operator == (const BadChannel & a) const{
    if(a.detid_ == detid_) return true;
    else return false;
  }

  uint32_t detid_;
  float   fraction_;

};

// check the overlap between two files giving a complate list of strip information
void mergeStripBadList(string file1, string file2, string outputDIR, string outputFileName){

  vector<BadChannel> vec1;
  vector<BadChannel> vec2;

  ifstream infile1 (file1.c_str());
  ifstream infile2 (file2.c_str());

  BadChannel* module = new BadChannel();
  if(infile1.is_open()){
    while(!infile1.eof()){
      infile1 >> module->detid_ >> module->fraction_ ;
      vec1.push_back(*module);
    }
  }
  infile1.close();

  if(infile2.is_open()){
    while(!infile2.eof()){
      infile2 >> module->detid_ >> module->fraction_ ;
      vec2.push_back(*module);
    }
  }
  infile2.close();

  long int total_vec1 = vec1.size();
  long int total_vec2 = vec2.size();

  // merging in a new vector 
  vector<BadChannel> vec3;
  for(auto element : vec1){
    auto element2 = std::find(vec2.begin(),vec2.end(),element);
    if(element2 != vec2.end()){
      module->detid_ = element.detid_;
      module->fraction_ = element.fraction_+element2->fraction_;
      vec3.push_back(*module);
    }
    else{
      vec3.push_back(element);
    }
  }
  
  for(auto element : vec2){
    auto element2 = std::find(vec1.begin(),vec1.end(),element);
    if(element2 == vec1.end()){
      vec3.push_back(element);
    }
  }

  // write a new otuput file
  ofstream mergedFile ((outputDIR+"/"+outputFileName).c_str());
  for(auto element : vec3){
    mergedFile<< element.detid_ <<" "<<element.fraction_<<" \n";
  }
  mergedFile.close();

}

