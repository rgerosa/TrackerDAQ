#include <algorithm>
#include "../BadTrackerStrip.h"

// check the overlap between two files giving a complate list of strip information
void checkOverlapWithOffline(string file1, string file2){

  vector<BadTrackerStrip> vec1;
  vector<BadTrackerStrip> vec2;

  ifstream infile1 (file1.c_str());
  ifstream infile2 (file2.c_str());

  BadTrackerStrip strip;
  cout<<"Reading the first file "<<endl;
  if(infile1.is_open()){
    while(!infile1.eof()){      
      infile1 >> strip.detid_ >> strip.lldCh_ >> strip.apvid_ >> strip.stripid_ ;
      vec1.push_back(strip);
    }
  }
  infile1.close();

  cout<<"Reading the second file "<<endl;
  if(infile2.is_open()){
    while(!infile2.eof()){
      infile2 >> strip.detid_ >> strip.lldCh_ >> strip.apvid_ >> strip.stripid_ ;
      vec2.push_back(strip);
    }
  }
  infile2.close();

  long int total_vec1 = vec1.size();
  long int total_vec2 = vec2.size();
  long int common_strip = 0;

  cout<<"Find overlaps "<<endl;
  for(auto element : vec1){    
    auto element2 = std::find(vec2.begin(),vec2.end(),element);
    if(element2 != vec2.end())
      common_strip++;
  }
  
  cout<<"Total number of strips from file 1 "<<total_vec1<<endl;
  cout<<"Total number of strips from file 2 "<<total_vec2<<endl;
  cout<<"Common tagged strips "<<common_strip<<endl;
  cout<<"Fraction wrt file 1 "<<double(common_strip)/total_vec1<<endl;
  cout<<"Fraction wrt file 2 "<<double(common_strip)/total_vec2<<endl;

}

