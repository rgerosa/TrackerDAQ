#include <algorithm>
#include "TrackerStrip.h"

// check the overlap between two files giving a complate list of strip information
void checkOverlap(string file1, string file2){

  vector<TrackerStrip> vec1;
  vector<TrackerStrip> vec2;

  ifstream infile1 (file1.c_str());
  ifstream infile2 (file2.c_str());

  TrackerStrip* strip = new TrackerStrip();
  if(infile1.is_open()){
    while(!infile1.eof()){
      infile1 >> strip->fecCrate_ >> strip->fecSlot_ >> strip->fecRing_ >> strip->ccuAdd_ >> strip->ccuCh_ >> strip->fedKey_ >> strip->lldCh_ >> strip->apvid_ >> strip->stripid_ ;
      vec1.push_back(*strip);
    }
  }
  infile1.close();

  if(infile2.is_open()){
    while(!infile2.eof()){
      infile2 >> strip->fecCrate_ >> strip->fecSlot_ >> strip->fecRing_ >> strip->ccuAdd_ >> strip->ccuCh_ >> strip->fedKey_ >> strip->lldCh_ >> strip->apvid_ >> strip->stripid_ ;
      vec2.push_back(*strip);
    }
  }
  infile2.close();

  long int total_vec1 = vec1.size();
  long int total_vec2 = vec2.size();
  long int common_strip;

  for(auto element : vec1){
    if(std::find(vec2.begin(),vec2.end(),element) != vec2.end())
      common_strip++;
  }
  
  cout<<"Total number of strips from file 1 "<<total_vec1<<endl;
  cout<<"Total number of strips from file 2 "<<total_vec2<<endl;
  cout<<"Common tagged strips "<<common_strip<<endl;
  cout<<"Fraction wrt file 1 "<<double(common_strip)/total_vec1<<endl;
  cout<<"Fraction wrt file 2 "<<double(common_strip)/total_vec2<<endl;

}

