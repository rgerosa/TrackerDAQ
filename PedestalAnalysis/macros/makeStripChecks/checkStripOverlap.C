#include <algorithm>
#include "../TrackerStrip.h"

// check the overlap between two files giving a complate list of strip information
void checkStripOverlap(string file1, string file2, bool fedKey_file1 = true, bool fedKey_file2 = true){

  vector<TrackerStrip> vec1_tecp;
  vector<TrackerStrip> vec1_tecm;
  vector<TrackerStrip> vec1_tob;
  vector<TrackerStrip> vec1_tib;
  vector<TrackerStrip> vec2_tecp;
  vector<TrackerStrip> vec2_tecm;
  vector<TrackerStrip> vec2_tob;
  vector<TrackerStrip> vec2_tib;

  ifstream infile1 (file1.c_str());
  ifstream infile2 (file2.c_str());

  TrackerStrip* strip = new TrackerStrip();
  cout<<"Loading file 1 "<<endl;
  if(infile1.is_open()){
    while(!infile1.eof()){
      if(not fedKey_file1)
	infile1 >> strip->fecCrate_ >> strip->fecSlot_ >> strip->fecRing_ >> strip->ccuAdd_ >> strip->ccuCh_ >> strip->lldCh_ >> strip->apvid_ >> strip->stripid_ ;
      else
	infile1 >> strip->fecCrate_ >> strip->fecSlot_ >> strip->fecRing_ >> strip->ccuAdd_ >> strip->ccuCh_ >> strip->fedKey_ >> strip->lldCh_ >> strip->apvid_ >> strip->stripid_ ;

      if(strip->fecCrate_ == 1)
	vec1_tib.push_back(*strip);
      else if(strip->fecCrate_ == 2)
	vec1_tecp.push_back(*strip);
      else if(strip->fecCrate_ == 3)
	vec1_tecm.push_back(*strip);
      else if(strip->fecCrate_ == 4)
	vec1_tob.push_back(*strip);
    }
  }
  infile1.close();

  cout<<"Loading file 2 "<<endl;
  if(infile2.is_open()){
    while(!infile2.eof()){
      if(not fedKey_file2)
	infile2 >> strip->fecCrate_ >> strip->fecSlot_ >> strip->fecRing_ >> strip->ccuAdd_ >> strip->ccuCh_ >> strip->lldCh_ >> strip->apvid_ >> strip->stripid_ ;
      else
	infile2 >> strip->fecCrate_ >> strip->fecSlot_ >> strip->fecRing_ >> strip->ccuAdd_ >> strip->ccuCh_ >> strip->fedKey_ >> strip->lldCh_ >> strip->apvid_ >> strip->stripid_ ;
      if(strip->fecCrate_ == 1)
	vec2_tib.push_back(*strip);
      else if(strip->fecCrate_ == 2)
	vec2_tecp.push_back(*strip);
      else if(strip->fecCrate_ == 3)
	vec2_tecm.push_back(*strip);
      else if(strip->fecCrate_ == 4)
	vec2_tob.push_back(*strip);
    }
  }
  infile2.close();
  
  ////////////
  cout<<"Comparison in TIB "<<endl;
  long int total_vec1_tib = vec1_tib.size();
  long int total_vec2_tib = vec2_tib.size();
  long int common_strip_tib = 0;

  for(auto element : vec1_tib){
    if(std::find(vec2_tib.begin(),vec2_tib.end(),element) != vec2_tib.end())
      common_strip_tib++;
  }
  
  cout<<"TIB: Total number of strips from file 1 "<<total_vec1_tib<<endl;
  cout<<"TIB: Total number of strips from file 2 "<<total_vec2_tib<<endl;
  cout<<"Common tagged strips "<<common_strip_tib<<endl;
  cout<<"Fraction wrt file 1 "<<double(common_strip_tib)/total_vec1_tib<<endl;
  cout<<"Fraction wrt file 2 "<<double(common_strip_tib)/total_vec2_tib<<endl;

  ////////////
  cout<<"Comparison in TOB "<<endl;
  long int total_vec1_tob = vec1_tob.size();
  long int total_vec2_tob = vec2_tob.size();
  long int common_strip_tob = 0;

  for(auto element : vec1_tob){
    if(std::find(vec2_tob.begin(),vec2_tob.end(),element) != vec2_tob.end())
      common_strip_tob++;
  }
  
  cout<<"TOB: Total number of strips from file 1 "<<total_vec1_tob<<endl;
  cout<<"TOB: Total number of strips from file 2 "<<total_vec2_tob<<endl;
  cout<<"Common tagged strips "<<common_strip_tob<<endl;
  cout<<"Fraction wrt file 1 "<<double(common_strip_tob)/total_vec1_tob<<endl;
  cout<<"Fraction wrt file 2 "<<double(common_strip_tob)/total_vec2_tob<<endl;


  ////////////
  cout<<"Comparison in TECP "<<endl;
  long int total_vec1_tecp = vec1_tecp.size();
  long int total_vec2_tecp = vec2_tecp.size();
  long int common_strip_tecp = 0;

  for(auto element : vec1_tecp){
    if(std::find(vec2_tecp.begin(),vec2_tecp.end(),element) != vec2_tecp.end())
      common_strip_tecp++;
  }
  
  cout<<"TECP: Total number of strips from file 1 "<<total_vec1_tecp<<endl;
  cout<<"TECP: Total number of strips from file 2 "<<total_vec2_tecp<<endl;
  cout<<"Common tagged strips "<<common_strip_tecp<<endl;
  cout<<"Fraction wrt file 1 "<<double(common_strip_tecp)/total_vec1_tecp<<endl;
  cout<<"Fraction wrt file 2 "<<double(common_strip_tecp)/total_vec2_tecp<<endl;

  ////////////
  cout<<"Comparison in TECM "<<endl;
  long int total_vec1_tecm = vec1_tecm.size();
  long int total_vec2_tecm = vec2_tecm.size();
  long int common_strip_tecm = 0;

  for(auto element : vec1_tecm){
    if(std::find(vec2_tecm.begin(),vec2_tecm.end(),element) != vec2_tecm.end())
      common_strip_tecm++;    
  }
  
  cout<<"TECM: Total number of strips from file 1 "<<total_vec1_tecm<<endl;
  cout<<"TECM: Total number of strips from file 2 "<<total_vec2_tecm<<endl;
  cout<<"Common tagged strips "<<common_strip_tecm<<endl;
  cout<<"Fraction wrt file 1 "<<double(common_strip_tecm)/total_vec1_tecm<<endl;
  cout<<"Fraction wrt file 2 "<<double(common_strip_tecm)/total_vec2_tecm<<endl;

}

