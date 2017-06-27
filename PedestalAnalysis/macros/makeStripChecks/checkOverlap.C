/// check overlap across two bad channel files in the format (detid - number of bad strips)
void checkOverlap(string file1, string file2){

  map<string,string> map1; // det-id, number of strips
  map<string,string> map2;

  ifstream infile1 (file1.c_str());
  ifstream infile2 (file2.c_str());

  if(infile1.is_open()){
    while(!infile1.eof()){
      string detid, nstrip;
      infile1 >> detid >> nstrip;
      map1[detid] = nstrip;
    }
  }
  infile1.close();

  if(infile2.is_open()){
    while(!infile2.eof()){
      string detid, nstrip;
      infile2 >> detid >> nstrip;
      map2[detid] = nstrip;
    }
  }
  infile2.close();

  long int total = 0;
  long int common_detid = 0;
  long int common_detid_nstrip = 0;

  for(auto element : map1){
    total++;
    if(map2.find(element.first) != map2.end()){
      common_detid++;
      if(element.second == map2[element.first])
	common_detid_nstrip++;
    }
  }
  
  cout<<"Total number of strips "<<total<<" common det id "<<common_detid<<" "<<100*float(common_detid)/total<<" % "<<endl;
  cout<<"Total number of strips "<<total<<" common det id and nstrip "<<common_detid_nstrip<<" "<<100*float(common_detid_nstrip)/total<<" % "<<endl;

}

