// Merge all the trees beloging to a delay setting: destination = merged file, reference is a single to take the delay map, sources are files to be merged
void TreeMerge(string destination, string reference, string inputPath, bool cancelInputFiles = false) {

  system(("rm "+destination).c_str());

  // make list of all inputFIles before creting the merge oen
  system(("ls "+inputPath+" | grep root > file.temp").c_str());
  ifstream infile;
  string line;
  vector<string> fileList;
  infile.open("file.temp",ifstream::in);
  if(infile.is_open()){
    while(!infile.eof()){
      getline(infile,line);
      if(line != "" and TString(line).Contains(".root"))
	fileList.push_back(line);
    }
  }
  system("rm file.temp");

  // make hadd from a command line
  gSystem->Exec(Form("hadd -f %s %s/*root",destination.c_str(),inputPath.c_str()));

  std::cout << "Reindexing the element of the merged tree..."  << std::endl;
  TFile* f = TFile::Open(destination.c_str(),"update");

  TFile* f2 = TFile::Open(reference.c_str());
  TTree* tree = (TTree*)f2->FindObjectAny("psumap");
  if(tree == 0 or tree == NULL){
    cout<<"[TreeMerge] no psumap file --> stop here"<<endl;
    return;
  }

  f->cd();
  TTree *newtree = tree->CloneTree();
  newtree->BuildIndex("dcuId"); // deti-id index
  newtree->Write("psumap",TObject::kOverwrite);

  // Look at the redoutMap
  tree = (TTree*)f2->FindObjectAny("readoutMap");
  if(tree == 0 or tree == NULL){
    cout<<"[TreeMerge] no readoutMap file --> stop here"<<endl;
    return;
  }
  
  f->cd();
  newtree = tree->CloneTree();
  newtree->BuildIndex("detid"); // det-id index
  newtree->Write("readoutMap",TObject::kOverwrite);
  f2->Close();
  
  //close
  std::cout << "Saving merged file." << std::endl;
  f->Write();
  f->Close();

  //cancel all the single files, keeping only the unmerged one
  if(cancelInputFiles){
    for(auto fileName : fileList){
      cout<<"rm "+inputPath+"/"+fileName<<endl;
      system(("rm "+inputPath+"/"+fileName).c_str());
    }
  }
}
