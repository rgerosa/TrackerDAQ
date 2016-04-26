void TreeMerge(string destination, string reference, string sources ) {
  
  // make hadd from a command line
  gSystem->Exec(Form("hadd -f %s %s",destination.c_str(),sources.c_str()));

  // then recompute the indices
  std::cout << "Reindexing the element of the merged tree..."  << std::endl;
  TFile* f = TFile::Open(destination.c_str(),"update");

  //clusters tree is the default one --> other added as friend in the analysis
  // build an index on events tree according to runId and eventid
  TTree* tree = (TTree*)f->Get("analysis/trackerDPG/events");
  if(tree && tree->GetEntries())  
    tree->BuildIndex("runid","eventid"); // make a index based on runid:eventid

  // build an index on vertex tree according to vertexid and eventid
  tree = (TTree*)f->Get("analysis/trackerDPG/vertices");
  if(tree && tree->GetEntries())  tree->BuildIndex("vertexid","eventid"); // make a index based on runid:eventid
  // build an index on vertex tree according to trackid and eventid 
  int i=0;
  while((tree = (TTree*)f->Get(Form("analysis/trackerDPG/tracks%d",i)))) {
    if(tree->GetEntries()) 
      tree->BuildIndex(Form("trackid%d",i),"eventid"); // index based on trackid and eventid
    i++;
  }

  // copy them from a referece
  std::cout << "patching..." << std::endl;
  TFile* f2 = TFile::Open(reference.c_str());
  tree = (TTree*)f2->Get("analysis/trackerDPG/psumap");
  f->cd("analysis/trackerDPG");
  TTree *newtree = tree->CloneTree();
  // add an index per dcuId
  newtree->BuildIndex("dcuId"); // deti-id index
  newtree->Write("psumap",TObject::kOverwrite);

  tree = (TTree*)f2->Get("analysis/trackerDPG/readoutMap");
  f->cd("analysis/trackerDPG");
  newtree = tree->CloneTree();
  // add an index per detid
  newtree->BuildIndex("detid"); // det-id index
  newtree->Write("readoutMap",TObject::kOverwrite);
  f2->Close();
  
  //close
  std::cout << "Saving merged file." << std::endl;
  f->Write();
  f->Close();
}
