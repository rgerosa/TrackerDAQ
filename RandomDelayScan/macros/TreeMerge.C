//TODO: retest write on castor
void TreeMerge(const char* destination, const char* reference, const char* sources ) 
{
  gSystem->Exec(Form("hadd -f %s %s",destination,sources));

  // then recompute the indices
  std::cout << "Reindexing..."  << std::endl;
  //TFile* f = TFile::Open("trackerDPG_PFGskim3-BeamBackground_merged.root","update");
  TFile* f = TFile::Open(destination,"update");
  TTree* tree;
  tree = (TTree*)f->Get("analysis/trackerDPG/events");
  if(tree && tree->GetEntries()) tree->BuildIndex("runid","eventid");
  tree = (TTree*)f->Get("analysis/trackerDPG/vertices");
  if(tree && tree->GetEntries()) tree->BuildIndex("vertexid","eventid");
  int i=0;
  while((tree = (TTree*)f->Get(Form("analysis/trackerDPG/tracks%d",i)))) {
    if(tree->GetEntries()) tree->BuildIndex(Form("trackid%d",i),"eventid");
    i++;
  }

  // what to do for psumap and readoutmap ?
  // copy from any source file
  std::cout << "patching..." << std::endl;
  TFile* f2 = TFile::Open(reference);
  tree = (TTree*)f2->Get("analysis/trackerDPG/psumap");
  f->cd("analysis/trackerDPG");
  TTree *newtree = tree->CloneTree();
  newtree->Write("psumap",TObject::kOverwrite);
//  tree->Write("psumap",TObject::kOverwrite);
  tree = (TTree*)f2->Get("analysis/trackerDPG/readoutMap");
  f->cd("analysis/trackerDPG");
  newtree = tree->CloneTree();
  newtree->Write("readoutMap",TObject::kOverwrite);
//  tree->Write("readoutMap",TObject::kOverwrite);
  f2->Close();

  //close
  std::cout << "Saving merged file." << std::endl;
  f->Write();
  f->Close();
}
