#include <map>
#include <vector>
#include <TTree.h>
#include <TFile.h>
#include <stdint.h>
#include <iostream>

void fillDelayTree(const char* file1, const char* file2, const char* file3)
{
  TTree* tree = NULL;
  std::map<uint32_t, vector<float> > delayMap;
  // fill the map from the 3 input files
  TFile* file1_ = TFile::Open(file1);
  tree = (TTree*)file1_->Get("tree");
  Double_t DETID = 0.;
  Double_t DELAY = 0.;
  tree->SetBranchAddress("DETID", &DETID);
  tree->SetBranchAddress("DELAY", &DELAY);
  for(int i=0; i<tree->GetEntries();++i) {
    tree->GetEntry(i);
    delayMap[uint32_t(round(DETID+0.1))] = std::vector<float>(3);
    delayMap[uint32_t(round(DETID+0.1))][0] = DELAY;
  }
  file1_->Close();
  TFile* file2_ = TFile::Open(file2);
  tree = (TTree*)file2_->Get("tree");
  tree->SetBranchAddress("DETID", &DETID);
  tree->SetBranchAddress("DELAY", &DELAY);
  for(int i=0; i<tree->GetEntries();++i) {
    tree->GetEntry(i);
    if(delayMap.find(uint32_t(round(DETID+0.1)))==delayMap.end())
       delayMap[uint32_t(round(DETID+0.1))] = std::vector<float>(3);
    delayMap[uint32_t(round(DETID+0.1))][1] = DELAY;
  }
  file2_->Close();
  TFile* file3_ = TFile::Open(file3);
  tree = (TTree*)file3_->Get("tree");
  tree->SetBranchAddress("DETID", &DETID);
  tree->SetBranchAddress("DELAY", &DELAY);
  for(int i=0; i<tree->GetEntries();++i) {
    tree->GetEntry(i);
    if(delayMap.find(uint32_t(round(DETID+0.1)))==delayMap.end())
       delayMap[uint32_t(round(DETID+0.1))] = std::vector<float>(3);
    delayMap[uint32_t(round(DETID+0.1))][2] = DELAY;
  }
  file3_->Close();
  // create the output and fill the tree there
  uint32_t detid;
  float DELAY1, DELAY2, DELAY3;
  TFile* file0_ = TFile::Open("PLLdelays.root","RECREATE");
  tree = new TTree("PLLdelays","PLLdelays");
  tree->Branch("detid",&detid,"detid/i");
  tree->Branch("delay1",&DELAY1);
  tree->Branch("delay2",&DELAY2);
  tree->Branch("delay3",&DELAY3);
  for(std::map<uint32_t, vector<float> >::const_iterator it = delayMap.begin(); it!=delayMap.end(); ++it) {
    detid  = it->first;
    DELAY1 = it->second[0];
    DELAY2 = it->second[1];
    DELAY3 = it->second[2];
    tree->Fill();
  }
  tree->BuildIndex("detid");
  tree->Write();
  file0_->Write();
  file0_->Close();
}
