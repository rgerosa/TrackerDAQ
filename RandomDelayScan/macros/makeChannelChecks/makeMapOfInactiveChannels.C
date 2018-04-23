#include <string>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "../delayUtils.h"

void  makeMapOfInactiveChannels(string inputFileName, 
				string outputDIR,
				string inputPSUMap = "$CMSSW_BASE/src/TrackerDAQAnalysis/RandomDelayScan/data/readoutMapWithPSU.root"
				){

  setTDRStyle();
  gROOT->SetBatch(kTRUE);

  // open and create TTreeReader for the input tree
  std::shared_ptr<TFile> inputFile (TFile::Open(inputFileName.c_str()));
  std::shared_ptr<TTree> inputTree ((TTree*) inputFile->FindObjectAny("delayCorrection"));

  system(("mkdir -p "+outputDIR).c_str());

  uint32_t Detid_i;
  uint16_t fedCh_i, fedId_i, fecCrate_i, fecRing_i, ccuAdd_i, ccuChan_i, lldChannel_i, fecSlot_i;
  int notFound_i;

  inputTree->SetBranchAddress("Detid",&Detid_i);
  inputTree->SetBranchAddress("fedCh",&fedCh_i);
  inputTree->SetBranchAddress("fedId",&fedId_i);
  inputTree->SetBranchAddress("fecCrate",&fecCrate_i);
  inputTree->SetBranchAddress("fecSlot",&fecSlot_i);
  inputTree->SetBranchAddress("fecRing",&fecRing_i);
  inputTree->SetBranchAddress("ccuAdd",&ccuAdd_i);
  inputTree->SetBranchAddress("ccuChan",&ccuChan_i);
  inputTree->SetBranchAddress("lldChannel",&lldChannel_i);
  inputTree->SetBranchAddress("notFound",&notFound_i);

  TFile* inputPSUMapFile = TFile::Open(inputPSUMap.c_str(),"READ");
  TTree* psuTree         = (TTree*) inputPSUMapFile->Get("readoutMap");
  uint32_t detid, dcuId;
  string PSUName;
  psuTree->SetBranchAddress("PSUname",&PSUName);
  psuTree->SetBranchAddress("detid",&detid);
  psuTree->SetBranchAddress("dcuId",&dcuId);

  // output file and output tree structure
  std::cout<<"### Start loop "<<std::endl;

  ofstream inactiveChannelFile ((outputDIR+"/inactiveChannels.txt").c_str());
  ofstream inactiveChannelMap  ((outputDIR+"/inactiveChannelsMap.txt").c_str());
 
  for(int iEvent = 0; iEvent < inputTree->GetEntries(); iEvent++){
    inputTree->GetEntry(iEvent);
    if(notFound_i == 0) continue; // skip bad channels in the fed
    inactiveChannelMap <<Detid_i<<" 1 "<<"\n";
    psuTree->GetEntryWithIndex(Detid_i);    
    inactiveChannelFile <<" Detid_i "<<Detid_i<<" fedId "<<fedId_i<<" fedCh "<<fedCh_i<<" fecCrate "<<fecCrate_i<<" fecSlot "<<fecSlot_i<<" fecRing "<<fecRing_i<<" ccuAdd "<<ccuAdd_i<<" ccuChan "<<ccuChan_i<<" lldChannel "<<lldChannel_i<<" PSU "<<PSUName<<"\n";
    
  }
  std::cout<<"### Close files "<<std::endl;
  inactiveChannelFile.close();
  inactiveChannelMap.close();
  

}


