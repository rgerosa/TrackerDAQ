
// Generate a file with a simple readout map from a DPG tree
void generateCablingTree(string inputDPGTree, string ouputDIR, string ouputFileName){
  
  //Open file
  TFile* inputFile = TFile::Open(inputDPGTree.c_str(),"READ");
  inputFile->cd();
  TTree* readoutMap = (TTree*) inputFile->FindObjectAny("readoutMap");

  system(("mkdir -p "+ouputDIR).c_str());
  
  TFile* outputFile = new TFile((ouputDIR+"/"+ouputFileName).c_str(),"RECREATE");
  outputFile->cd();
  readoutMap->SetBranchStatus("*",kFALSE);
  readoutMap->SetBranchStatus("fecCrate",kTRUE);
  readoutMap->SetBranchStatus("fecSlot",kTRUE);
  readoutMap->SetBranchStatus("fecRing",kTRUE);
  readoutMap->SetBranchStatus("ccuAdd",kTRUE);
  readoutMap->SetBranchStatus("ccuChan",kTRUE);
  readoutMap->SetBranchStatus("lldChannel",kTRUE);
  readoutMap->SetBranchStatus("fedId",kTRUE);
  readoutMap->SetBranchStatus("fedCh",kTRUE);
  readoutMap->SetBranchStatus("detid",kTRUE);
  
  TTree* outputTree = (TTree*) readoutMap->CopyTree("");
  outputTree->SetName("readoutMap");
  outputTree->BuildIndex("detid");
  outputTree->Write(outputTree->GetName(),TObject::kOverwrite);
  outputFile->Close();
  inputFile->Close();
}
