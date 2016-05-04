#include <string>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"

// take as input the output root file produced by delayValidationPerModule (TTree with floating point correction for each detId).
// It creates a new file with a TTree with only Detid, fedChannel, delay in step of 1.04 ns
void  delayCorrectionPerModule(string fileName, string outputDIR, string outputName){

  // compute corrections
  std::cout<<"###############################"<<std::endl;
  std::cout<<"#### computing corrections ####"<<std::endl;
  std::cout<<"###############################"<<std::endl;

  // open and create TTreeReader for the input tree
  std::shared_ptr<TFile> inputFile (TFile::Open(fileName.c_str()));
  std::shared_ptr<TTree> inputTree ((TTree*) inputFile->FindObjectAny("delayCorrection"));

  system(("mkdir -p "+outputDIR).c_str());
  
  TTreeReader reader(inputTree.get());
  TTreeReaderValue<uint32_t> Detid_i    (reader,"Detid");
  TTreeReaderValue<uint16_t> fedCh_i    (reader,"fedCh");
  TTreeReaderValue<float>    delayCorr_i (reader,"delayCorr");

  // output file and output tree structure
  std::shared_ptr<TFile> outputFile (new TFile((outputDIR+"/"+outputName+".root").c_str(),"RECREATE"));
  outputFile->cd();
  std::shared_ptr<TTree> outputTree (new TTree("delayCorrection","delayCorrection"));
  uint32_t Detid;
  uint32_t fedCh;
  float    delayCorr;
  outputTree->Branch("Detid",&Detid,"Detid/I");
  outputTree->Branch("fedCh",&fedCh,"fedCh/I");
  outputTree->Branch("delayCorr",&delayCorr,"delayCorr/F");

  std::cout<<"### Start loop "<<std::endl;
  // start loop on the input tree
  while(reader.Next()){
    Detid = *Detid_i;
    fedCh = *fedCh_i;
    delayCorr = std::round(*delayCorr_i*24./25.)*25/24;// delay in unitis of 24/25    
    outputTree->Fill();
  }

  std::cout<<"### Loop finished "<<std::endl;

  // write output
  outputFile->cd();
  outputTree->BuildIndex("Detid");  
  outputTree->Write();
  return;
}


