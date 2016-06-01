#include <string>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"

#include "delayUtils.h"

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
  // to reconstruct the Gaussian fit vs delay for each module
  TTreeReaderValue<float>    measuredMeanAmplitude_i (reader,"measuredMeanAmplitude");
  TTreeReaderValue<float>    measuredSigma_i (reader,"measuredSigma");
  TTreeReaderValue<float>    measuredDelay_i (reader,"measuredDelay");

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
  
  std::map<std::string,std::string> rawDelayMap;
  std::map<std::string,std::string> delayMap;
  std::map<std::string,std::string> signalIncreaseVsRawDelayMap;
  std::map<std::string,std::string> signalIncreaseVsDelayMap;

  std::cout<<"### Start loop "<<std::endl;
  // start loop on the input tree
  vector<double> limits;
  setLimitsAndBinning("delay",limits);
  TF1 fitfunc ("fitfunc","[0]*TMath::Gaus(x,[1],[2])",limits.front(),limits.back());

  while(reader.Next()){
    Detid = *Detid_i;
    fedCh = *fedCh_i;
    delayCorr = std::round(*delayCorr_i*24./25.)*25/24;// delay in unitis of 24/25    
    outputTree->Fill();
    
    rawDelayMap[to_string(Detid)] = to_string(*delayCorr_i);
    delayMap[to_string(Detid)]    = to_string(delayCorr);

    // estimate the gain in signal amplitude
    if(*delayCorr_i == 0){ // this happens either when no data for the module are available or some quality cuts
      signalIncreaseVsRawDelayMap[to_string(Detid)] = to_string(1);
      signalIncreaseVsDelayMap[to_string(Detid)] = to_string(1);
    }
    else if(*measuredDelay_i == *delayCorr_i){ // we need to reconstruct the full Gaussian fit of the fedChannel (AOH channel)
      fitfunc.SetParameter(0,*measuredMeanAmplitude_i);
      fitfunc.SetParameter(1,*measuredDelay_i);
      fitfunc.SetParameter(2,*measuredSigma_i);
      signalIncreaseVsRawDelayMap[to_string(Detid)] = to_string(fitfunc.Eval(*measuredDelay_i)/fitfunc.Eval(0));
      signalIncreaseVsDelayMap[to_string(Detid)] = to_string(fitfunc.Eval(delayCorr)/fitfunc.Eval(0));
    }
    else{
      cerr<<"Huston we have a problem in the input delay tree --> this should never happen "<<endl;
    }    
  }

  std::cout<<"### Loop finished "<<std::endl;  
  std::cout<<"### Make output text file"<<std::endl;

  ofstream rawDelayFile ((outputDIR+"/rawDelayCorrection.txt").c_str());
  for(auto imap : rawDelayMap){
    rawDelayFile << imap.first << "   "<<imap.second<<"\n";
  }
  rawDelayFile.close();

  ofstream delayFile ((outputDIR+"/delayCorrection.txt").c_str());
  for(auto imap : delayMap){
    delayFile << imap.first << "   "<<imap.second<<"\n";
  }
  delayFile.close();

  ofstream signalGainRawDelay ((outputDIR+"/signalGainRawDelay.txt").c_str());
  for(auto imap : signalIncreaseVsRawDelayMap)
    signalGainRawDelay << imap.first << "   "<<imap.second<<"\n";
  signalGainRawDelay.close();

  ofstream signalGainDelay ((outputDIR+"/signalGainDelay.txt").c_str());
  for(auto imap : signalIncreaseVsDelayMap)
    signalGainDelay << imap.first << "   "<<imap.second<<"\n";
  signalGainDelay.close();
    
  // write output
  outputFile->cd();
  outputTree->BuildIndex("Detid");  
  outputTree->Write();
  return;
}


