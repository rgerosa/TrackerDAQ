#include <string>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"

#include "../delayUtils.h"

// take as input the output root file produced by delayValidationPerModule (TTree with floating point correction for each detId).
// It creates a new file with a TTree with only Detid, fedChannel, delay in step of 1.04 ns
void  delayCorrectionPerModule(string fileName, string outputDIR, string outputName, bool saveFits = false, float delayCutForPlotting = 4, string observable = "maxCharge"){

  setTDRStyle();
  gROOT->SetBatch(kTRUE);

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
  TTreeReaderValue<uint16_t> fedId_i    (reader,"fedId");
  TTreeReaderValue<uint16_t> fecRing_i  (reader,"fecRing");
  TTreeReaderValue<uint16_t> ccuAdd_i    (reader,"ccuAdd");
  TTreeReaderValue<uint16_t> ccuChan_i    (reader,"ccuChan");
  TTreeReaderValue<int>      notFound_i (reader,"notFound");
  TTreeReaderValue<float>    delayCorr_i (reader,"delayCorr");
  // to reconstruct the Gaussian fit vs delay for each module
  TTreeReaderValue<float>    measuredMeanAmplitude_i (reader,"measuredMeanAmplitude");
  TTreeReaderValue<float>    measuredSigma_i (reader,"measuredSigma");
  TTreeReaderValue<float>    measuredDelay_i (reader,"measuredDelay");
  TTreeReaderValue<vector<float> >  amplitude (reader,"amplitude");
  TTreeReaderValue<vector<float> >  amplitudeUnc (reader,"amplitudeUnc");

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
  std::map<std::string,std::string> largeDelayMap;

  std::cout<<"### Start loop "<<std::endl;
  // start loop on the input tree
  vector<double> limits;
  setLimitsAndBinning("delay",limits);
  TF1  fitfunc ("fitfunc","[0]*TMath::Gaus(x,[1],[2])",limits.front(),limits.back());
  TH1F profile ("profile","",limits.size()-1,&limits[0]);
  TCanvas canvas ("canvas","",600,600);
  canvas.cd();
  while(reader.Next()){
    Detid = *Detid_i;
    fedCh = *fedCh_i;
    delayCorr = std::round(*delayCorr_i*24./25.)*25/24;// delay in unitis of 24/25        
    if(delayCorr == -0) delayCorr = 0;
    outputTree->Fill();
    
    if(*notFound_i) continue;

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
      if(fabs(*delayCorr_i) > delayCutForPlotting){
	largeDelayMap[to_string(Detid)] = to_string(*fedId_i)+" "+to_string(*fecRing_i)+" "+to_string(*ccuAdd_i)+" "+to_string(*ccuChan_i);
	profile.Reset();
	if(saveFits){
	  for(int iBin = 1; iBin < profile.GetNbinsX(); iBin++){
	    if(amplitude->at(iBin) == 0) continue;
	    profile.SetBinContent(iBin,amplitude->at(iBin));
	    profile.SetBinError(iBin,amplitudeUnc->at(iBin));
	  }
	  profile.GetXaxis()->SetTitle("delay [ns]");
	  if(observable == "maxCharge")
	    profile.GetYaxis()->SetTitle("Amplitude [ADC]");	
	  else
	    profile.GetYaxis()->SetTitle("S/N");	
	  profile.SetMarkerColor(kBlack);
	  profile.SetMarkerSize(1);
	  profile.SetMarkerStyle(20);
	  profile.Draw("PE");
	  fitfunc.SetLineColor(kRed);
	  fitfunc.SetLineWidth(2);	  
	  fitfunc.Draw("Lsame");
	  profile.GetYaxis()->SetRangeUser(fitfunc.GetMinimum(limits.front(),limits.back())*0.75,fitfunc.GetMaximum(limits.front(),limits.back())*1.25);
	  CMS_lumi(&canvas,"");
	  
	  canvas.SaveAs((outputDIR+"/plot_detid_"+to_string(Detid)+".png").c_str(),"png");
	  canvas.SaveAs((outputDIR+"/plot_detid_"+to_string(Detid)+".pdf").c_str(),"pdf");
	}
      }
    }
    else{
      cerr<<"Huston we have a problem in the input delay tree --> this should never happen "<<endl;
    }    
  }

  std::cout<<"### Loop finished "<<std::endl;  
  std::cout<<"### Make output text file"<<std::endl;

  ofstream rawDelayFile ((outputDIR+"/rawDelayCorrection_"+observable+".txt").c_str());
  for(auto imap : rawDelayMap){
    rawDelayFile << imap.first << "   "<<imap.second<<"\n";
  }
  rawDelayFile.close();

  ofstream delayFile ((outputDIR+"/delayCorrection_"+observable+".txt").c_str());
  for(auto imap : delayMap){
    delayFile << imap.first << "   "<<imap.second<<"\n";
  }
  delayFile.close();

  ofstream signalGainRawDelay ((outputDIR+"/signalGainRawDelay_"+observable+".txt").c_str());
  for(auto imap : signalIncreaseVsRawDelayMap)
    signalGainRawDelay << imap.first << "   "<<imap.second<<"\n";
  signalGainRawDelay.close();

  //  ofstream signalGainDelay ((outputDIR+"/signalGainDelay_"+observable+".txt").c_str());
  //  for(auto imap : signalIncreaseVsDelayMap)
  //    signalGainDelay << imap.first << "   "<<imap.second<<"\n";
  //  signalGainDelay.close();

  //  ofstream largeDelayChannel((outputDIR+"/largeDelayChannel"+observable+".txt").c_str());
  //  for(auto imap : largeDelayMap)
  //    largeDelayChannel << imap.first << "   "<<imap.second<<"\n";
  //  largeDelayChannel.close();

  // write output
  outputFile->cd();
  outputTree->BuildIndex("Detid");  
  outputTree->Write();
  return;
}


