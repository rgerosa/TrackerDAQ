#include <string>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "../delayUtils.h"

void  makeCheckDistributionLargeDelay(string inputFileName, 
				      string outputDIR, 
				      float  delayThreshold,
				      string observable){

  setTDRStyle();
  gROOT->SetBatch(kTRUE);

  // compute corrections
  std::cout<<"###############################"<<std::endl;
  std::cout<<"#### computing corrections ####"<<std::endl;
  std::cout<<"###############################"<<std::endl;

  // open and create TTreeReader for the input tree
  std::shared_ptr<TFile> inputFile (TFile::Open(inputFileName.c_str()));
  std::shared_ptr<TTree> inputTree ((TTree*) inputFile->FindObjectAny("delayCorrection"));
  
  system(("mkdir -p "+outputDIR).c_str());
  
  TTreeReader reader(inputTree.get());
  TTreeReaderValue<uint32_t> Detid_i     (reader,"Detid");
  TTreeReaderValue<uint16_t> fedCh_i     (reader,"fedCh");
  TTreeReaderValue<uint16_t> fedId_i     (reader,"fedId");
  TTreeReaderValue<uint16_t> fecRing_i   (reader,"fecRing");
  TTreeReaderValue<uint16_t> fecSlot_i   (reader,"fecSlot");
  TTreeReaderValue<uint16_t> fecCrate_i   (reader,"fecCrate");
  TTreeReaderValue<uint16_t> ccuAdd_i    (reader,"ccuAdd");
  TTreeReaderValue<uint16_t> ccuChan_i   (reader,"ccuChan");
  TTreeReaderValue<uint16_t> lldChannel_i   (reader,"lldChannel");
  TTreeReaderValue<int>      notFound_i  (reader,"notFound");
  TTreeReaderValue<float>    delayCorr_i (reader,"delayCorr");
  // to reconstruct the Gaussian fit vs delay for each module
  TTreeReaderValue<float>    measuredMeanAmplitude_i (reader,"measuredMeanAmplitude");
  TTreeReaderValue<float>    measuredSigma_i (reader,"measuredSigma");
  TTreeReaderValue<float>    measuredDelay_i (reader,"measuredDelay");
  TTreeReaderValue<vector<float> >  amplitude (reader,"amplitude");
  TTreeReaderValue<vector<float> >  amplitudeUnc (reader,"amplitudeUnc");

  // output file and output tree structure
  std::cout<<"### Start loop "<<std::endl;

  // start loop on the input tree 
  vector<double> limits;
  setLimitsAndBinning("delay",limits);
  TCanvas canvas ("canvas","",600,600);
  
  map<uint32_t,TH1F*> largeDelayHistogramMap;
  map<uint32_t,TF1*>  largeDelayFunctionMap;

  TLegend leg (0.3,0.3,0.6,0.5);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);

  ofstream channelsLargeDelay ((outputDIR+"/channelsLargeDelay.txt").c_str());
  ofstream channelsLargeDelayMap ((outputDIR+"/channelsLargeDelayMap.txt").c_str());
  
  while(reader.Next()){

    if(*notFound_i) continue; // skip bad channels in the fed
    if(fabs(*delayCorr_i) < delayThreshold) continue;

    channelsLargeDelay<<" Detid_i "<<*Detid_i<<" fedId "<<*fedId_i<<" fedCh "<<*fedCh_i<<" fecCrate "<<*fecCrate_i<<" fecSlot "<<*fecSlot_i<<" fecRing "<<*fecRing_i<<" ccuAdd "<<*ccuAdd_i<<" ccuChan "<<*ccuChan_i<<" lldChannel "<<*lldChannel_i<<"\n";
    channelsLargeDelayMap<<*Detid_i<<" 1"<<"\n";

    if((largeDelayHistogramMap[*Detid_i] != 0 and largeDelayHistogramMap[*Detid_i] != NULL) and 
       (largeDelayFunctionMap[*Detid_i] != 0 and largeDelayFunctionMap[*Detid_i] != NULL)) continue; // means that the same detId have been already studied

    if(largeDelayHistogramMap[*Detid_i] == 0 or largeDelayHistogramMap[*Detid_i] == NULL)
      largeDelayHistogramMap[*Detid_i] = new TH1F(Form("histo_largeDelay_detid_%d",*Detid_i),"",limits.size()-1,&limits[0]);
    if(largeDelayFunctionMap[*Detid_i] == 0 or largeDelayFunctionMap[*Detid_i] == NULL)
      largeDelayFunctionMap[*Detid_i] = new TF1(Form("func_largeDelay_detid_%d",*Detid_i),"[0]*TMath::Gaus(x,[1],[2])",limits.front(),limits.back());

    // build the histogram and the function
    largeDelayFunctionMap[*Detid_i]->SetParameter(0,*measuredMeanAmplitude_i);
    largeDelayFunctionMap[*Detid_i]->SetParameter(1,*measuredDelay_i);
    largeDelayFunctionMap[*Detid_i]->SetParameter(2,*measuredSigma_i);
    largeDelayFunctionMap[*Detid_i]->SetLineColor(kRed);
    largeDelayFunctionMap[*Detid_i]->SetLineWidth(2);
    
    for(int iBin = 1; iBin < largeDelayHistogramMap[*Detid_i]->GetNbinsX(); iBin++){
      if(amplitude->at(iBin) == 0) continue;
      largeDelayHistogramMap[*Detid_i]->SetBinContent(iBin,amplitude->at(iBin));
      largeDelayHistogramMap[*Detid_i]->SetBinError(iBin,amplitudeUnc->at(iBin));
    }
    largeDelayHistogramMap[*Detid_i]->GetXaxis()->SetTitle("delay [ns]");
    if(observable == "maxCharge")
      largeDelayHistogramMap[*Detid_i]->GetYaxis()->SetTitle("corrected signal (ADC)");
    else if(observable == "")
      largeDelayHistogramMap[*Detid_i]->GetYaxis()->SetTitle("corrected S/N");
    largeDelayHistogramMap[*Detid_i]->SetMarkerColor(kBlack);
    largeDelayHistogramMap[*Detid_i]->SetMarkerSize(1);
    largeDelayHistogramMap[*Detid_i]->SetMarkerStyle(20);

    largeDelayHistogramMap[*Detid_i]->GetYaxis()->SetRangeUser(largeDelayFunctionMap[*Detid_i]->GetMinimum(limits.front(),limits.back())*0.75,largeDelayFunctionMap[*Detid_i]->GetMaximum(limits.front(),limits.back())*1.25);
    largeDelayHistogramMap[*Detid_i]->Draw("PE");
    largeDelayFunctionMap[*Detid_i]->Draw("Lsame");
    CMS_lumi(&canvas,"");

    uint32_t subdetid    = int((*Detid_i-0x10000000)/0x2000000);
    string postfix = "TIB";
    if(subdetid == 4)
      postfix = "TID";
    if(subdetid == 5)
      postfix = "TOB";
    if(subdetid == 6)
      postfix = "TEC";
    
    leg.Clear();
    leg.AddEntry((TObject*)(0),Form("%s, Fed = %d, Detid = %d",postfix.c_str(),*fedId_i,*Detid_i),"");
    leg.Draw("same");
    canvas.SaveAs((outputDIR+Form("/distribution_%s_detid_%d.png",postfix.c_str(),*Detid_i)).c_str(),"png");
    canvas.SaveAs((outputDIR+Form("/distribution_%s_detid_%d.pdf",postfix.c_str(),*Detid_i)).c_str(),"pdf");
  }

  channelsLargeDelay.close();
  channelsLargeDelayMap.close();
}


