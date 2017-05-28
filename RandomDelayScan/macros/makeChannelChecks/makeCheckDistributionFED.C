#include <string>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "../delayUtils.h"

void  makeCheckDistributionFED(string inputFileName, 
			       string outputDIR, 
			       int    fedNumber,
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
  std::cout<<"### Start loop "<<std::endl;

  // start loop on the input tree 
  vector<double> limits;
  setLimitsAndBinning("delay",limits);
  TCanvas canvas ("canvas","",600,600);
  
  map<uint32_t,TH1F*> fedHistogramMap;
  map<uint32_t,TF1*>  fedFunctionMap;

  TLegend leg (0.3,0.3,0.6,0.5);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  
  while(reader.Next()){

    if(*notFound_i) continue; // skip bad channels in the fed
    if(*fedId_i != fedNumber) continue; // skip channels not belonging to the interesting FED

    if((fedHistogramMap[*Detid_i] != 0 and fedHistogramMap[*Detid_i] != NULL) and 
       (fedFunctionMap[*Detid_i] != 0 and fedFunctionMap[*Detid_i] != NULL)) continue; // means that the same detId have been already studied

    if(fedHistogramMap[*Detid_i] == 0 or fedHistogramMap[*Detid_i] == NULL)
      fedHistogramMap[*Detid_i] = new TH1F(Form("histo_fed_%d_detid_%d",*fedId_i,*Detid_i),"",limits.size()-1,&limits[0]);
    if(fedFunctionMap[*Detid_i] == 0 or fedFunctionMap[*Detid_i] == NULL)
      fedFunctionMap[*Detid_i] = new TF1(Form("func_fed_%d_detid_%d",*fedId_i,*Detid_i),"[0]*TMath::Gaus(x,[1],[2])",limits.front(),limits.back());


    // build the histogram and the function
    fedFunctionMap[*Detid_i]->SetParameter(0,*measuredMeanAmplitude_i);
    fedFunctionMap[*Detid_i]->SetParameter(1,*measuredDelay_i);
    fedFunctionMap[*Detid_i]->SetParameter(2,*measuredSigma_i);
    fedFunctionMap[*Detid_i]->SetLineColor(kRed);
    fedFunctionMap[*Detid_i]->SetLineWidth(2);
    
    for(int iBin = 1; iBin < fedHistogramMap[*Detid_i]->GetNbinsX(); iBin++){
      if(amplitude->at(iBin) == 0) continue;
      fedHistogramMap[*Detid_i]->SetBinContent(iBin,amplitude->at(iBin));
      fedHistogramMap[*Detid_i]->SetBinError(iBin,amplitudeUnc->at(iBin));
    }
    fedHistogramMap[*Detid_i]->GetXaxis()->SetTitle("delay [ns]");
    if(observable == "maxCharge")
      fedHistogramMap[*Detid_i]->GetYaxis()->SetTitle("corrected signal (ADC)");
    else if(observable == "")
      fedHistogramMap[*Detid_i]->GetYaxis()->SetTitle("corrected S/N");
    fedHistogramMap[*Detid_i]->SetMarkerColor(kBlack);
    fedHistogramMap[*Detid_i]->SetMarkerSize(1);
    fedHistogramMap[*Detid_i]->SetMarkerStyle(20);

    fedHistogramMap[*Detid_i]->GetYaxis()->SetRangeUser(fedFunctionMap[*Detid_i]->GetMinimum(limits.front(),limits.back())*0.75,fedFunctionMap[*Detid_i]->GetMaximum(limits.front(),limits.back())*1.25);
    fedHistogramMap[*Detid_i]->Draw("PE");
    fedFunctionMap[*Detid_i]->Draw("Lsame");
    CMS_lumi(&canvas,"");

    leg.Clear();
    leg.AddEntry((TObject*)(0),Form("Fed = %d, Detid = %d",*fedId_i,*Detid_i),"");
    leg.Draw("same");
    canvas.SaveAs((outputDIR+Form("/distribution_fed_%d_detid_%d.png",*fedId_i,*Detid_i)).c_str(),"png");
    canvas.SaveAs((outputDIR+Form("/distribution_fed_%d_detid_%d.pdf",*fedId_i,*Detid_i)).c_str(),"pdf");
  }
}


