#include "../CMS_lumi.h"

static int  reductionFactor = 1;

void plotCanvas(TCanvas* canvas, vector<TH1F*> histo, string outputDIR){

  canvas->cd();  
  unsigned int i = 0;
  TLegend leg (0.6,0.6,0.9,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  for(auto hist : histo){
    hist->GetXaxis()->SetTitle("p-value");
    hist->GetYaxis()->SetTitle("Entries");
    hist->SetLineColor(i+1);
    hist->SetLineWidth(2);
    if(i == 0){
      hist->Draw("hist");
    }
    else hist->Draw("hist same");
    leg.AddEntry(hist,hist->GetTitle(),"L");
    i++;	       
  }
  leg.Draw("same");
  
  CMS_lumi(canvas,"",true);
  histo.at(0)->GetYaxis()->SetRangeUser(0,1000*1e3);
  canvas->SaveAs((outputDIR+"/testStatistics.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/testStatistics.pdf").c_str(),"pdf");
  histo.at(0)->GetYaxis()->SetRangeUser(10*1e2,100*1e6);
  canvas->SetLogy();
  canvas->SaveAs((outputDIR+"/testStatistics_log.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/testStatistics_log.pdf").c_str(),"pdf");

}

void plotTestStatistics(string inputFileName,string outputDIR){

  system(("mkdir -p "+outputDIR).c_str());

  setTDRStyle();
  gROOT->SetBatch(kTRUE);

  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  inputFile->cd();
  TTree* tree = (TTree*) inputFile->Get("pedestalFullNoise");

  TH1F* ksPValue       = new TH1F("ksPValue","Kolmogorov Smirnov",75,0,1);
  TH1F* adPValue = new TH1F("adPValue","Anderson Darling",  75,0,1);
  TH1F* jbPValue = new TH1F("jbPValue","Jacque Bera",75,0,1);
  TH1F* chi2PValue = new TH1F("chi2PValue","#chi^{2} test",75,0,1);

  ksPValue->Sumw2();
  adPValue->Sumw2();
  jbPValue->Sumw2();
  chi2PValue->Sumw2();

  float    fitChi2Probab, kSProbab, jBProbab, aDProbab;
  float    nBin, xMin, xMax;
  vector<float>* noiseDistribution = 0;

  tree->SetBranchStatus("*",kFALSE);
  tree->SetBranchStatus("fitChi2Probab",kTRUE);
  tree->SetBranchStatus("kSProbab",kTRUE);
  tree->SetBranchStatus("jBProbab",kTRUE);
  tree->SetBranchStatus("aDProbab",kTRUE);
  tree->SetBranchStatus("noiseDistribution",kTRUE);
  tree->SetBranchStatus("nBin",kTRUE);
  tree->SetBranchStatus("xMin",kTRUE);
  tree->SetBranchStatus("xMax",kTRUE);

  tree->SetBranchAddress("fitChi2Probab",&fitChi2Probab);
  tree->SetBranchAddress("kSProbab",&kSProbab);
  tree->SetBranchAddress("aDProbab",&aDProbab);
  tree->SetBranchAddress("jBProbab",&jBProbab);
  tree->SetBranchAddress("noiseDistribution",&noiseDistribution);
  tree->SetBranchAddress("nBin",&nBin);
  tree->SetBranchAddress("xMin",&xMin);
  tree->SetBranchAddress("xMax",&xMax);
 

  TCanvas* canvas = new TCanvas("canvas","canvas",600,650);
  TH1F* noiseHist = NULL;

  for(long int iChannel = 0; iChannel < tree->GetEntries(); iChannel++){
    tree->GetEntry(iChannel);
    cout.flush();
    if(iChannel %10000 == 0) cout<<"\r"<<"iChannel "<<100*double(iChannel)/(tree->GetEntries()/reductionFactor)<<" % ";
    if(iChannel > double(tree->GetEntries())/reductionFactor) break;

    if(noiseHist == NULL){
      noiseHist = new TH1F ("noiseHist","",nBin,xMin,xMax);
      noiseHist->Sumw2();
    }
    noiseHist->Reset();

    // create the noise distribution for the given strip                                                                                                                                               
    for(int iBin = 0; iBin < noiseDistribution->size(); iBin++){
      noiseHist->SetBinContent(iBin+1,noiseDistribution->at(iBin));
    }
    
    //////// basic quality cuts
    if(noiseHist->Integral() == 0) continue;
    if(noiseHist->GetRMS() < 1.5) continue;
    if((noiseHist->GetBinContent(1)+noiseHist->GetBinContent(noiseHist->GetNbinsX()))/noiseHist->Integral() >= 0.3) continue;
   
    ksPValue->Fill(kSProbab);
    chi2PValue->Fill(fitChi2Probab);
    adPValue->Fill(aDProbab);
    jbPValue->Fill(jBProbab);
  }

  vector<TH1F*> vecHisto; 
  vecHisto.push_back(ksPValue);
  vecHisto.push_back(adPValue);
  vecHisto.push_back(jbPValue);
  vecHisto.push_back(chi2PValue);

  plotCanvas(canvas,vecHisto,outputDIR);
}
