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

  TH1F* ksPValue       = new TH1F("ksPValue","Kolmogorov Smirnov",500,0,1);
  TH1F* adPValue = new TH1F("adPValue","Anderson Darling",  500,0,1);
  TH1F* jbPValue = new TH1F("jbPValue","Jacque Bera",500,0,1);
  TH1F* chi2PValue = new TH1F("chi2PValue","#chi^{2} test",500,0,1);

  ksPValue->Sumw2();
  adPValue->Sumw2();
  jbPValue->Sumw2();
  chi2PValue->Sumw2();

  float    fitChi2Probab, kSProbab, jBProbab, aDProbab;

  tree->SetBranchStatus("*",kFALSE);
  tree->SetBranchStatus("fitChi2Probab",kTRUE);
  tree->SetBranchStatus("kSProbab",kTRUE);
  tree->SetBranchStatus("jBProbab",kTRUE);
  tree->SetBranchStatus("aDProbab",kTRUE);

  tree->SetBranchAddress("fitChi2Probab",&fitChi2Probab);
  tree->SetBranchAddress("kSProbab",&kSProbab);
  tree->SetBranchAddress("aDProbab",&aDProbab);
  tree->SetBranchAddress("jBProbab",&jBProbab);

  TCanvas* canvas = new TCanvas("canvas","canvas",600,650);

  for(long int iChannel = 0; iChannel < tree->GetEntries(); iChannel++){
    tree->GetEntry(iChannel);
    cout.flush();
    if(iChannel %10000 == 0) cout<<"\r"<<"iChannel "<<100*double(iChannel)/(tree->GetEntries()/reductionFactor)<<" % ";
    if(iChannel > double(tree->GetEntries())/reductionFactor) break;

    ksPValue->Fill(kSProbab);
    chi2PValue->Fill(fitChi2Probab);
    adPValue->Fill(aDProbab);
    jbPValue->Fill(jBProbab);
  }

  vector<TH1F*> vecHisto; 
  vecHisto.push_back(ksPValue);
  //  vecHisto.push_back(adPValue);
  //  vecHisto.push_back(jbPValue);
  //  vecHisto.push_back(chi2PValue);

  plotCanvas(canvas,vecHisto,outputDIR);
}
