#include "../CMS_lumi.h"

void makeDelayVsPartition(string inputFileName, string outputDIR){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());

  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");

  TH1F* distribution_tib = (TH1F*) inputFile->Get("TIB_ring_1_mean");
  TH1F* distribution_tob = (TH1F*) inputFile->Get("TOB_ring_1_mean");
  TH1F* distribution_tid = (TH1F*) inputFile->Get("TID_ring_1_mean");
  TH1F* distribution_tecp_t = (TH1F*) inputFile->Get("TECPt_ring_1_mean");
  TH1F* distribution_tecp_T = (TH1F*) inputFile->Get("TECPT_ring_1_mean");
  TH1F* distribution_tecm_t = (TH1F*) inputFile->Get("TECMt_ring_1_mean");
  TH1F* distribution_tecm_T = (TH1F*) inputFile->Get("TECMT_ring_1_mean");

  TH1F* distribution_tec_T = (TH1F*) distribution_tecp_T->Clone("distribution_tec_T");
  TH1F* distribution_tec_t = (TH1F*) distribution_tecp_t->Clone("distribution_tec_t");
  distribution_tec_T->Reset();
  distribution_tec_t->Reset();
  
  for(int iBinX = 0; iBinX < distribution_tecp_t->GetNbinsX(); iBinX++){
    distribution_tec_t->SetBinContent(iBinX+1,(distribution_tecp_t->GetBinContent(iBinX+1)+distribution_tecm_t->GetBinContent(iBinX+1))/2);
    distribution_tec_t->SetBinError(iBinX+1,sqrt(distribution_tecp_t->GetBinError(iBinX+1)*distribution_tecp_t->GetBinError(iBinX+1)+distribution_tecm_t->GetBinError(iBinX+1)*distribution_tecm_t->GetBinError(iBinX+1))/sqrt(2));
    distribution_tec_T->SetBinContent(iBinX+1,(distribution_tecp_T->GetBinContent(iBinX+1)+distribution_tecm_T->GetBinContent(iBinX+1))/2);
    distribution_tec_T->SetBinError(iBinX+1,sqrt(distribution_tecp_T->GetBinError(iBinX+1)*distribution_tecp_T->GetBinError(iBinX+1)+distribution_tecm_T->GetBinError(iBinX+1)*distribution_tecm_T->GetBinError(iBinX+1))/sqrt(2));

  }

  TF1* gaus = new TF1("gaus","gaus(0)",-10,10);
  distribution_tec_t->Fit(gaus,"RMQ");
  distribution_tec_T->Fit(gaus,"RMQ");

  // plotting results
  TCanvas* canvas = new TCanvas("canvas","canvas",650,650);
  canvas->cd();

  TH1F* frame = (TH1F*) distribution_tib->Clone("frame");
  frame->Reset();
  frame->GetXaxis()->SetTitle("Delay (ns)");
  frame->GetXaxis()->SetTitleOffset(1);
  frame->GetXaxis()->SetTitleSize(0.045);
  
  frame->GetXaxis()->SetLabelSize(0.035);
  frame->GetYaxis()->SetTitle("Leading strip signal (ADC)");
  frame->GetYaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetLabelSize(0.035);
  frame->GetYaxis()->SetTitleOffset(1.2);
  frame->GetYaxis()->SetRangeUser(min(distribution_tib->GetMinimum(),min(distribution_tob->GetMinimum(),
				  min(distribution_tid->GetMinimum(),min(distribution_tec_T->GetMinimum(),
									 distribution_tec_t->GetMinimum()))))*0.5,		 
				  max(distribution_tib->GetMaximum(),max(distribution_tob->GetMaximum(),
				  max(distribution_tid->GetMaximum(),max(distribution_tec_T->GetMaximum(),
									 distribution_tec_t->GetMaximum()))))*1.5);
  frame->GetXaxis()->SetNdivisions(511);
  frame->Draw();

  CMS_lumi(canvas,"",false,false);  
    
  distribution_tib->SetMarkerColor(kBlack);
  distribution_tib->SetLineColor(kBlack);
  ((TF1*)distribution_tib->GetListOfFunctions()->At(0))->SetLineColor(kBlack);
  distribution_tib->SetMarkerSize(1);
  distribution_tib->SetMarkerStyle(20);
  distribution_tob->SetMarkerColor(TColor::GetColor("#CF3721"));
  distribution_tob->SetLineColor(TColor::GetColor("#CF3721"));
  ((TF1*)distribution_tob->GetListOfFunctions()->At(0))->SetLineColor(TColor::GetColor("#CF3721"));
  distribution_tob->SetMarkerSize(1);
  distribution_tob->SetMarkerStyle(20);
  distribution_tid->SetMarkerColor(kBlue);
  distribution_tid->SetLineColor(kBlue);
  ((TF1*)distribution_tid->GetListOfFunctions()->At(0))->SetLineColor(kBlue);
  distribution_tid->SetMarkerSize(1);
  distribution_tid->SetMarkerStyle(20);
  distribution_tec_T->SetLineColor(TColor::GetColor("#4D975D"));
  distribution_tec_T->SetMarkerColor(TColor::GetColor("#4D975D"));
  ((TF1*)distribution_tec_T->GetListOfFunctions()->At(0))->SetLineColor(TColor::GetColor("#4D975D"));
  distribution_tec_T->SetMarkerSize(1);
  distribution_tec_T->SetMarkerStyle(20);
  distribution_tec_t->SetMarkerColor(TColor::GetColor("#FAAF08"));
  distribution_tec_t->SetLineColor(TColor::GetColor("#FAAF08"));
  ((TF1*)distribution_tec_t->GetListOfFunctions()->At(0))->SetLineColor(TColor::GetColor("#FAAF08"));
  distribution_tec_t->SetMarkerSize(1);
  distribution_tec_t->SetMarkerStyle(20);

  distribution_tib->Draw("EPsame");
  distribution_tob->Draw("EPsame");
  distribution_tid->Draw("EPsame");
  distribution_tec_T->Draw("EPsame");
  distribution_tec_t->Draw("EPsame");

  TLegend leg (0.655,0.64,0.924,0.90);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry((TObject*)(0),"2017 Data","");
  leg.AddEntry(distribution_tib,"TIB","EP");
  leg.AddEntry(distribution_tob,"TOB","EP");
  leg.AddEntry(distribution_tid,"TID","EP");
  leg.AddEntry(distribution_tec_T,"TEC thick","EP");
  leg.AddEntry(distribution_tec_t,"TEC thin","EP");
  leg.Draw("same");

  canvas->SaveAs((outputDIR+"/clusterCharge_vs_delay_perPartition.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/clusterCharge_vs_delay_perPartition.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/clusterCharge_vs_delay_perPartition.root").c_str(),"root");
  
}
