#include "../CMS_lumi.h"

void makeDelayVsLayer(string inputFileName, string outputDIR){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());

  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");

  vector<TH1F*> distribution_layers;
  distribution_layers.push_back((TH1F*) inputFile->Get("TIB_layer_1_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TIB_layer_2_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TIB_layer_3_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TIB_layer_4_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TOB_layer_1_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TOB_layer_2_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TOB_layer_3_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TOB_layer_4_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TOB_layer_5_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TOB_layer_6_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TID_layer_1_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TID_layer_2_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TID_layer_3_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECPT_layer_1_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECPT_layer_2_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECPT_layer_3_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECPT_layer_4_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECPT_layer_5_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECPT_layer_6_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECPT_layer_7_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECPT_layer_8_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECPT_layer_9_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECPt_layer_1_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECPt_layer_2_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECPt_layer_3_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECPt_layer_4_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECPt_layer_5_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECPt_layer_6_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECPt_layer_7_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECPt_layer_8_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECPt_layer_9_mean"));

  distribution_layers.push_back((TH1F*) inputFile->Get("TECMT_layer_1_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECMT_layer_2_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECMT_layer_3_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECMT_layer_4_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECMT_layer_5_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECMT_layer_6_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECMT_layer_7_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECMT_layer_8_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECMT_layer_9_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECMt_layer_1_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECMt_layer_2_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECMt_layer_3_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECMt_layer_4_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECMt_layer_5_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECMt_layer_6_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECMt_layer_7_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECMt_layer_8_mean"));
  distribution_layers.push_back((TH1F*) inputFile->Get("TECMt_layer_9_mean"));

  TH1F* totalhisto = new TH1F("totalHisto","",distribution_layers.size(),0,distribution_layers.size()+1);
  for(unsigned int ihisto = 0; ihisto < distribution_layers.size(); ihisto++){
    totalhisto->SetBinContent(ihisto+1,((TF1*)(distribution_layers.at(ihisto)->GetListOfFunctions()->At(0)))->GetParameter(1));
    totalhisto->SetBinError(ihisto+1,((TF1*)(distribution_layers.at(ihisto)->GetListOfFunctions()->At(0)))->GetParError(1));
  }

  // plotting results
  TCanvas* canvas = new TCanvas("canvas","canvas",900,650);
  canvas->cd();

  gPad->SetBottomMargin(0.11);
  gPad->SetLeftMargin(0.1);

  TH1F* frame = (TH1F*) totalhisto->Clone("frame");
  frame->Reset();
  frame->GetXaxis()->SetTitle("Strip Tracker Layer");
  frame->GetXaxis()->SetTitleOffset(0.75);
  frame->GetYaxis()->SetTitle("Delay (ns)");
  frame->GetYaxis()->SetTitleOffset(0.8);
  frame->GetYaxis()->SetRangeUser(-3,3);
  frame->GetXaxis()->SetLabelSize(0);
  frame->Draw();

  TF1* line = new TF1("line","0",-0.5,distribution_layers.size()+2);
  line->SetLineColor(kRed);
  line->SetLineWidth(2);

  totalhisto->SetFillColor(TColor::GetColor("#4D975D"));
  totalhisto->SetLineColor(kBlack);
  totalhisto->SetMarkerColor(kBlack);
  totalhisto->SetMarkerStyle(20);
  totalhisto->SetMarkerSize(1);
  totalhisto->Draw("hist same");
  line->Draw("Lsame");
  totalhisto->Draw("EPsame");

  TLine* line_tib = new TLine(totalhisto->GetBinLowEdge(5),-2.5,totalhisto->GetBinLowEdge(5),2);
  line_tib->SetLineColor(kBlack);
  line_tib->SetLineStyle(7);
  line_tib->SetLineWidth(2);
  line_tib->Draw("same");
  TLine* line_tob = new TLine(totalhisto->GetBinLowEdge(11),-2.5,totalhisto->GetBinLowEdge(11),2);
  line_tob->SetLineColor(kBlack);
  line_tob->SetLineStyle(7);
  line_tob->SetLineWidth(2);
  line_tob->Draw("same");
  TLine* line_tid = new TLine(totalhisto->GetBinLowEdge(14),-2.5,totalhisto->GetBinLowEdge(14),2);
  line_tid->SetLineColor(kBlack);
  line_tid->SetLineStyle(7);
  line_tid->SetLineWidth(2);
  line_tid->Draw("same");
  TLine* line_tecp_T = new TLine(totalhisto->GetBinLowEdge(23),-2.5,totalhisto->GetBinLowEdge(23),2);
  line_tecp_T->SetLineColor(kBlack);
  line_tecp_T->SetLineStyle(7);
  line_tecp_T->SetLineWidth(2);
  line_tecp_T->Draw("same");
  TLine* line_tecp_t = new TLine(totalhisto->GetBinLowEdge(32),-2.5,totalhisto->GetBinLowEdge(32),2);
  line_tecp_t->SetLineColor(kBlack);
  line_tecp_t->SetLineStyle(7);
  line_tecp_t->SetLineWidth(2);
  line_tecp_t->Draw("same");
  TLine* line_tecm_T = new TLine(totalhisto->GetBinLowEdge(41),-2.5,totalhisto->GetBinLowEdge(41),2);
  line_tecm_T->SetLineColor(kBlack);
  line_tecm_T->SetLineStyle(7);
  line_tecm_T->SetLineWidth(2);
  line_tecm_T->Draw("same");

  TLatex *text_tib = new TLatex();
  text_tib->SetTextFont(62);
  text_tib->SetNDC();
  text_tib->SetTextSize(0.04);
  text_tib->DrawLatex(0.12,0.27,"TIB");

  TLatex *text_tob = new TLatex();
  text_tob->SetNDC();
  text_tob->SetTextSize(0.04);
  text_tob->SetTextFont(62);
  text_tob->DrawLatex(0.195,0.27,"TOB");

  TLatex *text_tid = new TLatex();
  text_tid->SetNDC();
  text_tid->SetTextSize(0.04);
  text_tid->SetTextFont(62);
  text_tid->DrawLatex(0.278,0.27,"TID");

  TLatex *text_tecp_T = new TLatex();
  text_tecp_T->SetNDC();
  text_tecp_T->SetTextSize(0.04);
  text_tecp_T->SetTextFont(62);
  text_tecp_T->DrawLatex(0.34,0.27,"TEC+ thick");

  TLatex *text_tecp_t = new TLatex();
  text_tecp_t->SetNDC();
  text_tecp_t->SetTextSize(0.04);
  text_tecp_t->SetTextFont(62);
  text_tecp_t->DrawLatex(0.50,0.27,"TEC+ thin");

  TLatex *text_tecm_T = new TLatex();
  text_tecm_T->SetNDC();
  text_tecm_T->SetTextSize(0.04);
  text_tecm_T->SetTextFont(62);
  text_tecm_T->DrawLatex(0.655,0.27,"TEC- thick");

  TLatex *text_tecm_t = new TLatex();
  text_tecm_t->SetNDC();
  text_tecm_t->SetTextSize(0.04);
  text_tecm_t->SetTextFont(62);
  text_tecm_t->DrawLatex(0.81,0.27,"TEC- thin");

  TLatex* cmsBanner = new TLatex();
  cmsBanner->SetTextSize(0.05);
  cmsBanner->SetNDC();
  cmsBanner->SetTextFont(62);
  cmsBanner->SetTextAlign(11);
  cmsBanner->DrawLatex(0.16,0.85,"CMS");

  cmsBanner->SetTextSize(0.045);
  cmsBanner->SetTextFont(52);
  cmsBanner->SetTextAlign(11);
  cmsBanner->DrawLatex(0.245,0.85,"Preliminary 2017");  

  cmsBanner->SetTextSize(0.045);
  cmsBanner->SetTextFont(42);
  cmsBanner->SetTextAlign(11);
  cmsBanner->DrawLatex(0.86,0.94,"13 TeV");

  canvas->RedrawAxis("sameaxis");

  canvas->SaveAs((outputDIR+"/clusterCharge_vs_delay_perLayer.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/clusterCharge_vs_delay_perLayer.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/clusterCharge_vs_delay_perLayer.root").c_str(),"root");
  
}
