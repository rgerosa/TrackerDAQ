#include "../CMS_lumi.h"

int apvInModule(const uint16_t & lldChan, const uint16_t & apvId){

  if(lldChan == 1 and apvId == 1) return 1;
  else if(lldChan == 1 and apvId == 2) return 2;
  else if(lldChan == 2 and apvId == 1) return 3;
  else if(lldChan == 2 and apvId == 2) return 4;
  else if(lldChan == 3 and apvId == 1) return 5;
  else if(lldChan == 3 and apvId == 2) return 6;
  return 0;

}

///////////////////////////////
void plotBadStripsCumulative(string inputFileName, string outputPlotDIR){

  setTDRStyle();
  system(("mkdir -p "+outputPlotDIR).c_str());
  gROOT->SetBatch(kTRUE);

  gStyle->SetLabelSize(0.04, "XY");
  gStyle->SetLabelSize(0.03, "Z");
  gStyle->SetPadGridX(1);
  
  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  TTree* tree = (TTree*) inputFile->Get("badStripTree");

  TTreeReader reader(tree);
  TTreeReaderValue<uint32_t> detid    (reader,"detid");
  TTreeReaderValue<uint16_t> fecCrate (reader,"fecCrate");
  TTreeReaderValue<uint16_t> fecSlot  (reader,"fecSlot");
  TTreeReaderValue<uint16_t> fecRing  (reader,"fecRing");
  TTreeReaderValue<uint16_t> ccuAdd   (reader,"ccuAdd");
  TTreeReaderValue<uint16_t> ccuChan  (reader,"ccuChan");
  TTreeReaderValue<uint16_t> lldChannel  (reader,"lldChannel");
  TTreeReaderValue<uint16_t> fedId  (reader,"fedId");
  TTreeReaderValue<uint16_t> fedCh  (reader,"fedCh");
  TTreeReaderValue<uint16_t> apvId  (reader,"apvId");
  TTreeReaderValue<uint16_t> stripId  (reader,"stripId");
  TTreeReaderValue<uint16_t> badStrip  (reader,"badStrip");
  
  TH2F* badStrip_perAPV = new TH2F("badStrip_perAPV","",128,1,129,3,0,3); // position in the APV
  TH2F* badStrip_perLLDChannel = new TH2F("badStrip_perLLDChannel","",3,1,4,3,0,3); // fiber within a module
  TH2F* badStrip_perModule = new TH2F("badStrip_perModule","",6,1,7,3,0,3); // APV chip within a module
  badStrip_perAPV->Sumw2();
  badStrip_perLLDChannel->Sumw2();
  badStrip_perModule->Sumw2();

  long int numberOfBadStrips = 0;

  /////-------------------
  while(reader.Next()){

    if(*badStrip  == 0) continue;
    numberOfBadStrips++;
    
    badStrip_perAPV->SetBinContent(badStrip_perAPV->GetXaxis()->FindBin(*stripId),2,
				   badStrip_perAPV->GetBinContent(badStrip_perAPV->GetXaxis()->FindBin(*stripId),2)+1);

    badStrip_perLLDChannel->SetBinContent(badStrip_perLLDChannel->GetXaxis()->FindBin(*lldChannel),2,
					  badStrip_perLLDChannel->GetBinContent(badStrip_perLLDChannel->GetXaxis()->FindBin(*lldChannel),2)+1);    

    badStrip_perModule->SetBinContent(apvInModule(*lldChannel,*apvId),2,
				      badStrip_perModule->GetBinContent(apvInModule(*lldChannel,*apvId),2)+1);
    
  }
  
  badStrip_perAPV->Scale(1./numberOfBadStrips);
  badStrip_perLLDChannel->Scale(1./numberOfBadStrips);
  badStrip_perModule->Scale(1./numberOfBadStrips);

  for(int iBinX = 0; iBinX < badStrip_perAPV->GetNbinsX(); iBinX++)
    badStrip_perAPV->GetXaxis()->SetBinLabel(iBinX+1,Form("%d",iBinX+1));

  for(int iBinX = 0; iBinX < badStrip_perLLDChannel->GetNbinsX(); iBinX++)
    badStrip_perLLDChannel->GetXaxis()->SetBinLabel(iBinX+1,Form("Fiber %d",iBinX+1));

  for(int iBinX = 0; iBinX < badStrip_perModule->GetNbinsX(); iBinX++)
    badStrip_perModule->GetXaxis()->SetBinLabel(iBinX+1,Form("APV %d",iBinX+1));
  
  badStrip_perAPV->GetXaxis()->SetTitle("");
  badStrip_perLLDChannel->GetXaxis()->SetTitle("");
  badStrip_perModule->GetXaxis()->SetTitle("");
  badStrip_perAPV->GetYaxis()->SetTitle("");
  badStrip_perLLDChannel->GetYaxis()->SetTitle("");
  badStrip_perModule->GetYaxis()->SetTitle("");
  badStrip_perAPV->GetYaxis()->SetLabelSize(0);
  badStrip_perLLDChannel->GetYaxis()->SetLabelSize(0);
  badStrip_perModule->GetYaxis()->SetLabelSize(0);

  badStrip_perAPV->GetZaxis()->SetTitle("fraction of bad strips");
  badStrip_perLLDChannel->GetZaxis()->SetTitle("fraction of bad strips");
  badStrip_perModule->GetZaxis()->SetTitle("fraction of bad strips");
  

  TCanvas* canvas = new TCanvas("canvas","canvas",2000,1000);
  canvas->SetLeftMargin(0.03);
  canvas->SetRightMargin(0.12);
  badStrip_perAPV->Draw("colz");
  CMS_lumi(canvas,"",false,false,false);

  canvas->SaveAs((outputPlotDIR+"/badStrip_perAPV.png").c_str(),"png");
  canvas->SaveAs((outputPlotDIR+"/badStrip_perAPV.pdf").c_str(),"pdf");

  gStyle->SetLabelSize(0.06, "XY");

  TCanvas* canvas2 = new TCanvas("canvas2","canvas2",800,600);
  canvas2->SetLeftMargin(0.03);
  canvas2->SetRightMargin(0.15);
  badStrip_perLLDChannel->Draw("colz");
  CMS_lumi(canvas2,"",false,false,false);

  canvas2->SaveAs((outputPlotDIR+"/badStrip_perFiber.png").c_str(),"png");
  canvas2->SaveAs((outputPlotDIR+"/badStrip_perFiber.pdf").c_str(),"pdf");

  badStrip_perModule->Draw("colz");
  CMS_lumi(canvas2,"",false,false,false);

  canvas2->SaveAs((outputPlotDIR+"/badStrip_perModule.png").c_str(),"png");
  canvas2->SaveAs((outputPlotDIR+"/badStrip_perModule.pdf").c_str(),"pdf");

}
