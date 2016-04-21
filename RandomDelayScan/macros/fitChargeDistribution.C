#include <vector>
#include <iostream>
#include <TFile.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TEventList.h>
#include <TCut.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TTreeReader.h>

#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooDataHist.h"
#include "RooLandau.h"
#include "RooGaussian.h"
#include "RooFFTConvPdf.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"

#include "TkPulseShape.h"

const float xMin = 20;
const float xMax = 250;
const int   nBin = 40;
const bool  verbosity = false;

class chargeShape{

public:
  chargeShape(){};
  ~chargeShape(){};
  
  RooRealVar*    charge;
  RooDataHist*   chargeDistribution;
  RooRealVar*    mean_landau;
  RooRealVar*    sigma_landau;
  RooLandau*     landau;
  RooRealVar*    mean_gauss;
  RooRealVar*    sigma_gauss;
  RooGaussian*   gauss;
  RooFFTConvPdf* landauXgauss;
  RooRealVar*    normalization;
  RooExtendPdf*  totalPdf;
  RooFitResult*  fitResult;

};
 
void fitChargeDistribution(const char* file0, string outputDirectory){

  system(("mkdir -p "+outputDirectory).c_str());

  gStyle->SetOptStat(0);
  gROOT->SetBatch(kTRUE);
  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;

  // open the file and prepare the cluster tree
  cout<<"########### fitChargeDistribution analysis ##############"<<endl;
  TFile *_file0 = TFile::Open(file0);
  TTree* clusters = (TTree*)_file0->Get("analysis/trackerDPG/clusters");

  clusters->SetAlias("subdetid","int((detid-0x10000000)/0x2000000)");
  clusters->SetAlias("barrellayer","int((detid%33554432)/0x4000)");
  clusters->SetAlias("TIDlayer","int((detid%33554432)/0x800)%4");
  clusters->SetAlias("TECPlayer","int((detid%33554432)/0x4000)-32");
  clusters->SetAlias("TECMlayer","int((detid%33554432)/0x4000)-16");
  clusters->SetAlias("R","sqrt(clglobalX**2+clglobalY**2+clglobalZ**2)");
  
  TFile* outputFile = new TFile((outputDirectory+"/outputCharge.root").c_str(),"RECREATE");
  outputFile->cd();

  // apply common preselection cuts on events, track and cluster quantities
  // make a index as a funcion of det id
  std::map<uint32_t,TH1F*> chargeDistributionMap;

  TTreeReader reader(clusters);
  TTreeReaderValue<uint32_t> detid  (reader,"detid");
  TTreeReaderValue<float> maxCharge (reader,"maxCharge");
  TTreeReaderValue<bool>  onTrack   (reader,"onTrack");
  TTreeReaderValue<float> angle     (reader,"angle");
  
  uint32_t iCluster = 0;
  cout<<"### loop on clusters tree: nClusters "<<clusters->GetEntries()<<endl;
  while(reader.Next()){
    cout.flush();
    if(iCluster % 1000000 == 0) cout<<"\r"<<"iCluster "<<100*double(iCluster)/clusters->GetEntries()<<" % ";
    iCluster++;
    if(chargeDistributionMap[*detid] == 0)
      chargeDistributionMap[*detid] = new TH1F(Form("chargeDistribution_detid_%d",*detid),"",nBin,xMin,xMax);
    if(not *onTrack or *angle < 0 or *maxCharge >= 254) continue;
    chargeDistributionMap[*detid]->Fill(*maxCharge);    
  }
  cout<<endl;
  cout<<"#### Map size = "<<chargeDistributionMap.size()<<endl;
  uint32_t detId = 0;
  uint32_t invalidFits = 0;
  
  std::map<uint32_t,chargeShape*> fitMap;
  chargeShape* chargeDist = NULL;
  cout<<"#### Start charge shape fits: nFits = "<<chargeDistributionMap.size()<<endl;
  for(auto ihist : chargeDistributionMap){
    detId++;
    if(detId % 50 == 0) cout<<"\r"<<"iFit "<<100*double(detId)/chargeDistributionMap.size()<<" % ";
    
    chargeDist = new chargeShape();
    chargeDist->charge = new RooRealVar(Form("charge_detid_%d",ihist.first),"",xMin,xMax);
    RooArgList vars(*chargeDist->charge);
    chargeDist->chargeDistribution = new RooDataHist(ihist.second->GetName(),"",vars,ihist.second);
    chargeDist->mean_landau  = new RooRealVar(Form("mean_landau_detid_%d",ihist.first),"",ihist.second->GetMean(),xMin,xMax);
    chargeDist->sigma_landau = new RooRealVar(Form("sigma_landau_detid_%d",ihist.first),"",ihist.second->GetRMS(),1.,100);
    chargeDist->landau       =  new RooLandau(Form("landau_detid_%d",ihist.first),"",*chargeDist->charge,*chargeDist->mean_landau,*chargeDist->sigma_landau) ;
    chargeDist->mean_gauss   = new RooRealVar(Form("mean_gauss_detid_%d",ihist.first),"",0.,-100,100);
    chargeDist->sigma_gauss  = new RooRealVar(Form("sigma_gauss_detid_%d",ihist.first),"",10,1.,50);
    chargeDist->gauss        = new RooGaussian(Form("gauss_detid_%d",ihist.first),"",*chargeDist->charge,*chargeDist->mean_gauss,*chargeDist->sigma_gauss);
    chargeDist->landauXgauss = new RooFFTConvPdf(Form("landauXgauss_detid_%d",ihist.first),"",*chargeDist->charge,*chargeDist->landau,*chargeDist->gauss);
    chargeDist->normalization = new RooRealVar(Form("normalization_detid_%d",ihist.first),"",ihist.second->Integral(),ihist.second->Integral()/10,ihist.second->Integral()*10);
    chargeDist->totalPdf     = new RooExtendPdf(Form("totalPdf_detid_%d",ihist.first),"",*chargeDist->landauXgauss,*chargeDist->normalization);  
    chargeDist->fitResult = chargeDist->totalPdf->fitTo(*chargeDist->chargeDistribution,RooFit::Range(xMin,xMax),RooFit::Extended(kTRUE),RooFit::NumCPU(4),RooFit::Save(kTRUE));
    fitMap[ihist.first]  = chargeDist;
    if(verbosity){
      chargeDist->fitResult->Print();
      cout<<"Fit status "<<chargeDist->fitResult->status()<<" covQual "<<chargeDist->fitResult->covQual()<<" numInvalidNLL "<<chargeDist->fitResult->numInvalidNLL()<<" edm "<<chargeDist->fitResult->edm()<<" minNll "<<chargeDist->fitResult->minNll()<<endl;
    }

    if(chargeDist->fitResult->status() != 0){
      invalidFits++;
      cerr<<"Problem with fit for detId "<<ihist.first<<" status not zero ! : "<<chargeDist->fitResult->status()<<endl;
      chargeDist->fitResult->Print();
    }
    if(chargeDist->fitResult->covQual() != 3 and chargeDist->fitResult->status() == 0){
      invalidFits++;
      cerr<<"Problem with fit for detId "<<ihist.first<<" covariance matrix not status 3 ! .. status = "<<chargeDist->fitResult->covQual()<<endl;
      chargeDist->fitResult->Print();
    }
  }

  cout<<"#### end of fit stage: nFits "<<detId<<" Invalid Fits "<<invalidFits<<" : "<<100*double(invalidFits)/detId<<" % "<<endl;
  outputFile->cd();
  outputFile->mkdir("RooPlots");
  outputFile->cd("RooPlots");
  RooPlot* frame = NULL;
  cout<<"#### start storing plots "<<endl;
  TCanvas* canvas = new TCanvas("canvas","",600,600);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->cd();
  for(auto ifit : fitMap){    
    
    RooPlot* frame = ifit.second->charge->frame();
    frame->SetTitle("");
    frame->GetYaxis()->SetTitle("Events");
    frame->GetXaxis()->SetTitle("Cluster maxCharge");
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetLabelSize(0.038);
    frame->GetXaxis()->SetLabelSize(0.038);
    ifit.second->chargeDistribution->plotOn(frame,RooFit::MarkerColor(kBlack),RooFit::MarkerSize(1),RooFit::MarkerStyle(20),RooFit::LineColor(kBlack),RooFit::LineWidth(2),RooFit::DrawOption("EP"));
    
    RooArgList list = ifit.second->fitResult->floatParsFinal();
    ifit.second->totalPdf->paramOn(frame,RooFit::Parameters(list),RooFit::Layout(0.58,0.85,0.60));    
    ifit.second->totalPdf->plotOn(frame,RooFit::LineColor(kRed));        
    float chi2 = frame->chiSquare(list.getSize());
    TPaveLabel *t1 = new TPaveLabel(0.7,0.8,0.9,0.92, Form("#chi^{2}/ndf = %f",chi2),"BRNDC"); 
    t1->SetFillColor(0);
    t1->SetFillStyle(0);
    t1->SetBorderSize(0);
    t1->Draw("same");
    frame->addObject(t1) ; 
    frame->Write(ifit.second->charge->GetName(),TObject::kOverwrite);    
    //    for canvas
  }
  outputFile->cd();
  for(auto ifit : fitMap){    
    
    RooPlot* frame2 = ifit.second->charge->frame();
    frame2->SetTitle("");
    frame2->GetYaxis()->SetTitle("Events");
    frame2->GetXaxis()->SetTitle("Cluster maxCharge");
    frame2->GetYaxis()->SetTitleSize(0.045);
    frame2->GetXaxis()->SetTitleSize(0.045);
    frame2->GetYaxis()->SetLabelSize(0.038);
    frame2->GetXaxis()->SetLabelSize(0.038);
    ifit.second->chargeDistribution->plotOn(frame2,RooFit::MarkerColor(kBlack),RooFit::MarkerSize(1),RooFit::MarkerStyle(20),RooFit::LineColor(kBlack),RooFit::LineWidth(2),RooFit::DrawOption("EP"));
    RooArgList list = ifit.second->fitResult->floatParsFinal();
    ifit.second->totalPdf->paramOn(frame2,RooFit::Parameters(list),RooFit::Layout(0.58,0.85,0.70));    
    frame2->getAttText()->SetTextSize(0.015);
    ifit.second->totalPdf->plotOn(frame2,RooFit::LineColor(kRed));        
    float chi2 = frame2->chiSquare(list.getSize());
    TPaveLabel *t1 = new TPaveLabel(0.7,0.8,0.88,0.94, Form("#chi^{2}/ndf = %f",chi2),"BRNDC"); 
    t1->SetFillColor(0);
    t1->SetFillStyle(0);
    t1->SetBorderSize(0);
    t1->Draw("same");
    frame2->addObject(t1) ; 
    frame2->Draw();    
    canvas->SaveAs((outputDirectory+"/"+string(ifit.second->charge->GetName())+".png").c_str(),"png");
  }
  
  cout<<"#### start output text file "<<endl;
  ofstream dump((outputDirectory+"/dumpMapPeak.txt").c_str());
  for(auto ifit : fitMap){
    TF1* funz = ifit.second->totalPdf->asTF( RooArgList(*ifit.second->charge));
    dump<<ifit.first<<"   "<<funz->GetMaximumX(xMin,xMax)<<"\n";
  }
  dump.close();

  outputFile->Close();


}

