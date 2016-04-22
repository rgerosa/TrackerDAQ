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
#include "RooPlot.h"
#include "TPaveLabel.h"

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
 
void fitChargeDistribution(const char* file0, string outputDirectory, bool saveCanvas = false){

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
  std::map<uint32_t,std::shared_ptr<TH1F> > chargeDistributionMap;

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
    if(chargeDistributionMap[*detid].get() == 0 or chargeDistributionMap[*detid].get() == NULL)
      chargeDistributionMap[*detid] = std::shared_ptr<TH1F>(new TH1F(Form("chargeDistribution_detid_%d",*detid),"",nBin,xMin,xMax));
    if(not *onTrack or *angle < 0 or *maxCharge >= 254) continue;
    chargeDistributionMap[*detid]->Fill(*maxCharge);    
  }
  cout<<endl;
  cout<<"#### Map size = "<<chargeDistributionMap.size()<<endl;
  uint32_t detId = 0;
  uint32_t invalidFits = 0;
  
  cout<<"#### Start charge shape fits: nFits = "<<chargeDistributionMap.size()<<endl;
  outputFile->cd();
  outputFile->mkdir("RooPlots");
  outputFile->cd("RooPlots");

  map<string,string> mapCharge;

  for(auto ihist : chargeDistributionMap){
    cout.flush();
    if(detId % 100 == 0) cout<<"\r"<<"iFit "<<100*double(detId)/double(chargeDistributionMap.size())<<" % ";    
    detId++;

    RooRealVar  charge(Form("charge_detid_%d",ihist.first),"",xMin,xMax);
    RooArgList  vars(charge);
    RooDataHist chargeDistribution(ihist.second->GetName(),"",vars,ihist.second.get()); 
    RooRealVar  mean_landau(Form("mean_landau_detid_%d",ihist.first),"",ihist.second->GetMean(),xMin,xMax);
    RooRealVar  sigma_landau(Form("sigma_landau_detid_%d",ihist.first),"",ihist.second->GetRMS(),1.,100);
    RooLandau   landau(Form("landau_detid_%d",ihist.first),"",charge,mean_landau,sigma_landau) ;
    RooRealVar  mean_gauss (Form("mean_gauss_detid_%d",ihist.first),"",0.,-150,150);
    RooRealVar  sigma_gauss(Form("sigma_gauss_detid_%d",ihist.first),"",10,1.,50);
    RooGaussian gauss(Form("gauss_detid_%d",ihist.first),"",charge,mean_gauss,sigma_gauss);
    RooFFTConvPdf landauXgauss(Form("landauXgauss_detid_%d",ihist.first),"",charge,landau,gauss);
    RooRealVar normalization(Form("normalization_detid_%d",ihist.first),"",ihist.second->Integral(),ihist.second->Integral()/5,ihist.second->Integral()*5);
    RooExtendPdf totalPdf   (Form("totalPdf_detid_%d",ihist.first),"",landauXgauss,normalization);  

    RooFitResult* fitResult = totalPdf.fitTo(chargeDistribution,RooFit::Range(xMin,xMax),RooFit::Extended(kTRUE),RooFit::Save(kTRUE));
    if(verbosity){
      fitResult->Print();
      cout<<"Fit status "<<fitResult->status()<<" covQual "<<fitResult->covQual()<<" numInvalidNLL "<<fitResult->numInvalidNLL()<<" edm "<<fitResult->edm()<<" minNll "<<fitResult->minNll()<<endl;
    }

    if(fitResult->status() != 0){
      invalidFits++;
      cerr<<"Problem with fit for detId "<<ihist.first<<" status not zero ! : "<<fitResult->status()<<endl;
      fitResult->Print();
    }
    if(fitResult->covQual() <= 1 and fitResult->status() == 0){
      invalidFits++;
      cerr<<"Problem with fit for detId "<<ihist.first<<" covariance matrix not status 3 ! .. status = "<<fitResult->covQual()<<endl;
      fitResult->Print();
    }

    // plot to store in a root file
    RooPlot* frame = charge.frame();
    frame->SetTitle("");
    frame->GetYaxis()->SetTitle("Events");
    frame->GetXaxis()->SetTitle("Cluster maxCharge");
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetLabelSize(0.038);
    frame->GetXaxis()->SetLabelSize(0.038);
    chargeDistribution.plotOn(frame,RooFit::MarkerColor(kBlack),RooFit::MarkerSize(1),RooFit::MarkerStyle(20),RooFit::LineColor(kBlack),RooFit::LineWidth(2),RooFit::DrawOption("EP"));
    
    RooArgList parlist = fitResult->floatParsFinal();
    totalPdf.paramOn(frame,RooFit::Parameters(parlist),RooFit::Layout(0.58,0.85,0.60));    
    totalPdf.plotOn(frame,RooFit::LineColor(kRed));        
    float chi2 = frame->chiSquare(parlist.getSize());
    TPaveLabel *t1 = new TPaveLabel(0.7,0.8,0.9,0.92, Form("#chi^{2}/ndf = %f",chi2),"BRNDC"); 
    t1->SetFillColor(0);
    t1->SetFillStyle(0);
    t1->SetBorderSize(0);
    t1->Draw("same");
    frame->addObject(t1) ; 
    frame->Write(charge.GetName(),TObject::kOverwrite);    

    // plot to store in canvas
    if(saveCanvas){

      TCanvas* canvas = new TCanvas("canvas","",600,600);
      canvas->SetTickx();
      canvas->SetTicky();

      RooPlot* frame2 = charge.frame();
      frame2->SetTitle("");
      frame2->GetYaxis()->SetTitle("Events");
      frame2->GetXaxis()->SetTitle("Cluster maxCharge");
      frame2->GetYaxis()->SetTitleSize(0.045);
      frame2->GetXaxis()->SetTitleSize(0.045);
      frame2->GetYaxis()->SetLabelSize(0.038);
      frame2->GetXaxis()->SetLabelSize(0.038);
      chargeDistribution.plotOn(frame2,RooFit::MarkerColor(kBlack),RooFit::MarkerSize(1),RooFit::MarkerStyle(20),RooFit::LineColor(kBlack),RooFit::LineWidth(2),RooFit::DrawOption("EP"));
      //      totalPdf.paramOn(frame2,RooFit::Parameters(parlist),RooFit::Layout(0.58,0.85,0.70));    
      frame2->getAttText()->SetTextSize(0.015);
      totalPdf.plotOn(frame2,RooFit::LineColor(kRed));        
      TPaveLabel *t2 = new TPaveLabel(0.7,0.8,0.88,0.94, Form("#chi^{2}/ndf = %f",chi2),"BRNDC"); 
      t2->SetFillColor(0);
      t2->SetFillStyle(0);
      t2->SetBorderSize(0);
      t2->Draw("same");
      frame2->addObject(t2) ; 
      frame2->Draw();    
      canvas->SaveAs((outputDirectory+"/"+string(charge.GetName())+".png").c_str(),"png");
      if(t2) delete t2;
      if(canvas) delete canvas;
    }
    
    TF1* funz = totalPdf.asTF(vars);
    mapCharge[to_string(ihist.first)] = to_string(funz->GetMaximumX(xMin,xMax));
    if(funz) delete funz;
    if(t1) delete t1;
    if(fitResult) delete fitResult;
    
  }

  cout<<"#### end of fit stage: nFits "<<detId<<" Invalid Fits "<<invalidFits<<" : "<<100*double(invalidFits)/detId<<" % "<<endl;
  ofstream dump((outputDirectory+"/dumpMapPeak.txt").c_str());
  for(auto imap : mapCharge)
    dump<<imap.first<<"   "<<imap.second<<"\n";
  dump.close();
  outputFile->Close();
  
}

