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


void setLimitsAndBinning(string & observable, float & xMin, float & xMax, int & nBin){

  if(observable == "maxCharge"){
    xMin = 20;
    xMax = 250;
    nBin = 40;
  }
  else if(observable == "clSignalOverNoise" or observable == "clCorrectedSignalOverNoise"){
    xMin = 0;
    xMax = 80;
    nBin = 40;
  }
  else{
    xMin = 20;
    xMax = 250;
    nBin = 40;
  }
}
 
void fitChargeDistribution(string file0, 
			   string file1,
			   string outputDirectory, 
			   string observable = "maxCharge", 
			   float  delayMin = 0,
			   float  delayMax = 10,
			   bool   saveCanvas = false){
  
  system(("mkdir -p "+outputDirectory).c_str());

  gStyle->SetOptStat(0);
  gROOT->SetBatch(kTRUE);
  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;

  // open the file and prepare the cluster tree
  cout<<"########### fitChargeDistribution analysis ##############"<<endl;
  std::unique_ptr<TFile> _file0 (TFile::Open(file0.c_str()));
  std::unique_ptr<TFile> _file1 (TFile::Open(file1.c_str()));
  std::unique_ptr<TTree> clusters ((TTree*)_file0->Get("analysis/trackerDPG/clusters"));
  clusters->AddFriend((TTree*)_file1->Get("analysis/trackerDPG/readoutMap"));
  
  // apply common preselection cuts on events, track and cluster quantities
  // make a index as a funcion of det id
  std::map<uint32_t,std::shared_ptr<TH1F> > chargeDistributionMap;

  uint32_t detid;
  float    obs, maxCharge, angle, delay;
  bool     onTrack;

  // set the event selection
  std::auto_ptr<TEventList> selectedEvents (new TEventList("selectedEvents","selectedEvents"));
  clusters->Draw(">> selectedEvents",Form("onTrack && angle > 0 && maxCharge < 254 && abs(delay) >= %f && abs(delay) < %f",delayMin,delayMax));
  selectedEvents->Print();
  clusters->SetEventList(selectedEvents.get());
  
  // set only some branches
  clusters->SetBranchStatus("*",kFALSE);
  clusters->SetBranchStatus("detid",kTRUE);
  clusters->SetBranchStatus(observable.c_str(),kTRUE);
  clusters->SetBranchStatus("maxCharge",kTRUE);
  clusters->SetBranchStatus("onTrack",kTRUE);
  clusters->SetBranchStatus("angle",kTRUE);

  clusters->SetBranchAddress("detid",&detid);
  clusters->SetBranchAddress(observable.c_str(),&obs);
  clusters->SetBranchAddress("maxCharge",&maxCharge);
  clusters->SetBranchAddress("onTrack",&onTrack);
  clusters->SetBranchAddress("angle",&angle);
  
  cout<<"### loop on clusters tree: nClusters "<<clusters->GetEntries()<<" surviving selections "<<selectedEvents->GetN()<<endl;
  float xMin = 0., xMax = 0.;
  int   nBin = 0;
  // take ranges and binning asaf of the observable name
  setLimitsAndBinning(observable,xMin,xMax,nBin);
  
  // loop on the selected events to fill the histogra map per dei id
  long int iEntry = 0; 
  for(long int iCluster = 0; iCluster < selectedEvents->GetN(); iCluster++){
    iEntry = selectedEvents->GetEntry(iCluster);
    if(iEntry < 0) continue;
    clusters->GetEntry(iEntry);
    counter++;
    cout.flush();    
    if(iCluster % 1000000 == 0) cout<<"\r"<<"iCluster "<<100*double(iCluster)/selectedEvents->GetN()<<" % ";
    if(chargeDistributionMap[detid].get() == 0 or chargeDistributionMap[detid].get() == NULL)
      chargeDistributionMap[detid] = std::shared_ptr<TH1F>(new TH1F(Form("chargeDistribution_detid_%d",detid),"",nBin,xMin,xMax));
    chargeDistributionMap[detid]->SetDirectory(0); // detached from TFile
    chargeDistributionMap[detid]->Fill(obs);    
  }
  
  cout<<endl;
  cout<<"#### Map size = "<<chargeDistributionMap.size()<<endl;
  uint32_t detId = 0;
  uint32_t invalidFits = 0;
  cout<<"#### Start charge shape fits: nFits = "<<chargeDistributionMap.size()<<endl;
  std::unique_ptr<TFile> outputFile (new TFile((outputDirectory+"/output"+observable+".root").c_str(),"RECREATE"));

  map<string,string> mapCharge;
  std::auto_ptr<TCanvas> canvas(new TCanvas("canvas","",600,600));
  canvas->SetTickx();
  canvas->SetTicky();
  // start fit loop with landau convoluted with Gaussian
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
    RooRealVar   normalization(Form("normalization_detid_%d",ihist.first),"",ihist.second->Integral(),ihist.second->Integral()/5,ihist.second->Integral()*5);
    RooExtendPdf totalPdf (Form("totalPdf_detid_%d",ihist.first),"",landauXgauss,normalization);

    std::auto_ptr<RooFitResult> fitResult(totalPdf.fitTo(chargeDistribution,RooFit::Range(xMin,xMax),RooFit::Extended(kTRUE),RooFit::Save(kTRUE)));
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
    std::shared_ptr<RooPlot> frame (charge.frame());
    frame->SetName(Form("frame_%s",charge.GetName()));
    frame->SetTitle("");
    frame->GetYaxis()->SetTitle("Events");
    frame->GetXaxis()->SetTitle(("Cluster "+observable).c_str());
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetLabelSize(0.038);
    frame->GetXaxis()->SetLabelSize(0.038);
    chargeDistribution.plotOn(frame.get(),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(1),RooFit::MarkerStyle(20),RooFit::LineColor(kBlack),RooFit::LineWidth(2),RooFit::DrawOption("EP"));    
    RooArgList parlist = fitResult->floatParsFinal();
    totalPdf.paramOn(frame.get(),RooFit::Parameters(parlist),RooFit::Layout(0.58,0.85,0.60));    
    frame->getAttText()->SetTextSize(0.015);
    totalPdf.plotOn(frame.get(),RooFit::LineColor(kRed));        
    float chi2 = frame->chiSquare(parlist.getSize());
    std::auto_ptr<TPaveLabel> t1 (new TPaveLabel(0.7,0.8,0.88,0.94, Form("#chi^{2}/ndf = %f",chi2),"BRNDC")); 
    t1->SetFillColor(0);
    t1->SetFillStyle(0);
    t1->SetBorderSize(0);    
    canvas->cd();
    frame->Draw();
    t1->Draw("same");
    canvas->Write(frame->GetName()); 

    // plot to store in canvas
    if(saveCanvas){
      canvas->SaveAs((outputDirectory+"/"+string(charge.GetName())+".png").c_str(),"png");
      canvas->SaveAs((outputDirectory+"/"+string(charge.GetName())+".pdf").c_str(),"pdf");
    }
    
    RooArgList funzParam = RooArgList(*totalPdf.getParameters(chargeDistribution));
    std::auto_ptr<TF1> funz(totalPdf.asTF(charge,funzParam,charge));    // normalized to one
    /*
    // in case compare peak height: define integration region with the histogram bin width
    charge.setRange("max",
		    funz->GetMaximumX(xMin,xMax)-ihist.second->GetBinWidth(ihist.second->FindBin(funz->GetMaximumX(xMin,xMax)))/2,
		    funz->GetMaximumX(xMin,xMax)+ihist.second->GetBinWidth(ihist.second->FindBin(funz->GetMaximumX(xMin,xMax)))/2);
    // define an integral in that range
    RooRealVar* iVal = (RooRealVar*) totalPdf.createIntegral(charge, RooFit::NormSet(charge), RooFit::Range("max"));

    // multiply integral value time the normalization
    cout<<" hist max "<<ihist.second->GetMaximum()<<" funz "<<iVal->getVal()<<" "<<iVal->getVal()*normalization.getVal()<<endl;
    */
    mapCharge[to_string(ihist.first)] = to_string(funz->GetMaximumX(xMin,xMax));    
    
  }
    
  cout<<"#### end of fit stage: nFits "<<detId<<" Invalid Fits "<<invalidFits<<" : "<<100*double(invalidFits)/detId<<" % "<<endl;
  // make the detid:peak map to be plotted on  the tracker map through http://test-stripdbmonitor.web.cern.ch/test-stripdbmonitor/PrintTrackerMap/print_TrackerMap.php
  cout<<"#### Dump channel map"<<endl;
  ofstream dump((outputDirectory+"/dumpMapPeak"+observable+".txt").c_str());
  for(auto imap : mapCharge)
    dump<<imap.first<<"   "<<imap.second<<"\n";
  dump.close();
  
  cout<<"#### Close output Root file"<<endl;
  outputFile->Close("R");
  
}

