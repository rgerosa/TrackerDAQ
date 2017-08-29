#include <vector>
#include <iostream>

#include "TFile.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TEventList.h"
#include "TCut.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TTreeReader.h"

#include "../delayUtils.h"

const bool  verbosity = false;

// reduce the number of events by                                                                                                                                                                     
static int  reductionFactor = 1;

void makeChargeDistributionPerModule(string file0, 
				     string outputDirectory, 
				     string observable = "maxCharge", 			   
				     float  delayMin = 0,
				     float  delayMax = 10,
				     bool   applyCorrection = true,
				     bool   saveCanvas      = true,
				     bool   plotCanvas      = false){
  
  system(("mkdir -p "+outputDirectory).c_str());

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;

  // open the file and prepare the cluster tree
  cout<<"########### fitChargeDistribution analysis ##############"<<endl;
  std::unique_ptr<TFile> _file0 (TFile::Open(file0.c_str()));
  std::unique_ptr<TTree> clusters   ((TTree*)_file0->FindObjectAny("clusters"));
  std::unique_ptr<TTree> readoutMap ((TTree*)_file0->FindObjectAny("readoutMap"));
  
  // apply common preselection cuts on events, track and cluster quantities
  // make a index as a funcion of det id
  std::map<uint32_t,TH1F*> chargeDistributionMap;
  
  // set only some branches
  uint32_t detid;
  float    obs, clSignalOverNoise, clCorrectedSignalOverNoise;
  clusters->SetBranchStatus("*",kFALSE);
  clusters->SetBranchStatus("detid",kTRUE);
  clusters->SetBranchStatus("clCorrectedSignalOverNoise",kTRUE);
  clusters->SetBranchStatus("clSignalOverNoise",kTRUE);
  clusters->SetBranchStatus(observable.c_str(),kTRUE);
  clusters->SetBranchAddress("detid",&detid);
  clusters->SetBranchAddress("clSignalOverNoise",&clSignalOverNoise);
  clusters->SetBranchAddress("clCorrectedSignalOverNoise",&clCorrectedSignalOverNoise);
  clusters->SetBranchAddress(observable.c_str(),&obs);

  float delay;
  readoutMap->SetBranchStatus("*",kFALSE);
  readoutMap->SetBranchStatus("detid",kTRUE);
  readoutMap->SetBranchStatus("delay",kTRUE);
  readoutMap->SetBranchAddress("delay",&delay);

  float xMin = 0., xMax = 0.;
  int   nBin = 0;
  // take ranges and binning asaf of the observable name
  setLimitsAndBinning(observable,xMin,xMax,nBin);
  
  // loop on the selected events to fill the histogra map per dei id
  long int selectedEvents = 0; 
  for(long int iCluster = 0; iCluster < clusters->GetEntries()/reductionFactor; iCluster++){
    
    cout.flush();    
    if(iCluster % 100000 == 0) cout<<"\r"<<"iCluster "<<100*double(iCluster)/(clusters->GetEntries()/reductionFactor)<<" % ";

    // apply cluster selections
    clusters->GetEntry(iCluster);
    //take the related event in the readOutmap
    readoutMap->GetEntryWithIndex(detid);
    if(fabs(delay) < delayMin or fabs(delay) > delayMax) continue;
    
    selectedEvents++;
    // fill histograms
    if(chargeDistributionMap[detid] == 0 or chargeDistributionMap[detid] == NULL)
      chargeDistributionMap[detid] = new TH1F(Form("chargeDistribution_detid_%d",detid),"",nBin,xMin,xMax);
    chargeDistributionMap[detid]->SetDirectory(0); // detached from TFile
    if(not applyCorrection)
      chargeDistributionMap[detid]->Fill(obs);    
    else
      chargeDistributionMap[detid]->Fill(obs*clCorrectedSignalOverNoise/clSignalOverNoise);          
  }
    
  cout<<"#### Total events: "<<clusters->GetEntries()<<" selected events : "<<selectedEvents<<endl;
  cout<<"#### Map size = "<<chargeDistributionMap.size()<<endl;
  uint32_t detId = 0;
  uint32_t invalidFits = 0;
  cout<<"#### Start charge shape fits: nFits = "<<chargeDistributionMap.size()<<endl;
  std::unique_ptr<TFile> outputFile (new TFile((outputDirectory+"/output"+observable+".root").c_str(),"RECREATE"));

  map<string,string> mapPeakCharge;
  map<string,string> mapMeanCharge;
  TCanvas* canvas = new TCanvas("canvas","",625,650);
  canvas->SetTickx();
  canvas->SetTicky();
  
  // start fit loop with landau convoluted with Gaussian
  Double_t params [5]; // checked the order in the pdf
  float yMin   = 0, yMax = 0;
  int   nBinsY = 0;

  TPaveText leg (0.55,0.55,0.9,0.80,"NDC");
  leg.SetTextAlign(11);
  leg.SetTextFont(42);
  for(auto ihist : chargeDistributionMap){
    cout.flush();
    if(detId % 100 == 0) cout<<"\r"<<"iFit "<<100*double(detId)/double(chargeDistributionMap.size())<<" % ";    
    detId++;

    setLimitsAndBinning(observable,yMin,yMax,nBinsY);

    // make the observable
    RooRealVar* charge = new RooRealVar(Form("charge_detid_%d",ihist.first),"",xMin,xMax);
    RooArgList*  vars  = new RooArgList(*charge);
    // covert histogram in a RooDataHist
    RooDataHist* chargeDistribution = new RooDataHist(ihist.second->GetName(),"",*vars,ihist.second); 
    // Build a Landau PDF
    RooRealVar*  mean_landau  = new RooRealVar(Form("mean_landau_detid_%d",ihist.first),"",ihist.second->GetMean(),xMin,xMax);
    RooRealVar*  sigma_landau = new RooRealVar(Form("sigma_landau_detid_%d",ihist.first),"",ihist.second->GetRMS(),1.,100);
    RooLandau*   landau       = new RooLandau(Form("landau_detid_%d",ihist.first),"",*charge,*mean_landau,*sigma_landau) ;
    // Build a Gaussian PDF
    RooRealVar*  mean_gauss  = new RooRealVar(Form("mean_gauss_detid_%d",ihist.first),"",0.,-150,150);
    RooRealVar*  sigma_gauss = new RooRealVar(Form("sigma_gauss_detid_%d",ihist.first),"",10,1.,50);
    RooGaussian* gauss       = new RooGaussian(Form("gauss_detid_%d",ihist.first),"",*charge,*mean_gauss,*sigma_gauss);
    // Make a convolution
    RooFFTConvPdf* landauXgauss =  new RooFFTConvPdf(Form("landauXgauss_detid_%d",ihist.first),"",*charge,*landau,*gauss);
    // Make an extended PDF
    RooRealVar*   normalization = new RooRealVar(Form("normalization_detid_%d",ihist.first),"",ihist.second->Integral(),ihist.second->Integral()/5,ihist.second->Integral()*5);
    RooExtendPdf* totalPdf      = new RooExtendPdf(Form("totalPdf_detid_%d",ihist.first),"",*landauXgauss,*normalization);
    
    // Perform the FIT
    RooFitResult* fitResult = totalPdf->fitTo(*chargeDistribution,RooFit::Range(xMin,xMax),RooFit::Extended(kTRUE),RooFit::Save(kTRUE));
    if(verbosity)
      cout<<"Fit status "<<fitResult->status()<<" covQual "<<fitResult->covQual()<<" numInvalidNLL "<<fitResult->numInvalidNLL()<<" edm "<<fitResult->edm()<<" minNll "<<fitResult->minNll()<<endl;
    if(fitResult->status() != 0){
      invalidFits++;
      cerr<<"Problem with fit for detId "<<ihist.first<<" status not zero ! : "<<fitResult->status()<<endl;
      continue; // skip the channel
    }
    if(fitResult->covQual() <= 1 and fitResult->status() == 0){
      invalidFits++;
      cerr<<"Problem with fit for detId "<<ihist.first<<" covariance matrix not status 3 ! .. status = "<<fitResult->covQual()<<endl;
      continue; // skip the channel
    }

    // plot to store in a root file
    TH1* frame = (TH1*) ihist.second->Clone(Form("frame_%s",ihist.second->GetName()));
    frame->Reset();
    frame->GetYaxis()->SetTitle("a.u.");
    if(observable == "maxCharge")
      frame->GetXaxis()->SetTitle("Leading strip charge (ADC)");
    else
      frame->GetXaxis()->SetTitle("Signal over Noise");

    frame->GetXaxis()->SetTitleOffset(1.1);
    frame->GetYaxis()->SetTitleOffset(1.1);
    frame->GetYaxis()->SetTitleSize(0.040);
    frame->GetXaxis()->SetTitleSize(0.040);
    frame->GetYaxis()->SetLabelSize(0.035);
    frame->GetXaxis()->SetLabelSize(0.035);
    frame->Draw();

    ihist.second->SetMarkerColor(kBlack);
    ihist.second->SetMarkerSize(1);
    ihist.second->SetMarkerStyle(20);
    ihist.second->SetLineColor(kBlack);
    // convert the pdf as a TF1
    RooArgList pars = fitResult->floatParsFinal();
    TF1* fitfunc = totalPdf->asTF(*charge,pars,*charge);    
    //Divide the histogram by the bin width
    ihist.second->Scale(1./ihist.second->Integral(),"width");
    
    fitfunc->SetLineColor(kRed);
    fitfunc->SetLineWidth(2);
    fitfunc->Draw("Lsame");
    ihist.second->Draw("EPsame");
    frame->GetYaxis()->SetRangeUser(0,ihist.second->GetMaximum()*1.5);
    CMS_lumi(canvas,"",false);

    leg.Clear();
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.AddText("");
    leg.AddText(Form("Landau Mean  = %.1f #pm %.1f",mean_landau->getVal(),mean_landau->getError()));
    leg.AddText(Form("Landau Width = %.1f #pm %.1f",sigma_landau->getVal(),sigma_landau->getError()));
    leg.AddText(Form("Gauss Mean  = %.1f #pm %.1f",mean_gauss->getVal(),mean_gauss->getError()));
    leg.AddText(Form("Gauss Width = %.1f #pm %.1f",sigma_gauss->getVal(),sigma_gauss->getError()));
    leg.AddText(Form("Integral = %.1f #pm %.1f",normalization->getVal(),normalization->getError()));
    // to calculate the chi2
    RooPlot* rooframe = charge->frame();
    chargeDistribution->plotOn(rooframe);
    totalPdf->plotOn(rooframe);
    float chi2 = rooframe->chiSquare(pars.getSize());
    leg.AddText(Form("#chi2/ndf = %.3f",chi2));
    leg.Draw("same");
    if(detId % savePlotEvery == 0 and saveCanvas){
      canvas->Write(frame->GetName()); 
      // plot to store in canvas
      if(plotCanvas){
	canvas->SaveAs((outputDirectory+"/"+string(ihist.second->GetName())+".png").c_str(),"png");
	canvas->SaveAs((outputDirectory+"/"+string(ihist.second->GetName())+".pdf").c_str(),"pdf");
      }
    }
    
    mapPeakCharge[to_string(ihist.first)] = to_string(fitfunc->GetMaximumX(xMin,xMax));        
    mapMeanCharge[to_string(ihist.first)] = to_string(fitfunc->Mean(xMin,xMax,fitfunc->GetParameters()));

    if(totalPdf) delete totalPdf;
    if(landauXgauss) delete landauXgauss;
    if(gauss) delete gauss;
    if(landau) delete landau;
    if(chargeDistribution) delete chargeDistribution;
    if(normalization) delete normalization;
    if(sigma_gauss) delete sigma_gauss;
    if(mean_gauss) delete mean_gauss;
    if(mean_landau) delete mean_landau;
    if(sigma_landau) delete sigma_landau;
    if(charge) delete charge;
    if(fitResult) delete fitResult;
    if(vars) delete vars;
    if(fitfunc) delete fitfunc;
    if(rooframe) delete rooframe;
  }
    
  cout<<"#### end of fit stage: nFits "<<detId<<" Invalid Fits "<<invalidFits<<" : "<<100*double(invalidFits)/detId<<" % "<<endl;
  // make the detid:peak map to be plotted on  the tracker map through http://test-stripdbmonitor.web.cern.ch/test-stripdbmonitor/PrintTrackerMap/print_TrackerMap.php
  cout<<"#### Dump peak channel map"<<endl;
  ofstream dumpPeak((outputDirectory+"/dumpMapPeak"+observable+".txt").c_str());
  for(auto imap : mapPeakCharge)
    dumpPeak<<imap.first<<"   "<<imap.second<<"\n";
  dumpPeak.close();

  cout<<"#### Dump mean channel map"<<endl;
  ofstream dumpMean((outputDirectory+"/dumpMapMean"+observable+".txt").c_str());
  for(auto imap : mapMeanCharge)
    dumpMean<<imap.first<<"   "<<imap.second<<"\n";
  dumpMean.close();
  
  cout<<"#### Close output Root file"<<endl;
  outputFile->Close("R");
  
}

