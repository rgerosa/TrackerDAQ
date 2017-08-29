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
#include "../delayUtils.h"

using namespace std; 

// reduce the number of events by
static int  reductionFactor = 1;

// create the plots in R slices 
static std::map<int32_t, TH1F* > TIBMaxCharge;
static std::map<int32_t, TH1F* > TIDMaxCharge;
static std::map<int32_t, TH1F* > TOBMaxCharge;
static std::map<int32_t, TH1F* > TECTMaxCharge;
static std::map<int32_t, TH1F* > TECtMaxCharge;

static std::map<int32_t, TF1* > fitTIBMaxCharge;
static std::map<int32_t, TF1* > fitTIDMaxCharge;
static std::map<int32_t, TF1* > fitTOBMaxCharge;
static std::map<int32_t, TF1* > fitTECTMaxCharge;
static std::map<int32_t, TF1* > fitTECtMaxCharge;

static std::map<int32_t, RooAbsPdf* > pdfTIBMaxCharge;
static std::map<int32_t, RooAbsPdf* > pdfTIDMaxCharge;
static std::map<int32_t, RooAbsPdf* > pdfTOBMaxCharge;
static std::map<int32_t, RooAbsPdf* > pdfTECTMaxCharge;
static std::map<int32_t, RooAbsPdf* > pdfTECtMaxCharge;

TF1* makeLandauGausFit(TH1F* histoToFit, int & status, string subdetector, const float & delay, const string & observable, pair<float,RooAbsPdf*> & pair){

  float xMin   = 0, xMax = 0;
  int   nBinsX = 0;
  setLimitsAndBinning(observable,xMin,xMax,nBinsX);

  // make the observable                                                                                                                                                                            
  RooRealVar* charge = new RooRealVar(Form("charge_%s_delay_%f",subdetector.c_str(),delay),"",xMin,xMax);
  RooArgList*  vars  = new RooArgList(*charge);
  // covert histogram in a RooDataHist                                                                                                                                                               
  RooDataHist* chargeDistribution = new RooDataHist(histoToFit->GetName(),"",*vars,histoToFit);
  // Build a Landau PDF                                                                                                                                                                              
  RooRealVar*  mean_landau  = new RooRealVar(Form("mean_landau_%s_delay_%f",subdetector.c_str(),delay),"",histoToFit->GetMean(),xMin,xMax);
  RooRealVar*  sigma_landau = new RooRealVar(Form("sigma_landau_%s_delay_%f",subdetector.c_str(),delay),"",histoToFit->GetRMS(),1.,100);
  RooLandau*   landau       = new RooLandau(Form("landau_delay_%s_%f",subdetector.c_str(),delay),"",*charge,*mean_landau,*sigma_landau) ;
  // Build a Gaussian PDF                                                                                                                                                                            
  RooRealVar*  mean_gauss  = new RooRealVar(Form("mean_gauss_%s_delay_%f",subdetector.c_str(),delay),"",0.,-150,150);
  RooRealVar*  sigma_gauss = new RooRealVar(Form("sigma_gauss_%s_delay_%f",subdetector.c_str(),delay),"",10,1.,50);
  RooGaussian* gauss       = new RooGaussian(Form("gauss_%s_delay_%f",subdetector.c_str(),delay),"",*charge,*mean_gauss,*sigma_gauss);
  // Make a convolution                                                                                                                                                                              
  RooFFTConvPdf* landauXgauss =  new RooFFTConvPdf(Form("landauXgauss_%s_delay_%f",subdetector.c_str(),delay),"",*charge,*landau,*gauss);
  // Make an extended PDF                                                                                                                                                                            
  RooRealVar*   normalization = new RooRealVar(Form("normalization_%s_delay_%f",subdetector.c_str(),delay),"",histoToFit->Integral(),histoToFit->Integral()/5,histoToFit->Integral()*5);
  RooExtendPdf* totalPdf      = new RooExtendPdf(Form("totalPdf_%s_delay_%f",subdetector.c_str(),delay),"",*landauXgauss,*normalization);

  pair.first  = delay;
  pair.second = totalPdf;

  // Perform the FIT                                                                                                                                                                                
  RooFitResult* fitResult = totalPdf->fitTo(*chargeDistribution,RooFit::Range(xMin,xMax),RooFit::Extended(kTRUE),RooFit::Save(kTRUE));

  // convert the pdf as a TF1                                                                                                                                                                        
  RooArgList pars = fitResult->floatParsFinal();
  status = fitResult->status();

  return totalPdf->asTF(*charge,pars,*charge);

}


//// Per ring analysis
void paritionPlots(TTree* tree, 
		   TTree* map,
		   const std::string & observable,
		   const std::string & outputDIR){



  std::cout<<"######################################################"<<std::endl;
  std::cout << "Preparing Ring plot for the different partitions    "<< std::endl; 
  std::cout<<"######################################################"<<std::endl;

  cout<<"tree set branch status "<<endl;
  // TTree Reader appear not to be working with addFriend and EventList
  uint32_t detid;
  float    maxCharge, clCorrectedSignalOverNoise, clSignalOverNoise, clglobalX, clglobalY, clglobalZ, thickness, obs;
  tree->SetBranchStatus("*",kFALSE);
  tree->SetBranchStatus("detid",kTRUE);
  tree->SetBranchStatus("clCorrectedSignalOverNoise",kTRUE);
  tree->SetBranchStatus("clSignalOverNoise",kTRUE);
  tree->SetBranchStatus("clglobalX",kTRUE);
  tree->SetBranchStatus("clglobalY",kTRUE);
  tree->SetBranchStatus("clglobalZ",kTRUE);
  tree->SetBranchStatus("thickness",kTRUE);
  tree->SetBranchStatus(observable.c_str(),kTRUE);
  tree->SetBranchAddress("detid",&detid);
  tree->SetBranchAddress("clCorrectedSignalOverNoise",&clCorrectedSignalOverNoise);
  tree->SetBranchAddress("clSignalOverNoise",&clSignalOverNoise);
  tree->SetBranchAddress("clglobalX",&clglobalX);
  tree->SetBranchAddress("clglobalY",&clglobalY);
  tree->SetBranchAddress("clglobalZ",&clglobalZ);
  tree->SetBranchAddress("thickness",&thickness);
  tree->SetBranchAddress(observable.c_str(),&obs);

  float delay;
  map->SetBranchStatus("*",kFALSE);
  map->SetBranchStatus("detid",kTRUE);
  map->SetBranchStatus("delay",kTRUE);
  map->SetBranchAddress("delay",&delay);

  // create vectors for the different Profiles
  float yMin   = 0, yMax = 0;
  int   nBinsY = 0;
  setLimitsAndBinning(observable,yMin,yMax,nBinsY);
  
  cout<<"#### Run over the events"<<endl;
  // create vectors  std::cout<<"Tree with nEntries "<<tree->GetEntries()<<std::endl;
  long int iEvent = 0;
  for( ; iEvent < tree->GetEntries()/reductionFactor; iEvent++){    

    tree->GetEntry(iEvent);
    // take the map delay from the detid
    map->GetEntryWithIndex(detid);

    cout.flush();
    if(iEvent % 100000 == 0) cout<<"\r"<<"iEvent "<<100*double(iEvent)/(tree->GetEntries()/reductionFactor)<<" % ";
    
    uint32_t subdetid    = int((detid-0x10000000)/0x2000000);
    float    R           = sqrt(clglobalX*clglobalX+clglobalY*clglobalY+clglobalZ*clglobalZ);
    float    value       = 0;
    if(observable == "maxCharge")
      value = obs*(clCorrectedSignalOverNoise)/(clSignalOverNoise);
    else
      value = obs;
    
    if(subdetid == 3){ // TIB
      if(TIBMaxCharge[round(delay*10/10)] == 0 or TIBMaxCharge[round(delay*10/10)] == NULL)
	TIBMaxCharge[round(delay*10/10)] = new TH1F(Form("TIB_delay_%.1f",delay),"",nBinsY,yMin,yMax);
      TIBMaxCharge[round(delay*10/10)]->Fill(value);
    }
    else if(subdetid == 5){ // TOB
      if(TOBMaxCharge[round(delay*10/10)] == 0 or TOBMaxCharge[round(delay*10/10)] == NULL)
	TOBMaxCharge[round(delay*10/10)] = new TH1F(Form("TOB_delay_%.1f",delay),"",nBinsY,yMin,yMax);
      TOBMaxCharge[round(delay*10/10)]->Fill(value);
    }
    else if(subdetid == 4){//TID
      if(TIDMaxCharge[round(delay*10/10)] == 0 or TIDMaxCharge[round(delay*10/10)] == NULL)
	TIDMaxCharge[round(delay*10/10)] = new TH1F(Form("TID_delay_%.1f",delay),"",nBinsY,yMin,yMax);
      TIDMaxCharge[round(delay*10/10)]->Fill(value);
    }

    else if((subdetid == 6  and clglobalZ > 0 and thickness > 400) or (subdetid == 6  and clglobalZ < 0 and thickness > 400)){ //TEC thick
      if(TECTMaxCharge[round(delay*10/10)] == 0 or TECTMaxCharge[round(delay*10/10)] == NULL)
	TECTMaxCharge[round(delay*10/10)] = new TH1F(Form("TECT_delay_%.1f",delay),"",nBinsY,yMin,yMax);
      TECTMaxCharge[round(delay*10/10)]->Fill(value);
    }

    else if((subdetid == 6  and clglobalZ > 0 and thickness < 400) or (subdetid == 6  and clglobalZ < 0 and thickness < 400)){ //TEC thick
      if(TECtMaxCharge[round(delay*10/10)] == 0 or TECtMaxCharge[round(delay*10/10)] == NULL)
	TECtMaxCharge[round(delay*10/10)] = new TH1F(Form("TECt_delay_%.1f",delay),"",nBinsY,yMin,yMax);
      TECtMaxCharge[round(delay*10/10)]->Fill(value);
    }
  }
    
  std::cout<<std::endl;
  std::cout<<"Loop on events terminated"<<std::endl;

  std::cout<<"Make a fit of the cluster charge or S/N for TIB "<<endl;
  long int iBadChannelFit = 0;
  RooAbsPdf* pdf = NULL;
  TF1* fitfunc = NULL;
  pair<float,RooAbsPdf*> pair; 
  for(auto imap : TIBMaxCharge){ // loop on the different 
    int status = -1;
    fitfunc = makeLandauGausFit(imap.second,status,"TIB",imap.first,observable,pair);
    pdfTIBMaxCharge[imap.first] = pair.second;
    fitTIBMaxCharge[imap.first] = fitfunc;    
    if(status != 0) iBadChannelFit++;
  }  
  std::cout<<"Bad channel fit from Landau+Gaus fit in TIB for each delay "<<iBadChannelFit<<" over "<<TIBMaxCharge.size()<<std::endl;

  std::cout<<"Make a fit of the cluster charge or S/N for TOB "<<endl;
  iBadChannelFit = 0;
  for(auto imap : TOBMaxCharge){ // loop on the different 
    int status = -1;
    fitfunc = makeLandauGausFit(imap.second,status,"TOB",imap.first,observable,pair);
    pdfTOBMaxCharge[imap.first] = pair.second;
    fitTOBMaxCharge[imap.first] = fitfunc;    
    if(status != 0) iBadChannelFit++;
  }  
  std::cout<<"Bad channel fit from Landau+Gaus fit in TOB for each delay "<<iBadChannelFit<<" over "<<TOBMaxCharge.size()<<std::endl;

  std::cout<<"Make a fit of the cluster charge or S/N for TID "<<endl;
  iBadChannelFit = 0;
  for(auto imap : TIDMaxCharge){ // loop on the different 
    int status = -1;
    fitfunc = makeLandauGausFit(imap.second,status,"TID",imap.first,observable,pair);
    pdfTIDMaxCharge[imap.first] = pair.second;
    fitTIDMaxCharge[imap.first] = fitfunc;    
    if(status != 0) iBadChannelFit++;
  }  
  std::cout<<"Bad channel fit from Landau+Gaus fit in TID for each delay "<<iBadChannelFit<<" over "<<TIDMaxCharge.size()<<std::endl;

  std::cout<<"Make a fit of the cluster charge or S/N for TECT "<<endl;
  iBadChannelFit = 0;
  for(auto imap : TECTMaxCharge){ // loop on the different 
    int status = -1;
    fitfunc = makeLandauGausFit(imap.second,status,"TECT",imap.first,observable,pair);
    pdfTECTMaxCharge[imap.first] = pair.second;
    fitTECTMaxCharge[imap.first] = fitfunc;    
    if(status != 0) iBadChannelFit++;
  }  
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECT for each delay "<<iBadChannelFit<<" over "<<TECTMaxCharge.size()<<std::endl;

  std::cout<<"Make a fit of the cluster charge or S/N for TECt "<<endl;
  iBadChannelFit = 0;
  for(auto imap : TECtMaxCharge){ // loop on the different 
    int status = -1;
    fitfunc = makeLandauGausFit(imap.second,status,"TECT",imap.first,observable,pair);
    pdfTECtMaxCharge[imap.first] = pair.second;
    fitTECtMaxCharge[imap.first] = fitfunc;    
    if(status != 0) iBadChannelFit++;
  }  
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECt for each delay "<<iBadChannelFit<<" over "<<TECtMaxCharge.size()<<std::endl;
  
  // correct profiles and fit them with a gaussian
  std::cout<<"#################################"<<std::endl;
  std::cout<<"##### End of Ring Analysis #####"<<std::endl;
  std::cout<<"#################################"<<std::endl;
  
}

////
void plotDistributions(TCanvas* canvas, const string & outputDIR){

  TH1F* frame = (TH1F*) TIBMaxCharge[0]->Clone("frame");
  frame->Reset();

  for(auto imap : TIBMaxCharge){

    // normalize as a pdf
    TIBMaxCharge[imap.first]->Scale(1./TIBMaxCharge[imap.first]->Integral(),"width");
    TOBMaxCharge[imap.first]->Scale(1./TOBMaxCharge[imap.first]->Integral(),"width");
    TIDMaxCharge[imap.first]->Scale(1./TIDMaxCharge[imap.first]->Integral(),"width");
    TECTMaxCharge[imap.first]->Scale(1./TECTMaxCharge[imap.first]->Integral(),"width");
    TECtMaxCharge[imap.first]->Scale(1./TECtMaxCharge[imap.first]->Integral(),"width");
    

    frame->GetXaxis()->SetTitle("leading strip charge (ADC)");
    frame->GetXaxis()->SetTitleOffset(1.1);
    frame->GetYaxis()->SetTitle("a.u.");
    frame->GetYaxis()->SetTitleOffset(1.35);
    frame->GetYaxis()->SetRangeUser(0,max(TIBMaxCharge[imap.first]->GetMaximum(),max(TOBMaxCharge[imap.first]->GetMaximum(),max(TIDMaxCharge[imap.first]->GetMaximum(),max(TECTMaxCharge[imap.first]->GetMaximum(),TECtMaxCharge[imap.first]->GetMaximum()))))*1.2);
    frame->Draw();

    CMS_lumi(canvas,"",false,false,0.4);
    
    fitTIBMaxCharge[imap.first]->SetLineColor(kBlack);
    fitTIBMaxCharge[imap.first]->SetLineWidth(2);
    fitTOBMaxCharge[imap.first]->SetLineColor(TColor::GetColor("#CF3721"));
    fitTOBMaxCharge[imap.first]->SetLineWidth(2);
    fitTIDMaxCharge[imap.first]->SetLineColor(kBlue);
    fitTIDMaxCharge[imap.first]->SetLineWidth(2);
    fitTECTMaxCharge[imap.first]->SetLineColor(TColor::GetColor("#3A8C4C"));
    fitTECTMaxCharge[imap.first]->SetLineWidth(2);
    fitTECtMaxCharge[imap.first]->SetLineColor(TColor::GetColor("#FAAF08"));
    fitTECtMaxCharge[imap.first]->SetLineWidth(2);

    TIBMaxCharge[imap.first]->SetMarkerColor(kBlack);
    TIBMaxCharge[imap.first]->SetLineColor(kBlack);
    TIBMaxCharge[imap.first]->SetMarkerSize(1);
    TIBMaxCharge[imap.first]->SetMarkerStyle(20);
    TOBMaxCharge[imap.first]->SetMarkerColor(TColor::GetColor("#CF3721"));
    TOBMaxCharge[imap.first]->SetLineColor(TColor::GetColor("#CF3721"));
    TOBMaxCharge[imap.first]->SetMarkerSize(1);
    TOBMaxCharge[imap.first]->SetMarkerStyle(20);
    TIDMaxCharge[imap.first]->SetMarkerColor(kBlue);
    TIDMaxCharge[imap.first]->SetLineColor(kBlue);
    TIDMaxCharge[imap.first]->SetMarkerSize(1);
    TIDMaxCharge[imap.first]->SetMarkerStyle(20);
    TECTMaxCharge[imap.first]->SetLineColor(TColor::GetColor("#3A8C4C"));
    TECTMaxCharge[imap.first]->SetMarkerColor(TColor::GetColor("#3A8C4C"));
    TECTMaxCharge[imap.first]->SetMarkerSize(1);
    TECTMaxCharge[imap.first]->SetMarkerStyle(20);
    TECtMaxCharge[imap.first]->SetMarkerColor(TColor::GetColor("#FAAF08"));
    TECtMaxCharge[imap.first]->SetLineColor(TColor::GetColor("#FAAF08"));
    TECtMaxCharge[imap.first]->SetMarkerSize(1);
    TECtMaxCharge[imap.first]->SetMarkerStyle(20);

    fitTIBMaxCharge[imap.first]->Draw("Lsame");
    fitTOBMaxCharge[imap.first]->Draw("Lsame");
    fitTIDMaxCharge[imap.first]->Draw("Lsame");
    fitTECTMaxCharge[imap.first]->Draw("Lsame");
    fitTECtMaxCharge[imap.first]->Draw("Lsame");
    TIBMaxCharge[imap.first]->Draw("EPsame");
    TOBMaxCharge[imap.first]->Draw("EPsame");
    TIDMaxCharge[imap.first]->Draw("EPsame");
    TECTMaxCharge[imap.first]->Draw("EPsame");
    TECtMaxCharge[imap.first]->Draw("EPsame");

    TLegend leg (0.55,0.56,0.82,0.82);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.AddEntry((TObject*)(0),"2017 Data","");
    leg.AddEntry(TIBMaxCharge[imap.first],"TIB","EP");
    leg.AddEntry(TOBMaxCharge[imap.first],"TOB","EP");
    leg.AddEntry(TIDMaxCharge[imap.first],"TID","EP");
    leg.AddEntry(TECTMaxCharge[imap.first],"TEC thick","EP");
    leg.AddEntry(TECtMaxCharge[imap.first],"TEC thin","EP");
    leg.Draw("same");

    canvas->SaveAs((outputDIR+"/clusterCharge_delay_"+to_string(imap.first)+".png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/clusterCharge_delay_"+to_string(imap.first)+".pdf").c_str(),"pdf");
  }

}


/// main function that run the analysis --> produce one plot per delay value showing the 4 partitions on the same canvas
void makeChargeDistributionPerPartition(string file0,  // inputfile
					string observable   = "maxCharge",   // observable to be considered: maxCharge, S/N ..etc
					string outputDIR    = "distributionMaxChargePerPartition" // output directory name
					){
  gStyle->SetOptStat(0);
  gROOT->SetBatch(kTRUE);
  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;

  // prepare style and load macros
  setTDRStyle();
  // not dump stat and fit info
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gROOT->SetBatch(kTRUE);
  gStyle->SetPadLeftMargin(0.14);

  system(("mkdir -p "+outputDIR).c_str());

  // open the file and prepare the cluster tree, adding the other trees as frined --> memory consuming
  std::cout<<"################################################"<<std::endl;
  std::cout<<"###### makeChargeDistributionPerPartition ######"<<std::endl;
  std::cout<<"################################################"<<std::endl;

  std::cout<<"Open Input Files"<<std::endl;
  TFile* _file0 (TFile::Open(file0.c_str()));
  TTree* clusters   ((TTree*)_file0->FindObjectAny("clusters"));
  TTree* readoutMap ((TTree*)_file0->FindObjectAny("readoutMap"));
  clusters->SetEventList(0);  
  
  // cumulate all rings
  paritionPlots(clusters,readoutMap,observable,outputDIR);
  TCanvas* canvas = new TCanvas("canvas","",600,650);
  canvas->cd();
  plotDistributions(canvas,outputDIR);
  
}

