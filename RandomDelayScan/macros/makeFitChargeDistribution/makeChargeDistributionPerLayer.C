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

// parametrize profile with a gaussian shape
static bool isGaussian = true;
// reduce the number of events by
static int  reductionFactor = 1;

static std::map<uint32_t,std::map<float,std::shared_ptr<TH1F> > > TIBlayers; // map layer:delay:distribution
static std::map<uint32_t,std::map<float,std::shared_ptr<TH1F> > > TOBlayers;
static std::map<uint32_t,std::map<float,std::shared_ptr<TH1F> > > TIDlayers;
static std::map<uint32_t,std::map<float,std::shared_ptr<TH1F> > > TECPTlayers;
static std::map<uint32_t,std::map<float,std::shared_ptr<TH1F> > > TECPtlayers;
static std::map<uint32_t,std::map<float,std::shared_ptr<TH1F> > > TECMTlayers;
static std::map<uint32_t,std::map<float,std::shared_ptr<TH1F> > > TECMtlayers;

static std::map<uint32_t,std::shared_ptr<TH1F> > TIBlayersMean; // map delay:distribution
static std::map<uint32_t,std::shared_ptr<TH1F> > TOBlayersMean;
static std::map<uint32_t,std::shared_ptr<TH1F> > TIDlayersMean;
static std::map<uint32_t,std::shared_ptr<TH1F> > TECPTlayersMean;
static std::map<uint32_t,std::shared_ptr<TH1F> > TECPtlayersMean;
static std::map<uint32_t,std::shared_ptr<TH1F> > TECMTlayersMean;
static std::map<uint32_t,std::shared_ptr<TH1F> > TECMtlayersMean;

static std::map<uint32_t,std::shared_ptr<TH1F> > TIBlayersMPV; // map delay:distribution
static std::map<uint32_t,std::shared_ptr<TH1F> > TOBlayersMPV;
static std::map<uint32_t,std::shared_ptr<TH1F> > TIDlayersMPV;
static std::map<uint32_t,std::shared_ptr<TH1F> > TECPTlayersMPV;
static std::map<uint32_t,std::shared_ptr<TH1F> > TECPtlayersMPV;
static std::map<uint32_t,std::shared_ptr<TH1F> > TECMTlayersMPV;
static std::map<uint32_t,std::shared_ptr<TH1F> > TECMtlayersMPV;

static TFile*   outputFitFile = NULL;
static TCanvas* outputCanvasFit = NULL;

int makeLandauGausFit(TH1F* histoToFit, TH1F* histoToFill, string subdetector, const float & delay, const string & observable, const string & outputDIR, const string & postfix = ""){

  float yMin   = 0, yMax = 0;
  int   nBinsY = 0;
  setLimitsAndBinning(observable,yMin,yMax,nBinsY);

  if(outputCanvasFit == NULL or outputCanvasFit == 0)
    outputCanvasFit = new TCanvas("canvas","",600,625);

  if(outputFitFile == NULL or outputFitFile == 0)
    outputFitFile = new TFile((outputDIR+"/outputFitCanvases_"+observable+"_"+postfix+".root").c_str(),"RECREATE");
  
  outputCanvasFit->SetTickx();
  outputCanvasFit->SetTicky();
  outputFitFile->cd();

  Double_t parameters[4];
  Double_t parametersHigh[4];
  Double_t parametersLow[4];
  const Double_t *fit_parameters;
  const Double_t *fit_parameters_error;
  Double_t chi2;
  Int_t    ndf;
  // Width of the landau distribution                                                                                                                                      
  parameters[0] = 2.; parametersHigh[0] = 50.; parametersLow[0] = 0.001;
  // MPV of landau peak                                                                                                                                                    
  parameters[1] = histoToFit->GetMean(); parametersHigh[1] = yMax; parametersLow[1] = yMin;
  // Total area                                                                                                                                                            
  parameters[2] = histoToFit->Integral(); parametersHigh[2] = histoToFit->Integral()*5; parametersLow[2] = histoToFit->Integral()/5;
  // width of gaussian                                                                                                                                                     
  parameters[3] = 10; parametersHigh[3] = 40; parametersLow[3] = 1.;
  
      // create the function                                                                                                                                                   
  TF1 *    fitfunc = new TF1(Form("fit_%s",histoToFit->GetName()),langaufun,yMin,yMax,4);
  fitfunc->SetParameters(parameters);
  fitfunc->SetParNames("Width","MPV","Area","GSigma");
  fitfunc->SetParLimits(0,parametersLow[0],parametersHigh[0]);
  fitfunc->SetParLimits(0,parametersLow[1],parametersHigh[1]);
  fitfunc->SetParLimits(0,parametersLow[2],parametersHigh[2]);
  fitfunc->SetParLimits(0,parametersLow[3],parametersHigh[3]);

  // define range for a better description of the peak                                                                                                                     
  if(observable == "maxCharge"){
    if(TString(subdetector).Contains("TIB")){ // TIB                                                                                                                        
      fitfunc->SetRange(histoToFit->GetBinCenter(histoToFit->GetMaximumBin())-1.2*histoToFit->GetRMS(),histoToFit->GetBinCenter(histoToFit->GetMaximumBin())+2.5*histoToFit->GetRMS());
    }
    else if(TString(subdetector).Contains("TID")) // TID                                                                                                                  
      fitfunc->SetRange(histoToFit->GetBinCenter(histoToFit->GetMaximumBin())-1.5*histoToFit->GetRMS(),histoToFit->GetBinCenter(histoToFit->GetMaximumBin())+2.5*histoToFit->GetRMS());
    else if(TString(subdetector).Contains("TOB")){ //TOB                                                                                                   
      fitfunc->SetRange(histoToFit->GetBinCenter(histoToFit->GetMaximumBin())-2.0*histoToFit->GetRMS(),histoToFit->GetBinCenter(histoToFit->GetMaximumBin())+3.0*histoToFit->GetRMS());
    }
    else if(TString(subdetector).Contains("TECInner"))
      fitfunc->SetRange(histoToFit->GetBinCenter(histoToFit->GetMaximumBin())-1.5*histoToFit->GetRMS(),histoToFit->GetBinCenter(histoToFit->GetMaximumBin())+2.5*histoToFit->GetRMS());
    else if(TString(subdetector).Contains("TECOuter"))
      fitfunc->SetRange(histoToFit->GetBinCenter(histoToFit->GetMaximumBin())-2.0*histoToFit->GetRMS(),histoToFit->GetBinCenter(histoToFit->GetMaximumBin())+3.0*histoToFit->GetRMS());
  }
  else
    fitfunc->SetRange(histoToFit->GetBinCenter(histoToFit->GetMaximumBin())-1.5*histoToFit->GetRMS(),histoToFit->GetBinCenter(histoToFit->GetMaximumBin())+histoToFit->GetRMS()*3.0);
  
  // make fit and get parameters                                                                                                                                           
  TFitResultPtr fitResult = histoToFit->Fit(fitfunc,"RSQN");
  if(not fitResult.Get() or fitResult->Status() != 0 or fitResult->CovMatrixStatus() <= 1){
    histoToFill->SetBinContent(histoToFill->FindBin(delay),histoToFit->GetBinCenter(histoToFit->GetMaximumBin()));
    histoToFill->SetBinError(histoToFill->FindBin(delay),histoToFit->GetMeanError());
    return -1;
  }
  else{
    fit_parameters = fitResult->GetParams();
    fit_parameters_error = fitResult->GetErrors();
    chi2   = fitResult->Chi2();
    ndf    = fitResult->Ndf();
    histoToFill->SetBinContent(histoToFill->FindBin(delay),fitfunc->GetMaximumX(yMin,yMax));
    histoToFill->SetBinError(histoToFill->FindBin(delay),fit_parameters_error[1]);
  }
  //draw fit
  // plot to store in a root file                                                                                                                                      
  outputCanvasFit->cd();
  histoToFit->GetYaxis()->SetTitle("Events");
  histoToFit->GetXaxis()->SetTitle(("Cluster "+observable).c_str());
  histoToFit->GetYaxis()->SetTitleSize(0.045);
  histoToFit->GetXaxis()->SetTitleSize(0.045);
  histoToFit->GetYaxis()->SetLabelSize(0.038);
  histoToFit->GetXaxis()->SetLabelSize(0.038);
  histoToFit->SetMarkerColor(kBlack);
  histoToFit->SetMarkerStyle(20);
  histoToFit->Draw("EP");
  fitfunc->SetLineWidth(2);
  fitfunc->SetLineColor(kRed);
  fitfunc->Draw("same");
  histoToFit->Draw("EPsame");
  TLegend leg (0.6,0.6,0.95,0.92);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry((TObject*)0,Form("#chi^{2}/ndf = %.2f",chi2/ndf),"");
  leg.AddEntry((TObject*)0,Form("Pdf Max = %.2f",fitfunc->GetMaximumX(yMin,yMax)),"");
  leg.AddEntry((TObject*)0,Form("Landau MPV = %.2f #pm %.2f",fit_parameters[1],fit_parameters_error[1]),"");
  leg.AddEntry((TObject*)0,Form("Landau Width = %.2f #pm %.2f",fit_parameters[0],fit_parameters_error[0]),"");
  leg.AddEntry((TObject*)0,Form("Gaussian #sigma = %.2f #pm %.2f",fit_parameters[3],fit_parameters_error[3]),"");
  leg.AddEntry((TObject*)0,Form("Normalization = %.2f #pm %.2f",fit_parameters[2],fit_parameters_error[2]),"");
  leg.Draw("same");
  outputCanvasFit->Modified();
  CMS_lumi(outputCanvasFit,"");
  outputCanvasFit->Write(histoToFit->GetName());
    
  return 1;
}



/// function that runs on the evnet and produce profiles for layers
void LayerPlots(const std::shared_ptr<TTree> & tree, 
		const std::shared_ptr<TTree> & map,
		const std::shared_ptr<TTree> & corrections,
		const string & observable,
		const string & outputDIR) {

  std::cout<<"######################################################"<<std::endl;
  std::cout << "Preparing Layer plot for the fifferent subdetectors "<< std::endl; 
  std::cout<<"######################################################"<<std::endl;
  
  // set branches for the cluster, readoutmap and no corrections trees
  uint32_t detid;
  float    clCorrectedSignalOverNoise, clSignalOverNoise, clglobalX, clglobalY, clglobalZ, thickness, obs;
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

  int   correction;
  corrections->SetBranchStatus("*",kFALSE);
  corrections->SetBranchStatus("detid",kTRUE);
  corrections->SetBranchStatus("correction",kTRUE);
  corrections->SetBranchAddress("correction",&correction);
  
  // create vectors for the different Profiles
  float yMin   = 0, yMax = 0;
  int   nBinsY = 0;
  vector<double> delayBins;
  setLimitsAndBinning("delay",delayBins);
  setLimitsAndBinning(observable,yMin,yMax,nBinsY);

  std::cout<<"Tree with nEntries "<<tree->GetEntries()<<std::endl;

  long int iEvent  = 0;
  for( ; iEvent < tree->GetEntries()/reductionFactor; iEvent++){    
    // take the event
    tree->GetEntry(iEvent);
    // take the map delay from the detid
    map->GetEntryWithIndex(detid);
    // take the correction from the detid
    corrections->GetEntryWithIndex(detid);

    cout.flush();
    if(iEvent % 100000 == 0) cout<<"\r"<<"iEvent "<<100*double(iEvent)/(tree->GetEntries()/reductionFactor)<<" % ";
    
    uint32_t subdetid    = int((detid-0x10000000)/0x2000000);
    uint32_t barrellayer = int((detid%33554432)/0x4000);
    uint32_t TIDlayer    = int((detid%33554432)/0x800)%4;
    uint32_t TECPlayer   = int((detid%33554432)/0x4000)-32;
    uint32_t TECMlayer   = int((detid%33554432)/0x4000)-16;
    float    R           = sqrt(clglobalX*clglobalX+clglobalY*clglobalY+clglobalZ*clglobalZ);

    float value = 0;
    if(observable == "maxCharge")
      value = obs*(clCorrectedSignalOverNoise)/(clSignalOverNoise);
    else 
      value = obs;

    // fill the maps
    if(subdetid == 3){
      if(TIBlayers[barrellayer][delay-correction].get() == 0 or TIBlayers[barrellayer][delay-correction].get() == NULL) 
	TIBlayers[barrellayer][delay-correction] = std::shared_ptr<TH1F> (new TH1F(Form("TIB_layer_%d_delay_%.1f",barrellayer,delay-correction),"",nBinsY,yMin,yMax));
      TH1::AddDirectory(kFALSE);
      TIBlayers[barrellayer][delay-correction]->Fill(value);
    }
    else if(subdetid == 5){
      if(TOBlayers[barrellayer][delay-correction].get() == 0 or TOBlayers[barrellayer][delay-correction].get() == NULL) 
	TOBlayers[barrellayer][delay-correction] = std::shared_ptr<TH1F> (new TH1F(Form("TOB_layer_%d_delay_%.1f",barrellayer,delay-correction),"",nBinsY,yMin,yMax));
      TH1::AddDirectory(kFALSE);
      TOBlayers[barrellayer][delay-correction]->Fill(value);
    }
    else if(subdetid == 4){
      if(TIDlayers[TIDlayer][delay-correction].get() == 0 or TIDlayers[TIDlayer][delay-correction].get() == NULL) 
	TIDlayers[TIDlayer][delay-correction] = std::shared_ptr<TH1F> (new TH1F(Form("TID_layer_%d_delay_%.1f",TIDlayer,delay-correction),"",nBinsY,yMin,yMax));
      TH1::AddDirectory(kFALSE);
      TIDlayers[TIDlayer][delay-correction]->Fill(value);
    }
    else if(subdetid == 6  and clglobalZ > 0 and thickness > 400){
      if(TECPTlayers[TECPlayer][delay-correction].get() == 0 or TECPTlayers[TECPlayer][delay-correction].get() == NULL) 
	TECPTlayers[TECPlayer][delay-correction] = std::shared_ptr<TH1F> (new TH1F(Form("TECPT_layer_%d_delay_%.1f",TECPlayer,delay-correction),"",nBinsY,yMin,yMax));
      TH1::AddDirectory(kFALSE);
      TECPTlayers[TECPlayer][delay-correction]->Fill(value);
    }
    else if(subdetid == 6  and clglobalZ > 0 and thickness < 400){
      if(TECPtlayers[TECPlayer][delay-correction].get() == 0 or TECPtlayers[TECPlayer][delay-correction].get() == NULL) 
	TECPtlayers[TECPlayer][delay-correction] = std::shared_ptr<TH1F> (new TH1F(Form("TECPt_layer_%d_delay_%.1f",TECPlayer,delay-correction),"",nBinsY,yMin,yMax));
      TH1::AddDirectory(kFALSE);
      TECPtlayers[TECPlayer][delay-correction]->Fill(value);
    }
    else if(subdetid == 6  and clglobalZ < 0 and thickness > 400){
      if(TECMTlayers[TECMlayer][delay-correction].get() == 0 or TECMTlayers[TECMlayer][delay-correction].get() == NULL) 
	TECMTlayers[TECMlayer][delay-correction] = std::shared_ptr<TH1F> (new TH1F(Form("TECMT_layer_%d_delay_%.1f",TECMlayer,delay-correction),"",nBinsY,yMin,yMax));
      TH1::AddDirectory(kFALSE);
      TECMTlayers[TECMlayer][delay-correction]->Fill(value);
    }
    else if(subdetid == 6  and clglobalZ < 0 and thickness < 400){
      if(TECMtlayers[TECMlayer][delay-correction].get() == 0 or TECMtlayers[TECMlayer][delay-correction].get() == NULL) 
	TECMtlayers[TECMlayer][delay-correction] = std::shared_ptr<TH1F> (new TH1F(Form("TECMt_layer_%d_delay_%.1f",TECMlayer,delay-correction),"",nBinsY,yMin,yMax));
      TH1::AddDirectory(kFALSE);
      TECMtlayers[TECMlayer][delay-correction]->Fill(value);
    }
  }
  
  std::cout<<std::endl;
  std::cout<<"Loop on events terminated"<<std::endl;

  std::cout<<"Build the mean and MPV distributions asaf of delay "<<endl;
  long int iBadChannelFit = 0;
  for(auto imap : TIBlayers){ // loop on the different layer
    if(TIBlayersMean[imap.first].get() == 0 or TIBlayersMean[imap.first].get() == NULL)
      TIBlayersMean[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TIB_layer_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    if(TIBlayersMPV[imap.first].get() == 0 or TIBlayersMPV[imap.first].get() == NULL)
      TIBlayersMPV[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TIB_layer_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    
    for(auto idelay : imap.second){
      // mean value
      TIBlayersMean[imap.first]->SetBinContent(TIBlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TIBlayersMean[imap.first]->SetBinError(TIBlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());
      
      int status = makeLandauGausFit(idelay.second.get(),TIBlayersMPV[imap.first].get(),"TIB",idelay.first,observable,outputDIR);
      if(status != 1)
	iBadChannelFit++;            
      
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TIB layers "<<iBadChannelFit<<" over "<<TIBlayers.size()*TIBlayers[1].size()<<std::endl;  
  iBadChannelFit = 0;
  for(auto imap : TOBlayers){ // loop on the different layer
    if(TOBlayersMean[imap.first].get() == 0 or TOBlayersMean[imap.first].get() == NULL)
      TOBlayersMean[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TOB_layer_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    if(TOBlayersMPV[imap.first].get() == 0 or TOBlayersMPV[imap.first].get() == NULL)
      TOBlayersMPV[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TOB_layer_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    
    for(auto idelay : imap.second){
      // mean value
      TOBlayersMean[imap.first]->SetBinContent(TOBlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TOBlayersMean[imap.first]->SetBinError(TOBlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());

      int status = makeLandauGausFit(idelay.second.get(),TOBlayersMPV[imap.first].get(),"TOB",idelay.first,observable,outputDIR);
      if(status != 1)
	iBadChannelFit++;      
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TOB layers "<<iBadChannelFit<<" over "<<TOBlayers.size()*TOBlayers[1].size()<<std::endl;  
  iBadChannelFit = 0;
  for(auto imap : TIDlayers){ // loop on the different layer
    if(TIDlayersMean[imap.first].get() == 0 or TIDlayersMean[imap.first].get() == NULL)
      TIDlayersMean[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TID_layer_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    if(TIDlayersMPV[imap.first].get() == 0 or TIDlayersMPV[imap.first].get() == NULL)
      TIDlayersMPV[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TID_layer_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    
    for(auto idelay : imap.second){
      // mean value
      TIDlayersMean[imap.first]->SetBinContent(TIDlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TIDlayersMean[imap.first]->SetBinError(TIDlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());

      int status = makeLandauGausFit(idelay.second.get(),TIDlayersMPV[imap.first].get(),"TID",idelay.first,observable,outputDIR);
      if(status != 1)
	iBadChannelFit++;      
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TID layers "<<iBadChannelFit<<" over "<<TIDlayers.size()*TIDlayers[1].size()<<std::endl;  

  iBadChannelFit = 0;
  for(auto imap : TECPTlayers){ // loop on the different layer
    if(TECPTlayersMean[imap.first].get() == 0 or TECPTlayersMean[imap.first].get() == NULL)
      TECPTlayersMean[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TECPT_layer_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    if(TECPTlayersMPV[imap.first].get() == 0 or TECPTlayersMPV[imap.first].get() == NULL)
      TECPTlayersMPV[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TECPT_layer_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    
    for(auto idelay : imap.second){
      // mean value
      TECPTlayersMean[imap.first]->SetBinContent(TECPTlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TECPTlayersMean[imap.first]->SetBinError(TECPTlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());

      int status = makeLandauGausFit(idelay.second.get(),TECPTlayersMPV[imap.first].get(),"TECOuter",idelay.first,observable,outputDIR);
      if(status != 1)
	iBadChannelFit++;      
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECPT layers "<<iBadChannelFit<<" over "<<TECPTlayers.size()*TECPTlayers[1].size()<<std::endl;  

  iBadChannelFit = 0;
  for(auto imap : TECPtlayers){ // loop on the different layer
    if(TECPtlayersMean[imap.first].get() == 0 or TECPtlayersMean[imap.first].get() == NULL)
      TECPtlayersMean[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TECPt_layer_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    if(TECPtlayersMPV[imap.first].get() == 0 or TECPtlayersMPV[imap.first].get() == NULL)
      TECPtlayersMPV[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TECPt_layer_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    
    for(auto idelay : imap.second){
      // mean value
      TECPtlayersMean[imap.first]->SetBinContent(TECPtlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TECPtlayersMean[imap.first]->SetBinError(TECPtlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());

      int status = makeLandauGausFit(idelay.second.get(),TECPtlayersMPV[imap.first].get(),"TECInner",idelay.first,observable,outputDIR);
      if(status != 1)
	iBadChannelFit++;      
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECPt layers "<<iBadChannelFit<<" over "<<TECPtlayers.size()*TECPtlayers[1].size()<<std::endl;  

  iBadChannelFit = 0;
  for(auto imap : TECMTlayers){ // loop on the different layer
    if(TECMTlayersMean[imap.first].get() == 0 or TECMTlayersMean[imap.first].get() == NULL)
      TECMTlayersMean[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TECMT_layer_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    if(TECMTlayersMPV[imap.first].get() == 0 or TECMTlayersMPV[imap.first].get() == NULL)
      TECMTlayersMPV[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TECMT_layer_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    
    for(auto idelay : imap.second){
      // mean value
      TECMTlayersMean[imap.first]->SetBinContent(TECMTlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TECMTlayersMean[imap.first]->SetBinError(TECMTlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());

      int status = makeLandauGausFit(idelay.second.get(),TECMTlayersMPV[imap.first].get(),"TECOuter",idelay.first,observable,outputDIR);
      if(status != 1)
	iBadChannelFit++;      
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECMT layers "<<iBadChannelFit<<" over "<<TECMTlayers.size()*TECMTlayers[1].size()<<std::endl;  

  iBadChannelFit = 0;
  for(auto imap : TECMtlayers){ // loop on the different layer
    if(TECMtlayersMean[imap.first].get() == 0 or TECMtlayersMean[imap.first].get() == NULL)
      TECMtlayersMean[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TECMt_layer_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    if(TECMtlayersMPV[imap.first].get() == 0 or TECMtlayersMPV[imap.first].get() == NULL)
      TECMtlayersMPV[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TECMt_layer_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    
    for(auto idelay : imap.second){
      // mean value
      TECMtlayersMean[imap.first]->SetBinContent(TECMtlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TECMtlayersMean[imap.first]->SetBinError(TECMtlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());

      int status = makeLandauGausFit(idelay.second.get(),TECMtlayersMPV[imap.first].get(),"TECInner",idelay.first,observable,outputDIR);
      if(status != 1)
	iBadChannelFit++;      
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECMt layers "<<iBadChannelFit<<" over "<<TECMtlayers.size()*TECMtlayers[1].size()<<std::endl;  

  // correct profiles and fit them with a gaussian
  std::cout<<"Analyze profiles TIB mean"<<std::endl;
  long int badFits = 0;
  for(auto iprof : TIBlayersMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TIB Mean bad gaussian fits "<<badFits<<std::endl;
  std::cout<<"Analyze TOB profiles mean"<<std::endl;
  badFits = 0;
  for(auto iprof : TOBlayersMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TOB Mean bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"Analyze TID profiles mean"<<std::endl;
  badFits = 0;
  for(auto iprof : TIDlayersMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TID Mean bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"Analyze TECPT profiles mean"<<std::endl;
  badFits = 0;
  for(auto iprof : TECPTlayersMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TECPT Mean bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"Analyze TECPt profiles mean"<<std::endl;
  badFits = 0;
  for(auto iprof : TECPtlayersMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TECPt Mean bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"Analyze TECMT profiles mean"<<std::endl;
  badFits = 0;
  for(auto iprof : TECMTlayersMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TECMT Mean bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"Analyze TECMt profiles mean"<<std::endl;
  badFits = 0;
  for(auto iprof : TECMtlayersMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TECMt Mean bad gaussian fits "<<badFits<<std::endl;

  // correct profiles and fit them with a gaussian
  std::cout<<"Analyze TIB profiles mpv"<<std::endl;
  badFits = 0;
  for(auto iprof : TIBlayersMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TIB MPV bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"Analyze TOB profiles mpv"<<std::endl;
  badFits = 0;
  for(auto iprof : TOBlayersMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TOB MPV bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"Analyze TID profiles mpv"<<std::endl;
  badFits = 0;
  for(auto iprof : TIDlayersMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TID MPV bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"Analyze TECPT profiles mpv"<<std::endl;
  badFits = 0;
  for(auto iprof : TECPTlayersMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TECPT MPV bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"Analyze TECPt profiles mpv"<<std::endl;
  badFits = 0;
  for(auto iprof : TECPtlayersMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TECPt MPV bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"Analyze TECMT profiles mpv"<<std::endl;
  badFits = 0;
  for(auto iprof : TECMTlayersMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TECMT MPV bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"Analyze TECMt profiles mpv"<<std::endl;
  badFits = 0;
  for(auto iprof : TECMtlayersMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TECMt MPV bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"#################################"<<std::endl;
  std::cout<<"##### End of Layer Analysis #####"<<std::endl;
  std::cout<<"#################################"<<std::endl;

}

// create the plots in R slices 
static std::map<uint32_t, std::map<float,std::shared_ptr<TH1F> > > TIBrs;
static std::map<uint32_t, std::map<float,std::shared_ptr<TH1F> > > TIDrs;
static std::map<uint32_t, std::map<float,std::shared_ptr<TH1F> > > TOBrs;
static std::map<uint32_t, std::map<float,std::shared_ptr<TH1F> > > TECPTrs;
static std::map<uint32_t, std::map<float,std::shared_ptr<TH1F> > > TECPtrs;
static std::map<uint32_t, std::map<float,std::shared_ptr<TH1F> > > TECMTrs;
static std::map<uint32_t, std::map<float,std::shared_ptr<TH1F> > > TECMtrs;
static std::map<uint32_t, std::map<float,std::shared_ptr<TH1F> > > TIB;
static std::map<uint32_t, std::map<float,std::shared_ptr<TH1F> > > TID;
static std::map<uint32_t, std::map<float,std::shared_ptr<TH1F> > > TOB;
static std::map<uint32_t, std::map<float,std::shared_ptr<TH1F> > > TECPT;
static std::map<uint32_t, std::map<float,std::shared_ptr<TH1F> > > TECPt;
static std::map<uint32_t, std::map<float,std::shared_ptr<TH1F> > > TECMT;
static std::map<uint32_t, std::map<float,std::shared_ptr<TH1F> > > TECMt;

static std::map<uint32_t, std::shared_ptr<TH1F> > TIBrsMean;
static std::map<uint32_t, std::shared_ptr<TH1F> > TIDrsMean;
static std::map<uint32_t, std::shared_ptr<TH1F> > TOBrsMean;
static std::map<uint32_t, std::shared_ptr<TH1F> > TECPTrsMean;
static std::map<uint32_t, std::shared_ptr<TH1F> > TECPtrsMean;
static std::map<uint32_t, std::shared_ptr<TH1F> > TECMTrsMean;
static std::map<uint32_t, std::shared_ptr<TH1F> > TECMtrsMean;

static std::map<uint32_t, std::shared_ptr<TH1F> > TIBMean;
static std::map<uint32_t, std::shared_ptr<TH1F> > TIDMean;
static std::map<uint32_t, std::shared_ptr<TH1F> > TOBMean;
static std::map<uint32_t, std::shared_ptr<TH1F> > TECPTMean;
static std::map<uint32_t, std::shared_ptr<TH1F> > TECPtMean;
static std::map<uint32_t, std::shared_ptr<TH1F> > TECMTMean;
static std::map<uint32_t, std::shared_ptr<TH1F> > TECMtMean;

static std::map<uint32_t, std::shared_ptr<TH1F> > TIBrsMPV;
static std::map<uint32_t, std::shared_ptr<TH1F> > TIDrsMPV;
static std::map<uint32_t, std::shared_ptr<TH1F> > TOBrsMPV;
static std::map<uint32_t, std::shared_ptr<TH1F> > TECPTrsMPV;
static std::map<uint32_t, std::shared_ptr<TH1F> > TECPtrsMPV;
static std::map<uint32_t, std::shared_ptr<TH1F> > TECMTrsMPV;
static std::map<uint32_t, std::shared_ptr<TH1F> > TECMtrsMPV;
static std::map<uint32_t, std::shared_ptr<TH1F> > TIBMPV;
static std::map<uint32_t, std::shared_ptr<TH1F> > TIDMPV;
static std::map<uint32_t, std::shared_ptr<TH1F> > TOBMPV;
static std::map<uint32_t, std::shared_ptr<TH1F> > TECPTMPV;
static std::map<uint32_t, std::shared_ptr<TH1F> > TECPtMPV;
static std::map<uint32_t, std::shared_ptr<TH1F> > TECMTMPV;
static std::map<uint32_t, std::shared_ptr<TH1F> > TECMtMPV;


//// Per ring analysis
void RPlots(const std::shared_ptr<TTree> & tree, 
	    const std::shared_ptr<TTree> & map,
	    const std::shared_ptr<TTree> & corrections,
	    const std::string & observable,
	    const std::string & outputDIR){



  std::cout<<"######################################################"<<std::endl;
  std::cout << "Preparing Ring plot for the fifferent subdetectors "<< std::endl; 
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

  int   correction;
  corrections->SetBranchStatus("*",kFALSE);
  corrections->SetBranchStatus("detid",kTRUE);
  corrections->SetBranchStatus("correction",kTRUE);
  corrections->SetBranchAddress("correction",&correction);

  // create vectors for the different Profiles
  float yMin = 0, yMax = 0;
  int   nBinsY = 0;
  vector<double> delayBins;
  setLimitsAndBinning("delay",delayBins);
  setLimitsAndBinning(observable,yMin,yMax,nBinsY);
  
  cout<<"create profiles  "<<endl;

  // create vectors  std::cout<<"Tree with nEntries "<<tree->GetEntries()<<std::endl;
  long int iEvent = 0;
  for( ; iEvent < tree->GetEntries()/reductionFactor; iEvent++){    
    tree->GetEntry(iEvent);
    // take the map delay from the detid
    map->GetEntryWithIndex(detid);
    // take the correction from the detid
    corrections->GetEntryWithIndex(detid);

    cout.flush();
    if(iEvent % 100000 == 0) cout<<"\r"<<"iEvent "<<100*double(iEvent)/(tree->GetEntries()/reductionFactor)<<" % ";
    
    uint32_t subdetid    = int((detid-0x10000000)/0x2000000);
    float    R           = sqrt(clglobalX*clglobalX+clglobalY*clglobalY+clglobalZ*clglobalZ);
    float    value       = 0;
    if(observable == "maxCharge")
      value = obs*(clCorrectedSignalOverNoise)/(clSignalOverNoise);
    else
      value = obs;

    if(subdetid == 3){
      int ring = int((R-TIBRing.rMin)/((TIBRing.rMax-TIBRing.rMin)/TIBRing.nDivision))+1;   
      if(TIBrs[ring][round((delay-correction)*10)/10].get() == 0 or TIBrs[ring][round((delay-correction)*10)/10].get() == NULL)      {
	TIBrs[ring][round((delay-correction)*10)/10]= std::shared_ptr<TH1F>(new TH1F(Form("TIB_ring_%d_delay_%.1f",ring,delay-correction),"",nBinsY,yMin,yMax));      
      }
      TIBrs[ring][round((delay-correction)*10)/10]->Fill(value);
    }    
    else if(subdetid == 5){
      int ring = int((R-TOBRing.rMin)/((TOBRing.rMax-TOBRing.rMin)/TOBRing.nDivision))+1;
      if(TOBrs[ring][round((delay-correction)*10)/10].get() == 0 or TOBrs[ring][round((delay-correction)*10)/10].get() == NULL)      
	TOBrs[ring][round((delay-correction)*10)/10]= std::shared_ptr<TH1F>(new TH1F(Form("TOB_ring_%d_delay_%.1f",ring,delay-correction),"",nBinsY,yMin,yMax));      
      TOBrs[ring][round((delay-correction)*10)/10]->Fill(value);
      
    }
    else if(subdetid == 4){
      int ring = int((R-TOBRing.rMin)/((TOBRing.rMax-TOBRing.rMin)/TOBRing.nDivision))+1;
      if(TIDrs[ring][round((delay-correction)*10)/10].get() == 0 or TIDrs[ring][round((delay-correction)*10)/10].get() == NULL)      
	TIDrs[ring][round((delay-correction)*10)/10]= std::shared_ptr<TH1F>(new TH1F(Form("TID_ring_%d_delay_%.1f",ring,delay-correction),"",nBinsY,yMin,yMax));      
      TIDrs[ring][round((delay-correction)*10)/10]->Fill(value);
      
    }

    else if(subdetid == 6  and clglobalZ > 0 and thickness > 400){
      int ring = int((R-TECRing.rMin)/((TECRing.rMax-TECRing.rMin)/TECRing.nDivision))+1;
      if(TECPTrs[ring][round((delay-correction)*10)/10].get() == 0 or TECPTrs[ring][round((delay-correction)*10)/10].get() == NULL)      
	TECPTrs[ring][round((delay-correction)*10)/10]= std::shared_ptr<TH1F>(new TH1F(Form("TECPT_ring_%d_delay_%.1f",ring,delay-correction),"",nBinsY,yMin,yMax));      
      TECPTrs[ring][round((delay-correction)*10)/10]->Fill(value);

    }
    else if(subdetid == 6  and clglobalZ > 0 and thickness < 400){
      int ring = int((R-TECRing.rMin)/((TECRing.rMax-TECRing.rMin)/TECRing.nDivision))+1;
      if(TECPtrs[ring][round((delay-correction)*10)/10].get() == 0 or TECPtrs[ring][round((delay-correction)*10)/10].get() == NULL)      
	TECPtrs[ring][round((delay-correction)*10)/10]= std::shared_ptr<TH1F>(new TH1F(Form("TECPt_ring_%d_delay_%.1f",ring,delay-correction),"",nBinsY,yMin,yMax));      
      TECPtrs[ring][round((delay-correction)*10)/10]->Fill(value);

    }
    else if(subdetid == 6  and clglobalZ < 0 and thickness > 400){
      int ring = int((R-TECRing.rMin)/((TECRing.rMax-TECRing.rMin)/TECRing.nDivision))+1;
      if(TECMTrs[ring][round((delay-correction)*10)/10].get() == 0 or TECMTrs[ring][round((delay-correction)*10)/10].get() == NULL)      
	TECMTrs[ring][round((delay-correction)*10)/10]= std::shared_ptr<TH1F>(new TH1F(Form("TECMT_ring_%d_delay_%.1f",ring,delay-correction),"",nBinsY,yMin,yMax));      
      TECMTrs[ring][round((delay-correction)*10)/10]->Fill(value);


    }
    else if(subdetid == 6  and clglobalZ < 0 and thickness < 400){
      int ring = int((R-TECRing.rMin)/((TECRing.rMax-TECRing.rMin)/TECRing.nDivision))+1;
      if(TECMtrs[ring][round((delay-correction)*10)/10].get() == 0 or TECMtrs[ring][round((delay-correction)*10)/10].get() == NULL)      
	TECMtrs[ring][round((delay-correction)*10)/10]= std::shared_ptr<TH1F>(new TH1F(Form("TECMt_ring_%d_delay_%.1f",ring,delay-correction),"",nBinsY,yMin,yMax));      
      TECMtrs[ring][round((delay-correction)*10)/10]->Fill(value);

    }
  }

  std::cout<<"Build the mean and MPV distributions asaf of delay "<<endl;
  long int iBadChannelFit = 0;
  for(auto imap : TIBrs){ // loop on the different layer                                                                                                                     
    if(TIBrsMean[imap.first].get() == 0 or TIBrsMean[imap.first].get() == NULL)
      TIBrsMean[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TIB_ring_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    if(TIBrsMPV[imap.first].get() == 0 or TIBrsMPV[imap.first].get() == NULL)
      TIBrsMPV[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TIB_ring_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    
    for(auto idelay : imap.second){
      // mean value                                                                                                                                                             
      TIBrsMean[imap.first]->SetBinContent(TIBrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TIBrsMean[imap.first]->SetBinError(TIBrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());

      int status = makeLandauGausFit(idelay.second.get(),TIBrsMPV[imap.first].get(),"TIBrs",idelay.first,observable,outputDIR,"partitions");
      if(status != 1)
        iBadChannelFit++;
    }
  }

  std::cout<<"Bad channel fit from Landau+Gaus fit in TIB rings "<<iBadChannelFit<<" over "<<TIBrs.size()*TIBrs[1].size()<<std::endl;

  iBadChannelFit = 0;
  for(auto imap : TOBrs){ // loop on the different layer                                                                                                                     
    if(TOBrsMean[imap.first].get() == 0 or TOBrsMean[imap.first].get() == NULL)
      TOBrsMean[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TOB_ring_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    if(TOBrsMPV[imap.first].get() == 0 or TOBrsMPV[imap.first].get() == NULL)
      TOBrsMPV[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TOB_ring_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    
    for(auto idelay : imap.second){
      // mean value                                                                                                                                                             
      TOBrsMean[imap.first]->SetBinContent(TOBrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TOBrsMean[imap.first]->SetBinError(TOBrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());

      int status = makeLandauGausFit(idelay.second.get(),TOBrsMPV[imap.first].get(),"TOBrs",idelay.first,observable,outputDIR,"partitions");
      if(status != 1)
        iBadChannelFit++;
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TOB rings "<<iBadChannelFit<<" over "<<TOBrs.size()*TOBrs[1].size()<<std::endl;

  iBadChannelFit = 0;
  for(auto imap : TIDrs){ // loop on the different layer                                                                                                                     
    if(TIDrsMean[imap.first].get() == 0 or TIDrsMean[imap.first].get() == NULL)
      TIDrsMean[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TID_ring_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    if(TIDrsMPV[imap.first].get() == 0 or TIDrsMPV[imap.first].get() == NULL)
      TIDrsMPV[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TID_ring_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    
    for(auto idelay : imap.second){
      // mean value                                                                                                                                                             
      TIDrsMean[imap.first]->SetBinContent(TIDrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TIDrsMean[imap.first]->SetBinError(TIDrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());

      int status = makeLandauGausFit(idelay.second.get(),TIDrsMPV[imap.first].get(),"TIDrs",idelay.first,observable,outputDIR);
      if(status != 1)
        iBadChannelFit++;
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TID rings "<<iBadChannelFit<<" over "<<TIDrs.size()*TIDrs[1].size()<<std::endl;

  iBadChannelFit = 0;
  for(auto imap : TECPTrs){ // loop on the different layer                                                                                                                     
    if(TECPTrsMean[imap.first].get() == 0 or TECPTrsMean[imap.first].get() == NULL)
      TECPTrsMean[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TECPT_ring_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    if(TECPTrsMPV[imap.first].get() == 0 or TECPTrsMPV[imap.first].get() == NULL)
      TECPTrsMPV[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TECPT_ring_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    
    for(auto idelay : imap.second){
      // mean value                                                                                                                                                             
      TECPTrsMean[imap.first]->SetBinContent(TECPTrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TECPTrsMean[imap.first]->SetBinError(TECPTrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());

      int status = makeLandauGausFit(idelay.second.get(),TECPTrsMPV[imap.first].get(),"TECPTrs",idelay.first,observable,outputDIR,"partitions");
      if(status != 1)
        iBadChannelFit++;
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECPT rings "<<iBadChannelFit<<" over "<<TECPTrs.size()*TECPTrs[1].size()<<std::endl;

  iBadChannelFit = 0;
  for(auto imap : TECPtrs){ // loop on the different layer                                                                                                                     
    if(TECPtrsMean[imap.first].get() == 0 or TECPtrsMean[imap.first].get() == NULL)
      TECPtrsMean[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TECPt_ring_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    if(TECPtrsMPV[imap.first].get() == 0 or TECPtrsMPV[imap.first].get() == NULL)
      TECPtrsMPV[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TECPt_ring_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    
    for(auto idelay : imap.second){
      // mean value                                                                                                                                                             
      TECPtrsMean[imap.first]->SetBinContent(TECPtrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TECPtrsMean[imap.first]->SetBinError(TECPtrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());

      int status = makeLandauGausFit(idelay.second.get(),TECPtrsMPV[imap.first].get(),"TECPtrs",idelay.first,observable,outputDIR,"partitions");
      if(status != 1)
        iBadChannelFit++;
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECPt rings "<<iBadChannelFit<<" over "<<TECPtrs.size()*TECPtrs[1].size()<<std::endl;

  iBadChannelFit = 0;
  for(auto imap : TECMTrs){ // loop on the different layer                                                                                                                     
    if(TECMTrsMean[imap.first].get() == 0 or TECMTrsMean[imap.first].get() == NULL)
      TECMTrsMean[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TECMT_ring_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    if(TECMTrsMPV[imap.first].get() == 0 or TECMTrsMPV[imap.first].get() == NULL)
      TECMTrsMPV[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TECMT_ring_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    
    for(auto idelay : imap.second){
      // mean value                                                                                                                                                             
      TECMTrsMean[imap.first]->SetBinContent(TECMTrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TECMTrsMean[imap.first]->SetBinError(TECMTrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());

      int status = makeLandauGausFit(idelay.second.get(),TECMTrsMPV[imap.first].get(),"TECMTrs",idelay.first,observable,outputDIR,"partitions");
      if(status != 1)
        iBadChannelFit++;
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECMT rings "<<iBadChannelFit<<" over "<<TECMTrs.size()*TECMTrs[1].size()<<std::endl;

  iBadChannelFit = 0;
  for(auto imap : TECMtrs){ // loop on the different layer                                                                                                                     
    if(TECMtrsMean[imap.first].get() == 0 or TECMtrsMean[imap.first].get() == NULL)
      TECMtrsMean[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TECMt_ring_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    if(TECMtrsMPV[imap.first].get() == 0 or TECMtrsMPV[imap.first].get() == NULL)
      TECMtrsMPV[imap.first] = std::shared_ptr<TH1F>(new TH1F(Form("TECMt_ring_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]));
    TH1::AddDirectory(kFALSE);
    
    for(auto idelay : imap.second){
      // mean value                                                                                                                                                             
      TECMtrsMean[imap.first]->SetBinContent(TECMtrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TECMtrsMean[imap.first]->SetBinError(TECMtrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());

      int status = makeLandauGausFit(idelay.second.get(),TECMtrsMPV[imap.first].get(),"TECMtrs",idelay.first,observable,outputDIR,"partitions");
      if(status != 1)
        iBadChannelFit++;
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECMt rings "<<iBadChannelFit<<" over "<<TECMtrs.size()*TECMtrs[1].size()<<std::endl;
 
  
  std::cout<<std::endl;
  std::cout<<"Loop on events terminated"<<std::endl;

  // correct profiles and fit them with a gaussian
  std::cout<<"Analyze profiles TIB mean layers"<<std::endl;
  long int badFits = 0;
  for(auto iprof : TIBrsMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TIB Mean bad gaussian fits "<<badFits<<std::endl;
  std::cout<<"Analyze TOB profiles mean layers"<<std::endl;
  badFits = 0;
  for(auto iprof : TOBrsMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TOB Mean bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"Analyze TID profiles mean layers"<<std::endl;
  badFits = 0;
  for(auto iprof : TIDrsMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TID Mean bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"Analyze TECPT profiles mean layers"<<std::endl;
  badFits = 0;
  for(auto iprof : TECPTrsMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TECPT Mean bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"Analyze TECPt profiles mean layers"<<std::endl;
  badFits = 0;
  for(auto iprof : TECPtrsMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TECPt Mean bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"Analyze TECMT profiles mean layers"<<std::endl;
  badFits = 0;
  for(auto iprof : TECMTrsMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TECMT Mean bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"Analyze TECMt profiles mean layers"<<std::endl;
  badFits = 0;
  for(auto iprof : TECMtrsMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TECMt Mean bad gaussian fits "<<badFits<<std::endl;

  // correct profiles and fit them with a gaussian
  std::cout<<"Analyze TIB profiles mpv layers"<<std::endl;
  badFits = 0;
  for(auto iprof : TIBrsMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TIB MPV bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"Analyze TOB profiles mpv layers"<<std::endl;
  badFits = 0;
  for(auto iprof : TOBrsMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TOB MPV bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"Analyze TID profiles mpv layers"<<std::endl;
  badFits = 0;
  for(auto iprof : TIDrsMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TID MPV bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"Analyze TECPT profiles mpv layers"<<std::endl;
  badFits = 0;
  for(auto iprof : TECPTrsMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TECPT MPV bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"Analyze TECPt profiles mpv layers"<<std::endl;
  badFits = 0;
  for(auto iprof : TECPtrsMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TECPt MPV bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"Analyze TECMT profiles mpv layers"<<std::endl;
  badFits = 0;
  for(auto iprof : TECMTrsMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TECMT MPV bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"Analyze TECMt profiles mpv layers"<<std::endl;
  badFits = 0;
  for(auto iprof : TECMtrsMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
        badFits++;
      }
    }
  }
  std::cout<<"TECMt MPV bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"#################################"<<std::endl;
  std::cout<<"##### End of Ring Analysis #####"<<std::endl;
  std::cout<<"#################################"<<std::endl;
}


/// main function that run the analysis
void delayValidation(string file0,  // inputfile
		     string file1 = "nocorrection.root",  // possible file with correction
		     string observable   = "maxCharge",   // observable to be considered: maxCharge, S/N ..etc
		     bool plotPartitions = true, // best delay setting in each partition 
		     bool plotLayer      = true, // best delay setting per layer
		     bool plotSlices     = false, // best delay setting per ring
		     string outputDIR    = "prompt" // output directory name
		     ){

  // prepare style and load macros
  setTDRStyle();
  // not dump stat and fit info
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gROOT->SetBatch(kTRUE);

  system(("mkdir -p "+outputDIR).c_str());

  // open the file and prepare the cluster tree, adding the other trees as frined --> memory consuming
  std::cout<<"#############################"<<std::endl;
  std::cout<<"###### delayValidation ######"<<std::endl;
  std::cout<<"#############################"<<std::endl;

  std::cout<<"Open Input Files"<<std::endl;
  std::shared_ptr<TFile> _file0 (TFile::Open(file0.c_str()));
  std::shared_ptr<TFile> _file1 (TFile::Open(file1.c_str()));
  std::shared_ptr<TTree> clusters   ((TTree*)_file0->FindObjectAny("clusters"));
  std::shared_ptr<TTree> readoutMap ((TTree*)_file0->FindObjectAny("readoutMap"));
  std::shared_ptr<TTree> delayCorrections ((TTree*)_file1->FindObjectAny("delayCorrections"));  
  clusters->SetEventList(0);  

  // create the plots per layer
  if(plotLayer){

    // run per layer analysis
    LayerPlots(clusters,readoutMap,delayCorrections,observable,outputDIR);

    // canvas for layer one mean observable result
    std::shared_ptr<TCanvas> c1_mean = prepareCanvas("TIB_layers_mean",observable);
    plotAll(c1_mean,TIBlayersMean);
    c1_mean->Print(Form("%s/TIB_layers_mean.root",outputDIR.c_str()));
    // canvas for layer one MPV observable result 
    std::shared_ptr<TCanvas> c1_mpv = prepareCanvas("TIB_layers_mpv",observable);
    plotAll(c1_mpv,TIBlayersMPV);
    c1_mpv->Print(Form("%s/TIB_layers_mpv.root",outputDIR.c_str()));
    
    /// Layer TID
    std::shared_ptr<TCanvas> c2_mean = prepareCanvas("TID_layers_mean",observable);
    plotAll(c2_mean,TIDlayersMean);
    c2_mean->Print(Form("%s/TID_layers_mean.root",outputDIR.c_str()));
    std::shared_ptr<TCanvas> c2_mpv = prepareCanvas("TID_layers_mpv",observable);
    plotAll(c2_mpv,TIDlayersMPV);
    c2_mpv->Print(Form("%s/TID_layers_mpv.root",outputDIR.c_str()));
    
    /// Layer TOB
    std::shared_ptr<TCanvas> c3_mean = prepareCanvas("TOB_layers_mean",observable);
    plotAll(c3_mean,TOBlayersMean);
    c3_mean->Print(Form("%s/TOB_layers_mean.root",outputDIR.c_str()));
    std::shared_ptr<TCanvas> c3_mpv = prepareCanvas("TOB_layers_mpv",observable);
    plotAll(c3_mpv,TOBlayersMPV);
    c3_mpv->Print(Form("%s/TOB_layers_mpv.root",outputDIR.c_str()));

    /// Layer TECP thin sensors
    std::shared_ptr<TCanvas> c4_mean = prepareCanvas("TECPt_layers_mean",observable);
    plotAll(c4_mean,TECPtlayersMean);
    c4_mean->Print(Form("%s/TECPt_layers_mean.root",outputDIR.c_str()));
    std::shared_ptr<TCanvas> c4_mpv = prepareCanvas("TECPt_layers_mpv",observable);
    plotAll(c4_mpv,TECPtlayersMPV);
    c4_mpv->Print(Form("%s/TECPt_layers_mpv.root",outputDIR.c_str()));
    
    /// Layer TECP thick sensors
    std::shared_ptr<TCanvas> c5_mean = prepareCanvas("TECPT_layers_mean",observable);
    plotAll(c5_mean,TECPTlayersMean);
    c5_mean->Print(Form("%s/TECPT_layers_mean.root",outputDIR.c_str()));
    std::shared_ptr<TCanvas> c5_mpv = prepareCanvas("TECPT_layers_mpv",observable);
    plotAll(c5_mpv,TECPTlayersMPV);
    c5_mpv->Print(Form("%s/TECPT_layers_mpv.root",outputDIR.c_str()));

    /// Layer TECM thin sensors
    std::shared_ptr<TCanvas> c6_mean = prepareCanvas("TECMt_layers_mean",observable);
    plotAll(c6_mean,TECMtlayersMean);
    c6_mean->Print(Form("%s/TECMt_layers_mean.root",outputDIR.c_str()));
    std::shared_ptr<TCanvas> c6_mpv = prepareCanvas("TECMt_layers_mpv",observable);
    plotAll(c6_mpv,TECMtlayersMPV);
    c6_mpv->Print(Form("%s/TECMt_layers_mpv.root",outputDIR.c_str()));

    /// Layer TECM thick sensors
    std::shared_ptr<TCanvas> c7_mean = prepareCanvas("TECMT_layers_mean",observable);
    plotAll(c7_mean,TECMTlayersMean);
    c7_mean->Print(Form("%s/TECMT_layers_mean.root",outputDIR.c_str()));
    std::shared_ptr<TCanvas> c7_mpv = prepareCanvas("TECMT_layers_mpv",observable);
    plotAll(c7_mpv,TECMTlayersMPV);
    c7_mpv->Print(Form("%s/TECMT_layers_mpv.root",outputDIR.c_str()));

    std::vector<std::shared_ptr<TH1F> > alllayersMean;
    for(auto tib : TIBlayersMean)
      alllayersMean.push_back(tib.second);
    for(auto tob : TOBlayersMean)
      alllayersMean.push_back(tob.second);
    for(auto tid : TIDlayersMean)
      alllayersMean.push_back(tid.second);
    for(auto tec : TECPTlayersMean)
      alllayersMean.push_back(tec.second);
    for(auto tec : TECPtlayersMean)
      alllayersMean.push_back(tec.second);
    for(auto tec : TECMTlayersMean)
      alllayersMean.push_back(tec.second);
    for(auto tec : TECMtlayersMean)
      alllayersMean.push_back(tec.second);

    std::vector<std::shared_ptr<TH1F> > alllayersMPV;
    for(auto tib : TIBlayersMPV)
      alllayersMPV.push_back(tib.second);
    for(auto tob : TOBlayersMPV)
      alllayersMPV.push_back(tob.second);
    for(auto tid : TIDlayersMPV)
      alllayersMPV.push_back(tid.second);
    for(auto tec : TECPTlayersMPV)
      alllayersMPV.push_back(tec.second);
    for(auto tec : TECPtlayersMPV)
      alllayersMPV.push_back(tec.second);
    for(auto tec : TECMTlayersMPV)
      alllayersMPV.push_back(tec.second);
    for(auto tec : TECMtlayersMPV)
      alllayersMPV.push_back(tec.second);

    // store all the different plots --> plot maximum value for each layer    
    std::shared_ptr<TCanvas> c8_mean (new TCanvas("c_layers_mean","",800,650));  
    plotMaxima(c8_mean,alllayersMean,outputDIR,"layers_mean");
    std::shared_ptr<TCanvas> c8_mpv (new TCanvas("c_layers_mpv","",800,650));  
    plotMaxima(c8_mpv,alllayersMPV,outputDIR,"layers_mpv");

    TIBlayers.clear();
    TOBlayers.clear();
    TIDlayers.clear();
    TECPTlayers.clear();
    TECPtlayers.clear();
    TECMTlayers.clear();
    TECMtlayers.clear();

    TIBlayersMean.clear();
    TOBlayersMean.clear();
    TIDlayersMean.clear();
    TECPTlayersMean.clear();
    TECPtlayersMean.clear();
    TECMTlayersMean.clear();
    TECMtlayersMean.clear();
    alllayersMean.clear();
    
    TIBlayersMPV.clear();
    TOBlayersMPV.clear();
    TIDlayersMPV.clear();
    TECPTlayersMPV.clear();
    TECPtlayersMPV.clear();
    TECMTlayersMPV.clear();
    TECMtlayersMPV.clear();
    alllayersMPV.clear();
    
  }

  // Per rings
  if(plotSlices){

    // make the hisograms
    RPlots(clusters,readoutMap,delayCorrections,observable,outputDIR);
    
    // Per ring in TIB
    std::shared_ptr<TCanvas> c1b_mean = prepareCanvas("TIB_distance_mean",observable);
    plotAll(c1b_mean,TIBrsMean);
    c1b_mean->Print(Form("%s/TIB_distance_mean.root",outputDIR.c_str()));
    std::shared_ptr<TCanvas> c1b_mpv = prepareCanvas("TIB_distance_mpv",observable);
    plotAll(c1b_mpv,TIBrsMPV);
    c1b_mpv->Print(Form("%s/TIB_distance_mpv.root",outputDIR.c_str()));

    // Per ring in TID
    std::shared_ptr<TCanvas> c2b_mean = prepareCanvas("TID_distance_mean",observable);
    plotAll(c2b_mean,TIDrsMean);
    c2b_mean->Print(Form("%s/TID_distance_mean.root",outputDIR.c_str()));
    std::shared_ptr<TCanvas> c2b_mpv = prepareCanvas("TID_distance_mpv",observable);
    plotAll(c2b_mpv,TIDrsMPV);
    c2b_mpv->Print(Form("%s/TID_distance_mpv.root",outputDIR.c_str()));

    // Per ring in TOB
    std::shared_ptr<TCanvas> c3b_mean = prepareCanvas("TOB_distance_mean",observable);
    plotAll(c3b_mean,TOBrsMean);
    c3b_mean->Print(Form("%s/TOB_distance_mean.root",outputDIR.c_str()));
    std::shared_ptr<TCanvas> c3b_mpv = prepareCanvas("TOB_distance_mpv",observable);
    plotAll(c3b_mpv,TOBrsMPV);
    c3b_mpv->Print(Form("%s/TOB_distance_mpv.root",outputDIR.c_str()));

    // Per ring in TECM
    std::shared_ptr<TCanvas> c4b_mean = prepareCanvas("TECMT_distance_mean",observable);
    plotAll(c4b_mean,TECMTrsMean);
    c4b_mean->Print(Form("%s/TECMT_distance_mean.root",outputDIR.c_str()));
    std::shared_ptr<TCanvas> c4b_mpv = prepareCanvas("TECMT_distance_mpv",observable);
    plotAll(c4b_mpv,TECMTrsMPV);
    c4b_mpv->Print(Form("%s/TECMT_distance_mpv.root",outputDIR.c_str()));
    std::shared_ptr<TCanvas> c5b_mean = prepareCanvas("TECMt_distance_mean",observable);
    plotAll(c5b_mean,TECMtrsMean);
    c5b_mean->Print(Form("%s/TECMt_distance_mean.root",outputDIR.c_str()));
    std::shared_ptr<TCanvas> c5b_mpv = prepareCanvas("TECMt_distance_mpv",observable);
    plotAll(c5b_mpv,TECMtrsMPV);
    c5b_mpv->Print(Form("%s/TECMt_distance_mpv.root",outputDIR.c_str()));

    // Per ring in TECM
    std::shared_ptr<TCanvas> c6b_mean = prepareCanvas("TECPT_distance_mean",observable);
    plotAll(c6b_mean,TECPTrsMean);
    c6b_mean->Print(Form("%s/TECPT_distance_mean.root",outputDIR.c_str()));
    std::shared_ptr<TCanvas> c6b_mpv = prepareCanvas("TECPT_distance_mpv",observable);
    plotAll(c6b_mpv,TECPTrsMPV);
    c6b_mpv->Print(Form("%s/TECPT_distance_mpv.root",outputDIR.c_str()));
    std::shared_ptr<TCanvas> c7b_mean = prepareCanvas("TECPt_distance_mean",observable);
    plotAll(c7b_mean,TECPtrsMean);
    c7b_mean->Print(Form("%s/TECPt_distance_mean.root",outputDIR.c_str()));
    std::shared_ptr<TCanvas> c7b_mpv = prepareCanvas("TECPt_distance_mpv",observable);
    plotAll(c7b_mpv,TECPtrsMPV);
    c7b_mpv->Print(Form("%s/TECPt_distance_mpv.root",outputDIR.c_str()));

    std::vector<std::shared_ptr<TH1F> > allrsMean;
    for(auto tib : TIBrsMean)
      allrsMean.push_back(tib.second);
    for(auto tob : TOBrsMean)
      allrsMean.push_back(tob.second);
    for(auto tid : TIDrsMean)
      allrsMean.push_back(tid.second);
    for(auto tec : TECPTrsMean)
      allrsMean.push_back(tec.second);
    for(auto tec : TECPtrsMean)
      allrsMean.push_back(tec.second);
    for(auto tec : TECMTrsMean)
      allrsMean.push_back(tec.second);
    for(auto tec : TECMtrsMean)
      allrsMean.push_back(tec.second);

    std::vector<std::shared_ptr<TH1F> > allrsMPV;
    for(auto tib : TIBrsMPV)
      allrsMPV.push_back(tib.second);
    for(auto tob : TOBrsMPV)
      allrsMPV.push_back(tob.second);
    for(auto tid : TIDrsMPV)
      allrsMPV.push_back(tid.second);
    for(auto tec : TECPTrsMPV)
      allrsMPV.push_back(tec.second);
    for(auto tec : TECPtrsMPV)
      allrsMPV.push_back(tec.second);
    for(auto tec : TECMTrsMPV)
      allrsMPV.push_back(tec.second);
    for(auto tec : TECMtrsMPV)
      allrsMPV.push_back(tec.second);

    std::shared_ptr<TCanvas> c8b_mean (new TCanvas("c_rings_mean","",800,650));  
    plotMaxima(c8b_mean,allrsMean,outputDIR,"rings");
    std::shared_ptr<TCanvas> c8b_mpv (new TCanvas("c_rings_mpv","",800,650));  
    plotMaxima(c8b_mpv,allrsMPV,outputDIR,"rings");
    
    TIBrs.clear();
    TIDrs.clear();
    TOBrs.clear();
    TECPTrs.clear();
    TECPtrs.clear();
    TECMTrs.clear();
    TECMtrs.clear();

    TIBrsMean.clear();
    TIDrsMean.clear();
    TOBrsMean.clear();
    TECPTrsMean.clear();
    TECPtrsMean.clear();
    TECMTrsMean.clear();
    TECMtrsMean.clear();
    allrsMean.clear();

    TIBrsMPV.clear();
    TIDrsMPV.clear();
    TOBrsMPV.clear();
    TECPTrsMPV.clear();
    TECPtrsMPV.clear();
    TECMTrsMPV.clear();
    TECMtrsMPV.clear();
    allrsMPV.clear();
    
  }

  // Plot per partition 
  if(plotPartitions){

    //change ring definition to collapse all of them
    TIBRing.nDivision = 1;
    TIDRing.nDivision = 1;
    TOBRing.nDivision = 1;
    TECRing.nDivision = 1;
    
    // cumulate all rings
    RPlots(clusters,readoutMap,delayCorrections,observable,outputDIR);
        
    // create the plots per partition
    std::vector<std::shared_ptr<TH1F> > allPartitionMean;
    for(auto tib : TIBrsMean)
      allPartitionMean.push_back(tib.second);
    for(auto tob : TOBrsMean)
      allPartitionMean.push_back(tob.second);
    for(auto tid : TIDrsMean)
      allPartitionMean.push_back(tid.second);
    for(auto tec : TECPTrsMean)
      allPartitionMean.push_back(tec.second);
    for(auto tec : TECPtrsMean)
      allPartitionMean.push_back(tec.second);
    for(auto tec : TECMTrsMean)
      allPartitionMean.push_back(tec.second);
    for(auto tec : TECMtrsMean)
      allPartitionMean.push_back(tec.second);

    std::vector<std::shared_ptr<TH1F> > allPartitionMPV;
    for(auto tib : TIBrsMPV)
      allPartitionMPV.push_back(tib.second);
    for(auto tob : TOBrsMPV)
      allPartitionMPV.push_back(tob.second);
    for(auto tid : TIDrsMPV)
      allPartitionMPV.push_back(tid.second);
    for(auto tec : TECPTrsMPV)
      allPartitionMPV.push_back(tec.second);
    for(auto tec : TECPtrsMPV)
      allPartitionMPV.push_back(tec.second);
    for(auto tec : TECMTrsMPV)
      allPartitionMPV.push_back(tec.second);
    for(auto tec : TECMtrsMPV)
      allPartitionMPV.push_back(tec.second);


    std::shared_ptr<TCanvas> c1_mean = prepareCanvas("Partitions_mean",observable);
    plotAll(c1_mean,allPartitionMean,"ring1");
    c1_mean->Print(Form("%s/Partitions_mean.root",outputDIR.c_str()));

    std::shared_ptr<TCanvas> c1_mpv = prepareCanvas("Partitions_mpv",observable);
    plotAll(c1_mpv,allPartitionMPV,"ring1");
    c1_mpv->Print(Form("%s/Partitions_mpv.root",outputDIR.c_str()));

    TIBrs.clear();
    TIDrs.clear();
    TOBrs.clear();
    TECPTrs.clear();
    TECPtrs.clear();
    TECMTrs.clear();
    TECMtrs.clear();

    TIBrsMean.clear();
    TIDrsMean.clear();
    TOBrsMean.clear();
    TECPTrsMean.clear();
    TECPtrsMean.clear();
    TECMTrsMean.clear();
    TECMtrsMean.clear();

    TIBrsMPV.clear();
    TIDrsMPV.clear();
    TOBrsMPV.clear();
    TECPTrsMPV.clear();
    TECPtrsMPV.clear();
    TECMTrsMPV.clear();
    TECMtrsMPV.clear();

    allPartitionMean.clear();
    allPartitionMPV.clear();
  } 
}

