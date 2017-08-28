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

int makeLandauGausFit(TH1F* histoToFit, string subdetector, const float & delay, const string & observable){

  float yMin   = 0, yMax = 0;
  int   nBinsY = 0;
  setLimitsAndBinning(observable,yMin,yMax,nBinsY);

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
  
  // make fit and get parameters                                                                                                                                           
  TFitResultPtr fitResult = histoToFit->Fit(fitfunc,"RSQ");
  if(not fitResult.Get() or fitResult->Status() != 0 or fitResult->CovMatrixStatus() <= 1){
    return -1;
  }
  else{
    fit_parameters = fitResult->GetParams();
    fit_parameters_error = fitResult->GetErrors();
    chi2   = fitResult->Chi2();
    ndf    = fitResult->Ndf();
    return fitResult->Status();
  }
}



/// function that runs on the evnet and produce profiles for layers
void LayerPlots(const std::shared_ptr<TTree> & tree, 
		const std::shared_ptr<TTree> & map,
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
      if(TIBlayers[barrellayer][delay].get() == 0 or TIBlayers[barrellayer][delay].get() == NULL) 
	TIBlayers[barrellayer][delay] = std::shared_ptr<TH1F> (new TH1F(Form("TIB_layer_%d_delay_%.1f",barrellayer,delay),"",nBinsY,yMin,yMax));
      TH1::AddDirectory(kFALSE);
      TIBlayers[barrellayer][delay]->Fill(value);
    }
    else if(subdetid == 5){
      if(TOBlayers[barrellayer][delay].get() == 0 or TOBlayers[barrellayer][delay].get() == NULL) 
	TOBlayers[barrellayer][delay] = std::shared_ptr<TH1F> (new TH1F(Form("TOB_layer_%d_delay_%.1f",barrellayer,delay),"",nBinsY,yMin,yMax));
      TH1::AddDirectory(kFALSE);
      TOBlayers[barrellayer][delay]->Fill(value);
    }
    else if(subdetid == 4){
      if(TIDlayers[TIDlayer][delay].get() == 0 or TIDlayers[TIDlayer][delay].get() == NULL) 
	TIDlayers[TIDlayer][delay] = std::shared_ptr<TH1F> (new TH1F(Form("TID_layer_%d_delay_%.1f",TIDlayer,delay),"",nBinsY,yMin,yMax));
      TH1::AddDirectory(kFALSE);
      TIDlayers[TIDlayer][delay]->Fill(value);
    }
    else if(subdetid == 6  and clglobalZ > 0 and thickness > 400){
      if(TECPTlayers[TECPlayer][delay].get() == 0 or TECPTlayers[TECPlayer][delay].get() == NULL) 
	TECPTlayers[TECPlayer][delay] = std::shared_ptr<TH1F> (new TH1F(Form("TECPT_layer_%d_delay_%.1f",TECPlayer,delay),"",nBinsY,yMin,yMax));
      TH1::AddDirectory(kFALSE);
      TECPTlayers[TECPlayer][delay]->Fill(value);
    }
    else if(subdetid == 6  and clglobalZ > 0 and thickness < 400){
      if(TECPtlayers[TECPlayer][delay].get() == 0 or TECPtlayers[TECPlayer][delay].get() == NULL) 
	TECPtlayers[TECPlayer][delay] = std::shared_ptr<TH1F> (new TH1F(Form("TECPt_layer_%d_delay_%.1f",TECPlayer,delay),"",nBinsY,yMin,yMax));
      TH1::AddDirectory(kFALSE);
      TECPtlayers[TECPlayer][delay]->Fill(value);
    }
    else if(subdetid == 6  and clglobalZ < 0 and thickness > 400){
      if(TECMTlayers[TECMlayer][delay].get() == 0 or TECMTlayers[TECMlayer][delay].get() == NULL) 
	TECMTlayers[TECMlayer][delay] = std::shared_ptr<TH1F> (new TH1F(Form("TECMT_layer_%d_delay_%.1f",TECMlayer,delay),"",nBinsY,yMin,yMax));
      TH1::AddDirectory(kFALSE);
      TECMTlayers[TECMlayer][delay]->Fill(value);
    }
    else if(subdetid == 6  and clglobalZ < 0 and thickness < 400){
      if(TECMtlayers[TECMlayer][delay].get() == 0 or TECMtlayers[TECMlayer][delay].get() == NULL) 
	TECMtlayers[TECMlayer][delay] = std::shared_ptr<TH1F> (new TH1F(Form("TECMt_layer_%d_delay_%.1f",TECMlayer,delay),"",nBinsY,yMin,yMax));
      TH1::AddDirectory(kFALSE);
      TECMtlayers[TECMlayer][delay]->Fill(value);
    }
  }
  
  std::cout<<std::endl;
  std::cout<<"Loop on events terminated"<<std::endl;

  std::cout<<"Build the mean and MPV distributions asaf of delay "<<endl;
  long int iBadChannelFit = 0;
  for(auto imap : TIBlayers){ // loop on the different layer
    for(auto idelay : imap.second){
      int status = makeLandauGausFit(idelay.second.get(),"TIB",idelay.first,observable);
      if(status != 0)
	iBadChannelFit++;                  
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TIB layers "<<iBadChannelFit<<" over "<<TIBlayers.size()*TIBlayers[1].size()<<std::endl;  

  iBadChannelFit = 0;
  for(auto imap : TOBlayers){ // loop on the different layer
    for(auto idelay : imap.second){
      int status = makeLandauGausFit(idelay.second.get(),"TOB",idelay.first,observable);
      if(status != 0)
	iBadChannelFit++;      
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TOB layers "<<iBadChannelFit<<" over "<<TOBlayers.size()*TOBlayers[1].size()<<std::endl;  

  iBadChannelFit = 0;
  for(auto imap : TIDlayers){ // loop on the different layer
    for(auto idelay : imap.second){
      int status = makeLandauGausFit(idelay.second.get(),"TID",idelay.first,observable);
      if(status != 0)
	iBadChannelFit++;      
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TID layers "<<iBadChannelFit<<" over "<<TIDlayers.size()*TIDlayers[1].size()<<std::endl;  

  iBadChannelFit = 0;
  for(auto imap : TECPTlayers){ // loop on the different layer
    for(auto idelay : imap.second){
      int status = makeLandauGausFit(idelay.second.get(),"TECOuter",idelay.first,observable);
      if(status != 0)
	iBadChannelFit++;      
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECPT layers "<<iBadChannelFit<<" over "<<TECPTlayers.size()*TECPTlayers[1].size()<<std::endl;  
  
  iBadChannelFit = 0;
  for(auto imap : TECPtlayers){ // loop on the different layer
    for(auto idelay : imap.second){
      int status = makeLandauGausFit(idelay.second.get(),"TECInner",idelay.first,observable);
      if(status != 0)
	iBadChannelFit++;      
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECPt layers "<<iBadChannelFit<<" over "<<TECPtlayers.size()*TECPtlayers[1].size()<<std::endl;  

  iBadChannelFit = 0;
  for(auto imap : TECMTlayers){ // loop on the different layer
    for(auto idelay : imap.second){
      int status = makeLandauGausFit(idelay.second.get(),"TECOuter",idelay.first,observable);
      if(status != 0)
	iBadChannelFit++;      
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECMT layers "<<iBadChannelFit<<" over "<<TECMTlayers.size()*TECMTlayers[1].size()<<std::endl;  

  iBadChannelFit = 0;
  for(auto imap : TECMtlayers){ // loop on the different layer
    for(auto idelay : imap.second){
      int status = makeLandauGausFit(idelay.second.get(),"TECInner",idelay.first,observable);
      if(status != 0)
	iBadChannelFit++;      
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECMt layers "<<iBadChannelFit<<" over "<<TECMtlayers.size()*TECMtlayers[1].size()<<std::endl;  

  std::cout<<"#################################"<<std::endl;
  std::cout<<"##### End of Layer Analysis #####"<<std::endl;
  std::cout<<"#################################"<<std::endl;

}

////                                                                                                                                                                                                   
void plotDistributions(TCanvas* canvas, const string & outputDIR){

  // compare all layers in TIB
  vector<float> delayVal;
  for(auto ientry : TIBlayers){
    for(auto imap : ientry.second)
      delayVal.push_back(imap.first);			
    break;  
  }
  
  TH1F* frame = NULL;    
  // make histo for each delay value
  for(size_t idelay = 0; idelay < delayVal.size(); idelay++){

    vector<TH1F*> histoToPlot;
    // loop on the TIB mape
    for(auto ientry : TIBlayers){
      for(auto imap : ientry.second){
	if(imap.first != delayVal.at(idelay)) continue;
	histoToPlot.push_back(imap.second.get());
      }
    }

    if(frame == 0 or frame == NULL){
      frame = (TH1F*) histoToPlot.at(0)->Clone("frame");
      frame->Reset();
      frame->GetXaxis()->SetTitle("leading strip charge (ADC)");
      frame->GetXaxis()->SetTitleOffset(1.1);
      frame->GetYaxis()->SetTitle("Number of clusters");
      frame->GetYaxis()->SetTitleOffset(1.35);
    }
    frame->Draw();

    CMS_lumi(canvas,"",false,false,0.4);

    TLegend leg (0.55,0.56,0.82,0.82);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.AddEntry((TObject*)(0),"2017 Data","");

    int icolor = 1;
    double max = 0;
    for(auto histo : histoToPlot){
      histo->SetLineColor(icolor);
      histo->SetMarkerColor(icolor);
      histo->SetMarkerSize(1);
      histo->SetMarkerStyle(20);
      histo->SetLineWidth(2);
      ((TH1F*) histo->GetListOfFunctions()->At(0))->SetLineColor(icolor);
      ((TH1F*) histo->GetListOfFunctions()->At(0))->SetLineWidth(2);
      histo->Draw("EPsame");
      if(histo->GetMaximum() > max) max = histo->GetMaximum();
      leg.AddEntry(histo,Form("TIB layer %d",icolor),"EP");
      icolor++;
    }
    frame->GetYaxis()->SetRangeUser(0,max*1.5);
    leg.Draw("same");
    canvas->SaveAs((outputDIR+"/clusterCharge_TIB_delay_"+to_string(delayVal.at(idelay))+".png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/clusterCharge_TIB_delay_"+to_string(delayVal.at(idelay))+".pdf").c_str(),"pdf");
  }   

  // make histo for each delay value
  for(size_t idelay = 0; idelay < delayVal.size(); idelay++){
    vector<TH1F*> histoToPlot;
    // loop on the TOB mape
    for(auto ientry : TOBlayers){
      frame->Draw();
      for(auto imap : ientry.second){
	if(imap.first != delayVal.at(idelay)) continue;
	histoToPlot.push_back(imap.second.get());
      }
    }

    CMS_lumi(canvas,"",false,false,0.4);

    TLegend leg (0.55,0.56,0.82,0.82);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.AddEntry((TObject*)(0),"2017 Data","");

    int icolor = 1;
    double max = 0;
    for(auto histo : histoToPlot){
      histo->SetLineColor(icolor);
      histo->SetMarkerColor(icolor);
      histo->SetLineWidth(2);
      histo->SetMarkerSize(1);
      histo->SetMarkerStyle(20);
      ((TH1F*) histo->GetListOfFunctions()->At(0))->SetLineColor(icolor);
      ((TH1F*) histo->GetListOfFunctions()->At(0))->SetLineWidth(2);
      histo->Draw("EPsame");
      if(histo->GetMaximum() > max) max = histo->GetMaximum();
      leg.AddEntry(histo,Form("TOB layer %d",icolor),"EP");
      icolor++;
    }
    leg.Draw("same");
    frame->GetYaxis()->SetRangeUser(0,max*1.5);
    canvas->SaveAs((outputDIR+"/clusterCharge_TOB_delay_"+to_string(delayVal.at(idelay))+".png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/clusterCharge_TOB_delay_"+to_string(delayVal.at(idelay))+".pdf").c_str(),"pdf");    
  }   


  // make histo for each delay value
  for(size_t idelay = 0; idelay < delayVal.size(); idelay++){
    vector<TH1F*> histoToPlot;
    // loop on the TIB mape
    for(auto ientry : TIDlayers){
      frame->Draw();
      for(auto imap : ientry.second){
	if(imap.first != delayVal.at(idelay)) continue;
	histoToPlot.push_back(imap.second.get());
      }
    }

    CMS_lumi(canvas,"",false,false,0.4);

    TLegend leg (0.55,0.56,0.82,0.82);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.AddEntry((TObject*)(0),"2017 Data","");

    int icolor = 1;
    double max = 0;
    for(auto histo : histoToPlot){
      histo->SetLineColor(icolor);
      histo->SetMarkerColor(icolor);
      histo->SetLineWidth(2);
      histo->SetMarkerSize(1);
      histo->SetMarkerStyle(20);
      ((TH1F*) histo->GetListOfFunctions()->At(0))->SetLineColor(icolor);
      ((TH1F*) histo->GetListOfFunctions()->At(0))->SetLineWidth(2);
      histo->Draw("EPsame");
      if(histo->GetMaximum() > max) max = histo->GetMaximum();
      leg.AddEntry(histo,Form("TID layer %d",icolor),"EP");
      icolor++;
    }
    leg.Draw("same");
    frame->GetYaxis()->SetRangeUser(0,max*1.5);
    canvas->SaveAs((outputDIR+"/clusterCharge_TID_delay_"+to_string(delayVal.at(idelay))+".png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/clusterCharge_TID_delay_"+to_string(delayVal.at(idelay))+".pdf").c_str(),"pdf");    
  }   
}


/// main function that run the analysis
void makeChargeDistributionPerLayer(string file0,  // inputfile
				    string observable   = "maxCharge",   // observable to be considered: maxCharge, S/N ..etc
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
  std::cout<<"############################################"<<std::endl;
  std::cout<<"###### makeChargeDistributionPerLayer ######"<<std::endl;
  std::cout<<"############################################"<<std::endl;

  std::cout<<"Open Input Files"<<std::endl;
  std::shared_ptr<TFile> _file0 (TFile::Open(file0.c_str()));
  std::shared_ptr<TTree> clusters   ((TTree*)_file0->FindObjectAny("clusters"));
  std::shared_ptr<TTree> readoutMap ((TTree*)_file0->FindObjectAny("readoutMap"));
  clusters->SetEventList(0);  

  // run per layer analysis
  LayerPlots(clusters,readoutMap,observable,outputDIR);
  TCanvas* canvas = new TCanvas("canvas","",600,650);
  canvas->cd();
  plotDistributions(canvas,outputDIR);
  
}

