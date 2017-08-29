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

static std::map<uint32_t,std::map<float,TH1F*> > TIBlayers; // map layer:delay:distribution
static std::map<uint32_t,std::map<float,TH1F*> > TOBlayers;
static std::map<uint32_t,std::map<float,TH1F*> > TIDlayers;
static std::map<uint32_t,std::map<float,TH1F*> > TECPTlayers;
static std::map<uint32_t,std::map<float,TH1F*> > TECPtlayers;
static std::map<uint32_t,std::map<float,TH1F*> > TECMTlayers;
static std::map<uint32_t,std::map<float,TH1F*> > TECMtlayers;

static std::map<uint32_t,std::map<float,TF1*> > fitTIBlayers; // map layer:delay:distribution
static std::map<uint32_t,std::map<float,TF1*> > fitTOBlayers;
static std::map<uint32_t,std::map<float,TF1*> > fitTIDlayers;
static std::map<uint32_t,std::map<float,TF1*> > fitTECPTlayers;
static std::map<uint32_t,std::map<float,TF1*> > fitTECPtlayers;
static std::map<uint32_t,std::map<float,TF1*> > fitTECMTlayers;
static std::map<uint32_t,std::map<float,TF1*> > fitTECMtlayers;

static std::map<uint32_t,std::map<float,RooAbsPdf*> > pdfTIBlayers; // map layer:delay:distribution
static std::map<uint32_t,std::map<float,RooAbsPdf*> > pdfTOBlayers;
static std::map<uint32_t,std::map<float,RooAbsPdf*> > pdfTIDlayers;
static std::map<uint32_t,std::map<float,RooAbsPdf*> > pdfTECPTlayers;
static std::map<uint32_t,std::map<float,RooAbsPdf*> > pdfTECPtlayers;
static std::map<uint32_t,std::map<float,RooAbsPdf*> > pdfTECMTlayers;
static std::map<uint32_t,std::map<float,RooAbsPdf*> > pdfTECMtlayers;

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



/// function that runs on the evnet and produce profiles for layers
void LayerPlots(TTree* tree, 
		TTree* map,
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
      if(TIBlayers[barrellayer][delay] == 0 or TIBlayers[barrellayer][delay] == NULL) 
	TIBlayers[barrellayer][delay] = new TH1F(Form("TIB_layer_%d_delay_%.1f",barrellayer,delay),"",nBinsY,yMin,yMax);
      TH1::AddDirectory(kFALSE);
      TIBlayers[barrellayer][delay]->Fill(value);
    }
    else if(subdetid == 5){
      if(TOBlayers[barrellayer][delay] == 0 or TOBlayers[barrellayer][delay] == NULL) 
	TOBlayers[barrellayer][delay] = new TH1F(Form("TOB_layer_%d_delay_%.1f",barrellayer,delay),"",nBinsY,yMin,yMax);
      TH1::AddDirectory(kFALSE);
      TOBlayers[barrellayer][delay]->Fill(value);
    }
    else if(subdetid == 4){
      if(TIDlayers[TIDlayer][delay] == 0 or TIDlayers[TIDlayer][delay] == NULL) 
	TIDlayers[TIDlayer][delay] = new TH1F(Form("TID_layer_%d_delay_%.1f",TIDlayer,delay),"",nBinsY,yMin,yMax);
      TH1::AddDirectory(kFALSE);
      TIDlayers[TIDlayer][delay]->Fill(value);
    }
    else if(subdetid == 6  and clglobalZ > 0 and thickness > 400){
      if(TECPTlayers[TECPlayer][delay] == 0 or TECPTlayers[TECPlayer][delay] == NULL) 
	TECPTlayers[TECPlayer][delay] = new TH1F(Form("TECPT_layer_%d_delay_%.1f",TECPlayer,delay),"",nBinsY,yMin,yMax);
      TH1::AddDirectory(kFALSE);
      TECPTlayers[TECPlayer][delay]->Fill(value);
    }
    else if(subdetid == 6  and clglobalZ > 0 and thickness < 400){
      if(TECPtlayers[TECPlayer][delay] == 0 or TECPtlayers[TECPlayer][delay] == NULL) 
	TECPtlayers[TECPlayer][delay] = new TH1F(Form("TECPt_layer_%d_delay_%.1f",TECPlayer,delay),"",nBinsY,yMin,yMax);
      TH1::AddDirectory(kFALSE);
      TECPtlayers[TECPlayer][delay]->Fill(value);
    }
    else if(subdetid == 6  and clglobalZ < 0 and thickness > 400){
      if(TECMTlayers[TECMlayer][delay] == 0 or TECMTlayers[TECMlayer][delay] == NULL) 
	TECMTlayers[TECMlayer][delay] = new TH1F(Form("TECMT_layer_%d_delay_%.1f",TECMlayer,delay),"",nBinsY,yMin,yMax);
      TH1::AddDirectory(kFALSE);
      TECMTlayers[TECMlayer][delay]->Fill(value);
    }
    else if(subdetid == 6  and clglobalZ < 0 and thickness < 400){
      if(TECMtlayers[TECMlayer][delay] == 0 or TECMtlayers[TECMlayer][delay] == NULL) 
	TECMtlayers[TECMlayer][delay] = new TH1F(Form("TECMt_layer_%d_delay_%.1f",TECMlayer,delay),"",nBinsY,yMin,yMax);
      TH1::AddDirectory(kFALSE);
      TECMtlayers[TECMlayer][delay]->Fill(value);
    }
  }

  std::cout<<std::endl;
  std::cout<<"Loop on events terminated"<<std::endl;

  std::cout<<"Build the mean and MPV distributions asaf of delay "<<endl;
  RooAbsPdf* pdf = NULL;
  TF1* fitfunc = NULL;
  pair<float,RooAbsPdf*> pair;
  long int iBadChannelFit = 0;
  for(auto imap : TIBlayers){ // loop on the different layer
    for(auto idelay : imap.second){
      int status = -1;
      fitfunc = makeLandauGausFit(idelay.second,status,"TIB_L"+to_string(imap.first),idelay.first,observable,pair);      
      pdfTIBlayers[imap.first][idelay.first] = pair.second;
      fitTIBlayers[imap.first][idelay.first] = fitfunc;      
      if(status != 0) iBadChannelFit++;                  
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TIB layers "<<iBadChannelFit<<" over "<<TIBlayers.size()*TIBlayers[1].size()<<std::endl;  

  iBadChannelFit = 0;
  for(auto imap : TOBlayers){ // loop on the different layer
    for(auto idelay : imap.second){
      int status = -1;
      fitfunc = makeLandauGausFit(idelay.second,status,"TOB_L"+to_string(imap.first),idelay.first,observable,pair);
      pdfTOBlayers[imap.first][idelay.first] = pair.second;
      fitTOBlayers[imap.first][idelay.first] = fitfunc;
      if(status != 0) iBadChannelFit++;      
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TOB layers "<<iBadChannelFit<<" over "<<TOBlayers.size()*TOBlayers[1].size()<<std::endl;  

  iBadChannelFit = 0;
  for(auto imap : TIDlayers){ // loop on the different layer
    for(auto idelay : imap.second){
      int status = -1;
      fitfunc = makeLandauGausFit(idelay.second,status,"TID_L"+to_string(imap.first),idelay.first,observable,pair);
      pdfTIDlayers[imap.first][idelay.first] = pair.second;
      fitTIDlayers[imap.first][idelay.first] = fitfunc;      
      if(status != 0) iBadChannelFit++;      
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TID layers "<<iBadChannelFit<<" over "<<TIDlayers.size()*TIDlayers[1].size()<<std::endl;  

  iBadChannelFit = 0;
  for(auto imap : TECPTlayers){ // loop on the different layer
    for(auto idelay : imap.second){
      int status = -1;
      fitfunc = makeLandauGausFit(idelay.second,status,"TECPThick_L"+to_string(imap.first),idelay.first,observable,pair);
      pdfTECPTlayers[imap.first][idelay.first] = pair.second;
      fitTECPTlayers[imap.first][idelay.first] = fitfunc;      
      if(status != 0) iBadChannelFit++;      
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECPT layers "<<iBadChannelFit<<" over "<<TECPTlayers.size()*TECPTlayers[1].size()<<std::endl;  
  
  iBadChannelFit = 0;
  for(auto imap : TECPtlayers){ // loop on the different layer
    for(auto idelay : imap.second){
      int status = -1;
      fitfunc = makeLandauGausFit(idelay.second,status,"TECPthin_L"+to_string(imap.first),idelay.first,observable,pair);
      pdfTECPtlayers[imap.first][idelay.first] = pair.second;
      fitTECPtlayers[imap.first][idelay.first] = fitfunc;      
      if(status != 0) iBadChannelFit++;      
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECPt layers "<<iBadChannelFit<<" over "<<TECPtlayers.size()*TECPtlayers[1].size()<<std::endl;  

  iBadChannelFit = 0;
  for(auto imap : TECMTlayers){ // loop on the different layer
    for(auto idelay : imap.second){
      int status = -1;
      fitfunc = makeLandauGausFit(idelay.second,status,"TECMThick_L"+to_string(imap.first),idelay.first,observable,pair);
      pdfTECMTlayers[imap.first][idelay.first] = pair.second;
      fitTECMTlayers[imap.first][idelay.first] = fitfunc;      
      if(status != 0) iBadChannelFit++;      
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECMT layers "<<iBadChannelFit<<" over "<<TECMTlayers.size()*TECMTlayers[1].size()<<std::endl;  

  iBadChannelFit = 0;
  for(auto imap : TECMtlayers){ // loop on the different layer
    for(auto idelay : imap.second){
      int status = -1;
      fitfunc = makeLandauGausFit(idelay.second,status,"TECMthin_L"+to_string(imap.first),idelay.first,observable,pair);
      pdfTECMtlayers[imap.first][idelay.first] = pair.second;
      fitTECMtlayers[imap.first][idelay.first] = fitfunc;      
      if(status != 0) iBadChannelFit++;      
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECMt layers "<<iBadChannelFit<<" over "<<TECMtlayers.size()*TECMtlayers[1].size()<<std::endl;  

  std::cout<<"#################################"<<std::endl;
  std::cout<<"##### End of Layer Analysis #####"<<std::endl;
  std::cout<<"#################################"<<std::endl;

}

//// ---------
void plotHisto(TCanvas* canvas, const vector<TH1F*> & histoToPlot, const vector<TF1*> & funcToPlot, const string & outputDIR, const float & delayVal, const string & postfix){

  TH1F* frame = (TH1F*) histoToPlot.at(0)->Clone("frame");    
  frame->Reset();
  frame->GetXaxis()->SetTitle("leading strip charge (ADC)");
  frame->GetXaxis()->SetTitleOffset(1.1);
  frame->GetYaxis()->SetTitle("a.u.");
  frame->GetYaxis()->SetTitleOffset(1.35);
  frame->Draw();
  CMS_lumi(canvas,"",false,false,0.4);

  TLegend leg (0.55,0.56,0.82,0.82);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry((TObject*)(0),"2017 Data","");
  
  int icolor = 1;
  for(auto func : funcToPlot){
    func->SetLineColor(icolor);
    func->SetLineWidth(2);
    func->Draw("Lsame");
    icolor++;
  }

  icolor = 1;
  double max = 0;

  for(auto histo : histoToPlot){
    histo->Scale(1./histo->Integral(),"width");
    histo->SetLineColor(icolor);
    histo->SetMarkerColor(icolor);
    histo->SetMarkerSize(1);
    histo->SetMarkerStyle(20);
    histo->SetLineWidth(2);
    histo->Draw("EPsame");
    if(histo->GetMaximum() > max) max = histo->GetMaximum();
    leg.AddEntry(histo,Form("TIB layer %d",icolor),"EP");
    icolor++;
  }
  

  frame->GetYaxis()->SetRangeUser(0,max*1.5);
  leg.Draw("same");

  canvas->SaveAs((outputDIR+"/clusterCharge_"+postfix+"_delay_"+to_string(delayVal)+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/clusterCharge_"+postfix+"_delay_"+to_string(delayVal)+".pdf").c_str(),"pdf");

  if(frame) delete frame;
}

////                                                                                                                                                                                                  
void plotDistributions(TCanvas* canvas, const string & outputDIR){

  cout<<"##### Plot the final distributions #####"<<endl;

  // compare all layers in TIB
  vector<float> delayVal;
  for(auto ientry : TIBlayers){
    for(auto imap : ientry.second)
      delayVal.push_back(imap.first);			
    break;  
  }
  
  // make histo for each delay value
  for(size_t idelay = 0; idelay < delayVal.size(); idelay++){
    vector<TH1F*> histoToPlot;
    vector<TF1*> funcToPlot;
    // loop on the TIB mape
    for(auto ientry : TIBlayers){
      for(auto imap : ientry.second){
	if(imap.first != delayVal.at(idelay)) continue;
	histoToPlot.push_back(imap.second);
	funcToPlot.push_back(fitTIBlayers[ientry.first][imap.first]);
      }
    }
    plotHisto(canvas,histoToPlot,funcToPlot,outputDIR,delayVal.at(idelay),"TIB");
  }

  // make histo for each delay value
  for(size_t idelay = 0; idelay < delayVal.size(); idelay++){
    vector<TH1F*> histoToPlot;
    vector<TF1*> funcToPlot;
    // loop on the TOB mape
    for(auto ientry : TOBlayers){
      for(auto imap : ientry.second){
	if(imap.first != delayVal.at(idelay)) continue;
	histoToPlot.push_back(imap.second);
	funcToPlot.push_back(fitTOBlayers[ientry.first][imap.first]);
      }
    }
    plotHisto(canvas,histoToPlot,funcToPlot,outputDIR,delayVal.at(idelay),"TOB");
  }

  // make histo for each delay value
  for(size_t idelay = 0; idelay < delayVal.size(); idelay++){
    vector<TH1F*> histoToPlot;
    vector<TF1*> funcToPlot;
    // loop on the TID mape
    for(auto ientry : TIDlayers){
      for(auto imap : ientry.second){
	if(imap.first != delayVal.at(idelay)) continue;
	histoToPlot.push_back(imap.second);
	funcToPlot.push_back(fitTIDlayers[ientry.first][imap.first]);
      }
    }
    plotHisto(canvas,histoToPlot,funcToPlot,outputDIR,delayVal.at(idelay),"TID");
  }

  // make histo for each delay value
  for(size_t idelay = 0; idelay < delayVal.size(); idelay++){
    vector<TH1F*> histoToPlot;
    vector<TF1*> funcToPlot;
    // loop on the TOB mape
    for(auto ientry : TECPTlayers){
      for(auto imap : ientry.second){
	if(imap.first != delayVal.at(idelay)) continue;
	histoToPlot.push_back(imap.second);
	funcToPlot.push_back(fitTECPTlayers[ientry.first][imap.first]);
      }
    }
    plotHisto(canvas,histoToPlot,funcToPlot,outputDIR,delayVal.at(idelay),"TECPT");
  }

  // make histo for each delay value
  for(size_t idelay = 0; idelay < delayVal.size(); idelay++){
    vector<TH1F*> histoToPlot;
    vector<TF1*> funcToPlot;
    // loop on the TOB mape
    for(auto ientry : TECPtlayers){
      for(auto imap : ientry.second){
	if(imap.first != delayVal.at(idelay)) continue;
	histoToPlot.push_back(imap.second);
	funcToPlot.push_back(fitTECPtlayers[ientry.first][imap.first]);
      }
    }
    plotHisto(canvas,histoToPlot,funcToPlot,outputDIR,delayVal.at(idelay),"TECPt");
  }

  // make histo for each delay value
  for(size_t idelay = 0; idelay < delayVal.size(); idelay++){
    vector<TH1F*> histoToPlot;
    vector<TF1*> funcToPlot;
    // loop on the TOB mape
    for(auto ientry : TECMTlayers){
      for(auto imap : ientry.second){
	if(imap.first != delayVal.at(idelay)) continue;
	histoToPlot.push_back(imap.second);
	funcToPlot.push_back(fitTECMTlayers[ientry.first][imap.first]);
      }
    }
    plotHisto(canvas,histoToPlot,funcToPlot,outputDIR,delayVal.at(idelay),"TECMT");
  }

  // make histo for each delay value
  for(size_t idelay = 0; idelay < delayVal.size(); idelay++){
    vector<TH1F*> histoToPlot;
    vector<TF1*> funcToPlot;
    // loop on the TOB mape
    for(auto ientry : TECMtlayers){
      for(auto imap : ientry.second){
	if(imap.first != delayVal.at(idelay)) continue;
	histoToPlot.push_back(imap.second);
	funcToPlot.push_back(fitTECMtlayers[ientry.first][imap.first]);
      }
    }
    plotHisto(canvas,histoToPlot,funcToPlot,outputDIR,delayVal.at(idelay),"TECMt");
  } 
}


/// main function that run the analysis
void makeChargeDistributionPerLayer(string file0,  // inputfile
				    string observable   = "maxCharge",   // observable to be considered: maxCharge, S/N ..etc
				    string outputDIR    = "prompt" // output directory name
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

  system(("mkdir -p "+outputDIR).c_str());

  // open the file and prepare the cluster tree, adding the other trees as frined --> memory consuming
  std::cout<<"############################################"<<std::endl;
  std::cout<<"###### makeChargeDistributionPerLayer ######"<<std::endl;
  std::cout<<"############################################"<<std::endl;

  std::cout<<"Open Input Files"<<std::endl;
  TFile* _file0  = TFile::Open(file0.c_str());
  TTree* clusters    = (TTree*)_file0->FindObjectAny("clusters");
  TTree* readoutMap  = (TTree*)_file0->FindObjectAny("readoutMap");
  clusters->SetEventList(0);  

  // run per layer analysis
  LayerPlots(clusters,readoutMap,observable,outputDIR);
  TCanvas* canvas = new TCanvas("canvas","",600,650);
  canvas->cd();
  plotDistributions(canvas,outputDIR);
  
}

