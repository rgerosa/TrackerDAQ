#include <vector>
#include <iostream>
#include <fstream>
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
#include <TChain.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TLeaf.h>

#include "../delayUtils.h"

#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooExtendPdf.h"
#include "RooDataHist.h"
#include "RooLandau.h"
#include "RooGaussian.h"
#include "RooFFTConvPdf.h"


using namespace RooFit;
 
// parametrize profile with a gaussian shape
static bool isGaussian = true;
// reduce the number of events by
static int reductionFactor = 1;
// min number of filled bins
static int minFilledBin    = 4;
// max allowed delay
static float peakBoundary  = 9.36;
// min amplitude
static float amplitudeMin  = 40.; 
// min signficance delay/sigma_delay
static float significance  = 2;
// min delay to apply significance cut
static float minDelayForSignificance = 2;
// max sigma allowed
static float maxSigma = 20;
// decide whether use histograms values or fitted ones
static bool doFitOfClusterShape = true;
// decide whether use histograms values or fitted ones
static bool storeFitOfClusterShape = true;

static std::vector<TFile* > files;
 

/// function that runs on the evnet and produce profiles for layers
void ChannnelPlots(const std::vector<TTree* > & tree, 
		   const std::vector<TTree* > & map,
		   TTree* corrections,		   
		   std::map<uint32_t,std::map<float,std::shared_ptr<TH1F> > > channelMap,
		   std::map<uint32_t,std::shared_ptr<TH1F> > & channelMapMean,
		   std::map<uint32_t,std::shared_ptr<TH1F> > & channelMapMPV,
		   std::map<uint32_t,int > & fitStatusMean,
		   std::map<uint32_t,int > & fitStatusMPV,
		   const string & observable,
		   const string & outputDIR) {

  std::cout<<"###############################################################"<<std::endl;
  std::cout << "Preparing Layer plot for the for all the different channels "<< std::endl; 
  std::cout<<"##############################################################"<<std::endl;

  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;

  if(observable == "clCorrectedSignalOverNoise" or observable == "clSignalOverNoise")
    amplitudeMin = 3.;

  // create vectors for the different Profiles
  float yMin   = 0, yMax = 0;
  int   nBinsY = 0;
  vector<double> delayBins;
  setLimitsAndBinning("delay",delayBins);
  setLimitsAndBinning(observable,yMin,yMax,nBinsY);

  int   correction;
  corrections->SetBranchStatus("*",kFALSE);
  corrections->SetBranchStatus("detid",kTRUE);
  corrections->SetBranchStatus("correction",kTRUE);
  corrections->SetBranchAddress("correction",&correction);

  if(tree.size() != map.size())
    std::cout<<"### Different size between delay and readoutMaps "<<std::endl;

  for(int iTree = 0; iTree < tree.size(); iTree++){
    std::cout<<"### Start analysis of Tree "<<iTree<<" of "<<tree.size()<<std::endl;
    // set branches for the cluster, readoutmap and no corrections trees
    uint32_t detid;
    float    clCorrectedSignalOverNoise, clSignalOverNoise, obs;
    tree.at(iTree)->SetBranchStatus("*",kFALSE);
    tree.at(iTree)->SetBranchStatus("detid",kTRUE);
    tree.at(iTree)->SetBranchStatus("clCorrectedSignalOverNoise",kTRUE);
    tree.at(iTree)->SetBranchStatus("clSignalOverNoise",kTRUE);
    tree.at(iTree)->SetBranchStatus(observable.c_str(),kTRUE);
    tree.at(iTree)->SetBranchAddress("detid",&detid);
    tree.at(iTree)->SetBranchAddress("clCorrectedSignalOverNoise",&clCorrectedSignalOverNoise);
    tree.at(iTree)->SetBranchAddress("clSignalOverNoise",&clSignalOverNoise);
    tree.at(iTree)->SetBranchAddress(observable.c_str(),&obs);

    float delay;
    map.at(iTree)->SetBranchStatus("*",kFALSE);
    map.at(iTree)->SetBranchStatus("detid",kTRUE);
    map.at(iTree)->SetBranchStatus("delay",kTRUE);
    map.at(iTree)->SetBranchAddress("delay",&delay);

    for(long int iEvent = 0; iEvent < tree.at(iTree)->GetEntries()/reductionFactor; iEvent++){
      cout.flush();
      if(iEvent % 100000 == 0) cout<<"\r"<<"iEvent "<<100*double(iEvent)/(tree.at(iTree)->GetEntries()/reductionFactor)<<" % ";

      // take the event                                                                                                                                                  
      tree.at(iTree)->GetEntry(iEvent);
      // take the map delay from the detid                                                                                                                                  
      map.at(iTree)->GetEntryWithIndex(detid);
      // take the correction from the detid
      corrections->GetEntryWithIndex(detid);

      if(channelMap[detid][delay-correction].get() == 0 or channelMap[detid][delay-correction].get() == NULL){
	TH1::AddDirectory(kFALSE);
	channelMap[detid][delay-correction] =  std::shared_ptr<TH1F> (new TH1F(Form("detid_%d_delay_%.2f",detid,delay-correction),"",nBinsY,yMin,yMax));
	channelMap[detid][delay-correction]->Sumw2();
      }
      float value = 0;
      if(observable == "maxCharge")
	value = obs*(clCorrectedSignalOverNoise)/(clSignalOverNoise);
      else 	
	value = obs;
      channelMap[detid][delay-correction]->Fill(value);
    }
    std::cout<<std::endl;
    if(iTree != 0){
      map.at(iTree)->SetDirectory(0);          
      files.at(iTree)->Close(); // to prevent high RAM occupancy    
    }
  }
  std::cout<<"Loop on events terminated"<<std::endl;  

  std::cout<<"Analyze each single channel --> MPV and Mean "<<std::endl;
  long int iChannel = 0;
  TFile* outfile = NULL;
  if(storeFitOfClusterShape){
    outfile = new TFile((outputDIR+"/distributions_"+observable+".root").c_str(),"RECREATE");
    outfile->cd();
  }
  std::auto_ptr<TCanvas> canvas(new TCanvas("canvas","",600,625));
  canvas->SetTickx();
  canvas->SetTicky();
  
  // make the convuluated fit as suggested in https://root.cern.ch/root/html/tutorials/fit/langaus.C.html
  Double_t parameters[4];
  Double_t parametersHigh[4];
  Double_t parametersLow[4];
  const Double_t *fit_parameters;
  const Double_t *fit_parameters_error;
  Double_t chi2;
  Int_t    ndf;  

  long int iBadChannelFit = 0;

  for(auto iMap : channelMap){
    iChannel++;
    cout.flush();
    if(iChannel % 100 == 0) cout<<"\r"<<"iChannel "<<100*double(iChannel)/(channelMap.size())<<" % ";
    TH1::AddDirectory(kFALSE);
    channelMapMPV[iMap.first]  = std::shared_ptr<TH1F> (new TH1F(Form("detid_%d_mpv",iMap.first),"",delayBins.size()-1,&delayBins[0]));    
    TH1::AddDirectory(kFALSE);
    channelMapMean[iMap.first] = std::shared_ptr<TH1F> (new TH1F(Form("detid_%d_mean",iMap.first),"",delayBins.size()-1,&delayBins[0]));    

    for(int iBin = 0; iBin < channelMapMPV[iMap.first]->GetNbinsX()+1; iBin++){
      channelMapMPV[iMap.first]->SetBinContent(iBin+1,0.);
      channelMapMPV[iMap.first]->SetBinError(iBin+1,0.);
      channelMapMean[iMap.first]->SetBinContent(iBin+1,0.);
      channelMapMean[iMap.first]->SetBinError(iBin+1,0.);
    }

    for(std::map<float,shared_ptr<TH1F> >::const_iterator itHisto = iMap.second.begin(); itHisto != iMap.second.end(); itHisto++){
      
      // to fill the Mean --> use always the histogram properties given its binning
      channelMapMean[iMap.first]->SetBinContent(channelMapMPV[iMap.first]->FindBin(itHisto->first),itHisto->second->GetMean());
      // as error use the error on the mean as a proxy
      channelMapMean[iMap.first]->SetBinError(channelMapMPV[iMap.first]->FindBin(itHisto->first),itHisto->second->GetMeanError());

      if(not doFitOfClusterShape){ // take the MPV as the maximu of the histograms and store histogram in case
	channelMapMPV[iMap.first]->SetBinContent(channelMapMPV[iMap.first]->FindBin(itHisto->first),itHisto->second->GetBinCenter(itHisto->second->GetMaximumBin()));
	// as error use the error on the mean as a proxy
	channelMapMPV[iMap.first]->SetBinError(channelMapMPV[iMap.first]->FindBin(itHisto->first),itHisto->second->GetMeanError());
	
	// to fill the MPV: delay found by means of FindBin, take the bin center of the maximum bin --> is ok as soon as the histogram has enough bins
	if(storeFitOfClusterShape and iChannel % 20 == 0){
	  if(not outfile->GetDirectory(Form("Delay_%.2f",itHisto->first)))
	    outfile->mkdir(Form("Delay_%.2f",itHisto->first));
	  outfile->cd(Form("Delay_%.2f",itHisto->first));
	  itHisto->second->Write();
	  outfile->cd();	            	  
	}
      }  
      else{
	
	// Width of the landau distribution
	parameters[0] = 2.; parametersHigh[0] = 50.; parametersLow[0] = 0.001;
	// MPV of landau peak
	parameters[1] = itHisto->second->GetMean(); parametersHigh[1] = yMax; parametersLow[1] = yMin;
	// Total area
	parameters[2] = itHisto->second->Integral(); parametersHigh[2] = itHisto->second->Integral()*5; parametersLow[2] = itHisto->second->Integral()/5;
	// width of gaussian
	parameters[3] = 10; parametersHigh[3] = 40; parametersLow[3] = 1.;
	// create the function
	TF1 *    fitfunc = new TF1(Form("fit_%s",itHisto->second->GetName()),langaufun,yMin,yMax,4);
	fitfunc->SetParameters(parameters);
	fitfunc->SetParNames("Width","MPV","Area","GSigma");
	fitfunc->SetParLimits(0,parametersLow[0],parametersHigh[0]);
	fitfunc->SetParLimits(1,parametersLow[1],parametersHigh[1]);
	fitfunc->SetParLimits(2,parametersLow[2],parametersHigh[2]);
	fitfunc->SetParLimits(3,parametersLow[3],parametersHigh[3]);

	uint32_t subdetid    = int((iMap.first-0x10000000)/0x2000000);
	uint32_t barrellayer = int((iMap.first%33554432)/0x4000);
	uint32_t TIDlayer    = int((iMap.first%33554432)/0x800)%4;
	uint32_t TECPlayer   = int((iMap.first%33554432)/0x4000)-32;
	uint32_t TECMlayer   = int((iMap.first%33554432)/0x4000)-16;
	
	// define range for a better description of the peak
	if(observable == "maxCharge"){	  
	  if(subdetid == 3) // TIB
	    fitfunc->SetRange(itHisto->second->GetBinCenter(itHisto->second->GetMaximumBin())-1.5*itHisto->second->GetRMS(),itHisto->second->GetBinCenter(itHisto->second->GetMaximumBin())+2.5*itHisto->second->GetRMS());
	  if(subdetid == 4) // TID
	    fitfunc->SetRange(itHisto->second->GetBinCenter(itHisto->second->GetMaximumBin())-1.5*itHisto->second->GetRMS(),itHisto->second->GetBinCenter(itHisto->second->GetMaximumBin())+2.5*itHisto->second->GetRMS());
	  else if(subdetid == 5) //TOB
	    fitfunc->SetRange(itHisto->second->GetBinCenter(itHisto->second->GetMaximumBin())-2.0*itHisto->second->GetRMS(),itHisto->second->GetBinCenter(itHisto->second->GetMaximumBin())+3.0*itHisto->second->GetRMS());
	  else if(subdetid == 6 and (TECPlayer < 6 or TECMlayer < 6))
	    fitfunc->SetRange(itHisto->second->GetBinCenter(itHisto->second->GetMaximumBin())-1.5*itHisto->second->GetRMS(),itHisto->second->GetBinCenter(itHisto->second->GetMaximumBin())+2.5*itHisto->second->GetRMS());
	  else if(subdetid == 6 and (TECMlayer >= 6 or TECMlayer >= 6))
	    fitfunc->SetRange(itHisto->second->GetBinCenter(itHisto->second->GetMaximumBin())-2.0*itHisto->second->GetRMS(),itHisto->second->GetBinCenter(itHisto->second->GetMaximumBin())+3.0*itHisto->second->GetRMS());
	}
	else
	  fitfunc->SetRange(itHisto->second->GetBinCenter(itHisto->second->GetMaximumBin())-1.5*itHisto->second->GetRMS(),itHisto->second->GetBinCenter(itHisto->second->GetMaximumBin())+itHisto->second->GetRMS()*3.0);

	// make fit and get parameters
	TFitResultPtr fitResult = itHisto->second->Fit(fitfunc,"RSQN");	
	if(not fitResult.Get() or fitResult->Status() != 0 or fitResult->CovMatrixStatus() <= 1){ 
	  iBadChannelFit++;	    
	  channelMapMPV[iMap.first]->SetBinContent(channelMapMPV[iMap.first]->FindBin(itHisto->first),itHisto->second->GetBinCenter(itHisto->second->GetMaximumBin()));
	  // as error use the error on the mean as a proxy
	  channelMapMPV[iMap.first]->SetBinError(channelMapMPV[iMap.first]->FindBin(itHisto->first),itHisto->second->GetMeanError());	  
	}
	else{

	  fit_parameters = fitResult->GetParams();
	  fit_parameters_error = fitResult->GetErrors();
	  chi2   = fitResult->Chi2();
	  ndf    = fitResult->Ndf();  
	  
 	  channelMapMPV[iMap.first]->SetBinContent(channelMapMPV[iMap.first]->FindBin(itHisto->first),fitfunc->GetMaximumX(yMin,yMax));
	  // as error use the error on Landau MPV as a proxy
	  channelMapMPV[iMap.first]->SetBinError(channelMapMPV[iMap.first]->FindBin(itHisto->first),fit_parameters_error[1]);
	  
	  if(storeFitOfClusterShape and iChannel %20 == 0){
	    if(not outfile->GetDirectory(Form("Delay_%.2f",itHisto->first)))
	      outfile->mkdir(Form("Delay_%.2f",itHisto->first));
	    outfile->cd(Form("Delay_%.2f",itHisto->first));
	    
	    // plot to store in a root file                                                                                                                                 
	    canvas->cd();
	    itHisto->second->GetYaxis()->SetTitle("Events");
	    itHisto->second->GetXaxis()->SetTitle(("Cluster "+observable).c_str());
	    itHisto->second->GetYaxis()->SetTitleSize(0.045);
	    itHisto->second->GetXaxis()->SetTitleSize(0.045);
	    itHisto->second->GetYaxis()->SetLabelSize(0.038);
	    itHisto->second->GetXaxis()->SetLabelSize(0.038);
	    itHisto->second->SetMarkerColor(kBlack);
	    itHisto->second->SetMarkerStyle(20);
	    itHisto->second->Draw("EP");
	    fitfunc->SetLineWidth(2);
	    fitfunc->SetLineColor(kRed);
	    fitfunc->Draw("same");
	    itHisto->second->Draw("EPsame");
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
	    canvas->Modified();
	    CMS_lumi(canvas.get(),"");
	    canvas->Write(itHisto->second->GetName());
	  }
	}
	delete fitfunc;	
      }      
    }    
  }
  std::cout<<std::endl;
  if(doFitOfClusterShape)
    std::cout<<"Bad channel fit from Landau+Gaus fit "<<iBadChannelFit<<std::endl;  
  
  // correct profiles and fit them with a gaussian
  std::cout<<"Analyze each single channel --> Profile "<<std::endl;
  long int iFit = 0;
  long int badFits = 0;
  long int noFitResult = 0;
  for(auto iprof : channelMapMean){
    cout.flush();
    if(iFit % 100 == 0) cout<<"\r"<<"iFit Profile "<<100*double(iFit)/(channelMapMean.size())<<" % ";
    if(iprof.second->Integral() == 0) {
      noFitResult++;
      continue;
    }
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q",false);
    if(result.Get()){
      fitStatusMean[iprof.first] = result->Status()+result->CovMatrixStatus() ;
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
    	badFits++;
      }
      iFit++;
    }
  }
 
  std::cout<<std::endl;
  std::cout<<"NFits performed on Mean "<<iFit<<" bad ones : "<<100*double(badFits)/iFit<<" % "<<" noFit result "<<100*double(noFitResult)/iFit<<std::endl;

  std::cout<<"Start making MPV fits "<<std::endl;
  iFit = 0;
  badFits = 0;
  noFitResult = 0;
  for(auto iprof : channelMapMPV){
    cout.flush();
    if(iFit % 100 == 0) cout<<"\r"<<"iFit Profile "<<100*double(iFit)/(channelMapMPV.size())<<" % ";
    if(iprof.second->Integral() == 0) {
      noFitResult++;
      continue;
    }
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      fitStatusMPV[iprof.first] = result->Status()+result->CovMatrixStatus() ;
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
	badFits++;
      }
      iFit++;
    }
  }
  std::cout<<std::endl;
  std::cout<<"NFits performed on MPV "<<iFit<<" bad ones : "<<100*double(badFits)/iFit<<" % "<<" noFit result "<<100*double(noFitResult)/iFit<<std::endl;
  
  std::cout<<"###################################"<<std::endl;
  std::cout<<"##### End of Channel Analysis #####"<<std::endl;
  std::cout<<"###################################"<<std::endl; 

}


void saveOutputTree(const vector<TTree* > & readoutMap, map<uint32_t,shared_ptr<TH1F> > & mapHist, map<uint32_t,int> & fitStatus, const string & outputDIR, const string & postfix, const string & observable){

  TFile*  outputFile;
  TTree*  outputTree;

  std::cout<<"### Dumpt output tree for tkCommissioner "<<std::endl;
  outputFile = new TFile((outputDIR+"/tree_delayCorrection_"+postfix+".root").c_str(),"RECREATE");
  outputFile->cd();
  // branches definition for the output tree
  readoutMap.at(0)->SetBranchStatus("*",kTRUE);
  readoutMap.at(0)->SetBranchStatus("moduleName",kFALSE);
  readoutMap.at(0)->SetBranchStatus("moduleId",kFALSE);
  readoutMap.at(0)->SetBranchStatus("delay",kFALSE);  
  outputTree = (TTree*) readoutMap.at(0)->CopyTree("");
  outputTree->SetName("delayCorrection");
  
  // set branch
  float    measuredDelay, measuredDelayUnc;
  float    measuredMeanAmplitude, measuredMeanAmplitudeUnc;
  float    measuredSigma, measuredSigmaUnc;
  float    delayCorr;
  int      fitRejected;
  int      nFilledBinRejected;
  int      sigmaRejected;
  int      amplitudeRejected;
  int      significanceRejected;
  int      notFound;
  int      peakOutBoundaryRejected;
  vector<float> amplitude;        
  vector<float> amplitudeUnc;        
  uint32_t detid;
  
  outputFile->cd();
  outputTree->SetBranchStatus("*",kTRUE);
  outputTree->GetBranch("detid")->SetName("Detid");
  outputTree->GetLeaf("detid")->SetName("Detid");
  
  TBranch* bmeasuredDelay    = outputTree->Branch("measuredDelay",&measuredDelay,"measuredDelay/F");
  TBranch* bmeasuredDelayUnc = outputTree->Branch("measuredDelayUnc",&measuredDelayUnc,"measuredDelayUnc/F");
  TBranch* bmeasuredMeanAmplitude    = outputTree->Branch("measuredMeanAmplitude",&measuredMeanAmplitude,"measuredMeanAmplitude/F");
  TBranch* bmeasuredMeanAmplitudeUnc = outputTree->Branch("measuredMeanAmplitudeUnc",&measuredMeanAmplitudeUnc,"measuredMeanAmplitudeUnc/F");
  TBranch* bmeasuredSigma    = outputTree->Branch("measuredSigma",&measuredSigma,"measuredSigma/F");
  TBranch* bmeasuredSigmaUnc = outputTree->Branch("measuredSigmaUnc",&measuredSigmaUnc,"measuredSigmaUnc/F");
  TBranch* bamplitude        = outputTree->Branch("amplitude","std::vector<float>",&amplitude);
  TBranch* bamplitudeUnc     = outputTree->Branch("amplitudeUnc","std::vector<float>",&amplitudeUnc);
  TBranch* bnotFound             = outputTree->Branch("notFound",&notFound,"notFound/I");
  TBranch* bsignificanceRejected = outputTree->Branch("significanceRejected",&significanceRejected,"significanceRejected/I");
  TBranch* bamplitudeRejected    = outputTree->Branch("amplitudeRejected",&amplitudeRejected,"amplitudeRejected/I");
  TBranch* bsigmaRejected        = outputTree->Branch("sigmaRejected",&sigmaRejected,"sigmaRejected/I");
  TBranch* bnFilledBinRejected   = outputTree->Branch("nFilledBinRejected",&nFilledBinRejected,"nFilledBinRejected/I");
  TBranch* bpeakOutBoundaryRejected   = outputTree->Branch("peakOutBoundaryRejected",&peakOutBoundaryRejected,"peakOutBoundaryRejected/I");
  TBranch* bfitRejected          = outputTree->Branch("fitRejected",&fitRejected,"fitRejected/I");
  TBranch* bdelayCorr        = outputTree->Branch("delayCorr",&delayCorr,"delayCorr/F");
  
  readoutMap.at(0)->SetBranchAddress("detid",&detid);
  outputTree->SetBranchAddress("Detid",&detid);
  
  long int notFoundChannels = 0;
  
  cout<<"### Start the loop on the channel map "<<endl;      
  for(long int iChannel = 0; iChannel < readoutMap.at(0)->GetEntries(); iChannel++){
    cout.flush();
    if(iChannel % 100 == 0) cout<<"\r"<<"iChannel "<<100*double(iChannel)/(readoutMap.at(0)->GetEntries())<<" % ";
    readoutMap.at(0)->GetEntry(iChannel);
        
    measuredDelay    = 0.;
    measuredDelayUnc = 0.;
    measuredMeanAmplitude    = 0.;
    measuredMeanAmplitudeUnc = 0.;
    measuredSigma    = 0.;
    measuredSigmaUnc = 0.;
    delayCorr        = 0.;
    notFound             = 0;
    significanceRejected = 0;
    amplitudeRejected    = 0;
    sigmaRejected        = 0;
    nFilledBinRejected   = 0;
    fitRejected          = 0;
    peakOutBoundaryRejected = 0;
    amplitude.clear();
    amplitudeUnc.clear();

    if(mapHist[detid] == 0 or mapHist[detid] == NULL or mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName())) == 0 or mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName())) == NULL){
      notFound = 1;       
      notFoundChannels++;	
    }    
    else{
      // fill the amplitude branch
      for(int iBin = 0; iBin < mapHist[detid]->GetNbinsX(); iBin++){
	amplitude.push_back(mapHist[detid]->GetBinContent(iBin));
	amplitudeUnc.push_back(mapHist[detid]->GetBinError(iBin));
      }
      
      measuredDelay    = mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParameter(1);
      measuredDelayUnc = mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParError(1);
      measuredMeanAmplitude    = mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParameter(0);
      measuredMeanAmplitudeUnc = mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParError(0);
      measuredSigma    = mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParameter(2);
      measuredSigmaUnc = mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParError(2);	
      
      /// apply selections for good fit --> check number of filled bins
      int nFiledBins = getFilledBins(mapHist[detid]);
      if(nFiledBins < minFilledBin){	
	nFilledBinRejected = 1;
	//store it in the output canvas file
	if(not outputFile->GetDirectory("binEntries"))
	  outputFile->mkdir("binEntries");
	outputFile->cd("binEntries");
	if(gDirectory->GetListOfKeys()->FindObject(Form("detid_%d",detid)) == 0){
	  std::shared_ptr<TCanvas> c1b = prepareCanvas(Form("detid_%d",detid),observable);
	  plotAll(c1b,mapHist[detid]);
	  c1b->Write();
	}
	outputFile->cd();	
      }
      
      // fit status: 0 for the fit, 3 for the covariance matrix
      if(fitStatus[detid] != 3){	  
	fitRejected = 1;
	//store it in the output canvas file
	if(not outputFile->GetDirectory("fitStatus"))
	  outputFile->mkdir("fitStatus");
	outputFile->cd("fitStatus");
	if(gDirectory->GetListOfKeys()->FindObject(Form("detid_%d",detid)) == 0){
	  std::shared_ptr<TCanvas> c1b = prepareCanvas(Form("detid_%d",detid),observable);
	  plotAll(c1b,mapHist[detid]);
	  c1b->Write();
	}
	outputFile->cd();	
      }
      
      // mean value outside boundary
      if(fabs(mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParameter(1)) > peakBoundary){	    
	peakOutBoundaryRejected = 1;	    
	//store it in the output canvas file
	if(not outputFile->GetDirectory("peakOutRange"))
	  outputFile->mkdir("peakOutRange");
	outputFile->cd("peakOutRange");
	if(gDirectory->GetListOfKeys()->FindObject(Form("detid_%d",detid)) == 0){
	  std::shared_ptr<TCanvas> c1b = prepareCanvas(Form("detid_%d",detid),observable);
	  plotAll(c1b,mapHist[detid]);
	  c1b->Write();
	}
	outputFile->cd();	
      } 
      
      // amplitude cut
      if(fabs(mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParameter(0)) < amplitudeMin){	    
	amplitudeRejected = 1;	    
	//store it in the output canvas file
	if(not outputFile->GetDirectory("smallAmplitude"))
	  outputFile->mkdir("smallAmplitude");
	outputFile->cd("smallAmplitude");
	if(gDirectory->GetListOfKeys()->FindObject(Form("detid_%d",detid)) == 0){
	  std::shared_ptr<TCanvas> c1b = prepareCanvas(Form("detid_%d",detid),observable);
	  plotAll(c1b,mapHist[detid]);
	  c1b->Write();
	}
	outputFile->cd();	
      } 
      
      // amplitude cut
      if(fabs(mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParameter(1))/fabs(mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParError(1)) < significance and fabs(mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParameter(1)) > minDelayForSignificance){	    
	significanceRejected = 1;	    
	//store it in the output canvas file
	if(not outputFile->GetDirectory("nonSignificantLargeShift"))
	  outputFile->mkdir("nonSignificantLargeShift");
	outputFile->cd("nonSignificantLargeShift");
	if(gDirectory->GetListOfKeys()->FindObject(Form("detid_%d",detid)) == 0){
	  std::shared_ptr<TCanvas> c1b = prepareCanvas(Form("detid_%d",detid),observable);
	  plotAll(c1b,mapHist[detid]);
	  c1b->Write();
	}
	outputFile->cd();	
      } 
      
      // sigma cut
      if(mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParameter(2) > maxSigma){	    
	sigmaRejected = 1;
	//store it in the output canvas file                                                                                                                               
	if(not outputFile->GetDirectory("largeGaussSigma"))
	  outputFile->mkdir("largeGaussSigma");
	outputFile->cd("largeGaussSigma");
	if(gDirectory->GetListOfKeys()->FindObject(Form("detid_%d",detid)) == 0){
	  std::shared_ptr<TCanvas> c1b = prepareCanvas(Form("detid_%d",detid),observable);
	  plotAll(c1b,mapHist[detid]);
	  c1b->Write();
	}
	outputFile->cd();
      }
      
      if(not notFound and not significanceRejected and not amplitudeRejected and not sigmaRejected and not nFilledBinRejected and not fitRejected and not peakOutBoundaryRejected)	
	delayCorr = mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParameter(1);
    }
    

    bmeasuredDelay->Fill();
    bmeasuredDelayUnc->Fill();
    bmeasuredMeanAmplitude->Fill();
    bmeasuredMeanAmplitudeUnc->Fill();
    bmeasuredSigma->Fill();
    bmeasuredSigmaUnc->Fill();
    bamplitude->Fill();      
    bamplitudeUnc->Fill();      
    bnotFound->Fill();
    bsignificanceRejected->Fill();
    bamplitudeRejected->Fill();
    bsigmaRejected->Fill();
    bnFilledBinRejected->Fill();
    bfitRejected->Fill();      
    bpeakOutBoundaryRejected->Fill();
    bdelayCorr->Fill();    
  }    

  cout<<endl;
  cout<<"### Write output tree for tkCommissioner: not found channels i.e. no clusters "<<100*float(notFoundChannels)/outputTree->GetEntries()<<" % "<<endl;    
  outputFile->cd();
  outputTree->BuildIndex("Detid");
  outputTree->Write(outputTree->GetName(),TObject::kOverwrite);
}


// by dedault some fit results are store in a root file, a text dump: detid fitted delay is produced, as well as deid delay uncertainty
void delayValidationPerModule(			      
			      string inputDIR, 
			      string file1    = "../data/nocorrection.root", // no correction file
			      string postfix      = "merged", // sub string to be used to find merged files
			      string observable   = "maxCharge",
			      string outputDIR    = "/home/rgerosa/TrackerDAQ/DELAY_SCAN/SkimmedTrees/DelayPerModule_test_maxCharge/",
			      bool   saveMeanCanvas     =  false,
			      bool   saveMPVCanvas      =  false,
			      bool   saveCorrectionTree =  true){

  // prepare style and load macros
  setTDRStyle();
  gROOT->SetBatch(kTRUE);

  system(("mkdir -p "+outputDIR).c_str());

  // open the file and prepare the cluster tree, adding the other trees as frined --> memory consuming
  std::cout<<"######################################"<<std::endl;
  std::cout<<"###### delayValidationPerModule ######"<<std::endl;
  std::cout<<"######################################"<<std::endl;


  std::cout<<"### Make input file list"<<std::endl;
  system(("find "+inputDIR+" -name \"*"+postfix+"*.root\" > file.temp").c_str());
  std::ifstream infile;
  string line;
  vector<string> fileList;
  infile.open("file.temp",ifstream::in);
  if(infile.is_open()){
    while(!infile.eof()){
      getline(infile,line);
      if(line != "" and TString(line).Contains(".root"))
        fileList.push_back(line);
    }
  }  
  system("rm file.temp");

  std::sort(fileList.begin(),fileList.end());

  std::cout<<"### Build the input TTree Vector"<<std::endl;
  std::vector< TTree*> clusters;
  std::vector< TTree*> readoutMap;
  for(auto ifile : fileList){    
    files.push_back(TFile::Open(ifile.c_str(),"READ"));    
    cout<<"file "<<ifile<<endl;
    clusters.push_back((TTree*) files.back()->FindObjectAny("clusters"));
    readoutMap.push_back((TTree*) files.back()->FindObjectAny("readoutMap"));
  }
  
  TFile* _file1 = TFile::Open(file1.c_str());
  TTree* delayCorrections = (TTree*)_file1->FindObjectAny("delayCorrections");
  
  //map with key the detId number, one profile associated to it
  std::map<uint32_t,std::map<float,std::shared_ptr<TH1F> > > channelMap; // for each delay value (i.e. run) gram for each detid storing the observable distribution
  std::map<uint32_t,std::shared_ptr<TH1F> > channelMapMean;
  std::map<uint32_t,std::shared_ptr<TH1F> > channelMapMPV;
  std::map<uint32_t,int > fitStatusMean;
  std::map<uint32_t,int > fitStatusMPV;
  ChannnelPlots(clusters,readoutMap,delayCorrections,channelMap,channelMapMean,channelMapMPV,fitStatusMean,fitStatusMPV,observable,outputDIR);  
  
  // dumpt in a text file to be displayed on the tracker map and uncertainty  
  // produce outputs
  if(not saveCorrectionTree){
    cout<<"### Dump peak in a text file for Mean"<<endl;
    ofstream dumpMean((outputDIR+"/dumpBestDelay_Mean_"+observable+".txt").c_str());
    int nullPointers = 0;
    for(auto imap : channelMapMean){    
      if(imap.second->GetFunction(Form("Gaus_%s",imap.second->GetName())))
	dumpMean<<imap.first<<"   "<<imap.second->GetFunction(Form("Gaus_%s",imap.second->GetName()))->GetParameter(1)<<"\n";
      else
	nullPointers++;
    }
    dumpMean.close();
    std::cout<<"### Null pointers "<<100*float(nullPointers)/channelMap.size()<<" %s "<<endl;

    cout<<"### Dump peak error in a text file for Mean"<<endl;
    ofstream dumpMeanError((outputDIR+"/dumpBestDelayError_Mean_"+observable+".txt").c_str());
    for(auto imap : channelMapMean){
      if(imap.second->GetFunction(Form("Gaus_%s",imap.second->GetName())))
	dumpMeanError<<imap.first<<"   "<<imap.second->GetFunction(Form("Gaus_%s",imap.second->GetName()))->GetParError(1)<<"\n";
    }  
    dumpMeanError.close();

    cout<<"### Dump peak in a text file for MPV"<<endl;
    ofstream dumpMPV((outputDIR+"/dumpBestDelay_MPV_"+observable+".txt").c_str());
    nullPointers = 0;
    for(auto imap : channelMapMPV){    
      if(imap.second->GetFunction(Form("Gaus_%s",imap.second->GetName())))
	dumpMPV<<imap.first<<"   "<<imap.second->GetFunction(Form("Gaus_%s",imap.second->GetName()))->GetParameter(1)<<"\n";
      else
	nullPointers++;
    }
    dumpMPV.close();
    std::cout<<"### Null pointers "<<100*float(nullPointers)/channelMap.size()<<" %s "<<endl;

    cout<<"### Dump peak error in a text file for MPV"<<endl;
    ofstream dumpMPVError((outputDIR+"/dumpBestDelayError_MPV_"+observable+".txt").c_str());
    for(auto imap : channelMapMPV){
      if(imap.second->GetFunction(Form("Gaus_%s",imap.second->GetName())))
	dumpMPVError<<imap.first<<"   "<<imap.second->GetFunction(Form("Gaus_%s",imap.second->GetName()))->GetParError(1)<<"\n";
    }  
    dumpMPVError.close();
  }
  
  if(saveMeanCanvas)
    saveAll(channelMapMean,outputDIR,observable,"Mean"); // save all the plots into a single root files
  if(saveMPVCanvas)
    saveAll(channelMapMPV,outputDIR,observable,"MPV");
  
  //// outptut tree with analysis result
  if(saveCorrectionTree){
    saveOutputTree(readoutMap,channelMapMean,fitStatusMean,outputDIR,"Mean",observable);
    saveOutputTree(readoutMap,channelMapMPV,fitStatusMPV,outputDIR,"MPV",observable);
  }
  cout<<"### Clear and Close "<<endl;  
  channelMap.clear();
  channelMapMPV.clear();
  channelMapMean.clear();
  readoutMap.clear();
  files.clear();
  
}

