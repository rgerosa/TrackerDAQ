#include "../CMS_lumi.h"
#include "../TrackerStrip.h"

///Gaussian quantiles
static float quantile       = 0.3173;
static float quantile1sigma = 0.317310507863;
static float quantile2sigma = 0.045500263896;
static float quantile3sigma = 0.002699796063;
static float quantile4sigma = 0.000063342484;
static float quantile5sigma = 0.000000573303;

/// skip particular FEDs
static vector<uint16_t> skipFEDid = {75};
static int reductionFactor        = 1;

// minimum RMS to be considered reasonable
static float minimumRMS          = 2;
static float maximumRMS          = 25;
static float maximumOverflow     = 0.3;
static float maximumSignificance = 10;

// double peaked strips
static float distance_cut   = 1;
static float ashman_cut     = 2;
static float chi2_cut       = 0.05;
static float amplitude_cut  = 0.85;
static float bimodality_cut = 0.55;

///////////////////
class noiseStrip {

public:

  noiseStrip(){
    nstrip = 0;
    noiseVal = 0;
    isDivided = false;
  };
  ~noiseStrip(){};

  int nstrip;
  float noiseVal;
  bool isDivided;

};

// store an output canvas for the test statistics
void storeOutputCanvas(TCanvas* canvas, TH1F * distribution, const TString & name, const string & outputDIR){

  canvas->cd();
  distribution->SetMarkerColor(kBlack);
  distribution->SetMarkerStyle(20);
  distribution->SetMarkerSize(1);
  distribution->GetXaxis()->SetTitle("Test statistics");
  distribution->GetYaxis()->SetTitle("Events");
  distribution->Draw("hist");
  CMS_lumi(canvas,"",true);

  canvas->SaveAs((outputDIR+"/"+string(name)+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+string(name)+".pdf").c_str(),"pdf");
    
} 

// Storing an output canvas with the noise distribution and two alternative fitting functions -> for double peaked channels
void storeOutputCanvas(TCanvas* canvas, TH1F * noiseDistribution,  TF1 * fitFunction, TF1 * fitFunction2, const TString & name){

  canvas->cd();
  noiseDistribution->SetMarkerColor(kBlack);
  noiseDistribution->SetMarkerStyle(20);
  noiseDistribution->SetMarkerSize(1);
  noiseDistribution->GetXaxis()->SetTitle("noise (ADC)");
  noiseDistribution->GetYaxis()->SetTitle("Events");
  noiseDistribution->Draw("EP");
  CMS_lumi(canvas,"",true);
  
  fitFunction->SetLineColor(kRed);
  fitFunction->SetLineWidth(2);
  fitFunction->Draw("Lsame");

  fitFunction2->SetLineColor(kBlue);
  fitFunction2->SetLineWidth(2);
  fitFunction2->Draw("Lsame");
  noiseDistribution->Draw("EPsame");

  canvas->Write(name);

}

// summary plot for a give strip
void storeOutputCanvas(TCanvas* canvas, TH1F * noiseDistribution,  TF1 * fitFunction, const TString & name, map<string,string> & parameters){

  canvas->cd();
  noiseDistribution->SetMarkerColor(kBlack);
  noiseDistribution->SetMarkerStyle(20);
  noiseDistribution->SetMarkerSize(1);
  noiseDistribution->GetXaxis()->SetTitle("noise (ADC)");
  noiseDistribution->GetYaxis()->SetTitle("Events");
  noiseDistribution->Draw("EP");
  CMS_lumi(canvas,"",true);
  
  fitFunction->SetLineColor(kRed);
  fitFunction->SetLineWidth(2);
  fitFunction->Draw("Lsame");
  noiseDistribution->Draw("EPsame");

  TLegend leg (0.6,0.55,0.95,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.SetBorderSize(0);
  leg.AddEntry((TObject*)0,Form("P(#chi^{2}) = %s",parameters["fitChi2Probab"].c_str()),"");
  leg.AddEntry((TObject*)0,Form("P(KS) = %s",parameters["kSProbab"].c_str()),"");
  leg.AddEntry((TObject*)0,Form("P(JB) = %s",parameters["jBProbab"].c_str()),"");
  leg.AddEntry((TObject*)0,Form("P(AD) = %s",parameters["aDProbab"].c_str()),"");
  leg.AddEntry((TObject*)0,Form("Mean = %s",parameters["fitGausMean"].c_str()),"");
  leg.AddEntry((TObject*)0,Form("Sigma = %s",parameters["fitGausSigma"].c_str()),"");
  leg.AddEntry((TObject*)0,Form("Skewness = %s",parameters["noiseSkewness"].c_str()),"");
  leg.AddEntry((TObject*)0,Form("Kurtosis = %s",parameters["noiseKurtosis"].c_str()),"");
  leg.Draw("same");
  canvas->Write(name);

}


////////// main function

void plotPedestalAnalysis(string inputFileName, string outputDIR, bool testDoubleGaussianChannels = true){

  system(("mkdir -p "+outputDIR).c_str());

  setTDRStyle();
  gROOT->SetBatch(kTRUE);

  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  inputFile->cd();

  // read the input tree
  TTree* tree = (TTree*) inputFile->Get("pedestalFullNoise");

  uint32_t detid,fedKey;
  uint16_t fecCrate,fecSlot, fecRing, ccuAdd, ccuChan, lldChannel, fedId, fedCh, apvId, stripId;
  float    fitChi2Probab, kSProbab, jBProbab, aDProbab, fitChi2;
  float    noiseSkewness, noiseKurtosis;
  float    fitGausMean, fitGausSigma, fitGausNormalization;
  float    fitGausMeanError, fitGausSigmaError, fitGausNormalizationError;
  vector<float>* noiseDistribution = 0;
  vector<float>* noiseDistributionError = 0;
  float    nBin, xMin, xMax;
  float    pedestal, noise, noiseRMS;

  tree->SetBranchStatus("*",kFALSE);
  tree->SetBranchStatus("detid",kTRUE);
  tree->SetBranchStatus("ccuChan",kTRUE);
  tree->SetBranchStatus("lldChannel",kTRUE);
  tree->SetBranchStatus("apvId",kTRUE);
  tree->SetBranchStatus("stripId",kTRUE);
  tree->SetBranchStatus("pedestal",kTRUE);
  tree->SetBranchStatus("noise",kTRUE);
  tree->SetBranchStatus("noiseRMS",kTRUE);

  tree->SetBranchAddress("detid",&detid);
  tree->SetBranchAddress("lldChannel",&lldChannel);
  tree->SetBranchAddress("apvId",&apvId);
  tree->SetBranchAddress("stripId",&stripId);
  tree->SetBranchAddress("pedestal",&pedestal);
  tree->SetBranchAddress("noise",&noise);
  tree->SetBranchAddress("noiseRMS",&noiseRMS);

  /// Basic info
  bool isfound = false;
  // loop on the pedestak analysis tree
  map<string,noiseStrip*> noiseSpreadAPV;
  map<string,noiseStrip*> noiseMeanAPV;

  // First loop to evaluate mean value of the noise per APV level
  cout<<"Loop to evaluate Mean noise across each APV"<<endl;
  for(long int iChannel = 0; iChannel < tree->GetEntries(); iChannel++){
    tree->GetEntry(iChannel);
    cout.flush();
    if(iChannel %10000 == 0) cout<<"\r"<<"iChannel "<<100*double(iChannel)/(tree->GetEntries()/reductionFactor)<<" % ";
    if(iChannel > double(tree->GetEntries())/reductionFactor) break;

    // skip problematic fed id
    isfound = false;
    for(auto skipfed : skipFEDid){
      if(fedId == skipfed) isfound = true;
    }
    if(isfound) continue;

    // to evaluate mean value and spread of noise
    string nameNoise = "Detid_"+to_string(detid)+"_lldCh"+to_string(lldChannel)+"_apv_"+to_string(apvId);
    
    if(noiseMeanAPV[nameNoise] == 0 or noiseMeanAPV[nameNoise] == NULL)
      noiseMeanAPV[nameNoise] = new noiseStrip();
    noiseMeanAPV[nameNoise]->nstrip += 1;
    noiseMeanAPV[nameNoise]->noiseVal += noise;
    noiseMeanAPV[nameNoise]->isDivided = false;    
  }

  std::cout<<std::endl;
  // Normalize it
  for(auto iapv : noiseMeanAPV){
    iapv.second->noiseVal /= float(iapv.second->nstrip);
    iapv.second->isDivided = true;
  }

  cout<<"Loop to evaluate Noise spread across each APV"<<endl;
  for(long int iChannel = 0; iChannel < tree->GetEntries(); iChannel++){
    tree->GetEntry(iChannel);
    cout.flush();
    if(iChannel %10000 == 0) cout<<"\r"<<"iChannel "<<100*double(iChannel)/(tree->GetEntries()/reductionFactor)<<" % ";
    if(iChannel > double(tree->GetEntries())/reductionFactor) break;

    // skip problematic fed id
    isfound = false;
    for(auto skipfed : skipFEDid){
      if(fedId == skipfed) isfound = true;
    }
    if(isfound) continue;

    // to evaluate mean value and spread of noise
    string nameNoise = "Detid_"+to_string(detid)+"_lldCh"+to_string(lldChannel)+"_apv_"+to_string(apvId);
    
    if(noiseSpreadAPV[nameNoise] == 0 or noiseSpreadAPV[nameNoise] == NULL)
      noiseSpreadAPV[nameNoise] = new noiseStrip();
    noiseSpreadAPV[nameNoise]->nstrip += 1;
    noiseSpreadAPV[nameNoise]->noiseVal += (noiseRMS-noiseMeanAPV[nameNoise]->noiseVal)*(noiseRMS-noiseMeanAPV[nameNoise]->noiseVal);
    noiseSpreadAPV[nameNoise]->isDivided = false;    
  }
  std::cout<<std::endl;

  for(auto iapv : noiseSpreadAPV){
    float val = iapv.second->noiseVal;
    iapv.second->noiseVal = sqrt(val/(iapv.second->nstrip-1));
    iapv.second->isDivided = true;
  }

  ////////////////// --------------------------------------------    
  tree->SetBranchStatus("fedId",kTRUE);
  tree->SetBranchStatus("fedKey",kTRUE);
  tree->SetBranchStatus("fecCrate",kTRUE);
  tree->SetBranchStatus("fecSlot",kTRUE);
  tree->SetBranchStatus("fecRing",kTRUE);
  tree->SetBranchStatus("ccuAdd",kTRUE);
  tree->SetBranchStatus("fedCh",kTRUE);
  tree->SetBranchStatus("fitChi2",kTRUE);
  tree->SetBranchStatus("fitChi2Probab",kTRUE);
  tree->SetBranchStatus("kSProbab",kTRUE);
  tree->SetBranchStatus("jBProbab",kTRUE);
  tree->SetBranchStatus("aDProbab",kTRUE);
  tree->SetBranchStatus("fitGausNormalization",kTRUE);
  tree->SetBranchStatus("fitGausMean",kTRUE);
  tree->SetBranchStatus("fitGausSigma",kTRUE);
  tree->SetBranchStatus("fitGausNormalizationError",kTRUE);
  tree->SetBranchStatus("fitGausMeanError",kTRUE);
  tree->SetBranchStatus("fitGausSigmaError",kTRUE);
  tree->SetBranchStatus("noiseSkewness",kTRUE);
  tree->SetBranchStatus("noiseKurtosis",kTRUE);
  tree->SetBranchStatus("noiseDistributionError",kTRUE);
  tree->SetBranchStatus("noiseDistribution",kTRUE);
  tree->SetBranchStatus("nBin",kTRUE);
  tree->SetBranchStatus("xMin",kTRUE);
  tree->SetBranchStatus("xMax",kTRUE);

  tree->SetBranchAddress("fedKey",&fedKey);
  tree->SetBranchAddress("fecCrate",&fecCrate);
  tree->SetBranchAddress("fecSlot",&fecSlot);
  tree->SetBranchAddress("fecRing",&fecRing);
  tree->SetBranchAddress("ccuAdd",&ccuAdd);
  tree->SetBranchAddress("ccuChan",&ccuChan);
  tree->SetBranchAddress("fedId",&fedId);
  tree->SetBranchAddress("fedCh",&fedCh);
  tree->SetBranchAddress("fitGausNormalization",&fitGausNormalization);
  tree->SetBranchAddress("fitGausMean",&fitGausMean);
  tree->SetBranchAddress("fitGausSigma",&fitGausSigma);
  tree->SetBranchAddress("fitGausNormalizationError",&fitGausNormalizationError);
  tree->SetBranchAddress("fitGausMeanError",&fitGausMeanError);
  tree->SetBranchAddress("fitGausSigmaError",&fitGausSigmaError);
  tree->SetBranchAddress("fitChi2",&fitChi2);
  tree->SetBranchAddress("fitChi2Probab",&fitChi2Probab);
  tree->SetBranchAddress("noiseSkewness",&noiseSkewness);
  tree->SetBranchAddress("noiseKurtosis",&noiseKurtosis);
  tree->SetBranchAddress("kSProbab",&kSProbab);
  tree->SetBranchAddress("aDProbab",&aDProbab);
  tree->SetBranchAddress("jBProbab",&jBProbab);
  tree->SetBranchAddress("noiseDistributionError",&noiseDistributionError);
  tree->SetBranchAddress("noiseDistribution",&noiseDistribution);
  tree->SetBranchAddress("nBin",&nBin);
  tree->SetBranchAddress("xMin",&xMin);
  tree->SetBranchAddress("xMax",&xMax);

  ///// Bad Strip list
  vector<TrackerStrip> badStrip;
  vector<TrackerStrip> badStripConservative;
  vector<TrackerStrip> badStripVsOffline;
  vector<TrackerStrip> badStripVsOfflineWithSignificance;
  ////
  TCanvas* canvas = new TCanvas("canvas","canvas",600,650);

  // Null integral
  TFile* badStripsNullIntegral = new TFile((outputDIR+"/badStripsNullIntegral.root").c_str(),"RECREATE");
  // small RMS
  TFile* badStripsSmallRMS = new TFile((outputDIR+"/badStripsSmallRMS.root").c_str(),"RECREATE");
  // overlfow
  TFile* badStripsOverflow = new TFile((outputDIR+"/badStripsOverflow.root").c_str(),"RECREATE");
  // overlfow
  TFile* badStripsNoiseSignificance = new TFile((outputDIR+"/badStripsNoiseSignificance.root").c_str(),"RECREATE");
  // rejected by Anderson Darling test
  TFile* badADTest   = new TFile((outputDIR+"/badStripsADTest.root").c_str(),"RECREATE");
  // rejected by KS test
  TFile* badKSTest   = new TFile((outputDIR+"/badStripsKSTest.root").c_str(),"RECREATE");
  // rejected by JB test
  TFile* badJBTest   = new TFile((outputDIR+"/badStripsJBTest.root").c_str(),"RECREATE");
  // rejected by Chi2 probability
  TFile* badChi2Test = new TFile((outputDIR+"/badStripsChi2Test.root").c_str(),"RECREATE");
  // selected bad strips using different methods
  TFile* badCombinedTest = new TFile((outputDIR+"/badStripsCombined.root").c_str(),"RECREATE");
  // selected bad strips using different methods
  TFile* badCombinedConservativeTest = new TFile((outputDIR+"/badStripsCombinedConservative.root").c_str(),"RECREATE");
  // Bad JB but good for AD and KS
  TFile* badJBNotADNotKSTest  = new TFile((outputDIR+"/badStripsJBNotADNotKS.root").c_str(),"RECREATE");
  // bad KS but good AD
  TFile* badKSNotADTest       = new TFile((outputDIR+"/badStripsKSNotAD.root").c_str(),"RECREATE");
  // bad Chi2 but good JB, good AD and KS
  TFile* badChi2NotKSandJBandADTest = new TFile((outputDIR+"/badStripsChi2NotKSandJBandAD.root").c_str(),"RECREATE");
  // strips passing conservative but not tight
  TFile* badCombinedPassingConservativeNotTight = new TFile((outputDIR+"/badCombinedPassingConservativeNotTight.root").c_str(),"RECREATE");
  // strips passing tight but not conservative
  TFile* badCombinedPassingTightNotConservative = new TFile((outputDIR+"/badCombinedPassingTightNotConservative.root").c_str(),"RECREATE");
  
  //// counters 
  long int nbadNullIntegral = 0;
  long int nbadSmallRMS = 0;
  long int nbadOverflow = 0;
  long int nbadNoiseSignificance = 0;
  long int nbadKSTest = 0;
  long int nbadJBTest = 0;
  long int nbadADTest = 0;
  long int nbadChi2Test = 0;
  long int nbadCombinedTest = 0;
  long int nbadCombinedConservativeTest = 0;
  long int nbadJBNotADNotKSTest = 0;
  long int nbadKSNotADTest  = 0;
  long int nbadChi2NotKSandJBandADTest = 0;
  long int nbadCombinedPassingTightNotConservative = 0;
  long int nbadCombinedPassingConservativeNotTight = 0;
  
  // map of bad channels according to different methods
  map<uint32_t,uint32_t> moduleDenominator;
  map<uint32_t,uint32_t> moduleNumerator;
  map<uint32_t,uint32_t> moduleNumeratorConservative;
  map<uint32_t,uint32_t> moduleNumeratorConservativeNotTight;
  map<uint32_t,uint32_t> moduleNumeratorTightNotConservative;
  map<uint32_t,uint32_t> moduleNumeratorNullIntegral;
  map<uint32_t,uint32_t> moduleNumeratorSmallRMS;
  map<uint32_t,uint32_t> moduleNumeratorOverflow;
  map<uint32_t,uint32_t> moduleNumeratorNoiseSignificance;
  map<uint32_t,uint32_t> moduleNumeratorKS;
  map<uint32_t,uint32_t> moduleNumeratorAD;
  map<uint32_t,uint32_t> moduleNumeratorJB;
  map<uint32_t,uint32_t> moduleNumeratorChi2;
  map<uint32_t,uint32_t> moduleNumeratorDoublePeak;
  
  ////// Double peaked channels
  long int nbadDoublePeakDistance   = 0;  // distance between two peak
  long int nbadDoublePeakAshman     = 0;  // Ashman distance
  long int nbadDoublePeakChi2       = 0;  // Chi2
  long int nbadDoublePeakAmplitude  = 0;  // Peak amplitude
  long int nbadDoublePeakBimodality = 0;  // Bimodality
  long int nbadDoublePeakCombined   = 0;  // Combined criteria


  TH1F* chi2Distance    = new TH1F("chi2Distance","",100,0,1);
  TH1F* peakDistance    = new TH1F("peakDistance","",100,0,3);
  TH1F* ashmanDistance  = new TH1F("ashmanDistance","",100,0,5);
  TH1F* bimodalityDistance = new TH1F("bimodalityDistance","",100,0,1);
  TH1F* amplitudeRatioDistance = new TH1F("amplitudeRatioDistance","",100,0,3);

  chi2Distance->Sumw2();
  peakDistance->Sumw2();
  ashmanDistance->Sumw2();
  bimodalityDistance->Sumw2();
  amplitudeRatioDistance->Sumw2();

  TFile* multiPeakChannelsChi2       = new TFile((outputDIR+"/multiPeakChannelsChi2.root").c_str(),"RECREATE");
  TFile* multiPeakChannelsDistance   = new TFile((outputDIR+"/multiPeakChannelsDistance.root").c_str(),"RECREATE");
  TFile* multiPeakChannelsAshman     = new TFile((outputDIR+"/multiPeakChannelsAshman.root").c_str(),"RECREATE");
  TFile* multiPeakChannelsAmplitude  = new TFile((outputDIR+"/multiPeakChannelsAmplitude.root").c_str(),"RECREATE");
  TFile* multiPeakChannelsBimodality = new TFile((outputDIR+"/multiPeakChannelsBimodality.root").c_str(),"RECREATE");
  TFile* multiPeakChannelsCombined   = new TFile((outputDIR+"/multiPeakChannelsCombined.root").c_str(),"RECREATE");

  int nonNullBins = 0;
  float chi2Ratio = 0;
  float distance  = 0;
  float ashman    = 0;
  float bimodality = 0;
  float amplitudeRatio = 0;

  
  // loop on the pedestal analysis tree
  string fedKeyStr ;
  TString name ;
  TH1F* noiseHist = NULL;
  TF1*  noiseFit       = NULL;
  TF1*  noiseFit2Gaus  = NULL;
  TFitResultPtr result;
  std::map<string,string> fitParam;

 
  cout<<"Loop to make test statistics"<<endl;
  for(long int iChannel = 0; iChannel < tree->GetEntries(); iChannel++){
    tree->GetEntry(iChannel);
    cout.flush();
    if(iChannel %10000 == 0) cout<<"\r"<<"iChannel "<<100*double(iChannel)/(tree->GetEntries()/reductionFactor)<<" % ";
    if(iChannel > double(tree->GetEntries())/reductionFactor) break;

    // skip problematic fed id
    isfound = false;
    for(auto skipfed : skipFEDid){
      if(fedId == skipfed) isfound = true;
    }
    if(isfound) continue;

    // make selections to identify bad noisy channels (not gaussian ones)
    std::stringstream stream;
    stream << std::hex << fedKey;
    fedKeyStr = stream.str();
    if(fedKeyStr.size() == 4)
      name = Form("fecCrate%d_fecSlot%d_fecRing%d_ccuAdd%d_ccuCh%d_fedKey0x0000%s_lldCh%d_apv%d_strip%d",fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,fedKeyStr.c_str(),lldChannel,apvId,stripId);
    else if(fedKeyStr.size() == 5)
      name = Form("fecCrate%d_fecSlot%d_fecRing%d_ccuAdd%d_ccuCh%d_fedKey0x000%s_lldCh%d_apv%d_strip%d",fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,fedKeyStr.c_str(),lldChannel,apvId,stripId);

    fitParam.clear();
    stringstream sMean;
    sMean << std::scientific << fitGausMean;
    fitParam["fitGausMean"]   = sMean.str();
    stringstream sSigma;
    sSigma << std::scientific << fitGausSigma;
    fitParam["fitGausSigma"]  = sSigma.str();
    stringstream sSkew;
    sSkew << std::scientific << noiseSkewness;
    fitParam["noiseSkewness"] = sSkew.str();
    stringstream sKurt;
    sKurt << std::scientific << noiseKurtosis;
    fitParam["noiseKurtosis"] = sKurt.str();
    stringstream sAD;
    sAD << std::scientific << aDProbab;
    fitParam["aDProbab"] = sAD.str();
    stringstream sKS;
    sKS << std::scientific << kSProbab;
    fitParam["kSProbab"] = sKS.str();
    stringstream sJB;
    sJB << std::scientific << jBProbab;
    fitParam["jBProbab"] = sJB.str();
    stringstream sChi2;
    sChi2 << std::scientific << fitChi2Probab;
    fitParam["fitChi2Probab"] = sChi2.str();

    moduleDenominator[detid] = moduleDenominator[detid]+1;

    if(noiseHist == NULL){
      noiseHist = new TH1F ("noiseHist","",nBin,xMin,xMax);
      noiseHist->Sumw2();
    }
    noiseHist->Reset();

    // create the noise distribution for the given strip
    for(int iBin = 0; iBin < noiseDistribution->size(); iBin++){
      noiseHist->SetBinContent(iBin+1,noiseDistribution->at(iBin));
      noiseHist->SetBinError(iBin+1,noiseDistributionError->at(iBin));
    }

    if(noiseFit == NULL)
      noiseFit = new TF1 ("noiseFist","gaus(0)",xMin,xMax);
    
    // set the parameters from the old fits
    noiseFit->SetRange(xMin,xMax);
    noiseFit->SetParameters(fitGausNormalization,fitGausMean,fitGausSigma);
    noiseFit->SetParError(0,fitGausNormalizationError);
    noiseFit->SetParError(1,fitGausMeanError);
    noiseFit->SetParError(2,fitGausSigmaError);

     // detect for each channel the number of non-null bins
    nonNullBins = 0;
    for(int iBin = 0; iBin < noiseHist->GetNbinsX(); iBin++){
      if(noiseHist->GetBinContent(iBin+1) != 0) nonNullBins++;
    }

    // set noise mean and RMS (spread)    
    string nameNoise = "Detid_"+to_string(detid)+"_lldCh"+to_string(lldChannel)+"_apv_"+to_string(apvId);
    if(noiseSpreadAPV[nameNoise]->isDivided == false and noiseMeanAPV[nameNoise]->isDivided == false){
      noiseSpreadAPV[nameNoise]->noiseVal = sqrt(noiseSpreadAPV[nameNoise]->noiseVal*noiseSpreadAPV[nameNoise]->noiseVal-noiseMeanAPV[nameNoise]->noiseVal*noiseMeanAPV[nameNoise]->noiseVal)/sqrt(float(noiseMeanAPV[nameNoise]->nstrip)-1);
      noiseMeanAPV[nameNoise]->noiseVal /= float(noiseMeanAPV[nameNoise]->nstrip);
      noiseSpreadAPV[nameNoise]->isDivided = true;
      noiseMeanAPV[nameNoise]->isDivided = true;
    }
    
    
    //Null integral
    if(noiseHist->Integral() <= 0){
       badStripsNullIntegral->cd();
       nbadNullIntegral++;
       moduleNumeratorNullIntegral[detid] +=1;
       continue;
    }
    
    /// very small RMS
    if(noiseHist->GetRMS() < minimumRMS){
      badStripsSmallRMS->cd();
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);
      nbadSmallRMS++;
      moduleNumeratorSmallRMS[detid] +=1;
      badStrip.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));
      badStripConservative.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));
      badStripVsOffline.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));  
      badStripVsOfflineWithSignificance.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));  
      continue;
    }

    /// overflow
    if((noiseHist->GetBinContent(1)+noiseHist->GetBinContent(noiseHist->GetNbinsX()))/noiseHist->Integral() >= maximumOverflow){
      badStripsOverflow->cd();
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);
      nbadOverflow++;
      moduleNumeratorOverflow[detid] +=1;
      badStrip.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));
      badStripConservative.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));
      badStripVsOffline.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));  
      badStripVsOfflineWithSignificance.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));  
      continue;
    }

    if(fabs(noiseHist->GetRMS()-noiseMeanAPV[nameNoise]->noiseVal)/noiseSpreadAPV[nameNoise]->noiseVal > maximumSignificance){
      badStripsNoiseSignificance->cd();
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);
      nbadNoiseSignificance++;
      moduleNumeratorNoiseSignificance[detid] +=1;
      badStrip.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));
      badStripConservative.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));
      badStripVsOfflineWithSignificance.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));  
      continue;
    }

    // probability for KS test smaller than a given CL --> three sigma means 1% of probability
    if(kSProbab < quantile3sigma){
      badKSTest->cd();
      nbadKSTest++;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);      
      moduleNumeratorKS[detid] += 1;
    }
    
    // probability for the JB test to be smaller than a given CL
    if(jBProbab < quantile5sigma){
      badJBTest->cd();
      nbadJBTest++;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);      
      moduleNumeratorJB[detid] += 1;
    }

    // Chi2 probability
    if(fitChi2Probab < quantile5sigma){
      badChi2Test->cd();
      nbadChi2Test++;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);      
      moduleNumeratorChi2[detid] += 1;
    }

    // Anderson Darling test
    if(aDProbab < quantile3sigma){
      badADTest->cd();
      nbadADTest++;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);
    }

    // relatively low anderson darling probability but KS less than 1%
    if(kSProbab < quantile3sigma and aDProbab > quantile3sigma and aDProbab < quantile){
      badKSNotADTest->cd();
      nbadKSNotADTest++;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);      
      moduleNumeratorKS[detid] += 1;
    }

    //// bad for JB but not for AD (AD sufficiently low [-CL,1sigma]
    if(jBProbab < quantile5sigma and aDProbab > quantile3sigma and aDProbab < quantile and kSProbab > quantile3sigma){
      badJBNotADNotKSTest->cd();
      nbadJBNotADNotKSTest++;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);      
      moduleNumeratorAD[detid] += 1;
    }


    if(fitChi2Probab < quantile5sigma and jBProbab > quantile5sigma and aDProbab > quantile3sigma and aDProbab < quantile and kSProbab > quantile3sigma){
      badChi2NotKSandJBandADTest->cd();
      nbadChi2NotKSandJBandADTest++;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);      
    }

    // combining all of them --> conservative approach
    bool passConservative = false;
    if(aDProbab < quantile3sigma and kSProbab < quantile3sigma and jBProbab < quantile3sigma and fitChi2Probab < quantile3sigma){
      passConservative = true;
      badCombinedConservativeTest->cd();
      nbadCombinedConservativeTest++;
      moduleNumeratorConservative[detid] = moduleNumeratorConservative[detid]+1;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);      
      // store it in order to dispaly on the tracker map
      badStripConservative.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));
    }

    //// combining all of them --> aggressive approach
    bool passTight = false;
    if(aDProbab < quantile3sigma or 
       (aDProbab > quantile3sigma and aDProbab < quantile and kSProbab < quantile3sigma) or 
       (aDProbab > quantile3sigma and aDProbab < quantile and kSProbab > quantile3sigma and jBProbab < quantile5sigma)){

      passTight = true;
      badCombinedTest->cd();
      nbadCombinedTest++;
      moduleNumerator[detid] = moduleNumerator[detid]+1;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);      

      name = Form("fecCrate%d_fecSlot%d_fecRing%d_ccuAdd%d_ccuCh%d_fedKey0x0000%s_lldCh%d_apv%d_strip%d",fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,fedKeyStr.c_str(),lldChannel,apvId,stripId);      
      // store it in order to dispaly on the tracker map
      badStrip.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));
      
      //try to identify double peaked channels between the ones marked as bad by the analysis 
      if(testDoubleGaussianChannels){
	if(noiseFit2Gaus == NULL)
	  // double gaussian in which the sigma is constrained to be the same --> identifing clear two peak channels
	  noiseFit2Gaus = new TF1("dgaus","[0]*exp(-((x-[1])*(x-[1]))/(2*[2]*[2]))+[3]*exp(-((x-[4])*(x-[4]))/(2*[5]*[5]))",xMin,xMax);
	
	noiseFit2Gaus->SetRange(xMin,xMax);
	noiseFit2Gaus->SetParameter(0,fitGausNormalization/2);
	noiseFit2Gaus->SetParameter(3,fitGausNormalization/2);
	noiseFit2Gaus->SetParameter(1,1.);
	noiseFit2Gaus->SetParameter(4,-1.);
	noiseFit2Gaus->SetParameter(2,fitGausSigma);
	noiseFit2Gaus->SetParameter(5,fitGausSigma);
	noiseFit2Gaus->SetParLimits(1,0.,xMax);
	noiseFit2Gaus->SetParLimits(4,xMin,0);
	result = noiseHist->Fit(noiseFit2Gaus,"QSR");

	chi2Ratio = 0;
	distance  = 0;
	ashman    = 0;
	bimodality = 0;
	amplitudeRatio = 0;
	////////////
	if(result.Get() and noiseHist->Integral() != 0){

	  //compute the chi2 ratio between 1 gauss and 2 gauss fits
	  chi2Ratio = 0.5*ROOT::Math::chisquared_cdf_c((fitChi2/(result->Ndf()+3))/(result->Chi2()/result->Ndf()),1);
	  chi2Distance->Fill(chi2Ratio);	
		
	  //distance between peaks
	  distance = fabs(noiseFit2Gaus->GetParameter(1)-noiseFit2Gaus->GetParameter(4))/(2*sqrt(noiseFit2Gaus->GetParameter(2)*noiseFit2Gaus->GetParameter(5)));
	  peakDistance->Fill(distance);

	  // ashman distance
	  ashman   = TMath::Power(2,0.5)*abs(noiseFit2Gaus->GetParameter(1)-noiseFit2Gaus->GetParameter(4))/(sqrt(pow(noiseFit2Gaus->GetParameter(2),2)+pow(noiseFit2Gaus->GetParameter(5),2)));
	  ashmanDistance->Fill(ashman);	 

	  // bimodality coefficient
	  if(nonNullBins > 3)
	    bimodality = (noiseHist->GetSkewness()*noiseHist->GetSkewness()+1)/(noiseHist->GetKurtosis()+3*(nonNullBins-1)*(nonNullBins-1)/((nonNullBins-2)*(nonNullBins-3)));
	  else
	    bimodality = (noiseHist->GetSkewness()*noiseHist->GetSkewness()+1)/(noiseHist->GetKurtosis());
	  bimodalityDistance->Fill(bimodality);	  

	  // ratio of amplitudes
	  amplitudeRatio = std::min(noiseFit2Gaus->GetParameter(0),noiseFit2Gaus->GetParameter(3))/std::max(noiseFit2Gaus->GetParameter(0),noiseFit2Gaus->GetParameter(3));
	  amplitudeRatioDistance->Fill(amplitudeRatio);	
	  
	  ///// --> flagged by a simple distance 
	  if(distance > distance_cut){
	    multiPeakChannelsDistance->cd();
	    storeOutputCanvas(canvas,noiseHist,noiseFit,noiseFit2Gaus,name);
	    nbadDoublePeakDistance++;
	  }

	  /// ashman coefficient
	  if(ashman > ashman_cut){
	    multiPeakChannelsAshman->cd();
	    storeOutputCanvas(canvas,noiseHist,noiseFit,noiseFit2Gaus,name);
	    nbadDoublePeakAshman++;
	  }
	  /// chi2 ratio
	  if(chi2Ratio < chi2_cut){
	    multiPeakChannelsChi2->cd();
	    storeOutputCanvas(canvas,noiseHist,noiseFit,noiseFit2Gaus,name);
	    nbadDoublePeakChi2++;
	  }
	  /// amplitude ratios
	  if(amplitudeRatio > amplitude_cut){
	    multiPeakChannelsAmplitude->cd();
	    storeOutputCanvas(canvas,noiseHist,noiseFit,noiseFit2Gaus,name);
	    nbadDoublePeakAmplitude++;
	  }
	  //// bimodality
	  if(bimodality > bimodality_cut){
	    multiPeakChannelsBimodality->cd();
	    storeOutputCanvas(canvas,noiseHist,noiseFit,noiseFit2Gaus,name);
	    nbadDoublePeakBimodality++;
	  }
	  /// combo of ashman and amplitude ratio
	  if(ashman > ashman_cut && amplitudeRatio > amplitude_cut){
	    multiPeakChannelsCombined->cd();
	    storeOutputCanvas(canvas,noiseHist,noiseFit,noiseFit2Gaus,name);
	    nbadDoublePeakCombined++;
	    moduleNumeratorDoublePeak[detid]++;
	  }
	}
      }
    }

    if(passConservative and not passTight){
      badCombinedPassingConservativeNotTight->cd();
      nbadCombinedPassingTightNotConservative++;
      moduleNumeratorConservativeNotTight[detid] = moduleNumeratorConservativeNotTight[detid]+1;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);
          
    }

    else if(not passConservative and passTight){
      badCombinedPassingTightNotConservative->cd();
      nbadCombinedPassingConservativeNotTight++;
      moduleNumeratorTightNotConservative[detid] = moduleNumeratorTightNotConservative[detid]+1;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);    
    }

  }
  
  // plot the chi2 and peak distance --> distributions for all these channels tested as possible candidates for 2 peak strips
  if(testDoubleGaussianChannels){
    storeOutputCanvas(canvas,chi2Distance,"chi2TestStatistics",outputDIR);
    storeOutputCanvas(canvas,peakDistance,"peakDistanceTestStatistics",outputDIR);
    storeOutputCanvas(canvas,ashmanDistance,"ashmanTestStatistics",outputDIR);
    storeOutputCanvas(canvas,amplitudeRatioDistance,"amplitudeRatioDistance",outputDIR);
    storeOutputCanvas(canvas,bimodalityDistance,"bimodalityDistance",outputDIR);
  }
  
  std::cout<<std::endl;

  // Closing files
  badStripsNullIntegral->Close();
  badStripsSmallRMS->Close();
  badStripsOverflow->Close();
  badStripsNoiseSignificance->Close();
  badKSTest->Close();
  badADTest->Close();
  badJBTest->Close();
  badChi2Test->Close();
  badCombinedConservativeTest->Close();
  badCombinedPassingTightNotConservative->Close();
  badCombinedPassingConservativeNotTight->Close();
  badCombinedTest->Close();
  badKSNotADTest->Close();
  badJBNotADNotKSTest->Close();
  badChi2NotKSandJBandADTest->Close();
  multiPeakChannelsCombined->Close();

  // Closing files
  multiPeakChannelsChi2->Close();
  multiPeakChannelsDistance->Close();
  multiPeakChannelsAshman->Close();
  multiPeakChannelsAmplitude->Close();
  multiPeakChannelsBimodality->Close();

  // output statistics
  cout<<"#### Bad Null Integral "<<nbadNullIntegral<<" --> "<<double(nbadNullIntegral)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad Small RMS    "<<nbadSmallRMS<<" --> "<<double(nbadSmallRMS)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad Overflow     "<<nbadOverflow<<" --> "<<double(nbadOverflow)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad Noise Significance "<<nbadNoiseSignificance<<" --> "<<double(nbadNoiseSignificance)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad AD Test Channels "<<nbadADTest<<" ---> "<<double(nbadADTest)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad KS Test Channels "<<nbadKSTest<<" ---> "<<double(nbadKSTest)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad JB Test Channels "<<nbadJBTest<<" ---> "<<double(nbadJBTest)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad Chi2 Test Channels "<<nbadChi2Test<<" ---> "<<double(nbadChi2Test)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad KS but not AD Test Channels "<<nbadKSNotADTest<<" ---> "<<double(nbadKSNotADTest)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad JB but not KS and not AD Test Channels "<<nbadJBNotADNotKSTest<<" ---> "<<double(nbadJBNotADNotKSTest)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad Chi2 but not KS and JB and AD Test Channels "<<nbadChi2NotKSandJBandADTest<<" ---> "<<double(nbadChi2NotKSandJBandADTest)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad Combined Test Channels "<<nbadCombinedTest<<" ---> "<<double(nbadCombinedTest)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad Combined Conservative Test Channels "<<nbadCombinedConservativeTest<<" ---> "<<double(nbadCombinedConservativeTest)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad Combined Conservative Not Tight Channels "<<nbadCombinedPassingConservativeNotTight<<" ---> "<<double(nbadCombinedPassingConservativeNotTight)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad Combined Tight Not Conservative Channels "<<nbadCombinedPassingTightNotConservative<<" ---> "<<double(nbadCombinedPassingTightNotConservative)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;

  if(testDoubleGaussianChannels){
    cout<<"###############################"<<endl;
    cout<<"#### Multiple peak finder ####"<<endl;
    cout<<"##############################"<<endl;
    cout<<"Two peak by Chi2 "<<nbadDoublePeakChi2<<" --> "<<double(nbadDoublePeakChi2)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
    cout<<"Two peak by Distance "<<nbadDoublePeakDistance<<" --> "<<double(nbadDoublePeakDistance)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
    cout<<"Two peak by Ashman "<<nbadDoublePeakAshman<<" --> "<<double(nbadDoublePeakAshman)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
    cout<<"Two peak by Amplitude "<<nbadDoublePeakAmplitude<<" --> "<<double(nbadDoublePeakAmplitude)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
    cout<<"Two peak by Bimodality "<<nbadDoublePeakBimodality<<" --> "<<double(nbadDoublePeakBimodality)/(tree->GetEntries()/reductionFactor)<<" % "<<endl;
    cout<<"Two peak by Combined "<<nbadDoublePeakCombined<<" ---> "<<double(nbadDoublePeakCombined)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  }

  //// Text file to plot on the tracker map

  /// -------> 
  ofstream nchannelMapNullIntegral ((outputDIR+"/numberBadChannelsNullIntegral.txt").c_str());
  for(auto module : moduleNumeratorNullIntegral)
    nchannelMapNullIntegral << module.first <<"  "<< double(moduleNumeratorNullIntegral[module.first])/double(moduleDenominator[module.first]) << "\n";
  nchannelMapNullIntegral.close();
  
  /// -------> 
  ofstream nchannelMapSmallRMS ((outputDIR+"/numberBadChannelsSmallRMS.txt").c_str());
  for(auto module : moduleNumeratorSmallRMS)
    nchannelMapSmallRMS << module.first <<"  "<< double(moduleNumeratorSmallRMS[module.first])/double(moduleDenominator[module.first]) << "\n";
  nchannelMapSmallRMS.close();


  /// -------> 
  ofstream nchannelMapOverflow ((outputDIR+"/numberBadChannelsOverflow.txt").c_str());
  for(auto module : moduleNumeratorOverflow)
    nchannelMapOverflow << module.first <<"  "<< double(moduleNumeratorOverflow[module.first])/double(moduleDenominator[module.first]) << "\n";
  nchannelMapOverflow.close();

  /// -------> 
  ofstream nchannelMapNoiseSignificance ((outputDIR+"/numberBadChannelsNoiseSignificance.txt").c_str());
  for(auto module : moduleNumeratorNoiseSignificance)
    nchannelMapNoiseSignificance << module.first <<"  "<< double(moduleNumeratorNoiseSignificance[module.first])/double(moduleDenominator[module.first]) << "\n";
  nchannelMapNoiseSignificance.close();

  /// -------> 
  ofstream nchannelMap ((outputDIR+"/numberBadChannels.txt").c_str());
  for(auto module : moduleNumerator)
    nchannelMap << module.first <<"  "<< double(moduleNumerator[module.first])/double(moduleDenominator[module.first]) << "\n";
  nchannelMap.close();

  /// -------> 
  ofstream nchannelMapConservative ((outputDIR+"/numberBadChannelsConservative.txt").c_str());
  for(auto module : moduleNumeratorConservative)
    nchannelMapConservative << module.first <<"  "<< double(moduleNumeratorConservative[module.first])/double(moduleDenominator[module.first]) << "\n";
  nchannelMapConservative.close();

  /// -------> 
  ofstream nchannelMapKS ((outputDIR+"/numberBadChannelsKS.txt").c_str());
  for(auto module : moduleNumeratorKS)
    nchannelMapKS << module.first <<"  "<< double(moduleNumeratorKS[module.first])/double(moduleDenominator[module.first]) << "\n";
  nchannelMapKS.close();

  /// -------> 
  ofstream nchannelMapJB ((outputDIR+"/numberBadChannelsJB.txt").c_str());
  for(auto module : moduleNumeratorJB)
    nchannelMapJB << module.first <<"  "<< double(moduleNumeratorJB[module.first])/double(moduleDenominator[module.first]) << "\n";
  nchannelMapJB.close();

  /// -------> 
  ofstream nchannelMapAD ((outputDIR+"/numberBadChannelsAD.txt").c_str());
  for(auto module : moduleNumeratorAD)
    nchannelMapAD << module.first <<"  "<< double(moduleNumeratorAD[module.first])/double(moduleDenominator[module.first]) << "\n";
  nchannelMapAD.close();

  /// -------> 
  ofstream nchannelMapDoublePeak ((outputDIR+"/numberBadChannelsDoublePeak.txt").c_str());
  for(auto module : moduleNumeratorDoublePeak)
    nchannelMapDoublePeak << module.first <<"  "<< double(moduleNumeratorDoublePeak[module.first])/double(moduleDenominator[module.first]) << "\n";
  nchannelMapDoublePeak.close();

  /// -------> 
  ofstream nchannelMapTightNotConservative ((outputDIR+"/numberBadChannelsPassingTightNotConservative.txt").c_str());
  for(auto module : moduleNumeratorTightNotConservative)
    nchannelMapTightNotConservative << module.first <<"  "<< double(moduleNumeratorTightNotConservative[module.first])/double(moduleDenominator[module.first]) << "\n";
  nchannelMapTightNotConservative.close();

  /// -------> 
  ofstream nchannelMapConservativeNotTight ((outputDIR+"/numberBadChannelsPassingConservativeNotTight.txt").c_str());
  for(auto module : moduleNumeratorConservativeNotTight)
    nchannelMapConservativeNotTight << module.first <<"  "<< double(moduleNumeratorConservativeNotTight[module.first])/double(moduleDenominator[module.first]) << "\n";
  nchannelMapConservativeNotTight.close();


  ////////// More detailed info
  
  // ------> detailed info of bad strips
  ofstream badStripDump ((outputDIR+"/badStripDump.txt").c_str());
  for(auto badstrip : badStrip){
    badStripDump<< badstrip.detid_ <<" "<<badstrip.lldCh_<<" "<<badstrip.apvid_<<" "<<badstrip.stripid_<<" \n";
  }

  badStripDump.close();

  // ------> detailed info of bad strips
  ofstream badStripDumpConservative ((outputDIR+"/badStripDumpConservative.txt").c_str());
  for(auto badstrip : badStripConservative){
    badStripDumpConservative<< badstrip.detid_ <<" "<<badstrip.lldCh_<<" "<<badstrip.apvid_<<" "<<badstrip.stripid_<<" \n";
  }

  badStripDumpConservative.close();

  // ------> detailed info of bad strips to be compared with offline
  ofstream badStripDumpVsOffline ((outputDIR+"/badStripDumpVsOffline.txt").c_str());
  for(auto badstrip : badStripVsOffline){
    badStripDumpVsOffline<< badstrip.detid_ <<" "<<badstrip.lldCh_<<" "<<badstrip.apvid_<<" "<<badstrip.stripid_<<" \n";
  }

  // ------> detailed info of bad strips to be compared with offline
  ofstream badStripDumpVsOfflineWithSignificance ((outputDIR+"/badStripDumpVsOfflineWithSignificance.txt").c_str());
  for(auto badstrip : badStripVsOfflineWithSignificance){
    badStripDumpVsOfflineWithSignificance<< badstrip.detid_ <<" "<<badstrip.lldCh_<<" "<<badstrip.apvid_<<" "<<badstrip.stripid_<<" \n";
  }
  
  badStripDumpVsOfflineWithSignificance.close();
}
