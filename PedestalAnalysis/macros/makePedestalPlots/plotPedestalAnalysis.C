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
static float maximumMean        = 20;
static float minimumRMS          = 2;
static float maximumRMS          = 30;
static float maximumSignificance = 10;

//to flag strips with long tail
static float minKurtosis = 2;
static float minIntegral5sigma = 0.0005;

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

  // read the input tree
  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  inputFile->cd();
  TTree* tree = (TTree*) inputFile->Get("pedestalFullNoise");

  TCanvas* canvas = new TCanvas("canvas","canvas",600,650);

  // input variables
  uint32_t detid,fedKey;
  uint16_t fecCrate,fecSlot, fecRing, ccuAdd, ccuChan, lldChannel, fedId, fedCh, apvId, stripId;
  int      isNullHisto, isNullFit;
  float    noiseMean, noiseRMS, noiseSkewness, noiseKurtosis, noiseSignificance, noiseIntegral5Sigma;
  float    fitGausMean, fitGausSigma, fitGausNormalization;
  float    fitChi2, fitChi2Probab, kSProbab, jBProbab, aDProbab;
  float    nBin, xMin, xMax;
  float    noise; 
  vector<float>* noiseDistribution = 0;
  vector<float>* noiseDistributionError = 0;

  tree->SetBranchStatus("*",kFALSE);
  tree->SetBranchStatus("detid",kTRUE);
  tree->SetBranchStatus("lldChannel",kTRUE);
  tree->SetBranchStatus("apvId",kTRUE);
  tree->SetBranchStatus("fedId",kTRUE);
  tree->SetBranchStatus("stripId",kTRUE);
  tree->SetBranchStatus("noise",kTRUE);
  tree->SetBranchStatus("fedKey",kTRUE);
  tree->SetBranchStatus("fecCrate",kTRUE);
  tree->SetBranchStatus("fecSlot",kTRUE);
  tree->SetBranchStatus("fecRing",kTRUE);
  tree->SetBranchStatus("ccuChan",kTRUE);
  tree->SetBranchStatus("ccuAdd",kTRUE);
  tree->SetBranchStatus("isNullHisto",kTRUE);
  tree->SetBranchStatus("isNullFit",kTRUE);
  tree->SetBranchStatus("fedCh",kTRUE);
  tree->SetBranchStatus("fitChi2",kTRUE);
  tree->SetBranchStatus("fitChi2Probab",kTRUE);
  tree->SetBranchStatus("kSProbab",kTRUE);
  tree->SetBranchStatus("jBProbab",kTRUE);
  tree->SetBranchStatus("aDProbab",kTRUE);
  tree->SetBranchStatus("fitGausNormalization",kTRUE);
  tree->SetBranchStatus("fitGausMean",kTRUE);
  tree->SetBranchStatus("fitGausSigma",kTRUE);
  tree->SetBranchStatus("noiseMean",kTRUE);
  tree->SetBranchStatus("noiseRMS",kTRUE);
  tree->SetBranchStatus("noiseSignificance",kTRUE);
  tree->SetBranchStatus("noiseSkewness",kTRUE);
  tree->SetBranchStatus("noiseKurtosis",kTRUE);
  tree->SetBranchStatus("noiseIntegral5Sigma",kTRUE);
  tree->SetBranchStatus("noiseDistributionError",kTRUE);
  tree->SetBranchStatus("noiseDistribution",kTRUE);
  tree->SetBranchStatus("nBin",kTRUE);
  tree->SetBranchStatus("xMin",kTRUE);
  tree->SetBranchStatus("xMax",kTRUE);

  tree->SetBranchAddress("detid",&detid);
  tree->SetBranchAddress("lldChannel",&lldChannel);
  tree->SetBranchAddress("apvId",&apvId);
  tree->SetBranchAddress("stripId",&stripId);
  tree->SetBranchAddress("fedId",&fedId);
  tree->SetBranchAddress("noise",&noise);
  tree->SetBranchAddress("fedKey",&fedKey);
  tree->SetBranchAddress("fecCrate",&fecCrate);
  tree->SetBranchAddress("fecSlot",&fecSlot);
  tree->SetBranchAddress("fecRing",&fecRing);
  tree->SetBranchAddress("ccuAdd",&ccuAdd);
  tree->SetBranchAddress("ccuChan",&ccuChan);
  tree->SetBranchAddress("fedId",&fedId);
  tree->SetBranchAddress("fedCh",&fedCh);
  tree->SetBranchAddress("isNullHisto",&isNullHisto);
  tree->SetBranchAddress("isNullFit",&isNullFit);
  tree->SetBranchAddress("fitGausNormalization",&fitGausNormalization);
  tree->SetBranchAddress("fitGausMean",&fitGausMean);
  tree->SetBranchAddress("fitGausSigma",&fitGausSigma);
  tree->SetBranchAddress("fitChi2",&fitChi2);
  tree->SetBranchAddress("fitChi2Probab",&fitChi2Probab);
  tree->SetBranchAddress("noiseMean",&noiseMean);
  tree->SetBranchAddress("noiseRMS",&noiseRMS);
  tree->SetBranchAddress("noiseSignificance",&noiseSignificance);
  tree->SetBranchAddress("noiseIntegral5Sigma",&noiseIntegral5Sigma);
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
  vector<TrackerStrip> badStripAggressive;
  vector<TrackerStrip> badStripConservative;
  vector<TrackerStrip> badStripFinal;
  vector<TrackerStrip> badStripVsOffline;
  vector<TrackerStrip> badStripVsOfflineWithSignificance;

  map<uint32_t,uint32_t> moduleDenominator;

  /// --- first set of bad strips
  // Null integral
  TFile* badStripsNullIntegral = new TFile((outputDIR+"/badStripsNullIntegral.root").c_str(),"RECREATE");
  // Large Mean
  TFile* badStripsLargeMean = new TFile((outputDIR+"/badStripsLargeMean.root").c_str(),"RECREATE");
  // Small RMS
  TFile* badStripsSmallRMS = new TFile((outputDIR+"/badStripsSmallRMS.root").c_str(),"RECREATE");
  // Large RMS
  TFile* badStripsLargeRMS = new TFile((outputDIR+"/badStripsLargeRMS.root").c_str(),"RECREATE");
  // Noise significance
  TFile* badStripsNoiseSignificance = new TFile((outputDIR+"/badStripsNoiseSignificance.root").c_str(),"RECREATE");

  long int nbadNullIntegral = 0;
  long int nbadLargeMean = 0;
  long int nbadSmallRMS  = 0;
  long int nbadLargeRMS  = 0;
  long int nbadNoiseSignificance = 0;

  map<uint32_t,uint32_t> moduleNumeratorNullIntegral;
  map<uint32_t,uint32_t> moduleNumeratorLargeMean;
  map<uint32_t,uint32_t> moduleNumeratorSmallRMS;
  map<uint32_t,uint32_t> moduleNumeratorLargeRMS;
  map<uint32_t,uint32_t> moduleNumeratorNoiseSignificance;

  // ----- Use test statistics  
  // rejected by Anderson Darling test
  TFile* badADTest   = new TFile((outputDIR+"/badStripsADTest.root").c_str(),"RECREATE");
  // rejected by KS test
  TFile* badKSTest   = new TFile((outputDIR+"/badStripsKSTest.root").c_str(),"RECREATE");
  // rejected by JB test
  TFile* badJBTest   = new TFile((outputDIR+"/badStripsJBTest.root").c_str(),"RECREATE");
  // rejected by Chi2 probability
  TFile* badChi2Test = new TFile((outputDIR+"/badStripsChi2Test.root").c_str(),"RECREATE");

  long int nbadKSTest = 0;
  long int nbadJBTest = 0;
  long int nbadADTest = 0;
  long int nbadChi2Test = 0;

  map<uint32_t,uint32_t> moduleNumeratorKS;
  map<uint32_t,uint32_t> moduleNumeratorAD;
  map<uint32_t,uint32_t> moduleNumeratorJB;
  map<uint32_t,uint32_t> moduleNumeratorChi2;

  // studying overlap between methods
  // bad KS but good AD
  TFile* badKSNotADTest       = new TFile((outputDIR+"/badStripsKSNotAD.root").c_str(),"RECREATE");
  // Bad JB but good for AD and KS
  TFile* badJBNotADNotKSTest  = new TFile((outputDIR+"/badStripsJBNotADNotKS.root").c_str(),"RECREATE");
  // bad Chi2 but good JB, good AD and KS
  TFile* badChi2NotKSandJBandADTest = new TFile((outputDIR+"/badStripsChi2NotKSandJBandAD.root").c_str(),"RECREATE");

  long int nbadJBNotADNotKSTest = 0;
  long int nbadKSNotADTest  = 0;
  long int nbadChi2NotKSandJBandADTest = 0;

  // Tagging tails
  TFile* badAsymmetricTails = new TFile((outputDIR+"/badStripsAsymmetricTails.root").c_str(),"RECREATE");

  long int nbadAsymmetricTails = 0;

  map<uint32_t,uint32_t> moduleNumeratorAsymmetricTails;

  // Combine properties ////
  // selected bad strips using different methods
  TFile* badCombinedAggressive = new TFile((outputDIR+"/badStripsCombinedAggressive.root").c_str(),"RECREATE");
  // selected bad strips using different methods
  TFile* badCombinedConservative = new TFile((outputDIR+"/badStripsCombinedConservative.root").c_str(),"RECREATE");
  // strips passing tight but not conservative
  TFile* badCombinedPassingTightNotConservative = new TFile((outputDIR+"/badCombinedPassingTightNotConservative.root").c_str(),"RECREATE");
  
  //// counters 
  long int nbadCombinedAggressive = 0;
  long int nbadCombinedConservative = 0;
  long int nbadCombinedPassingTightNotConservative = 0;
  long int nbadCombinedFinal = 0;
  
  // map of bad channels according to different methods
  map<uint32_t,uint32_t> moduleNumeratorAggressive;
  map<uint32_t,uint32_t> moduleNumeratorConservative;
  map<uint32_t,uint32_t> moduleNumeratorTightNotConservative;
  map<uint32_t,uint32_t> moduleNumeratorFinal;

  // Double peaked distributions
  map<uint32_t,uint32_t> moduleNumeratorDoublePeak;

  long int nbadDoublePeakDistance   = 0;  // distance between two peak
  long int nbadDoublePeakAshman     = 0;  // Ashman distance
  long int nbadDoublePeakChi2       = 0;  // Chi2
  long int nbadDoublePeakAmplitude  = 0;  // Peak amplitude
  long int nbadDoublePeakBimodality = 0;  // Bimodality
  long int nbadDoublePeakCombined   = 0;  // Combined criteria

  TFile* multiPeakChannelsChi2       = new TFile((outputDIR+"/multiPeakChannelsChi2.root").c_str(),"RECREATE");
  TFile* multiPeakChannelsDistance   = new TFile((outputDIR+"/multiPeakChannelsDistance.root").c_str(),"RECREATE");
  TFile* multiPeakChannelsAshman     = new TFile((outputDIR+"/multiPeakChannelsAshman.root").c_str(),"RECREATE");
  TFile* multiPeakChannelsAmplitude  = new TFile((outputDIR+"/multiPeakChannelsAmplitude.root").c_str(),"RECREATE");
  TFile* multiPeakChannelsBimodality = new TFile((outputDIR+"/multiPeakChannelsBimodality.root").c_str(),"RECREATE");
  TFile* multiPeakChannelsCombined   = new TFile((outputDIR+"/multiPeakChannelsCombined.root").c_str(),"RECREATE");


  ///// Useful things
  int   nonNullBins = 0;
  float chi2Ratio  = 0;
  float distance   = 0;
  float ashman     = 0;
  float bimodality = 0;
  float amplitudeRatio = 0;
  
  string fedKeyStr ;
  TString name ;
  TH1F* noiseHist   = NULL;
  TF1*  noiseFit    = NULL;
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
    bool isfound = false;
    for(auto skipfed : skipFEDid){
      if(fedId == skipfed) isfound = true;
    }
    if(isfound) continue;

    /// create noise histogram
    if(noiseHist == NULL){
      noiseHist = new TH1F ("noiseHist","",nBin,xMin,xMax);
      noiseHist->Sumw2();
    }
    noiseHist->Reset();

    for(int iBin = 0; iBin < noiseDistribution->size(); iBin++){
      noiseHist->SetBinContent(iBin+1,noiseDistribution->at(iBin));
      noiseHist->SetBinError(iBin+1,noiseDistributionError->at(iBin));
    }

    // detect for each channel the number of non-null bins --> useful for bimodality test
    nonNullBins = 0;
    for(int iBin = 0; iBin < noiseHist->GetNbinsX(); iBin++){
      if(noiseHist->GetBinContent(iBin+1) != 0) nonNullBins++;
    }

    // Create the gaussiaan fit for each strip
    if(noiseFit == NULL)
      noiseFit = new TF1 ("noiseFist","gaus(0)",xMin,xMax);
    
    // set the parameters from the old fits
    noiseFit->SetRange(xMin,xMax);
    noiseFit->SetParameters(fitGausNormalization,fitGausMean,fitGausSigma);

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

    bool passFinalSelection = false;
    
    //Null integral
    if(noiseHist->Integral() <= 0){
      if(isNullHisto != 1) 
	cerr<<"Histogram with null integral from bins but marked as empty in the analysis -> check "<<endl;
      badStripsNullIntegral->cd();
      nbadNullIntegral++;
      moduleNumeratorNullIntegral[detid] +=1;
      continue;
    }

    // Large mean
    if(fabs(noiseMean) > maximumMean and not passFinalSelection){

      badStripsLargeMean->cd();
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);
      nbadLargeMean++;
      moduleNumeratorLargeMean[detid] += 1;
      
      // strips always declared as bad                                                                                                                                                               
      badStripFinal.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));
      badStripAggressive.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));
      badStripConservative.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));
      badStripVsOffline.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));
      badStripVsOfflineWithSignificance.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));

      // final one
      nbadCombinedFinal++;
      moduleNumeratorFinal[detid] += 1;
      passFinalSelection = true;
    }

    /// very small RMS
    if(noiseRMS < minimumRMS and not passFinalSelection){

      badStripsSmallRMS->cd();
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);
      nbadSmallRMS++;
      moduleNumeratorSmallRMS[detid] +=1;
      
      // strips always declared as bad
      badStripFinal.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));
      badStripAggressive.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));
      badStripConservative.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));
      badStripVsOffline.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));
      badStripVsOfflineWithSignificance.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));

      // final one
      nbadCombinedFinal++;
      moduleNumeratorFinal[detid] += 1;
      passFinalSelection = true;
    }

    /// Large RMS
    if(noiseRMS > maximumRMS and not passFinalSelection){
      badStripsLargeRMS->cd();
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);
      nbadLargeRMS++;
      moduleNumeratorLargeRMS[detid] +=1;

      badStripFinal.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));
      badStripAggressive.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));
      badStripConservative.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));
      badStripVsOffline.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));
      badStripVsOfflineWithSignificance.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));

      // final one
      nbadCombinedFinal++;
      moduleNumeratorFinal[detid] += 1;
      passFinalSelection = true;
    }
    
    // Noise Significance
    if(fabs(noiseSignificance) > maximumSignificance and not passFinalSelection){
      badStripsNoiseSignificance->cd();
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);
      nbadNoiseSignificance++;
      moduleNumeratorNoiseSignificance[detid] +=1;

      badStripFinal.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));
      badStripAggressive.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));
      badStripConservative.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));
      badStripVsOfflineWithSignificance.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));

      // final one
      nbadCombinedFinal++;
      moduleNumeratorFinal[detid] += 1;
      passFinalSelection = true;      
    }

    ////// Looking at the normality tests
  
    // probability for KS test smaller than a given CL --> three sigma means 1% of probability
    if(kSProbab < quantile3sigma and not passFinalSelection){
      badKSTest->cd();
      nbadKSTest++;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);      
      moduleNumeratorKS[detid] += 1;
    }
    
    // probability for the JB test to be smaller than a given CL
    if(jBProbab < quantile5sigma and not passFinalSelection){
      badJBTest->cd();
      nbadJBTest++;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);      
      moduleNumeratorJB[detid] += 1;
    }

    // Chi2 probability
    if(fitChi2Probab < quantile5sigma and not passFinalSelection){
      badChi2Test->cd();
      nbadChi2Test++;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);      
      moduleNumeratorChi2[detid] += 1;
    }

    // Anderson Darling test
    if(aDProbab < quantile3sigma and not passFinalSelection){
      badADTest->cd();
      nbadADTest++;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);
    }

    /// Study complementarity

    // relatively low anderson darling probability but KS less than 1%
    if(kSProbab < quantile3sigma and 
       aDProbab > quantile3sigma and aDProbab < quantile and not passFinalSelection){
      badKSNotADTest->cd();
      nbadKSNotADTest++;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);      
      moduleNumeratorKS[detid] += 1;
    }

    //// bad for JB but not for AD (AD sufficiently low [-CL,1sigma]
    if(jBProbab < quantile5sigma and 
       aDProbab > quantile3sigma and aDProbab < quantile and 
       kSProbab > quantile3sigma and not passFinalSelection){
      badJBNotADNotKSTest->cd();
      nbadJBNotADNotKSTest++;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);      
      moduleNumeratorAD[detid] += 1;
    }


    if(fitChi2Probab < quantile5sigma and 
       jBProbab > quantile5sigma and 
       aDProbab > quantile3sigma and aDProbab < quantile and 
       kSProbab > quantile3sigma and not passFinalSelection){
      badChi2NotKSandJBandADTest->cd();
      nbadChi2NotKSandJBandADTest++;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);      
    }

    // combining all of them --> conservative approach --> all the methods should show a low p-value to flag a strips as bad
    bool passConservative = false;
    if(aDProbab < quantile3sigma and kSProbab < quantile3sigma and jBProbab < quantile3sigma and fitChi2Probab < quantile3sigma and not passFinalSelection){
      passConservative = true;
      badCombinedConservative->cd();
      nbadCombinedConservative++;
      moduleNumeratorConservative[detid] = moduleNumeratorConservative[detid]+1;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);      
      // store it in order to dispaly on the tracker map
      badStripConservative.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));
      // final one
      nbadCombinedFinal++;
      moduleNumeratorFinal[detid] += 1;
    }

    //// combining all of them --> aggressive approach
    bool passTight = false;
    if((aDProbab < quantile3sigma or 
	(aDProbab > quantile3sigma and aDProbab < quantile and kSProbab < quantile3sigma) or 
	(aDProbab > quantile3sigma and aDProbab < quantile and kSProbab > quantile3sigma and jBProbab < quantile5sigma)) and not passFinalSelection){
	 
      passTight = true;
      badCombinedAggressive->cd();
      nbadCombinedAggressive++;
      moduleNumeratorAggressive[detid] = moduleNumeratorAggressive[detid]+1;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);      

      // store it in order to dispaly on the tracker map
      badStripAggressive.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,detid,apvId,stripId));      
    }

    // passing tight but nont conservative
    if(not passConservative and passTight and not passFinalSelection){
      badCombinedPassingTightNotConservative->cd();
      nbadCombinedPassingTightNotConservative++;
      moduleNumeratorTightNotConservative[detid] = moduleNumeratorTightNotConservative[detid]+1;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);    
    }
    
    // update here the status
    if(passConservative) passFinalSelection = true;

    // identify asymm tails
    bool passTail = false;
    if(noiseKurtosis > minKurtosis and noiseIntegral5Sigma > minIntegral5sigma and not passFinalSelection){
      passTail = true;
      badAsymmetricTails->cd();
      nbadAsymmetricTails++;
      moduleNumeratorAsymmetricTails[detid] = moduleNumeratorAsymmetricTails[detid]+1;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);

      // final one
      nbadCombinedFinal++;
      moduleNumeratorFinal[detid] += 1;

    }
    
    if(passTail) passFinalSelection = true;
    
    // try to identify double peaked strips
    if(passFinalSelection and testDoubleGaussianChannels){
      
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
	
	//distance between peaks
	distance = fabs(noiseFit2Gaus->GetParameter(1)-noiseFit2Gaus->GetParameter(4))/(2*sqrt(noiseFit2Gaus->GetParameter(2)*noiseFit2Gaus->GetParameter(5)));
	
	// ashman distance
	ashman   = TMath::Power(2,0.5)*abs(noiseFit2Gaus->GetParameter(1)-noiseFit2Gaus->GetParameter(4))/(sqrt(pow(noiseFit2Gaus->GetParameter(2),2)+pow(noiseFit2Gaus->GetParameter(5),2)));

	// bimodality coefficient
	if(nonNullBins > 3)
	  bimodality = (noiseHist->GetSkewness()*noiseHist->GetSkewness()+1)/(noiseHist->GetKurtosis()+3*(nonNullBins-1)*(nonNullBins-1)/((nonNullBins-2)*(nonNullBins-3)));
	else
	    bimodality = (noiseHist->GetSkewness()*noiseHist->GetSkewness()+1)/(noiseHist->GetKurtosis());
	
	// ratio of amplitudes
	amplitudeRatio = std::min(noiseFit2Gaus->GetParameter(0),noiseFit2Gaus->GetParameter(3))/std::max(noiseFit2Gaus->GetParameter(0),noiseFit2Gaus->GetParameter(3));
	
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

  // plot the chi2 and peak distance --> distributions for all these channels tested as possible candidates for 2 peak strips
  std::cout<<std::endl;

  // Closing files
  badStripsNullIntegral->Close();
  badStripsSmallRMS->Close();
  badStripsLargeRMS->Close();
  badStripsNoiseSignificance->Close();
  badKSTest->Close();
  badADTest->Close();
  badJBTest->Close();
  badChi2Test->Close();
  badKSNotADTest->Close();
  badJBNotADNotKSTest->Close();
  badChi2NotKSandJBandADTest->Close();
  badCombinedConservative->Close();
  badCombinedAggressive->Close();
  badCombinedPassingTightNotConservative->Close();
  
  // Closing files
  multiPeakChannelsChi2->Close();
  multiPeakChannelsDistance->Close();
  multiPeakChannelsAshman->Close();
  multiPeakChannelsAmplitude->Close();
  multiPeakChannelsBimodality->Close();
  multiPeakChannelsCombined->Close();

  // output sstatistics
  cout<<"#### Bad Null Integral "<<nbadNullIntegral<<" --> "<<double(nbadNullIntegral)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad Large Mean    "<<nbadLargeMean<<" --> "<<double(nbadLargeMean)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad Small RMS    "<<nbadSmallRMS<<" --> "<<double(nbadSmallRMS)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad Large RMS    "<<nbadLargeRMS<<" --> "<<double(nbadLargeRMS)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad Noise Significance "<<nbadNoiseSignificance<<" --> "<<double(nbadNoiseSignificance)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad AD Test Channels "<<nbadADTest<<" ---> "<<double(nbadADTest)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad KS Test Channels "<<nbadKSTest<<" ---> "<<double(nbadKSTest)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad JB Test Channels "<<nbadJBTest<<" ---> "<<double(nbadJBTest)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad Chi2 Test Channels "<<nbadChi2Test<<" ---> "<<double(nbadChi2Test)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad KS but not AD Test Channels "<<nbadKSNotADTest<<" ---> "<<double(nbadKSNotADTest)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad JB but not KS and not AD Test Channels "<<nbadJBNotADNotKSTest<<" ---> "<<double(nbadJBNotADNotKSTest)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad Chi2 but not KS and JB and AD Test Channels "<<nbadChi2NotKSandJBandADTest<<" ---> "<<double(nbadChi2NotKSandJBandADTest)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad Combined Aggressive Channels "<<nbadCombinedAggressive<<" ---> "<<double(nbadCombinedAggressive)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad Combined Conservative Channels "<<nbadCombinedConservative<<" ---> "<<double(nbadCombinedConservative)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad Combined Tight Not Conservative Channels "<<nbadCombinedPassingTightNotConservative<<" ---> "<<double(nbadCombinedPassingTightNotConservative)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad Asymmetric Tails "<<nbadAsymmetricTails<<" ---> "<<double(nbadAsymmetricTails)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  cout<<"#### Bad Combined Final Channels "<<nbadCombinedFinal<<" ---> "<<double(nbadCombinedFinal)/(tree->GetEntries()/reductionFactor)*100<<" % "<<endl;
  

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
  ofstream nchannelMapLargeMean ((outputDIR+"/numberBadChannelsLargeMean.txt").c_str());
  for(auto module : moduleNumeratorLargeMean)
    nchannelMapLargeMean << module.first <<"  "<< double(moduleNumeratorLargeMean[module.first])/double(moduleDenominator[module.first]) << "\n";
  nchannelMapLargeMean.close();
  
  /// -------> 
  ofstream nchannelMapSmallRMS ((outputDIR+"/numberBadChannelsSmallRMS.txt").c_str());
  for(auto module : moduleNumeratorSmallRMS)
    nchannelMapSmallRMS << module.first <<"  "<< double(moduleNumeratorSmallRMS[module.first])/double(moduleDenominator[module.first]) << "\n";
  nchannelMapSmallRMS.close();

  /// -------> 
  ofstream nchannelMapLargeRMS ((outputDIR+"/numberBadChannelsLargeRMS.txt").c_str());
  for(auto module : moduleNumeratorLargeRMS)
    nchannelMapLargeRMS << module.first <<"  "<< double(moduleNumeratorLargeRMS[module.first])/double(moduleDenominator[module.first]) << "\n";
  nchannelMapLargeRMS.close();

  /// -------> 
  ofstream nchannelMapNoiseSignificance ((outputDIR+"/numberBadChannelsNoiseSignificance.txt").c_str());
  for(auto module : moduleNumeratorNoiseSignificance)
    nchannelMapNoiseSignificance << module.first <<"  "<< double(moduleNumeratorNoiseSignificance[module.first])/double(moduleDenominator[module.first]) << "\n";
  nchannelMapNoiseSignificance.close();

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
  ofstream nchannelMapAggressive ((outputDIR+"/numberBadChannelsAggressive.txt").c_str());
  for(auto module : moduleNumeratorAggressive)
    nchannelMapAggressive << module.first <<"  "<< double(moduleNumeratorAggressive[module.first])/double(moduleDenominator[module.first]) << "\n";
  nchannelMapAggressive.close();

  /// -------> 
  ofstream nchannelMapConservative ((outputDIR+"/numberBadChannelsConservative.txt").c_str());
  for(auto module : moduleNumeratorConservative)
    nchannelMapConservative << module.first <<"  "<< double(moduleNumeratorConservative[module.first])/double(moduleDenominator[module.first]) << "\n";
  nchannelMapConservative.close();

  /// -------> 
  ofstream nchannelMapTightNotConservative ((outputDIR+"/numberBadChannelsPassingTightNotConservative.txt").c_str());
  for(auto module : moduleNumeratorTightNotConservative)
    nchannelMapTightNotConservative << module.first <<"  "<< double(moduleNumeratorTightNotConservative[module.first])/double(moduleDenominator[module.first]) << "\n";
  nchannelMapTightNotConservative.close();


  /// -------> 
  ofstream nchannelMapDoublePeak ((outputDIR+"/numberBadChannelsDoublePeak.txt").c_str());
  for(auto module : moduleNumeratorDoublePeak)
    nchannelMapDoublePeak << module.first <<"  "<< double(moduleNumeratorDoublePeak[module.first])/double(moduleDenominator[module.first]) << "\n";
  nchannelMapDoublePeak.close();


  /// -------> 
  ofstream nchannelMapAsymmetricTails ((outputDIR+"/numberBadChannelsAsymmetricTails.txt").c_str());
  for(auto module : moduleNumeratorAsymmetricTails)
    nchannelMapAsymmetricTails << module.first <<"  "<< double(moduleNumeratorAsymmetricTails[module.first])/double(moduleDenominator[module.first]) << "\n";
  nchannelMapAsymmetricTails.close();

  /// -------> 
  ofstream nchannelMapFinal ((outputDIR+"/numberBadChannelsFinal.txt").c_str());
  for(auto module : moduleNumeratorFinal)
    nchannelMapFinal << module.first <<"  "<< double(moduleNumeratorFinal[module.first])/double(moduleDenominator[module.first]) << "\n";
  nchannelMapFinal.close();

  ////////// More detailed info  

  // ------> detailed info of bad strips
  ofstream badStripDumpFinal ((outputDIR+"/badStripDumpFinal.txt").c_str());
  for(auto badstrip : badStripFinal){
    badStripDumpFinal<< badstrip.detid_ <<" "<<badstrip.lldCh_<<" "<<badstrip.apvid_<<" "<<badstrip.stripid_<<" \n";
  }

  badStripDumpFinal.close();

  // ------> detailed info of bad strips
  ofstream badStripDumpConservative ((outputDIR+"/badStripDumpConservative.txt").c_str());
  for(auto badstrip : badStripConservative){
    badStripDumpConservative<< badstrip.detid_ <<" "<<badstrip.lldCh_<<" "<<badstrip.apvid_<<" "<<badstrip.stripid_<<" \n";
  }

  badStripDumpConservative.close();

  // ------> detailed info of bad strips
  ofstream badStripDumpAggressive ((outputDIR+"/badStripDumpAggressive.txt").c_str());
  for(auto badstrip : badStripAggressive){
    badStripDumpAggressive<< badstrip.detid_ <<" "<<badstrip.lldCh_<<" "<<badstrip.apvid_<<" "<<badstrip.stripid_<<" \n";
  }

  badStripDumpAggressive.close();

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
