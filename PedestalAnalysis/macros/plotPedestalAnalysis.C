#include "CMS_lumi.h"
#include "TrackerStrip.h"

static float quantile       = 0.3173;
static float quantile1sigma = 0.317310507863;
static float quantile2sigma = 0.045500263896;
static float quantile3sigma = 0.002699796063;
static float quantile4sigma = 0.000063342484;
static float quantile5sigma = 0.000000573303;

static vector<uint16_t> skipFEDid = {75};
static int reductionFactor     = 1;
static int nFilledBinSelection = 10;
static float minimumRMS        = 1.5;

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

void plotPedestalAnalysis(string inputFileName, string outputDIR, bool testDoubleGaussianChannels){

  system(("mkdir -p "+outputDIR).c_str());

  setTDRStyle();
  gROOT->SetBatch(kTRUE);

  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  inputFile->cd();
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

  tree->SetBranchStatus("*",kFALSE);
  tree->SetBranchStatus("detid",kTRUE);
  tree->SetBranchStatus("fedKey",kTRUE);
  tree->SetBranchStatus("fecCrate",kTRUE);
  tree->SetBranchStatus("fecSlot",kTRUE);
  tree->SetBranchStatus("fecRing",kTRUE);
  tree->SetBranchStatus("ccuAdd",kTRUE);
  tree->SetBranchStatus("ccuChan",kTRUE);
  tree->SetBranchStatus("lldChannel",kTRUE);
  tree->SetBranchStatus("fedId",kTRUE);
  tree->SetBranchStatus("fedCh",kTRUE);
  tree->SetBranchStatus("apvId",kTRUE);
  tree->SetBranchStatus("stripId",kTRUE);
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
  tree->SetBranchStatus("noiseDistribution",kTRUE);
  tree->SetBranchStatus("noiseDistributionError",kTRUE);
  tree->SetBranchStatus("nBin",kTRUE);
  tree->SetBranchStatus("xMin",kTRUE);
  tree->SetBranchStatus("xMax",kTRUE);

  tree->SetBranchAddress("detid",&detid);
  tree->SetBranchAddress("fedKey",&fedKey);
  tree->SetBranchAddress("fecCrate",&fecCrate);
  tree->SetBranchAddress("fecSlot",&fecSlot);
  tree->SetBranchAddress("fecRing",&fecRing);
  tree->SetBranchAddress("ccuAdd",&ccuAdd);
  tree->SetBranchAddress("ccuChan",&ccuChan);
  tree->SetBranchAddress("lldChannel",&lldChannel);
  tree->SetBranchAddress("fedId",&fedId);
  tree->SetBranchAddress("fedCh",&fedCh);
  tree->SetBranchAddress("apvId",&apvId);
  tree->SetBranchAddress("stripId",&stripId);
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
  tree->SetBranchAddress("noiseDistribution",&noiseDistribution);
  tree->SetBranchAddress("noiseDistributionError",&noiseDistributionError);
  tree->SetBranchAddress("nBin",&nBin);
  tree->SetBranchAddress("xMin",&xMin);
  tree->SetBranchAddress("xMax",&xMax);

  // file for bad identified strips using Anderson Darling test as benchmark
  TFile* badStripsNFilledBins = new TFile((outputDIR+"/badStripsNFilledBins.root").c_str(),"RECREATE");
  TFile* badADTest   = new TFile((outputDIR+"/badStripsADTest.root").c_str(),"RECREATE");
  TFile* badKSTest   = new TFile((outputDIR+"/badStripsKSTest.root").c_str(),"RECREATE");
  TFile* badJBTest   = new TFile((outputDIR+"/badStripsJBTest.root").c_str(),"RECREATE");
  TFile* badChi2Test = new TFile((outputDIR+"/badStripsChi2Test.root").c_str(),"RECREATE");
  TFile* badCombinedTest = new TFile((outputDIR+"/badStripsCombined.root").c_str(),"RECREATE");
  TFile* badJBNotADNotKSTest  = new TFile((outputDIR+"/badStripsJBNotADNotKS.root").c_str(),"RECREATE");
  TFile* badKSNotADTest       = new TFile((outputDIR+"/badStripKSNotAD.root").c_str(),"RECREATE");
  TFile* badChi2NotKSandJBandADTest = new TFile((outputDIR+"/badStripsChi2NotKSandJBandAD.root").c_str(),"RECREATE");
  ////
  long int nbadNFilledBins = 0;
  long int nbadKSTest = 0;
  long int nbadJBTest = 0;
  long int nbadADTest = 0;
  long int nbadChi2Test = 0;
  long int nbadCombinedTest = 0;
  long int nbadJBNotADNotKSTest = 0;
  long int nbadKSNotADTest  = 0;
  long int nbadChi2NotKSandJBandADTest = 0;

  TCanvas* canvas = new TCanvas("canvas","canvas",600,650);

  map<uint32_t,uint32_t> moduleDenominator;
  map<uint32_t,uint32_t> moduleNumerator;
  map<uint32_t,uint32_t> moduleNumeratorFilledBins;
  map<uint32_t,uint32_t> moduleNumeratorKS;
  map<uint32_t,uint32_t> moduleNumeratorAD;
  map<uint32_t,uint32_t> moduleNumeratorJB;
  map<uint32_t,uint32_t> moduleNumeratorChi2;
  map<uint32_t,uint32_t> moduleNumeratorDoublePeak;
  
  /////
  vector<TrackerStrip> badStrip;

  //////
  long int nbadDoublePeakDistance   = 0;
  long int nbadDoublePeakAshman     = 0;
  long int nbadDoublePeakChi2       = 0;
  long int nbadDoublePeakAmplitude  = 0;
  long int nbadDoublePeakBimodality = 0;
  long int nbadDoublePeakCombined   = 0;

  TH1F* noiseHist      = NULL;
  TF1*  noiseFit       = NULL;
  TF1*  noiseFit2Gaus  = NULL;
  TFitResultPtr result;

  TH1F* chi2Distance    = new TH1F("chi2Distance","",100,0,1);
  chi2Distance->Sumw2();
  TH1F* peakDistance    = new TH1F("peakDistance","",100,0,3);
  peakDistance->Sumw2();
  TH1F* ashmanDistance  = new TH1F("ashmanDistance","",100,0,5);
  ashmanDistance->Sumw2();
  TH1F* bimodalityDistance = new TH1F("bimodalityDistance","",100,0,1);
  bimodalityDistance->Sumw2();
  TH1F* amplitudeRatioDistance = new TH1F("amplitudeRatioDistance","",100,0,3);
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

  bool isfound = false;
  string fedKeyStr ;
  TString name ;
  std::map<string,string> fitParam;
 
  // loop on the pedestak analysis tree
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
    stringstream sKS;
    sKS << std::scientific << kSProbab;
    fitParam["kSProbab"] = sKS.str();
    stringstream sJB;
    sJB << std::scientific << jBProbab;
    fitParam["jBProbab"] = sJB.str();
    stringstream sChi2;
    sChi2 << std::scientific << fitChi2Probab;
    fitParam["fitChi2Probab"] = sChi2.str();
    stringstream sAD;
    sAD << std::scientific << aDProbab;
    fitParam["aDProbab"] = sAD.str();

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

    /// skip channels with too few filled bins or a extremely small RMS
    if(nonNullBins < nFilledBinSelection or noiseHist->GetRMS() < minimumRMS){
      badStripsNFilledBins->cd();
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);
      nbadNFilledBins++;
      moduleNumeratorFilledBins[detid] +=1;
      continue;
    }


    // probability for KS test smaller than a given CL
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
    if(fitChi2Probab < quantile4sigma){
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


    if(fitChi2Probab < quantile4sigma and jBProbab > quantile5sigma and aDProbab > quantile3sigma and aDProbab < quantile and kSProbab > quantile3sigma){
      badChi2NotKSandJBandADTest->cd();
      nbadChi2NotKSandJBandADTest++;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);      
    }


    //// combining all of them
    if(aDProbab < quantile3sigma or 
       (aDProbab > quantile3sigma and aDProbab < quantile and kSProbab < quantile3sigma) or 
       (aDProbab > quantile3sigma and aDProbab < quantile and kSProbab > quantile3sigma and jBProbab < quantile5sigma)  or 
       (aDProbab > quantile3sigma and aDProbab < quantile and kSProbab > quantile3sigma and jBProbab > quantile5sigma and fitChi2Probab < quantile4sigma)){

      badCombinedTest->cd();
      nbadCombinedTest++;
      moduleNumerator[detid] = moduleNumerator[detid]+1;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);      
      name = Form("fecCrate%d_fecSlot%d_fecRing%d_ccuAdd%d_ccuCh%d_fedKey0x0000%s_lldCh%d_apv%d_strip%d",fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,fedKeyStr.c_str(),lldChannel,apvId,stripId);      

      // store it in order to dispaly on the tracker map
      badStrip.push_back(TrackerStrip(fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,uint32_t(atoi(fedKeyStr.c_str())),lldChannel,apvId,stripId));

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
	if(result.Get() or noiseHist->Integral() == 0){

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
	  if(distance > 1){
	    multiPeakChannelsDistance->cd();
	    storeOutputCanvas(canvas,noiseHist,noiseFit,noiseFit2Gaus,name);
	    nbadDoublePeakDistance++;
	  }
	  /// ashman coefficient
	  if(ashman > 2){
	    multiPeakChannelsAshman->cd();
	    storeOutputCanvas(canvas,noiseHist,noiseFit,noiseFit2Gaus,name);
	    nbadDoublePeakAshman++;
	  }
	  /// chi2 ratio
	  if(chi2Ratio < 0.05){
	    multiPeakChannelsChi2->cd();
	    storeOutputCanvas(canvas,noiseHist,noiseFit,noiseFit2Gaus,name);
	    nbadDoublePeakChi2++;
	  }
	  /// amplitude ratios
	  if(amplitudeRatio > 0.85){
	    multiPeakChannelsAmplitude->cd();
	    storeOutputCanvas(canvas,noiseHist,noiseFit,noiseFit2Gaus,name);
	    nbadDoublePeakAmplitude++;
	  }
	  //// bimodality
	  if(bimodality > 0.55){
	    multiPeakChannelsBimodality->cd();
	    storeOutputCanvas(canvas,noiseHist,noiseFit,noiseFit2Gaus,name);
	    nbadDoublePeakBimodality++;
	  }
	  /// combo of ashman and amplitude ratio
	  if(ashman > 2 && amplitudeRatio > 0.85){
	    multiPeakChannelsCombined->cd();
	    storeOutputCanvas(canvas,noiseHist,noiseFit,noiseFit2Gaus,name);
	    nbadDoublePeakCombined++;
	    moduleNumeratorDoublePeak[detid]++;
	  }
	}
      }
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

  badStripsNFilledBins->Close();
  badKSTest->Close();
  badADTest->Close();
  badJBTest->Close();
  badChi2Test->Close();
  badCombinedTest->Close();
  badKSNotADTest->Close();
  badJBNotADNotKSTest->Close();
  badChi2NotKSandJBandADTest->Close();
  multiPeakChannelsCombined->Close();
  //////
  multiPeakChannelsChi2->Close();
  multiPeakChannelsDistance->Close();
  multiPeakChannelsAshman->Close();
  multiPeakChannelsAmplitude->Close();
  multiPeakChannelsBimodality->Close();

  cout<<"#### Bad Nfilled bins "<<nbadNFilledBins<<" --> "<<double(nbadNFilledBins)/(tree->GetEntries()/reductionFactor)<<" % "<<endl;
  cout<<"#### Bad AD Test Channels "<<nbadADTest<<" ---> "<<double(nbadADTest)/(tree->GetEntries()/reductionFactor)<<" % "<<endl;
  cout<<"#### Bad KS Test Channels "<<nbadKSTest<<" ---> "<<double(nbadKSTest)/(tree->GetEntries()/reductionFactor)<<" % "<<endl;
  cout<<"#### Bad JB Test Channels "<<nbadJBTest<<" ---> "<<double(nbadJBTest)/(tree->GetEntries()/reductionFactor)<<" % "<<endl;
  cout<<"#### Bad Chi2 Test Channels "<<nbadChi2Test<<" ---> "<<double(nbadChi2Test)/(tree->GetEntries()/reductionFactor)<<" % "<<endl;
  cout<<"#### Bad KS but not AD Test Channels "<<nbadKSNotADTest<<" ---> "<<double(nbadKSNotADTest)/(tree->GetEntries()/reductionFactor)<<" % "<<endl;
  cout<<"#### Bad JB but not KS and not AD Test Channels "<<nbadJBNotADNotKSTest<<" ---> "<<double(nbadJBNotADNotKSTest)/(tree->GetEntries()/reductionFactor)<<" % "<<endl;
  cout<<"#### Bad Chi2 but not KS and JB and AD Test Channels "<<nbadChi2NotKSandJBandADTest<<" ---> "<<double(nbadChi2NotKSandJBandADTest)/(tree->GetEntries()/reductionFactor)<<" % "<<endl;
  cout<<"#### Bad Combined Test Channels "<<nbadCombinedTest<<" ---> "<<double(nbadCombinedTest)/(tree->GetEntries()/reductionFactor)<<" % "<<endl;

  if(testDoubleGaussianChannels){
    cout<<"###############################"<<endl;
    cout<<"#### Multiple peak finder ####"<<endl;
    cout<<"##############################"<<endl;
    cout<<"Two peak by Chi2 "<<nbadDoublePeakChi2<<" --> "<<double(nbadDoublePeakChi2)/(tree->GetEntries()/reductionFactor)<<" % "<<endl;
    cout<<"Two peak by Distance "<<nbadDoublePeakDistance<<" --> "<<double(nbadDoublePeakDistance)/(tree->GetEntries()/reductionFactor)<<" % "<<endl;
    cout<<"Two peak by Ashman "<<nbadDoublePeakAshman<<" --> "<<double(nbadDoublePeakAshman)/(tree->GetEntries()/reductionFactor)<<" % "<<endl;
    cout<<"Two peak by Amplitude "<<nbadDoublePeakAmplitude<<" --> "<<double(nbadDoublePeakAmplitude)/(tree->GetEntries()/reductionFactor)<<" % "<<endl;
    cout<<"Two peak by Bimodality "<<nbadDoublePeakBimodality<<" --> "<<double(nbadDoublePeakBimodality)/(tree->GetEntries()/reductionFactor)<<" % "<<endl;
    cout<<"Two peak by Combibed "<<nbadDoublePeakCombined<<" ---> "<<double(nbadDoublePeakCombined)/(tree->GetEntries()/reductionFactor)<<" % "<<endl;
  }

  /// -------> 
  ofstream channelMap ((outputDIR+"/fractionOfGoodChannels.txt").c_str());
  for(auto module : moduleDenominator)
    channelMap << module.first <<"  "<< 1. - double(moduleNumerator[module.first])/double(moduleDenominator[module.first]) << "\n";
  channelMap.close();

  /// -------> 
  ofstream nchannelMapFilledBins ((outputDIR+"/numberBadChannelsFilledBins.txt").c_str());
  for(auto module : moduleNumeratorFilledBins)
    nchannelMapFilledBins << module.first <<"  "<< moduleNumeratorFilledBins[module.first] << "\n";
  nchannelMapFilledBins.close();

  /// -------> 
  ofstream nchannelMap ((outputDIR+"/numberBadChannels.txt").c_str());
  for(auto module : moduleNumerator)
    nchannelMap << module.first <<"  "<< moduleNumerator[module.first] << "\n";
  nchannelMap.close();

  /// -------> 
  //  ofstream nNonEdgechannelMap ((outputDIR+"/numberBadChannelsNonEdge.txt").c_str());
  //  for(auto module : moduleNonEdgeNumerator)
  //    nNonEdgechannelMap << module.first <<"  "<< moduleNonEdgeNumerator[module.first] << "\n";
  //  nNonEdgechannelMap.close();

  /// -------> 
  //  ofstream nAPVEdgechannelMap ((outputDIR+"/numberBadChannelsAPVEdge.txt").c_str());
  //  for(auto module : moduleAPVEdgeNumerator)
  //    nAPVEdgechannelMap << module.first <<"  "<< moduleAPVEdgeNumerator[module.first] << "\n";
  //  nAPVEdgechannelMap.close();

  /// -------> 
  ofstream nchannelMapKS ((outputDIR+"/numberBadChannelsKS.txt").c_str());
  for(auto module : moduleNumeratorKS)
    nchannelMapKS << module.first <<"  "<< moduleNumeratorKS[module.first] << "\n";
  nchannelMapKS.close();

  /// -------> 
  ofstream nchannelMapJB ((outputDIR+"/numberBadChannelsJB.txt").c_str());
  for(auto module : moduleNumeratorJB)
    nchannelMapJB << module.first <<"  "<< moduleNumeratorJB[module.first] << "\n";
  nchannelMapJB.close();

  /// -------> 
  ofstream nchannelMapAD ((outputDIR+"/numberBadChannelsAD.txt").c_str());
  for(auto module : moduleNumeratorAD)
    nchannelMapAD << module.first <<"  "<< moduleNumeratorAD[module.first] << "\n";
  nchannelMapAD.close();

  /// -------> 
  ofstream nchannelMapDoublePeak ((outputDIR+"/numberBadChannelsDoublePeak.txt").c_str());
  for(auto module : moduleNumeratorDoublePeak)
    nchannelMapDoublePeak << module.first <<"  "<< moduleNumeratorDoublePeak[module.first] << "\n";
  nchannelMapDoublePeak.close();
  
  // ------> detailed info of bad strips
  ofstream badStripDump ((outputDIR+"/badStripDump.txt").c_str());
  for(auto badstrip : badStrip){
    badStripDump<< badstrip.fecCrate_<<" "<<badstrip.fecSlot_<<" "<<badstrip.fecRing_<<" "<<badstrip.ccuAdd_<<" "<<badstrip.ccuCh_<<" "<<badstrip.fedKey_<<" "<<badstrip.lldCh_<<" "<<badstrip.apvid_<<" "<<badstrip.stripid_<<" \n";
  }
  badStripDump.close();

}
