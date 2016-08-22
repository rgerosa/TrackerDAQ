#include "CMS_lumi.h"

static float quantile       = 0.5;
static float quantile1sigma = 0.317310507863;
static float quantile2sigma = 0.045500263896;
static float quantile3sigma = 0.002699796063;
static float quantile4sigma = 0.000063342484;
static float quantile5sigma = 0.000000573303;

static int  reductionFactor = 1;

void storeOutputCanvas(TCanvas* canvas, TH1F & noiseDistribution,  TF1 & fitFunction, const TString & name, map<string,string> & parameters){

  canvas->cd();
  noiseDistribution.SetMarkerColor(kBlack);
  noiseDistribution.SetMarkerStyle(20);
  noiseDistribution.SetMarkerSize(1);
  noiseDistribution.GetXaxis()->SetTitle("noise (ADC)");
  noiseDistribution.GetYaxis()->SetTitle("Events");
  noiseDistribution.Draw("EP");
  CMS_lumi(canvas,"",true);
  
  fitFunction.SetLineColor(kRed);
  fitFunction.SetLineWidth(2);
  fitFunction.Draw("Lsame");
  noiseDistribution.Draw("EPsame");

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


void plotPedestalAnalysis(string inputFileName, string outputDIR){

  system(("mkdir -p "+outputDIR).c_str());

  setTDRStyle();
  gROOT->SetBatch(kTRUE);

  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  inputFile->cd();
  TTree* tree = (TTree*) inputFile->Get("pedestalFullNoise");

  uint32_t detid,fedKey;
  uint16_t fecCrate,fecSlot, fecRing, ccuAdd, ccuChan, lldChannel, fedId, fedCh, apvId, stripId;
  float    fitChi2Probab, kSProbab, jBProbab, aDProbab;
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

  TFile* badKsTest = new TFile((outputDIR+"/badStripsKsTest.root").c_str(),"RECREATE");
  TFile* badjBTest = new TFile((outputDIR+"/badStripsjBTest.root").c_str(),"RECREATE");
  TFile* badChi2Test = new TFile((outputDIR+"/badStripsChi2Test.root").c_str(),"RECREATE");
  TFile* badaDTest = new TFile((outputDIR+"/badStripsaDTest.root").c_str(),"RECREATE");
  TFile* badCombinedTest = new TFile((outputDIR+"/badStripsCombined.root").c_str(),"RECREATE");
  TFile* badJBNotKSTest = new TFile((outputDIR+"/badStripsjBNotKS.root").c_str(),"RECREATE");
  TFile* badaDNotKSandjBTest = new TFile((outputDIR+"/badStripaDNotKSNotjB.root").c_str(),"RECREATE");
  TFile* badChi2NotKSandjBandaDTest = new TFile((outputDIR+"/badStripsChi2NotKsandjBandaD.root").c_str(),"RECREATE");

  long int nbadKsTest = 0;
  long int nbadjBTest = 0;
  long int nbadaDTest = 0;
  long int nbadChi2Test = 0;
  long int nbadCombinedTest = 0;
  long int nbadJBNotKSTest = 0;
  long int nbadaDNotKSandjBTest = 0;
  long int nbadChi2NotKSandjBandaDTest = 0;

  TCanvas* canvas = new TCanvas("canvas","canvas",600,650);

  map<uint32_t,uint32_t> moduleDenominator;
  map<uint32_t,uint32_t> moduleNumerator;
 
  for(long int iChannel = 0; iChannel < tree->GetEntries(); iChannel++){
    tree->GetEntry(iChannel);
    cout.flush();
    if(iChannel %10000 == 0) cout<<"\r"<<"iChannel "<<100*double(iChannel)/(tree->GetEntries()/reductionFactor)<<" % ";
    if(iChannel > double(tree->GetEntries())/reductionFactor) break;

    // make selections to identify bad noisy channels (not gaussian ones)
    std::stringstream stream;
    stream << std::hex << fedKey;
    string fedKeyStr = stream.str();
    TString name ;
    if(fedKeyStr.size() == 4)
      name = Form("fecCrate%d_fecSlot%d_fecRing%d_ccuAdd%d_ccuCh%d_fedKey0x0000%s_lldCh%d_apv%d_strip%d",fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,fedKeyStr.c_str(),lldChannel,apvId,stripId);
    else if(fedKeyStr.size() == 5)
      name = Form("fecCrate%d_fecSlot%d_fecRing%d_ccuAdd%d_ccuCh%d_fedKey0x000%s_lldCh%d_apv%d_strip%d",fecCrate,fecSlot,fecRing,ccuAdd,ccuChan,fedKeyStr.c_str(),lldChannel,apvId,stripId);
      
    std::map<string,string> fitParam;
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

    TH1F noiseHist ("noiseHist","",nBin,xMin,xMax);
    for(int iBin = 0; iBin < noiseDistribution->size(); iBin++){
      noiseHist.SetBinContent(iBin+1,noiseDistribution->at(iBin));
      noiseHist.SetBinError(iBin+1,noiseDistributionError->at(iBin));
    }

    TF1  noiseFit  ("noiseFist","gaus(0)",xMin,xMax);
    noiseFit.SetParameters(fitGausNormalization,fitGausMean,fitGausSigma);
    noiseFit.SetParError(0,fitGausNormalizationError);
    noiseFit.SetParError(1,fitGausMeanError);
    noiseFit.SetParError(2,fitGausSigmaError);
        
    if(kSProbab < quantile2sigma){
      badKsTest->cd();
      nbadKsTest++;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);      
    }
    
    if(jBProbab < quantile5sigma){
      badjBTest->cd();
      nbadjBTest++;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);      
    }

    if(fitChi2Probab < quantile4sigma){
      badChi2Test->cd();
      nbadChi2Test++;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);      
    }

    if(aDProbab < quantile3sigma){
      badaDTest->cd();
      nbadaDTest++;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);
    }
    if(jBProbab < quantile5sigma and kSProbab > quantile2sigma and kSProbab < quantile){
      badJBNotKSTest->cd();
      nbadJBNotKSTest++;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);      
    }

    if(aDProbab < quantile3sigma and jBProbab > quantile5sigma and kSProbab > quantile2sigma and kSProbab < quantile){
      badaDNotKSandjBTest->cd();
      nbadaDNotKSandjBTest++;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);      
    }


    if(fitChi2Probab < quantile4sigma and jBProbab > quantile5sigma and kSProbab > quantile2sigma and kSProbab < quantile and aDProbab > quantile3sigma){
      badChi2NotKSandjBandaDTest->cd();
      nbadChi2NotKSandjBandaDTest++;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);      
    }

    if(kSProbab < quantile2sigma or (kSProbab > quantile2sigma and kSProbab < quantile and jBProbab < quantile5sigma) or (kSProbab > quantile2sigma and kSProbab < quantile and jBProbab > quantile5sigma and aDProbab < quantile3sigma)  or (kSProbab > quantile2sigma and kSProbab < quantile and jBProbab > quantile5sigma and aDProbab > quantile3sigma and fitChi2Probab < quantile4sigma)){
      badCombinedTest->cd();
      nbadCombinedTest++;
      moduleNumerator[detid] = moduleNumerator[detid]+1;
      storeOutputCanvas(canvas,noiseHist,noiseFit,name,fitParam);      
    }
  }
  std::cout<<std::endl;
  badKsTest->Close();
  badaDTest->Close();
  badjBTest->Close();
  badChi2Test->Close();
  badCombinedTest->Close();
  badJBNotKSTest->Close();
  badaDNotKSandjBTest->Close();
  badChi2NotKSandjBandaDTest->Close();

  cout<<"#### Bad KS Test Channels "<<nbadKsTest<<" ---> "<<double(nbadKsTest)/tree->GetEntries()<<" % "<<endl;
  cout<<"#### Bad JB Test Channels "<<nbadjBTest<<" ---> "<<double(nbadjBTest)/tree->GetEntries()<<" % "<<endl;
  cout<<"#### Bad AD Test Channels "<<nbadaDTest<<" ---> "<<double(nbadaDTest)/tree->GetEntries()<<" % "<<endl;
  cout<<"#### Bad Chi2 Test Channels "<<nbadChi2Test<<" ---> "<<double(nbadChi2Test)/tree->GetEntries()<<" % "<<endl;
  cout<<"#### Bad JB but not KS Test Channels "<<nbadJBNotKSTest<<" ---> "<<double(nbadJBNotKSTest)/tree->GetEntries()<<" % "<<endl;
  cout<<"#### Bad AD but not KS and not JB Test Channels "<<nbadaDNotKSandjBTest<<" ---> "<<double(nbadaDNotKSandjBTest)/tree->GetEntries()<<" % "<<endl;
  cout<<"#### Bad Chi2 but not KS and JB and AD Test Channels "<<nbadChi2NotKSandjBandaDTest<<" ---> "<<double(nbadChi2NotKSandjBandaDTest)/tree->GetEntries()<<" % "<<endl;
  cout<<"#### Bad Combined Test Channels "<<nbadCombinedTest<<" ---> "<<double(nbadCombinedTest)/tree->GetEntries()<<" % "<<endl;

  ofstream channelMap ((outputDIR+"/fractionOfGoodChannels.txt").c_str());
  for(auto module : moduleDenominator)
    channelMap << module.first <<"  "<< 1. - double(moduleNumerator[module.first])/double(moduleDenominator[module.first]) << "\n";
  channelMap.close();

  ofstream nchannelMap ((outputDIR+"/numberBadChannels.txt").c_str());
  for(auto module : moduleNumerator)
    nchannelMap << module.first <<"  "<< moduleNumerator[module.first] << "\n";
  nchannelMap.close();

}
