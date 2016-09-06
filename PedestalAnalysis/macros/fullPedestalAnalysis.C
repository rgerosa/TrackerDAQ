#include "../../../DataFormats/SiStripCommon/interface/SiStripFecKey.h"
#include "../../../DataFormats/SiStripCommon/interface/SiStripFedKey.h"

#include "CMS_lumi.h"

// reduce the number of channels by                                                                                                                                                                 
static int  reductionFactor = 1;
static bool generateRandomDistribution = false;

void fullPedestalAnalysis(string inputDIR, string outputDIR, string inputCablingMap, string outputFileName){

  gROOT->ProcessLine("gErrorIgnoreLevel = 1");
  
  // open the file and prepare the cluster tree, adding the other trees as frined --> memory consuming                                                                                                
  std::cout<<"##################################"<<std::endl;
  std::cout<<"###### fullPedestalAnalysis ######"<<std::endl;
  std::cout<<"##################################"<<std::endl;

  clock_t tStart = clock();

  // prepare style and load macros                                                                                                                                                                    
  setTDRStyle();
  gROOT->SetBatch(kTRUE);

  system(("mkdir -p "+outputDIR).c_str());
  ifstream file;

  std::cout<<"### Make input file list"<<std::endl;
  system(("find "+inputDIR+" -name \"*.root\" > file.temp").c_str());
  std::ifstream infile;
  string line;
  vector<string> fileList;
  infile.open("file.temp",ifstream::in);
  if(infile.is_open()){
    while(!infile.eof()){
      getline(infile,line);
      if(line != "" and TString(line).Contains(".root") and line !="\n"){
        fileList.push_back(line);
      }
    }
  }
  system("rm file.temp");
  std::sort(fileList.begin(),fileList.end());

  TFile* cablingFile = TFile::Open(inputCablingMap.c_str(),"READ");
  cablingFile->cd();
  TTree* readoutMap = (TTree*) cablingFile->FindObjectAny("readoutMap");
  TTreeReader reader(readoutMap);
  TTreeReaderValue<uint32_t> detid    (reader,"detid");
  TTreeReaderValue<uint16_t> fecCrate (reader,"fecCrate");
  TTreeReaderValue<uint16_t> fecSlot  (reader,"fecSlot");
  TTreeReaderValue<uint16_t> fecRing  (reader,"fecRing");
  TTreeReaderValue<uint16_t> ccuAdd   (reader,"ccuAdd");
  TTreeReaderValue<uint16_t> ccuChan  (reader,"ccuChan");
  TTreeReaderValue<uint16_t> lldChannel  (reader,"lldChannel");
  TTreeReaderValue<uint16_t> fedId  (reader,"fedId");
  TTreeReaderValue<uint16_t> fedCh  (reader,"fedCh");

  // output tree
  TFile* ouputTreeFile = new TFile((outputDIR+"/"+outputFileName).c_str(),"RECREATE");
  ouputTreeFile->cd();
  ouputTreeFile->SetCompressionLevel(0);
  TTree* outputTree = new TTree("pedestalFullNoise","pedestalFullNoise");
  
  // branches
  uint32_t detid_,fedKey_;
  uint16_t fecCrate_,fecSlot_, fecRing_, ccuAdd_, ccuChan_, lldChannel_, fedId_, fedCh_, apvId_, stripId_;
  float    noiseMean_,noiseRMS_, noiseSkewness_, noiseKurtosis_;
  float    fitChi2_, fitChi2Probab_, fitStatus_;
  float    fitGausMean_, fitGausSigma_, fitGausNormalization_;
  float    fitGausMeanError_, fitGausSigmaError_, fitGausNormalizationError_;
  float    noiseIntegral3Sigma_, noiseIntegral3SigmaFromFit_;
  float    noiseIntegral4Sigma_, noiseIntegral4SigmaFromFit_;
  float    noiseIntegral5Sigma_, noiseIntegral5SigmaFromFit_;
  float    kSValue_, kSProbab_, jBValue_, jBProbab_, aDValue_, aDProbab_;
  vector<float> noiseDistribution_, noiseDistributionError_;
  float xMin_, xMax_, nBin_ ;

  outputTree->Branch("detid",&detid_,"detid/i");
  outputTree->Branch("fedKey",&fedKey_,"fedKey/i");
  outputTree->Branch("fecCrate",&fecCrate_,"fecCrate/s");
  outputTree->Branch("fecSlot",&fecSlot_,"fecSlot/s");
  outputTree->Branch("fecRing",&fecRing_,"fecRing/s");
  outputTree->Branch("ccuAdd",&ccuAdd_,"ccuAdd/s");
  outputTree->Branch("ccuChan",&ccuChan_,"ccuChan/s");
  outputTree->Branch("lldChannel",&lldChannel_,"lldChannel/s");
  outputTree->Branch("fedId",&fedId_,"fedId/s");
  outputTree->Branch("fedCh",&fedCh_,"fedCh/s");
  outputTree->Branch("apvId",&apvId_,"apvId/s");
  outputTree->Branch("stripId",&stripId_,"stripId/s");

  outputTree->Branch("noiseMean",&noiseMean_,"noiseMean/F");
  outputTree->Branch("noiseRMS",&noiseRMS_,"noiseRMS/F");
  outputTree->Branch("noiseSkewness",&noiseSkewness_,"noiseSkewness/F");
  outputTree->Branch("noiseKurtosis",&noiseKurtosis_,"noiseKurtosis/F");
  outputTree->Branch("fitGausNormalization",&fitGausNormalization_,"fitGausNormalization/F");
  outputTree->Branch("fitGausMean",&fitGausMean_,"fitGausMean/F");
  outputTree->Branch("fitGausSigma",&fitGausSigma_,"fitGausSigma/F");
  outputTree->Branch("fitGausNormalizationError",&fitGausNormalizationError_,"fitGausNormalizationError/F");
  outputTree->Branch("fitGausMeanError",&fitGausMeanError_,"fitGausMeanError/F");
  outputTree->Branch("fitGausSigmaError",&fitGausSigmaError_,"fitGausSigmaError/F");
  outputTree->Branch("fitChi2",&fitChi2_,"fitChi2/F");
  outputTree->Branch("fitChi2Probab",&fitChi2Probab_,"fitChi2Probab/F");
  outputTree->Branch("fitStatus",&fitStatus_,"fitStatus_F");
  outputTree->Branch("noiseIntegral3Sigma",&noiseIntegral3Sigma_,"noiseIntegral3Sigma/F");
  outputTree->Branch("noiseIntegral3SigmaFromFit",&noiseIntegral3SigmaFromFit_,"noiseIntegral3SigmaFromFit/F");
  outputTree->Branch("noiseIntegral4Sigma",&noiseIntegral4Sigma_,"noiseIntegral4Sigma/F");
  outputTree->Branch("noiseIntegral4SigmaFromFit",&noiseIntegral4SigmaFromFit_,"noiseIntegral4SigmaFromFit/F");
  outputTree->Branch("noiseIntegral5Sigma",&noiseIntegral4Sigma_,"noiseIntegral5Sigma/F");
  outputTree->Branch("noiseIntegral5SigmaFromFit",&noiseIntegral4SigmaFromFit_,"noiseIntegral5SigmaFromFit/F");
  outputTree->Branch("kSValue",&kSValue_,"kSValue/F");
  outputTree->Branch("jBValue",&jBValue_,"jBValue/F");
  outputTree->Branch("aDValue",&aDValue_,"aDValue/F");
  outputTree->Branch("kSProbab",&kSProbab_,"kSProbab/F");
  outputTree->Branch("jBProbab",&jBProbab_,"jBProbab/F");
  outputTree->Branch("aDProbab",&aDProbab_,"aDProbab/F");
  outputTree->Branch("xMin",&xMin_,"xMin/F");
  outputTree->Branch("xMax",&xMax_,"xMax/F");
  outputTree->Branch("nBin",&nBin_,"nBin/F");

  bool histoBranches = false;

  // Loop on the file list to extract each histogram 2D DQM histo with full noise distribution  
  TH1F* histoNoiseStrip = NULL;
  TF1*  fitFunc = NULL;
  TH1F* randomHisto = NULL;
  TFitResultPtr result;
  for(auto file : fileList){
    cout<<"input file: "<<file<<endl;
    TFile* inputFile = TFile::Open(file.c_str(),"READ");
    inputFile->cd();
    // take into account that the DQM file structure for strips is always the same --> use cabling map to browse the histograms
    reader.SetEntry(0);
    TH2* histoNoise = NULL;
    long int iChannel = 0;
    int noFitResult = 0;
    while(reader.Next()){
      cout.flush();
      if(iChannel %10 == 0) cout<<"\r"<<"iChannel "<<100*double(iChannel)/(readoutMap->GetEntries()/reductionFactor)<<" % ";
      if(iChannel > double(readoutMap->GetEntries())/reductionFactor) break;
      iChannel++;
      TString objName;
      uint32_t fedKey =  SiStripFedKey(*fedId,SiStripFedKey::feUnit(*fedCh),SiStripFedKey::feChan(*fedCh)).key();
      std::stringstream stream;
      stream << std::hex << fedKey;
      string fedKeyStr = stream.str();
      if(fedKeyStr.size() == 4)
	objName = Form("DQMData/SiStrip/ControlView/FecCrate%d/FecSlot%d/FecRing%d/CcuAddr%d/CcuChan%d/ExpertHisto_PedsFullNoise_FedKey0x0000%s_LldChannel%d_Noise2D",*fecCrate,*fecSlot,*fecRing,*ccuAdd,*ccuChan,fedKeyStr.c_str(),*lldChannel);      
      else if(fedKeyStr.size() == 5)
	objName = Form("DQMData/SiStrip/ControlView/FecCrate%d/FecSlot%d/FecRing%d/CcuAddr%d/CcuChan%d/ExpertHisto_PedsFullNoise_FedKey0x000%s_LldChannel%d_Noise2D",*fecCrate,*fecSlot,*fecRing,*ccuAdd,*ccuChan,fedKeyStr.c_str(),*lldChannel);      
      else
	cerr<<"hex number to short "<<fedKeyStr<<" --> please check "<<endl;

      inputFile->GetObject(objName.Data(),histoNoise);
      // extract single strip noise histogram --> loop on the y-axis
      uint16_t apvID = 0;
      uint16_t stripID = 0;       
      if(histoNoiseStrip == 0 or histoNoiseStrip == NULL){
	histoNoiseStrip = new TH1F ("histoNoiseStrip","",histoNoise->GetNbinsX(),histoNoise->GetXaxis()->GetXmin(),histoNoise->GetXaxis()->GetXmax());
	histoNoiseStrip->Sumw2();
      }
      for(int iBinY = 0; iBinY < histoNoise->GetNbinsY(); iBinY++){
	histoNoiseStrip->Reset();
	histoNoiseStrip->SetDirectory(0);
	// two multiplexed APV per line
	if(iBinY < histoNoise->GetNbinsY()/2) apvID = 1;
	else apvID = 2;
	// strip id
	stripID++;
	if(stripID > 128) stripID = 1;
	// loop on x-axis bin
	for(int iBinX = 0; iBinX < histoNoise->GetNbinsX(); iBinX++){
	  histoNoiseStrip->SetBinContent(iBinX+1,histoNoise->GetBinContent(iBinX+1,iBinY+1));
	  histoNoiseStrip->SetBinError(iBinX+1,histoNoise->GetBinError(iBinX+1,iBinY+1));	    
	}
     	
	// to initialize branches
	detid_ = 0; fedKey_ = 0; fecCrate_ = 0; fecSlot_ = 0; fecRing_ = 0; ccuAdd_ = 0; ccuChan_ = 0; lldChannel_ = 0; fedId_ = 0; fedCh_ = 0; apvId_ = 0; stripId_ = 0; 
	noiseMean_ = 0.; noiseRMS_ =  0.; noiseSkewness_ = 0.; noiseKurtosis_ = 0.; 
	fitGausMean_ = 0.; fitGausSigma_ = 0.;fitGausNormalization_ = 0.;
	fitGausMeanError_ = 0.; fitGausSigmaError_ = 0.;fitGausNormalizationError_ = 0.;	  	  
	fitChi2_ = 0.; fitChi2Probab_ = 0.; fitStatus_ = -1.; 
	noiseIntegral3Sigma_ = 0.; noiseIntegral3SigmaFromFit_ = 0.; 
	noiseIntegral4Sigma_ = 0.; noiseIntegral4SigmaFromFit_ = 0.; 
	noiseIntegral5Sigma_ = 0.; noiseIntegral5SigmaFromFit_ = 0.; 
	kSProbab_ = 0.; jBProbab_ = 0.;
	kSValue_ = 0.; jBValue_ = 0.; 
	aDValue_= 0.; aDProbab_ = 0.;
	nBin_ = 0.; xMin_ = 0.; xMax_ = 0.;
	
	// basic info
	detid_ = *detid;
	fedKey_ = fedKey;
	fecCrate_ = *fecCrate;
	fecSlot_ = *fecSlot;
	fecRing_ = *fecRing;
	ccuAdd_  = *ccuAdd;
	ccuChan_ = *ccuChan;
	lldChannel_ = *lldChannel;
	fedId_   = *fedId;
	fedCh_   = *fedCh;
	apvId_   = apvID;
	stripId_ = stripID;
	
	// basic info of nioise distribution
	noiseMean_ = histoNoiseStrip->GetMean();
	noiseRMS_  = histoNoiseStrip->GetRMS();
	noiseSkewness_ = histoNoiseStrip->GetSkewness();
	noiseKurtosis_ = histoNoiseStrip->GetKurtosis();
	float integral = histoNoiseStrip->Integral();	
	noiseIntegral3Sigma_ = (histoNoiseStrip->Integral(histoNoiseStrip->FindBin(noiseMean_+noiseRMS_*3),histoNoiseStrip->GetNbinsX()+1) + histoNoiseStrip->Integral(0,histoNoiseStrip->FindBin(noiseMean_-noiseRMS_*3)))/integral;
	noiseIntegral4Sigma_ = (histoNoiseStrip->Integral(histoNoiseStrip->FindBin(noiseMean_+noiseRMS_*4),histoNoiseStrip->GetNbinsX()+1) + histoNoiseStrip->Integral(0,histoNoiseStrip->FindBin(noiseMean_-noiseRMS_*4)))/integral;
	noiseIntegral5Sigma_ = (histoNoiseStrip->Integral(histoNoiseStrip->FindBin(noiseMean_+noiseRMS_*5),histoNoiseStrip->GetNbinsX()+1) + histoNoiseStrip->Integral(0,histoNoiseStrip->FindBin(noiseMean_-noiseRMS_*5)))/integral;
	
	// make a gaussian fit	  	
	if(fitFunc == NULL or fitFunc == 0){
	  fitFunc = new TF1 ("fitFunc","gaus(0)",histoNoise->GetXaxis()->GetXmin(),histoNoise->GetXaxis()->GetXmax());
	}
	fitFunc->SetRange(histoNoise->GetXaxis()->GetXmin(),histoNoise->GetXaxis()->GetXmax());
	fitFunc->SetParameters(histoNoiseStrip->Integral(),noiseMean_,noiseRMS_);
	result = histoNoiseStrip->Fit(fitFunc,"QSR");

	if(result.Get()){
	    fitStatus_     = result->Status();
	    fitGausNormalization_  = fitFunc->GetParameter(0);
	    fitGausMean_   = fitFunc->GetParameter(1);
	    fitGausSigma_  = fitFunc->GetParameter(2);
	    fitGausNormalizationError_  = fitFunc->GetParError(0);
	    fitGausMeanError_  = fitFunc->GetParError(1);
	    fitGausSigmaError_ = fitFunc->GetParError(2);
	    fitChi2_           = result->Chi2();
	    fitChi2Probab_     = result->Prob();

	    noiseIntegral3SigmaFromFit_ = (histoNoiseStrip->Integral(histoNoiseStrip->FindBin(noiseMean_+fitGausSigma_*3),histoNoiseStrip->GetNbinsX()+1) + histoNoiseStrip->Integral(0,histoNoiseStrip->FindBin(noiseMean_-fitGausSigma_*3)))/histoNoiseStrip->Integral();
	    noiseIntegral4SigmaFromFit_ = (histoNoiseStrip->Integral(histoNoiseStrip->FindBin(noiseMean_+fitGausSigma_*4),histoNoiseStrip->GetNbinsX()+1) + histoNoiseStrip->Integral(0,histoNoiseStrip->FindBin(noiseMean_-fitGausSigma_*4)))/histoNoiseStrip->Integral();
	    noiseIntegral5SigmaFromFit_ = (histoNoiseStrip->Integral(histoNoiseStrip->FindBin(noiseMean_+fitGausSigma_*5),histoNoiseStrip->GetNbinsX()+1) + histoNoiseStrip->Integral(0,histoNoiseStrip->FindBin(noiseMean_-fitGausSigma_*5)))/histoNoiseStrip->Integral();
	    
	    jBValue_   = (histoNoiseStrip->Integral()/6)*(noiseSkewness_*noiseSkewness_+(noiseKurtosis_*noiseKurtosis_)/4);	  
	    jBProbab_  = ROOT::Math::chisquared_cdf_c(jBValue_,2);

	    if(randomHisto == 0 or randomHisto == NULL)
	      randomHisto = (TH1F*) histoNoiseStrip->Clone("randomHisto");	    	    
	    randomHisto->Reset();
	    randomHisto->SetDirectory(0);     
	
      
	    if(integral != 0){	      
	      if(generateRandomDistribution){
		randomHisto->FillRandom("fitFunc",histoNoiseStrip->Integral());	    
		kSValue_  = histoNoiseStrip->KolmogorovTest(randomHisto,"MN");
		kSProbab_ = histoNoiseStrip->KolmogorovTest(randomHisto,"N");	    
		aDValue_  = histoNoiseStrip->AndersonDarlingTest(randomHisto,"T");
		aDProbab_ = histoNoiseStrip->AndersonDarlingTest(randomHisto);
	      }
	      else{
		
		randomHisto->Add(fitFunc);		
		kSValue_  = histoNoiseStrip->KolmogorovTest(randomHisto,"MN"); 
		kSProbab_ = histoNoiseStrip->KolmogorovTest(randomHisto,"N");
		// AD test
		ROOT::Fit::BinData data1;
		ROOT::Fit::BinData data2;
		ROOT::Fit::FillData(data1,histoNoiseStrip,0);
		data2.Initialize(randomHisto->GetNbinsX()+1,1);
		for(int ibin = 0; ibin < randomHisto->GetNbinsX(); ibin++){ 
		  if(histoNoiseStrip->GetBinContent(ibin+1) != 0 or randomHisto->GetBinContent(ibin+1) >= 1)
		    data2.Add(randomHisto->GetBinCenter(ibin+1),randomHisto->GetBinContent(ibin+1),randomHisto->GetBinError(ibin+1));
		}
	  
		double probab;
		double value;
		ROOT::Math::GoFTest::AndersonDarling2SamplesTest(data1,data2,probab,value);
		aDValue_ = value;
		aDProbab_ = probab;
	      }
	    }
	}
	else
	  noFitResult++;
	
	if(not histoBranches){
	  noiseDistribution_.clear();
	  noiseDistributionError_.clear();
	  outputTree->Branch("noiseDistribution","vector<float>",&noiseDistribution_);
	  outputTree->Branch("noiseDistributionError","vector<float>",&noiseDistributionError_);
	  histoBranches = true;
	}
    
	// set histogram
	noiseDistribution_.clear();
	noiseDistributionError_.clear();
	for(int iBin = 0; iBin < histoNoiseStrip->GetNbinsX(); iBin++){
	  noiseDistribution_.push_back(histoNoiseStrip->GetBinContent(iBin+1));
	  noiseDistributionError_.push_back(histoNoiseStrip->GetBinError(iBin+1));	      
	}
    
	nBin_ = histoNoiseStrip->GetNbinsX();
	xMin_ = histoNoise->GetXaxis()->GetBinLowEdge(1);
	xMax_ = histoNoise->GetXaxis()->GetBinLowEdge(histoNoise->GetNbinsX()+1);

	// fill all branches for each strip
	ouputTreeFile->cd();
	outputTree->Fill();
      }
    }
    inputFile->Close();
    std::cout<<std::endl;
    cout<<"No fit results found for "<<100*double(noFitResult)/iChannel<<endl;
  }
  outputTree->BuildIndex("detid");
  outputTree->Write(outputTree->GetName(),TObject::kOverwrite);
  ouputTreeFile->Close();
  cablingFile->Close();

  /* Do your stuff here */
  cout<<"Time taken: "<<(double)(clock() - tStart)/CLOCKS_PER_SEC<<endl;  
}
