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
#include <TChain.h>
#include "delayUtils.h"
 
// parametrize profile with a gaussian shape
static bool isGaussian = true;
// reduce the number of events by
static int  reductionFactor = 1;
// min number of filled bins
static int minFilledBin = 11;
// max allowed delay
static float peakBoundary = 10.8;
// min amplitude
static float amplitudeMin = 40.;
// min signficance delay/sigma_delay
static float significance = 2;
// min delay to apply significance cut
static float minDelayForSignificance = 3;
// max sigma allowed
static float maxSigma = 20;

/// function that runs on the evnet and produce profiles for layers
void ChannnelPlots(const std::vector<std::shared_ptr<TTree> > & tree, 
		   const std::vector<std::shared_ptr<TTree> > & map,
		   const std::shared_ptr<TTree>  & corrections,
		   std::map<uint32_t,std::shared_ptr<TProfile> > & channelMap,
		   std::map<uint32_t,int > & fitStatus,
		   const string & observable) {

  std::cout<<"###############################################################"<<std::endl;
  std::cout << "Preparing Layer plot for the for all the different channels "<< std::endl; 
  std::cout<<"##############################################################"<<std::endl;

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

      // make the profile if not
      if(channelMap[detid].get() == 0 or channelMap[detid].get() == NULL)
	channelMap[detid] = std::shared_ptr<TProfile> (new TProfile(Form("detid_%d",detid),"",delayBins.size()-1,&delayBins[0],yMin,yMax));

      float value = 0;
      if(observable == "maxCharge")
	value = obs*(clCorrectedSignalOverNoise)/(clSignalOverNoise);
      else 
	value = obs;

      channelMap[detid]->Fill(delay-correction,value);      
    }
    std::cout<<std::endl;
  }

  // correct profiles and fit them with a gaussian
  std::cout<<"Analyze each single channel "<<std::endl;
  long int iFit = 0;
  long int badFits = 0;
  long int noFitResult = 0;
  for(auto iprof : channelMap){
    if(observable == "maxCharge")
      correctProfile(iprof.second);
    cout.flush();
    if(iFit % 100 == 0) cout<<"\r"<<"iFit "<<100*double(iFit)/(channelMap.size())<<" % ";
    if(iprof.second->Integral() == 0) {
      noFitResult++;
      continue;
    }
    TFitResultPtr result = fitProfile(iprof.second,isGaussian,"Q");
    fitStatus[iprof.first] = result->Status()+result->CovMatrixStatus() ;
    if(result->CovMatrixStatus() != 3 or result->Status() != 0){
      badFits++;
    }
    iFit++;
  }
  
  std::cout<<std::endl;
  std::cout<<"Loop on events terminated"<<std::endl;  
  std::cout<<"NFits performed "<<iFit<<" bad ones : "<<100*badFits/iFit<<" % "<<" noFit result "<<100*noFitResult/iFit<<std::endl;

  std::cout<<"###################################"<<std::endl;
  std::cout<<"##### End of Channel Analysis #####"<<std::endl;
  std::cout<<"###################################"<<std::endl;

}

// by dedault some fit results are store in a root file, a text dump: detid fitted delay is produced, as well as deid delay uncertainty
void delayValidationPerModule(string inputDIR,    // take the input directory where all merged files are located
			      string file1        = "nocorrection.root", // no correction file
			      string postfix      = "merged", // sub string to be used to find merged files
			      string observable   = "maxCharge",
			      string outputDIR    = "prompt",
			      bool   saveCanvas   =  false,
			      bool   saveCorrectionTree = true){

  // prepare style and load macros
  setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
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
  std::vector< shared_ptr<TTree> > clusters;
  std::vector< shared_ptr<TTree> > readoutMap;
  std::vector< shared_ptr<TFile> > files;
  for(auto ifile : fileList){
    files.push_back(shared_ptr<TFile>(TFile::Open(ifile.c_str())));
    clusters.push_back(shared_ptr<TTree>((TTree*) files.back()->FindObjectAny("clusters")));
    readoutMap.push_back(shared_ptr<TTree>((TTree*) files.back()->FindObjectAny("readoutMap")));
  }

  std::shared_ptr<TFile> _file1 (TFile::Open(file1.c_str()));
  std::shared_ptr<TTree> delayCorrections ((TTree*)_file1->FindObjectAny("delayCorrections"));

  //map with key the detId number, one profile associated to it
  std::map<uint32_t,std::shared_ptr<TProfile> > channelMap;
  std::map<uint32_t,int > fitStatus;
  ChannnelPlots(clusters,readoutMap,delayCorrections,channelMap,fitStatus,observable);  
  // dumpt in a text file to be displayed on the tracker map and uncertainty

  // produce outputs
  cout<<"### Dump peak in a text file"<<endl;
  ofstream dumpPeak((outputDIR+"/dumpBestDelay_"+observable+".txt").c_str());
  int nullPointers = 0;
  for(auto imap : channelMap){    
    if(imap.second->GetFunction(Form("Gaus_%s",imap.second->GetName())))
      dumpPeak<<imap.first<<"   "<<imap.second->GetFunction(Form("Gaus_%s",imap.second->GetName()))->GetParameter(1)<<"\n";
    else
      nullPointers++;
  }
  dumpPeak.close();

  std::cout<<"### Null pointers "<<100*float(nullPointers)/channelMap.size()<<" %s "<<endl;

  cout<<"### Dump peak error in a text file"<<endl;
  ofstream dumpMean((outputDIR+"/dumpBestDelayError_"+observable+".txt").c_str());
  for(auto imap : channelMap){
    if(imap.second->GetFunction(Form("Gaus_%s",imap.second->GetName())))
      dumpMean<<imap.first<<"   "<<imap.second->GetFunction(Form("Gaus_%s",imap.second->GetName()))->GetParError(1)<<"\n";
  }  
  dumpMean.close();
  
  if(saveCanvas)
    saveAll(channelMap,outputDIR,observable); // save all the plots into a single root files


  //// outptut tree with analysis result
  std::shared_ptr<TFile> outputFile;
  std::shared_ptr<TTree> outputTree;

  if(saveCorrectionTree){
 
    std::cout<<"### Dumpt output tree for tkCommissioner "<<std::endl;
    outputFile = std::shared_ptr<TFile>(new TFile((outputDIR+"/tree_delayCorrection.root").c_str(),"RECREATE"));

    // branches definition for the output tree
    readoutMap.at(0)->SetBranchStatus("*",kTRUE);
    readoutMap.at(0)->SetBranchStatus("moduleName",kFALSE);
    readoutMap.at(0)->SetBranchStatus("moduleId",kFALSE);
    readoutMap.at(0)->SetBranchStatus("delay",kFALSE);
    outputTree = std::shared_ptr<TTree>(readoutMap.at(0)->CopyTree(""));
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
    TBranch* bdelayCorr        = outputTree->Branch("delayCorr",&delayCorr,"delayCorr/F");

    TBranch* bnotFound             = outputTree->Branch("notFound",&notFound,"notFound/I");
    TBranch* bsignificanceRejected = outputTree->Branch("significanceRejected",&significanceRejected,"significanceRejected/I");
    TBranch* bamplitudeRejected    = outputTree->Branch("amplitudeRejected",&amplitudeRejected,"amplitudeRejected/I");
    TBranch* bsigmaRejected        = outputTree->Branch("sigmaRejected",&sigmaRejected,"sigmaRejected/I");
    TBranch* bnFilledBinRejected   = outputTree->Branch("nFilledBinRejected",&nFilledBinRejected,"nFilledBinRejected/I");
    TBranch* bpeakOutBoundaryRejected   = outputTree->Branch("peakOutBoundaryRejected",&peakOutBoundaryRejected,"peakOutBoundaryRejected/I");
    TBranch* bfitRejected          = outputTree->Branch("fitRejected",&fitRejected,"fitRejected/I");

    uint32_t detid;
    outputTree->SetBranchAddress("Detid",&detid);
    long int notFoundChannels = 0;
    
    cout<<"### Start the loop on the channel map "<<endl;    
    for(int iChannel = 0; iChannel < outputTree->GetEntries(); iChannel++){
      cout.flush();
      if(iChannel % 100 == 0) cout<<"\r"<<"iChannel "<<100*double(iChannel)/(outputTree->GetEntries())<<" % ";
      outputTree->GetEntry(iChannel);

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
      
      // not found detid --> i.e. no data there --> can be masked
      if(channelMap[detid] == 0 or channelMap[detid] == NULL or channelMap[detid]->GetFunction(Form("Gaus_%s",channelMap[detid]->GetName())) == 0 or channelMap[detid]->GetFunction(Form("Gaus_%s",channelMap[detid]->GetName())) == NULL){
	notFound = 1;       
	notFoundChannels++;	
      }

      else{
	// fill the amplitude branch
	for(int iBin = 0; iBin < channelMap[detid]->GetNbinsX(); iBin++){
	  amplitude.push_back(channelMap[detid]->GetBinContent(iBin));
	}
	
            
	measuredDelay    = channelMap[detid]->GetFunction(Form("Gaus_%s",channelMap[detid]->GetName()))->GetParameter(1);
	measuredDelayUnc = channelMap[detid]->GetFunction(Form("Gaus_%s",channelMap[detid]->GetName()))->GetParError(1);
	measuredMeanAmplitude    = channelMap[detid]->GetFunction(Form("Gaus_%s",channelMap[detid]->GetName()))->GetParameter(0);
	measuredMeanAmplitudeUnc = channelMap[detid]->GetFunction(Form("Gaus_%s",channelMap[detid]->GetName()))->GetParError(0);
	measuredSigma    = channelMap[detid]->GetFunction(Form("Gaus_%s",channelMap[detid]->GetName()))->GetParameter(2);
	measuredSigmaUnc = channelMap[detid]->GetFunction(Form("Gaus_%s",channelMap[detid]->GetName()))->GetParError(2);	

	bmeasuredDelay->Fill();
	bmeasuredDelayUnc->Fill();
	bmeasuredMeanAmplitude->Fill();
	bmeasuredMeanAmplitudeUnc->Fill();
	bmeasuredSigma->Fill();
	bmeasuredSigmaUnc->Fill();
	bamplitude->Fill();
      
	/// apply selections for good fit --> check number of filled bins
	int nFiledBins = getFilledBins(channelMap[detid]);
	if(nFiledBins < minFilledBin){	
	  nFilledBinRejected = 1;

	//store it in the output canvas file
	  if(not outputFile->GetDirectory("binEntries"))
	    outputFile->mkdir("binEntries");
	  outputFile->cd("binEntries");
	  if(gDirectory->GetListOfKeys()->FindObject(Form("detid_%d",detid)) == 0){
	    std::shared_ptr<TCanvas> c1b = prepareCanvas(Form("detid_%d",detid),observable);
	    plotAll(c1b,channelMap[detid]);
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
	    plotAll(c1b,channelMap[detid]);
	    c1b->Write();
	  }
	  outputFile->cd();	
	}
	
	// mean value outside boundary
	if(fabs(channelMap[detid]->GetFunction(Form("Gaus_%s",channelMap[detid]->GetName()))->GetParameter(1)) > peakBoundary){
	
	  peakOutBoundaryRejected = 1;

	  //store it in the output canvas file
	  if(not outputFile->GetDirectory("peakOutRange"))
	    outputFile->mkdir("peakOutRange");
	  outputFile->cd("peakOutRange");
	  if(gDirectory->GetListOfKeys()->FindObject(Form("detid_%d",detid)) == 0){
	    std::shared_ptr<TCanvas> c1b = prepareCanvas(Form("detid_%d",detid),observable);
	    plotAll(c1b,channelMap[detid]);
	    c1b->Write();
	  }
	  outputFile->cd();	
	} 
	
	// amplitude cut
	if(fabs(channelMap[detid]->GetFunction(Form("Gaus_%s",channelMap[detid]->GetName()))->GetParameter(0)) < amplitudeMin){
	  
	  amplitudeRejected = 1;

	  //store it in the output canvas file
	  if(not outputFile->GetDirectory("smallAmplitude"))
	    outputFile->mkdir("smallAmplitude");
	  outputFile->cd("smallAmplitude");
	  if(gDirectory->GetListOfKeys()->FindObject(Form("detid_%d",detid)) == 0){
	    std::shared_ptr<TCanvas> c1b = prepareCanvas(Form("detid_%d",detid),observable);
	    plotAll(c1b,channelMap[detid]);
	    c1b->Write();
	  }
	  outputFile->cd();	
	} 

	// amplitude cut
	if(fabs(channelMap[detid]->GetFunction(Form("Gaus_%s",channelMap[detid]->GetName()))->GetParameter(1))/fabs(channelMap[detid]->GetFunction(Form("Gaus_%s",channelMap[detid]->GetName()))->GetParError(1)) < significance and fabs(channelMap[detid]->GetFunction(Form("Gaus_%s",channelMap[detid]->GetName()))->GetParameter(1)) > minDelayForSignificance){
	
	  significanceRejected = 1;

	  //store it in the output canvas file
	  if(not outputFile->GetDirectory("nonSignificantLargeShift"))
	    outputFile->mkdir("nonSignificantLargeShift");
	  outputFile->cd("nonSignificantLargeShift");
	  if(gDirectory->GetListOfKeys()->FindObject(Form("detid_%d",detid)) == 0){
	    std::shared_ptr<TCanvas> c1b = prepareCanvas(Form("detid_%d",detid),observable);
	    plotAll(c1b,channelMap[detid]);
	    c1b->Write();
	  }
	  outputFile->cd();	
	} 

	// sigma cut
	if(channelMap[detid]->GetFunction(Form("Gaus_%s",channelMap[detid]->GetName()))->GetParameter(2) > maxSigma){
	  
	  sigmaRejected = 1;
	  //store it in the output canvas file                                                                                                                               
	  if(not outputFile->GetDirectory("largeGaussSigma"))
	    outputFile->mkdir("largeGaussSigma");
	  outputFile->cd("largeGaussSigma");
	  if(gDirectory->GetListOfKeys()->FindObject(Form("detid_%d",detid)) == 0){
	    std::shared_ptr<TCanvas> c1b = prepareCanvas(Form("detid_%d",detid),observable);
	    plotAll(c1b,channelMap[detid]);
	    c1b->Write();
	  }
	  outputFile->cd();
	}

	if(not notFound and not significanceRejected and not amplitudeRejected and not sigmaRejected and not nFilledBinRejected and not fitRejected and not peakOutBoundaryRejected)
	  
	  delayCorr = channelMap[detid]->GetFunction(Form("Gaus_%s",channelMap[detid]->GetName()))->GetParameter(1);
      }

      bnotFound->Fill();
      bsignificanceRejected->Fill();
      bamplitudeRejected->Fill();
      bsigmaRejected->Fill();
      bnFilledBinRejected->Fill();
      bfitRejected->Fill();      
      bdelayCorr->Fill();
      bpeakOutBoundaryRejected->Fill();
    }    
    cout<<endl;
    cout<<"### Write output tree for tkCommissioner: not found channels i.e. no clusters "<<100*float(notFoundChannels)/outputTree->GetEntries()<<" % "<<endl;    
    outputFile->cd();
    outputTree->BuildIndex("Detid");
    outputTree->Write(outputTree->GetName(),TObject::kOverwrite);
  }
  cout<<"### Clear and Close "<<endl;  
  clusters.clear();
  readoutMap.clear();
  files.clear();
}

