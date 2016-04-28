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
static int  reductionFactor = 100;

/// function that runs on the evnet and produce profiles for layers
void ChannnelPlots(const std::vector<std::shared_ptr<TTree> > & tree, 
		   const std::vector<std::shared_ptr<TTree> > & map,
		   const std::shared_ptr<TTree>  & corrections,
		   std::map<uint32_t,std::shared_ptr<TProfile> > & channelMap,
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


void delayValidationPerModule(string inputDIR,    // take the input directory where all merged files are located
			      string file1        = "nocorrection.root", // no correction file
			      string postfix      = "merged", // sub string to be used to find merged files
			      string observable   = "maxCharge",
			      string outputDIR    = "prompt"){

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
  ChannnelPlots(clusters,readoutMap,delayCorrections,channelMap,observable);  
  // dumpt in a text file to be displayed on the tracker map and uncertainty
  cout<<"#### Dump peak in a text file"<<endl;
  ofstream dumpPeak((outputDIR+"/dumpBestDelay_"+observable+".txt").c_str());
  int nullPointers = 0;
  for(auto imap : channelMap){    
    if(imap.second->GetFunction(Form("Gaus_%s",imap.second->GetName())))
      dumpPeak<<imap.first<<"   "<<imap.second->GetFunction(Form("Gaus_%s",imap.second->GetName()))->GetMaximumX(-10,10)<<"\n";
    else
      nullPointers++;
  }
  dumpPeak.close();

  cout<<"#### Dump peak error in a text file"<<endl;
  ofstream dumpMean((outputDIR+"/dumpBestDelayError_"+observable+".txt").c_str());
  for(auto imap : channelMap){
    if(imap.second->GetFunction(Form("Gaus_%s",imap.second->GetName())))
      dumpMean<<imap.first<<"   "<<imap.second->GetFunction(Form("Gaus_%s",imap.second->GetName()))->GetParError(1)<<"\n";
  }  
  dumpMean.close();

  std::cout<<"### Null pointers "<<100*float(nullPointers)/channelMap.size()<<" %s "<<endl;

  saveAll(channelMap,outputDIR,observable); // save all the plots into a single root files
 

  clusters.clear();
  readoutMap.clear();
  files.clear();


}

