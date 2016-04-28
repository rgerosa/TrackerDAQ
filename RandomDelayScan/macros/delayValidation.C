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

#include "delayUtils.h"
 
// parametrize profile with a gaussian shape
static bool isGaussian = true;
// reduce the number of events by
static int  reductionFactor = 1;

/// function that runs on the evnet and produce profiles for layers
void LayerPlots(const std::shared_ptr<TTree> & tree, 
		const std::shared_ptr<TTree> & map,
		const std::shared_ptr<TTree> & corrections,
		std::vector<std::shared_ptr<TProfile> > & TIBlayers,
		std::vector<std::shared_ptr<TProfile> > & TOBlayers,
		std::vector<std::shared_ptr<TProfile> > & TIDlayers,
		std::vector<std::shared_ptr<TProfile> > & TECPTlayers,
		std::vector<std::shared_ptr<TProfile> > & TECPtlayers,
		std::vector<std::shared_ptr<TProfile> > & TECMTlayers,
		std::vector<std::shared_ptr<TProfile> > & TECMtlayers,
		const string & observable) {

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

  int   correction;
  corrections->SetBranchStatus("*",kFALSE);
  corrections->SetBranchStatus("detid",kTRUE);
  corrections->SetBranchStatus("correction",kTRUE);
  corrections->SetBranchAddress("correction",&correction);
  
  // create vectors for the different Profiles
  float yMin = 0, yMax = 0, xMin = 0, xMax = 0;
  int   nBinsX = 0, nBinsY = 0;
  setLimitsAndBinning("delay",xMin,xMax,nBinsX);
  setLimitsAndBinning(observable,yMin,yMax,nBinsY);

  for(int ilayer = 0; ilayer < 4; ilayer++){
    std::shared_ptr<TProfile> temp (new TProfile(Form("TIB_layer%d",ilayer+1),"",nBinsX,xMin,xMax,yMin,yMax));
    TIBlayers.push_back(temp);
  }
  for(int ilayer = 0; ilayer < 6; ilayer++){
    std::shared_ptr<TProfile> temp (new TProfile(Form("TOB_layer%d",ilayer+1),"",nBinsX,xMin,xMax,yMin,yMax));
    TOBlayers.push_back(temp);
  }
  for(int ilayer = 0; ilayer < 3; ilayer++){
    std::shared_ptr<TProfile> temp (new TProfile(Form("TID_layer%d",ilayer+1),"",nBinsX,xMin,xMax,yMin,yMax));
    TIDlayers.push_back(temp);
  }
  for(int ilayer = 0; ilayer < 9; ilayer++){
    std::shared_ptr<TProfile> temp (new TProfile(Form("TECP_layer%d_T",ilayer+1),"",nBinsX,xMin,xMax,yMin,yMax));
    TECPTlayers.push_back(temp);
    temp = std::shared_ptr<TProfile>(new TProfile(Form("TECP_layer%d_t",ilayer+1),"",nBinsX,xMin,xMax,yMin,yMax));
    TECPtlayers.push_back(temp);
  }
  for(int ilayer = 0; ilayer < 9; ilayer++){
    std::shared_ptr<TProfile> temp (new TProfile(Form("TECM_layer%d_T",ilayer+1),"",nBinsX,xMin,xMax,yMin,yMax));
    TECMTlayers.push_back(temp);
    temp = std::shared_ptr<TProfile>(new TProfile(Form("TECM_layer%d_t",ilayer+1),"",nBinsX,xMin,xMax,yMin,yMax));
    TECMtlayers.push_back(temp);
  }

  std::cout<<"Tree with nEntries "<<tree->GetEntries()<<std::endl;

  long int iEvent  = 0;
  for( ; iEvent < tree->GetEntries()/reductionFactor; iEvent++){    
    // take the event
    tree->GetEntry(iEvent);
    // take the map delay from the detid
    map->GetEntryWithIndex(detid);
    // take the correction from the detid
    corrections->GetEntryWithIndex(detid);

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

    if(subdetid == 3)      
      TIBlayers.at(barrellayer-1)->Fill(delay-correction,value);
    else if(subdetid == 5)
      TOBlayers.at(barrellayer-1)->Fill(delay-correction,value);   
    else if(subdetid == 4)
      TIDlayers.at(TIDlayer-1)->Fill(delay-correction,value);    
    else if(subdetid == 6  and clglobalZ > 0 and thickness > 400)
      TECPTlayers.at(TECPlayer-1)->Fill(delay-correction,value);
    else if(subdetid == 6  and clglobalZ > 0 and thickness < 400)
      TECPtlayers.at(TECPlayer-1)->Fill(delay-correction,value);
    else if(subdetid == 6  and clglobalZ < 0 and thickness > 400)
      TECMTlayers.at(TECMlayer-1)->Fill(delay-correction,value);
    else if(subdetid == 6  and clglobalZ < 0 and thickness < 400)
      TECMtlayers.at(TECMlayer-1)->Fill(delay-correction,value);    
  }

  std::cout<<std::endl;
  std::cout<<"Loop on events terminated"<<std::endl;
  
  // correct profiles and fit them with a gaussian
  std::cout<<"Analyze TIB profiles"<<std::endl;
  for(auto iprof : TIBlayers){
    if(observable == "maxCharge")
      correctProfile(iprof);
    fitProfile(iprof,isGaussian);
  }
  std::cout<<"Analyze TOB profiles"<<std::endl;
  for(auto iprof : TOBlayers){
    if(observable == "maxCharge")
      correctProfile(iprof);
    fitProfile(iprof,isGaussian);
  }     
  std::cout<<"Analyze TID profiles"<<std::endl;
  for(auto iprof : TIDlayers){
    if(observable == "maxCharge")
      correctProfile(iprof);
    fitProfile(iprof,isGaussian);
  }
  std::cout<<"Analyze TECPT profiles"<<std::endl;
  for(auto iprof : TECPTlayers){
    if(observable == "maxCharge")
      correctProfile(iprof);
    fitProfile(iprof,isGaussian);
  }
  std::cout<<"Analyze TECPt profiles"<<std::endl;
  for(auto iprof : TECPtlayers){
    if(observable == "maxCharge")
      correctProfile(iprof);
    fitProfile(iprof,isGaussian);
  }
  std::cout<<"Analyze TECMT profiles"<<std::endl;
  for(auto iprof : TECMTlayers){
    if(observable == "maxCharge")
      correctProfile(iprof);
    fitProfile(iprof,isGaussian);
  }
  std::cout<<"Analyze TECMt profiles"<<std::endl;
  for(auto iprof : TECMtlayers){
    if(observable == "maxCharge")
      correctProfile(iprof);
    fitProfile(iprof,isGaussian);
  }  

  std::cout<<"#################################"<<std::endl;
  std::cout<<"##### End of Layer Analysis #####"<<std::endl;
  std::cout<<"#################################"<<std::endl;

}

//// Per ring analysis
void RPlots(const std::shared_ptr<TTree> & tree, 
	    const std::shared_ptr<TTree> & map,
	    const std::shared_ptr<TTree> & corrections,
	    std::vector<std::shared_ptr<TProfile> > & TIBrings, 
	    std::vector<std::shared_ptr<TProfile> > & TOBrings,
	    std::vector<std::shared_ptr<TProfile> > & TIDrings,
	    std::vector<std::shared_ptr<TProfile> > & TECPTrings,
	    std::vector<std::shared_ptr<TProfile> > & TECPtrings,
	    std::vector<std::shared_ptr<TProfile> > & TECMTrings,
	    std::vector<std::shared_ptr<TProfile> > & TECMtrings,
	    const std::string & observable){



  std::cout<<"######################################################"<<std::endl;
  std::cout << "Preparing Ring plot for the fifferent subdetectors "<< std::endl; 
  std::cout<<"######################################################"<<std::endl;

  cout<<"tree set branch status "<<endl;
  // TTree Reader appear not to be working with addFriend and EventList
  uint32_t detid;
  float    maxCharge, clCorrectedSignalOverNoise, clSignalOverNoise, clglobalX, clglobalY, clglobalZ, thickness, obs;
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

  int   correction;
  corrections->SetBranchStatus("*",kFALSE);
  corrections->SetBranchStatus("detid",kTRUE);
  corrections->SetBranchStatus("correction",kTRUE);
  corrections->SetBranchAddress("correction",&correction);

  // create vectors for the different Profiles
  float yMin = 0, yMax = 0, xMin = 0, xMax = 0;
  int   nBinsX = 0, nBinsY = 0;
  setLimitsAndBinning("delay",xMin,xMax,nBinsX);
  setLimitsAndBinning(observable,yMin,yMax,nBinsY);
  
  cout<<"create profiles  "<<endl;
  // create vectors
  for(int iring = 0; iring < TIBRing.nDivision; iring++){
    std::shared_ptr<TProfile> temp (new TProfile(Form("TIB_ring%d",iring+1),"",nBinsX,xMin,xMax,yMin,yMax));
    TIBrings.push_back(temp);
  }
  for(int iring = 0; iring < TOBRing.nDivision; iring++){
    std::shared_ptr<TProfile> temp (new TProfile(Form("TOB_ring%d",iring+1),"",nBinsX,xMin,xMax,yMin,yMax));
    TOBrings.push_back(temp);
  }
  for(int iring = 0; iring < TIDRing.nDivision; iring++){
    std::shared_ptr<TProfile> temp (new TProfile(Form("TID_ring%d",iring+1),"",nBinsX,xMin,xMax,yMin,yMax));
    TIDrings.push_back(temp);
  }
  for(int iring = 0; iring < TECRing.nDivision; iring++){
    std::shared_ptr<TProfile> temp (new TProfile(Form("TECP_ring%d_T",iring+1),"",nBinsX,xMin,xMax,yMin,yMax));
    TECPTrings.push_back(temp);
    temp = std::shared_ptr<TProfile>(new TProfile(Form("TECP_ring%d_t",iring+1),"",nBinsX,xMin,xMax,yMin,yMax));
    TECPtrings.push_back(temp);
  }
  for(int iring = 0; iring < TECRing.nDivision; iring++){
    std::shared_ptr<TProfile> temp (new TProfile(Form("TECM_ring%d_T",iring+1),"",nBinsX,xMin,xMax,yMin,yMax));
    TECMTrings.push_back(temp);
    temp = std::shared_ptr<TProfile>(new TProfile(Form("TECM_ring%d_t",iring+1),"",nBinsX,xMin,xMax,yMin,yMax));
    TECMtrings.push_back(temp);
  }

  std::cout<<"Tree with nEntries "<<tree->GetEntries()<<std::endl;

  long int iEvent = 0;
  for( ; iEvent < tree->GetEntries()/reductionFactor; iEvent++){    
    tree->GetEntry(iEvent);
    // take the map delay from the detid
    map->GetEntryWithIndex(detid);
    // take the correction from the detid
    corrections->GetEntryWithIndex(detid);

    cout.flush();
    if(iEvent % 100000 == 0) cout<<"\r"<<"iEvent "<<100*double(iEvent)/(tree->GetEntries()/reductionFactor)<<" % ";
    
    uint32_t subdetid    = int((detid-0x10000000)/0x2000000);
    float    R           = sqrt(clglobalX*clglobalX+clglobalY*clglobalY+clglobalZ*clglobalZ);
    float value = 0;
    if(observable == "maxCharge")
      value = obs*(clCorrectedSignalOverNoise)/(clSignalOverNoise);
    else
      value = obs;

    if(subdetid == 3)
      TIBrings.at(int((R-TIBRing.rMin)/((TIBRing.rMax-TIBRing.rMin)/TIBRing.nDivision)))->Fill(delay-correction,value);    
    else if(subdetid == 5)
      TOBrings.at(int((R-TOBRing.rMin)/((TOBRing.rMax-TOBRing.rMin)/TOBRing.nDivision)))->Fill(delay-correction,value);  
    else if(subdetid == 4)
      TIDrings.at(int((R-TIDRing.rMin)/((TIDRing.rMax-TIDRing.rMin)/TIDRing.nDivision)))->Fill(delay-correction,value);    
    else if(subdetid == 6  and clglobalZ > 0 and thickness > 400)
      TECPTrings.at(int((R-TECRing.rMin)/((TECRing.rMax-TECRing.rMin)/TECRing.nDivision)))->Fill(delay-correction,value);
    else if(subdetid == 6  and clglobalZ > 0 and thickness < 400)
      TECPtrings.at(int((R-TECRing.rMin)/((TECRing.rMax-TECRing.rMin)/TECRing.nDivision)))->Fill(delay-correction,value);
    else if(subdetid == 6  and clglobalZ < 0 and thickness > 400)
      TECMTrings.at(int((R-TECRing.rMin)/((TECRing.rMax-TECRing.rMin)/TECRing.nDivision)))->Fill(delay-correction,value);
    else if(subdetid == 6  and clglobalZ < 0 and thickness < 400)
      TECMtrings.at(int((R-TECRing.rMin)/((TECRing.rMax-TECRing.rMin)/TECRing.nDivision)))->Fill(delay-correction,value);
  }

  std::cout<<std::endl;
  std::cout<<"Loop on events terminated"<<std::endl;
  
  // correct profiles and fit them with a gaussian
  std::cout<<"Analyze TIB profiles"<<std::endl;
  for(auto iprof : TIBrings){
    if(observable == "maxCharge")
      correctProfile(iprof);
    fitProfile(iprof,isGaussian);
  }
  std::cout<<"Analyze TOB profiles"<<std::endl;
  for(auto iprof : TOBrings){
    if(observable == "maxCharge")
      correctProfile(iprof);
    fitProfile(iprof,isGaussian);
  }     
  std::cout<<"Analyze TID profiles"<<std::endl;
  for(auto iprof : TIDrings){
    if(observable == "maxCharge")
      correctProfile(iprof);
    fitProfile(iprof,isGaussian);
  }
  std::cout<<"Analyze TECPT profiles"<<std::endl;
  for(auto iprof : TECPTrings){
    if(observable == "maxCharge")
      correctProfile(iprof);
    fitProfile(iprof,isGaussian);
  }
  std::cout<<"Analyze TECPt profiles"<<std::endl;
  for(auto iprof : TECPtrings){
    if(observable == "maxCharge")
      correctProfile(iprof);
    fitProfile(iprof,isGaussian);
  }
  std::cout<<"Analyze TECMT profiles"<<std::endl;
  for(auto iprof : TECMTrings){
    if(observable == "maxCharge")
      correctProfile(iprof);
    fitProfile(iprof,isGaussian);
  }
  std::cout<<"Analyze TECMt profiles"<<std::endl;
  for(auto iprof : TECMtrings){
    if(observable == "maxCharge")
      correctProfile(iprof);
    fitProfile(iprof,isGaussian);
  }  

  std::cout<<"#################################"<<std::endl;
  std::cout<<"##### End of Ring Analysis #####"<<std::endl;
  std::cout<<"#################################"<<std::endl;

}

/// main function that run the analysis
void delayValidation(string file0, 
		     string file1 = "nocorrection.root", 
		     string observable   = "maxCharge",
		     bool plotPartitions = true, 
		     bool plotLayer      = true, 
		     bool plotSlices     = false, 
		     string outputDIR    = "prompt"){

  // prepare style and load macros
  setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gROOT->SetBatch(kTRUE);

  system(("mkdir -p "+outputDIR).c_str());

  // open the file and prepare the cluster tree, adding the other trees as frined --> memory consuming
  std::cout<<"#############################"<<std::endl;
  std::cout<<"###### delayValidation ######"<<std::endl;
  std::cout<<"#############################"<<std::endl;

  std::cout<<"Open Input Files"<<std::endl;
  std::auto_ptr<TFile> _file0 (TFile::Open(file0.c_str()));
  std::auto_ptr<TFile> _file1 (TFile::Open(file1.c_str()));
  std::shared_ptr<TTree> clusters   ((TTree*)_file0->FindObjectAny("clusters"));
  std::shared_ptr<TTree> readoutMap ((TTree*)_file0->FindObjectAny("readoutMap"));
  std::shared_ptr<TTree> delayCorrections ((TTree*)_file1->FindObjectAny("delayCorrections"));  
  clusters->SetEventList(0);  

  // create the plots per layer
  if(plotLayer){

    std::vector<std::shared_ptr<TProfile> > TIBlayers;
    std::vector<std::shared_ptr<TProfile> > TOBlayers;
    std::vector<std::shared_ptr<TProfile> > TIDlayers;
    std::vector<std::shared_ptr<TProfile> > TECPTlayers;
    std::vector<std::shared_ptr<TProfile> > TECPtlayers;
    std::vector<std::shared_ptr<TProfile> > TECMTlayers;
    std::vector<std::shared_ptr<TProfile> > TECMtlayers;
    LayerPlots(clusters,readoutMap,delayCorrections,TIBlayers,TOBlayers,TIDlayers,TECPTlayers,TECPtlayers,TECMTlayers,TECMtlayers,observable);
  
    std::shared_ptr<TCanvas> c1 = prepareCanvas("TIB_layers",observable);
    plotAll(c1,TIBlayers);
    c1->Print(Form("%s/TIB_layers.root",outputDIR.c_str()));

    std::shared_ptr<TCanvas> c2 = prepareCanvas("TID_layers",observable);
    plotAll(c2,TIDlayers);
    c2->Print(Form("%s/TID_layers.root",outputDIR.c_str()));

    std::shared_ptr<TCanvas> c3 = prepareCanvas("TOB_layers",observable);
    plotAll(c3,TOBlayers);
    c3->Print(Form("%s/TOB_layers.root",outputDIR.c_str()));

    std::shared_ptr<TCanvas> c4 = prepareCanvas("TECPt_layers",observable);
    plotAll(c4,TECPtlayers);
    c4->Print(Form("%s/TECPt_layers.root",outputDIR.c_str()));

    std::shared_ptr<TCanvas> c5 = prepareCanvas("TECPT_layers",observable);
    plotAll(c5,TECPTlayers);
    c5->Print(Form("%s/TECPT_layers.root",outputDIR.c_str()));

    std::shared_ptr<TCanvas> c6 = prepareCanvas("TECMt_layers",observable);
    plotAll(c6,TECMtlayers);
    c6->Print(Form("%s/TECMt_layers.root",outputDIR.c_str()));

    std::shared_ptr<TCanvas> c7 = prepareCanvas("TECMT_layers",observable);
    plotAll(c7,TECMTlayers);
    c7->Print(Form("%s/TECMT_layers.root",outputDIR.c_str()));

    std::vector<std::shared_ptr<TProfile> > alllayers;
    alllayers.reserve(TIBlayers.size()+TIDlayers.size()+TOBlayers.size()+TECPtlayers.size()+TECPTlayers.size()+TECMtlayers.size()+TECMTlayers.size());
    alllayers.insert(alllayers.end(),TIBlayers.begin(),TIBlayers.end());
    alllayers.insert(alllayers.end(),TIDlayers.begin(),TIDlayers.end());
    alllayers.insert(alllayers.end(),TOBlayers.begin(),TOBlayers.end());
    alllayers.insert(alllayers.end(),TECPtlayers.begin(),TECPtlayers.end());
    alllayers.insert(alllayers.end(),TECPTlayers.begin(),TECPTlayers.end());
    alllayers.insert(alllayers.end(),TECMtlayers.begin(),TECMtlayers.end());
    alllayers.insert(alllayers.end(),TECMTlayers.begin(),TECMTlayers.end());
    // store all the different plots
    std::shared_ptr<TCanvas> c8 (new TCanvas("c_layers","",800,650));  
    plotMaxima(c8,alllayers,outputDIR,"layers");
    
  }
  
  if(plotSlices){

    // create the plots in R slices 
    std::vector<std::shared_ptr<TProfile> > TIBrs;
    std::vector<std::shared_ptr<TProfile> > TIDrs;
    std::vector<std::shared_ptr<TProfile> > TOBrs;
    std::vector<std::shared_ptr<TProfile> > TECPTrs;
    std::vector<std::shared_ptr<TProfile> > TECPtrs;
    std::vector<std::shared_ptr<TProfile> > TECMTrs;
    std::vector<std::shared_ptr<TProfile> > TECMtrs;

    RPlots(clusters,readoutMap,delayCorrections,TIBrs,TOBrs,TIDrs,TECPTrs,TECPtrs,TECMTrs,TECMtrs,observable);
    
    std::shared_ptr<TCanvas> c1b = prepareCanvas("TIB_distance",observable);
    plotAll(c1b,TIBrs);
    c1b->Print(Form("%s/TIB_distance.root",outputDIR.c_str()));

    std::shared_ptr<TCanvas> c2b = prepareCanvas("TID_distance",observable);
    plotAll(c2b,TIDrs);
    c2b->Print(Form("%s/TID_distance.root",outputDIR.c_str()));
    
    std::shared_ptr<TCanvas> c3b = prepareCanvas("TOB_distance",observable);
    plotAll(c3b,TOBrs);
    c3b->Print(Form("%s/TOB_distance.root",outputDIR.c_str()));

    std::shared_ptr<TCanvas> c4b = prepareCanvas("TECPt_distance",observable);
    plotAll(c4b,TECPtrs);
    c4b->Print(Form("%s/TECPt_distance.root",outputDIR.c_str()));

    std::shared_ptr<TCanvas> c5b = prepareCanvas("TECPT_distance",observable);
    plotAll(c5b,TECPTrs);
    c5b->Print(Form("%s/TECPT_distance.root",outputDIR.c_str()));

    std::shared_ptr<TCanvas> c6b = prepareCanvas("TECMt_distance",observable);
    plotAll(c6b,TECMtrs);
    c6b->Print(Form("%s/TECMt_distance.root",outputDIR.c_str()));

    std::shared_ptr<TCanvas> c7b = prepareCanvas("TECMT_distance",observable);
    plotAll(c7b,TECMTrs);
    c7b->Print(Form("%s/TECMT_distance.root",outputDIR.c_str()));


    std::vector<std::shared_ptr<TProfile> > allrs;
    allrs.reserve(TIBrs.size()+TIDrs.size()+TOBrs.size()+TECMTrs.size()+TECMtrs.size()+TECPTrs.size()+TECPtrs.size());
    allrs.insert(allrs.end(),TIBrs.begin(),TIBrs.end());
    allrs.insert(allrs.end(),TIDrs.begin(),TIDrs.end());
    allrs.insert(allrs.end(),TOBrs.begin(),TOBrs.end());
    allrs.insert(allrs.end(),TECMTrs.begin(),TECMTrs.end());
    allrs.insert(allrs.end(),TECMtrs.begin(),TECMtrs.end());
    allrs.insert(allrs.end(),TECPTrs.begin(),TECPTrs.end());
    allrs.insert(allrs.end(),TECPtrs.begin(),TECPtrs.end());

    std::shared_ptr<TCanvas> c8b (new TCanvas("c_rings","",800,650));  
    plotMaxima(c8b,allrs,outputDIR,"rings");
  }
  
  if(plotPartitions){

    std::vector<std::shared_ptr<TProfile> > TIB;
    std::vector<std::shared_ptr<TProfile> > TID;
    std::vector<std::shared_ptr<TProfile> > TOB;  
    std::vector<std::shared_ptr<TProfile> > TECPT; 
    std::vector<std::shared_ptr<TProfile> > TECPt; 
    std::vector<std::shared_ptr<TProfile> > TECMT;
    std::vector<std::shared_ptr<TProfile> > TECMt;

    //change ring definition to collapse all of them
    TIBRing.nDivision = 1;
    TIDRing.nDivision = 1;
    TOBRing.nDivision = 1;
    TECRing.nDivision = 1;
    
    RPlots(clusters,readoutMap,delayCorrections,TIB,TOB,TID,TECPT,TECPt,TECMT,TECMt,observable);
        
    // create the plots per partition
    std::vector<std::shared_ptr<TProfile> > allPartitions;
    allPartitions.reserve(TIB.size()+TOB.size()+TID.size()+TECPT.size()+TECPt.size()+TECMT.size()+TECMt.size());
    allPartitions.insert(allPartitions.end(),TIB.begin(),TIB.end());
    allPartitions.insert(allPartitions.end(),TID.begin(),TID.end());
    allPartitions.insert(allPartitions.end(),TOB.begin(),TOB.end());
    allPartitions.insert(allPartitions.end(),TECMT.begin(),TECMT.end());
    allPartitions.insert(allPartitions.end(),TECMt.begin(),TECMt.end());
    allPartitions.insert(allPartitions.end(),TECPT.begin(),TECPT.end());
    allPartitions.insert(allPartitions.end(),TECPt.begin(),TECPt.end());

    std::shared_ptr<TCanvas> c1 = prepareCanvas("Partitions",observable);
    plotAll(c1,allPartitions);
    c1->Print(Form("%s/Partitions.root",outputDIR.c_str()));
  }
}

