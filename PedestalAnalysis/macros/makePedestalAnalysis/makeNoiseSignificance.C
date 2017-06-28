#include "../CMS_lumi.h"
#include "../TrackerStrip.h"

static float minimumRMS = 2;
static float maximumRMS = 30;
static int reductionFactor  = 1;

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

void makeNoiseSignificance(string inputFileName, string outputDIR, string outputName){

  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();
  gROOT->SetBatch(kTRUE);

  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");  // read the input tree                                                                                                           
  TTree* tree = (TTree*) inputFile->Get("pedestalFullNoise");

  // Set only useful branches
  uint32_t detid,fedKey;
  uint16_t fecCrate,fecSlot, fecRing, ccuAdd, ccuChan, lldChannel, fedId, fedCh, apvId, stripId;
  float    noise, noiseRMS;
  
  tree->SetBranchStatus("*",kFALSE);
  tree->SetBranchStatus("detid",kTRUE);
  tree->SetBranchStatus("lldChannel",kTRUE);
  tree->SetBranchStatus("apvId",kTRUE);
  tree->SetBranchStatus("stripId",kTRUE);
  tree->SetBranchStatus("noise",kTRUE);
  tree->SetBranchStatus("noiseRMS",kTRUE);

  tree->SetBranchAddress("detid",&detid);
  tree->SetBranchAddress("lldChannel",&lldChannel);
  tree->SetBranchAddress("apvId",&apvId);
  tree->SetBranchAddress("stripId",&stripId);
  tree->SetBranchAddress("noise",&noise);
  tree->SetBranchAddress("noiseRMS",&noiseRMS);

  // loop on the pedestak analysis tree                                                                                                                                                             
  map<string,noiseStrip*> noiseSpread;
  map<string,noiseStrip*> noiseMean;

  // First loop to evaluate mean value of the noise per APV level                                                                                                                                   
  cout<<"Loop to evaluate Mean noise across each APV"<<endl;

  for(long int iChannel = 0; iChannel < tree->GetEntries(); iChannel++){
    tree->GetEntry(iChannel);
    cout.flush();
    if(iChannel %10000 == 0) cout<<"\r"<<"iChannel "<<100*double(iChannel)/(tree->GetEntries()/reductionFactor)<<" % ";
    if(iChannel > double(tree->GetEntries())/reductionFactor) break;
    
    string nameNoise = "Detid_"+to_string(detid)+"_lldCh"+to_string(lldChannel)+"_apv_"+to_string(apvId);
    
    if(noiseMean[nameNoise] == 0 or noiseMean[nameNoise] == NULL)
      noiseMean[nameNoise] = new noiseStrip();
    noiseMean[nameNoise]->nstrip += 1;
    noiseMean[nameNoise]->noiseVal += noise;
    noiseMean[nameNoise]->isDivided = false;
  }

  std::cout<<std::endl;
  for(auto iapv : noiseMean){
    iapv.second->noiseVal /= float(iapv.second->nstrip);
    iapv.second->isDivided = true;
  }

  cout<<"Loop to evaluate Noise spread across each APV"<<endl;
  for(long int iChannel = 0; iChannel < tree->GetEntries(); iChannel++){
    tree->GetEntry(iChannel);
    cout.flush();
    if(iChannel %10000 == 0) cout<<"\r"<<"iChannel "<<100*double(iChannel)/(tree->GetEntries()/reductionFactor)<<" % ";
    if(iChannel > double(tree->GetEntries())/reductionFactor) break;

    // to evaluate mean value and spread of noise                                                                                                                                                      
    string nameNoise = "Detid_"+to_string(detid)+"_lldCh"+to_string(lldChannel)+"_apv_"+to_string(apvId);

    if(noiseSpread[nameNoise] == 0 or noiseSpread[nameNoise] == NULL)
      noiseSpread[nameNoise] = new noiseStrip();
    noiseSpread[nameNoise]->nstrip += 1;
    noiseSpread[nameNoise]->noiseVal += (noise-noiseMean[nameNoise]->noiseVal)*(noise-noiseMean[nameNoise]->noiseVal);
    noiseSpread[nameNoise]->isDivided = false;
  }
  std::cout<<std::endl;

  for(auto iapv : noiseSpread){
    float val = iapv.second->noiseVal;
    iapv.second->noiseVal = sqrt(val)/sqrt(iapv.second->nstrip-1);
    iapv.second->isDivided = true;
  }

  
  TFile* output = new TFile((outputDIR+"/"+outputName).c_str(),"RECREATE");
  output->cd();

  // make a clone of the tree --> clone all branches  
  tree->SetBranchStatus("*",kTRUE);
  TTree* outtree = tree->CloneTree(0);
  
  // new branches in the output tree
  float noiseSignificance, noiseMeanAPV, noiseSpreadAPV;
  TBranch* b_noiseMeanAPV =  outtree->Branch("noiseMeanAPV",&noiseMeanAPV,"noiseMeanAPV/F");
  TBranch* b_noiseSpreadAPV = outtree->Branch("noiseSpreadAPV",&noiseSpreadAPV,"noiseSpreadAPV/F");
  TBranch* b_noiseSignificance = outtree->Branch("noiseSignificance",&noiseSignificance,"noiseSignificance/F");

  cout<<"Loop to fill new branches"<<endl;
  for(long int iChannel = 0; iChannel < tree->GetEntries(); iChannel++){
    tree->GetEntry(iChannel);
    cout.flush();
    if(iChannel %10000 == 0) cout<<"\r"<<"iChannel "<<100*double(iChannel)/(tree->GetEntries()/reductionFactor)<<" % ";
    if(iChannel > double(tree->GetEntries())/reductionFactor) break;

    string nameNoise = "Detid_"+to_string(detid)+"_lldCh"+to_string(lldChannel)+"_apv_"+to_string(apvId);
    
    noiseMeanAPV = noiseMean[nameNoise]->noiseVal;
    noiseSpreadAPV = noiseSpread[nameNoise]->noiseVal;
    noiseSignificance = (noiseRMS-noiseMeanAPV)/noiseSpreadAPV;
    outtree->Fill();
    
  }
  outtree->Write(tree->GetName());
  output->Close();
}
