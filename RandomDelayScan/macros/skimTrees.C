#include <string>
#include "TFile.h"
#include "TTree.h"

// codes that run on TrackerDPGAnalysis output, applies selections according to event/tracks/vertex/cluster properities
// creates an output file where only the cluster tree for the selected events is stored, together with the PSU and delay maps

void skimTrees(string inputFileName, string outputFileName, bool isBOn = true) {

  cout<<"################################"<<endl;
  cout<<"#### Skim TrackerDPG Trees #####"<<endl;
  cout<<"################################"<<endl;
  
  TFile* inputFile = TFile::Open(inputFileName.c_str());

  cout<<"### Load clusters tree "<<endl;
  TTree* clustersTree = (TTree*) inputFile->FindObjectAny("clusters");
  if(clustersTree == 0 or clustersTree == NULL){
    cout<<"[skimTrees] no cluster tree found --> problem "<<endl;
    return;
  }

  // event tree
  cout<<"### Load and index the event tree "<<endl;
  TTree* eventTree = (TTree*) inputFile->FindObjectAny("events");
  if(eventTree == 0 or eventTree == NULL){
    cout<<"[skimTrees] no event tree found --> problem "<<endl;
    return;
  }
  eventTree->BuildIndex("runid","eventid");

  // track tree
  cout<<"### Load and index the tracks0 tree "<<endl;
  TTree* trackTree = (TTree*) inputFile->FindObjectAny("tracks0");
  if(trackTree == 0 or trackTree == NULL){
    cout<<"[skimTrees] no track tree found --> problem "<<endl;
    return;
  }
  trackTree->BuildIndex("trackid0","eventid");

  // vertex tree
  cout<<"### Load and index the vertex tree "<<endl;
  TTree* vertexTree = (TTree*) inputFile->FindObjectAny("vertices");
  if(vertexTree == 0 or vertexTree == NULL){
    cout<<"[skimTrees] no track tree found --> problem "<<endl;
    return;
  }
  vertexTree->BuildIndex("vertexid","eventid");
  

  //loop on clustersTree entry and
  cout<<"### Create outputFile, define selection and copyTree "<<endl;  
  TFile* outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  outputFile->cd();
  // in order to apply selections
  clustersTree->AddFriend(eventTree);
  clustersTree->AddFriend(trackTree);
  clustersTree->AddFriend(vertexTree);
  // apply selections
  string eventSelection;
  string trackSelection;
  string vertexSelection;
  string clusterSelection;
  if(isBOn){
    eventSelection   = "lowPixelProbabilityFraction[0] < 0.4 && lowPixelProbabilityFraction[0] >-0.5 && nVertices > 0 && ";
    trackSelection   = "pt > 1 && quality>2 && pterr/pt < 0.2 && dedx1 < 5 && ";
    vertexSelection  = "";
    clusterSelection = "onTrack && angle > 0 && maxCharge < 254";
  }
  else{
    eventSelection   = "nVertices>=0 && ";
    trackSelection   = "pt>= 0 && ";
    vertexSelection  = "";
    clusterSelection = "onTrack && angle > 0 && maxCharge < 254";
  }

  cout<<"### eventSelection:   "<<eventSelection<<endl;
  cout<<"### trackSelection:   "<<trackSelection<<endl;
  cout<<"### vertexSelection:  "<<vertexSelection<<endl;
  cout<<"### clusterSelection: "<<clusterSelection<<endl;
  cout<<"### totalSelection:   "<<eventSelection+trackSelection+vertexSelection+clusterSelection<<endl;
  TTree* outputTree = clustersTree->CopyTree((eventSelection+trackSelection+vertexSelection+clusterSelection).c_str());
  outputTree->Write("clusters",TObject::kOverwrite);
  
  cout<<"### Copy the PSU map in the output map "<<endl;
  //copy the PSU map
  TTree* PSUmap = (TTree*) inputFile->FindObjectAny("psumap");
  if(PSUmap == 0 or PSUmap == NULL){
    cout<<"[skimTrees] no PSU map found --> problem "<<endl;
    return;
  }      
  PSUmap->CloneTree()->Write("psumap",TObject::kOverwrite);

  cout<<"### Copy the readout map in the output map "<<endl;
  TTree* readoutMap = (TTree*) inputFile->FindObjectAny("readoutMap");
  if(readoutMap == 0 or readoutMap == NULL){
    cout<<"[skimTrees] no readoutMap found --> problem "<<endl;
    return;
  }      
  readoutMap->CloneTree()->Write("readoutMap",TObject::kOverwrite);

  std::cout << "Saving merged file." << std::endl;
  outputFile->Close();
  inputFile->Close();
}
