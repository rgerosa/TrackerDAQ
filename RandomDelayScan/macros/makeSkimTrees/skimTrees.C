#include <string>
#include "TFile.h"
#include "TTree.h"

// codes that run on TrackerDPGAnalysis output, applies selections according to event/tracks/vertex/cluster properities
// creates an output file where only the cluster tree for the selected events is stored, together with the PSU and delay maps

void skimTrees(string inputFileName, string outputFileName, bool isBOn = true) {


  cout<<"################################"<<endl;
  cout<<"#### Skim TrackerDPG Trees #####"<<endl;
  cout<<"################################"<<endl;

  vector<TFile*> inputFile;
  vector<TTree*> clustersTree;
  vector<TTree*> eventTree;
  vector<TTree*> trackTree;
  vector<TTree*> vertexTree;

  TString name (inputFileName.c_str());
  if(name.Contains(".root")){// single file
    inputFile.push_back(TFile::Open(inputFileName.c_str()));
    
    clustersTree.push_back((TTree*) inputFile.back()->FindObjectAny("clusters"));
    if(clustersTree.back() == 0 or clustersTree.back() == NULL){
      cout<<"[skimTrees] no cluster tree found --> problem "<<endl;
      return;
    }

    eventTree.push_back((TTree*) inputFile.back()->FindObjectAny("events"));
    if(eventTree.back() == 0 or eventTree.back() == NULL){
      cout<<"[skimTrees] no event tree found --> problem "<<endl;
      return;
    }    
    eventTree.back()->BuildIndex("runid","eventid");

    trackTree.push_back((TTree*) inputFile.back()->FindObjectAny("tracks0"));
    if(trackTree.back() == 0 or trackTree.back() == NULL){
      cout<<"[skimTrees] no track tree found --> problem "<<endl;
      return;
    }
    trackTree.back()->BuildIndex("trackid0","eventid");

    vertexTree.push_back((TTree*) inputFile.back()->FindObjectAny("vertices"));
    if(vertexTree.back() == 0 or vertexTree.back() == NULL){
      cout<<"[skimTrees] no track tree found --> problem "<<endl;
      return;
    }
    vertexTree.back()->BuildIndex("vertexid","eventid");
  }
  else{ // list of file in a text  one

    ifstream infile (inputFileName.c_str());
    if(infile.is_open()){
      string filename;
      while(!infile.eof()){
	getline(infile,filename);
	if(filename == "" or not TString(filename.c_str()).Contains("root")) continue;
	inputFile.push_back(TFile::Open(filename.c_str()));

	clustersTree.push_back((TTree*) inputFile.back()->FindObjectAny("clusters"));
	if(clustersTree.back() == 0 or clustersTree.back() == NULL){
	  cout<<"[skimTrees] no cluster tree found --> problem --> skip "<<endl;
	  continue;
	}
	
	eventTree.push_back((TTree*) inputFile.back()->FindObjectAny("events"));
	if(eventTree.back() == 0 or eventTree.back() == NULL){
	  cout<<"[skimTrees] no event tree found --> problem --> skip "<<endl;
	  continue;
	}    
	eventTree.back()->BuildIndex("runid","eventid");
	
	trackTree.push_back((TTree*) inputFile.back()->FindObjectAny("tracks0"));
	if(trackTree.back() == 0 or trackTree.back() == NULL){
	  cout<<"[skimTrees] no track tree found --> problem --> skip "<<endl;
	  continue;
	}
	trackTree.back()->BuildIndex("trackid0","eventid");

	vertexTree.push_back((TTree*) inputFile.back()->FindObjectAny("vertices"));
	if(vertexTree.back() == 0 or vertexTree.back() == NULL){
	  cout<<"[skimTrees] no track tree found --> problem --> skip "<<endl;
	  continue;
	}
	vertexTree.back()->BuildIndex("vertexid","eventid");
      }
    }
    infile.close();
  }
  
  //loop on clustersTree entry and
  cout<<"### Create outputFile, define selection and copyTree "<<endl;  
  TFile* outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  outputFile->cd();

  // apply selections
  string eventSelection;
  string trackSelection;
  string vertexSelection;
  string clusterSelection;
  if(isBOn){
    eventSelection   = "nVertices > 0 && ";
    trackSelection   = "pt > 0.9 && quality>2 && pterr/pt < 0.2 && ";
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

  vector<TTree*> outputTree;
  TList* tree_list = new TList();
  for(size_t itree = 0; itree < clustersTree.size(); itree++){

    // in order to apply selections
    clustersTree.at(itree)->AddFriend(eventTree.at(itree));
    clustersTree.at(itree)->AddFriend(trackTree.at(itree));
    clustersTree.at(itree)->AddFriend(vertexTree.at(itree));
    
    outputTree.push_back(clustersTree.at(itree)->CopyTree((eventSelection+trackSelection+vertexSelection+clusterSelection).c_str()));
    tree_list->Add(outputTree.back());
  }
  
  TTree* final_output_tree = TTree::MergeTrees(tree_list); 
  final_output_tree->Write("clusters",TObject::kOverwrite);

  /// Just use the first file here
  cout<<"### Copy the PSU map in the output map "<<endl;
  //copy the PSU map
  TTree* PSUmap = (TTree*) inputFile.front()->FindObjectAny("psumap");
  if(PSUmap == 0 or PSUmap == NULL){
    cout<<"[skimTrees] no PSU map found --> problem "<<endl;
    return;
  }      
  PSUmap->CloneTree()->Write("psumap",TObject::kOverwrite);

  cout<<"### Copy the readout map in the output map "<<endl;
  TTree* readoutMap = (TTree*) inputFile.front()->FindObjectAny("readoutMap");
  if(readoutMap == 0 or readoutMap == NULL){
    cout<<"[skimTrees] no readoutMap found --> problem "<<endl;
    return;
  }      
  readoutMap->CloneTree()->Write("readoutMap",TObject::kOverwrite);

  std::cout << "Saving merged file." << std::endl;
  outputFile->Close();
  for(auto file : inputFile)
    file->Close();
}
