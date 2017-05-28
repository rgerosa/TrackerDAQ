
void makePSUDetIdMap(string inputFileWithReadOutMap = "$CMSSW_BASE/src/TrackerDAQAnalysis/RandomDelayScan/data/readoutMap.root",
		     string PSUDCUMap      = "../../data/PSUmapping.csv",
		     string outputFileName = "$CMSSW_BASE/src/TrackerDAQAnalysis/RandomDelayScan/data/readoutMapWithPSU.root"){


  TFile* inputReadoutMap = TFile::Open(inputFileWithReadOutMap.c_str(),"READ");
  TTree* tree = (TTree*) inputReadoutMap->Get("readoutMap");

  // set branches
  uint32_t Detid_i, dcuId_i;
  tree->SetBranchAddress("dcuId",&dcuId_i);
  
  cout<<"##### Read PSU map "<<endl;
  map<uint32_t,string> PSUDUCmapping;
  ifstream mapFile (PSUDCUMap.c_str());
  string line;
  if(mapFile.is_open()){
    char buffer[1024];
    while(!mapFile.eof()){
      mapFile.getline(buffer,1024);      
      std::istringstream line(buffer);
      std::string name;
      std::string dcuNumber;
      // one line contains the PSU name + all dcuids connected to it.                                                                                                                                  
      line >> name;
      if(name == "") continue;
      while(!line.eof()) {
	line >> dcuNumber;
	if(PSUDUCmapping.find(stoi(dcuNumber)) != PSUDUCmapping.end() and PSUDUCmapping[stoi(dcuNumber)] != name){
	  cerr<<"Problem .. this dcu connected to more then one PG : dcu "<<dcuNumber<<" PG "<<name<<" PG2 "<<PSUDUCmapping[stoi(dcuNumber)]<<endl;
	  continue;
	}
	PSUDUCmapping[stoi(dcuNumber)] = name;
      }
    }
  }

  //create output tree and outputFile
  cout<<"##### Clone tree and loop on events"<<endl;
  TFile* outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  TTree* treeOut = (TTree*) tree->CloneTree(0);
  string PSUname;
  treeOut->Branch("PSUname",&PSUname);
 
  for(int iEvent = 0; iEvent < tree->GetEntries(); iEvent++){
    tree->GetEntry(iEvent);
    PSUname = PSUDUCmapping[dcuId_i];
    treeOut->Fill();
  }
  
  treeOut->BuildIndex("detid");
  treeOut->Write("readoutMap");
  outputFile->Close();
}
