void computeCorrections(const TGraph* delays, float* corrections)
{
  Double_t x,y;
  for(int i=0; i<delays->GetN(); ++i) {
    delays->GetPoint(i,x,y);
    corrections[int(x)] = int(y*24/25.+0.5)*25/24.;
    //corrections[int(x)] = y;
  }
}

int findBin(unsigned int detid, float globalX, float globalY, float globalZ) 
{
  // the definitions used are:
  // "subdetid","int((detid-0x10000000)/0x2000000)"
  // "R","sqrt(globalX**2+globalY**2+globalZ**2)"

  int subdetid = int((detid-0x10000000)/0x2000000);
  float R = sqrt(globalX*globalX+globalY*globalY+globalZ*globalZ);

  if(subdetid == 3) {
    int bin = 1+int((R-20)/10.);
    bin = bin>0 ? bin : 1;
    bin = bin<7 ? bin : 6;
    return bin-1;
  } else if (subdetid == 4) {
    int bin = 7+int((R-80)/10.);
    bin = bin>6 ? bin : 7;
    bin = bin<11 ? bin : 10;
    return bin-1;
  } else if (subdetid == 5) {
    int bin = 11+int((R-60)/5.);
    bin = bin>10 ? bin : 11;
    bin = bin<27 ? bin : 26;
    return bin-1;
  } else if (subdetid == 6) {
    int bin = 27+int((R-120)/20.);
    bin = bin>26 ? bin : 27;
    bin = bin<36 ? bin : 35;
    if(sqrt(globalX*globalX+globalY*globalY)<60) bin += 9;
    if(globalZ>0) bin += 18;
    return bin-1;
  } else {
    return 0;
  }

}

TTree* delayCorrection(TGraph* input, TTree* readoutmap, const char* output = "computedCorrections.root")
{
  // this routine takes as input a TGraph with best delays for various subsets of modules
  // and the readoutmap tree.
  // It creates a new readoutmap friend tree with the correction for each module
  // and puts it in output file.
  
  // warning: the findBin method must be adapted to the binning of the reference graph.

  // the code is organized as follow:
  // * first analyze the input to fill a map bin->correction
  // * then prepare the friend and loop on the map. For each entry one calls findBin(detid,R).

  // compute corrections
  std::cout << "computing corrections" << std::endl;
  float corrections[62] = {0.}; // array of corrections
  float theCorrection = 0; // correction 
  computeCorrections(input,corrections);

  // prepare reading of the readoutmap (taking only relevant variables)
  unsigned int detid; 
  float globalX, globalY, globalZ, delay;
  readoutmap->SetBranchAddress("detid",&detid);
  readoutmap->SetBranchAddress("globalX",&globalX);
  readoutmap->SetBranchAddress("globalY",&globalY);
  readoutmap->SetBranchAddress("globalZ",&globalZ);
  
  // create a file for output
  TFile outputFile(output,"recreate");
  TTree *correctionTree = new TTree("delayCorrections","delayCorrections");
  correctionTree->Branch("correction",&theCorrection,"correction/F");
  correctionTree->Branch("detid",&detid,"detid/i");

  // loop over map and fill correction tree
  std::vector<unsigned int> corrids;

  std::cout << "filling the tree" << std::endl;
  Int_t nentries = (Int_t)readoutmap->GetEntries();
  for (Int_t i = 0; i < nentries; i++){
    bool existsAlready = false;
    readoutmap->GetEntry(i);
    theCorrection = corrections[findBin(detid,globalX,globalY,globalZ)];
    for (std::size_t k = 0; k < corrids.size(); k++) {
        if (detid == corrids[k]) existsAlready = true;
    }
    
    if (!existsAlready) {
        corrids.push_back(detid);
        correctionTree->Fill();
    }
  }
  
  // write output
  outputFile.cd();
  correctionTree->BuildIndex("detid");
  correctionTree->Write();

  return correctionTree;
}


