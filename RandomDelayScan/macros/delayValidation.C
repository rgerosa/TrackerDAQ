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

#include "CMS_lumi.h"
#include "TkPulseShape.h"
 
// basic profile for maxCharge
static TProfile frame("frame","",20,-10,10,0,256);
// parametrize profile with a gaussian shape
static bool isGaussian = true;
// reduce the number of events by
static int  reductionFactor = 20;

// struct to handle ring definition
class trackerRing{

public:
  trackerRing(){};
  ~trackerRing(){};
  trackerRing(float rMin, float rMax, int nDivision):
    rMin(rMin),
    rMax(rMax),
    nDivision(nDivision){};

  float rMin;
  float rMax;
  int   nDivision;
};

static trackerRing TIBRing(20,80,6);
static trackerRing TIDRing(80,120,4);
static trackerRing TOBRing(60,150,18);
static trackerRing TECRing(120,300,9);

////
float limit(float SoNcut){
  return 3.814567e+00+8.336601e+00*SoNcut-1.511334e-01*pow(SoNcut,2);
}

float correctMeasurement(float mean, float SoNcut){
  if(mean>limit(SoNcut))
      return -8.124872e+00+9.860108e-01*mean-3.618158e-03*pow(mean,2)+2.037263e-05*pow(mean,3);
  else return 0.;
}

void correctProfile(const std::shared_ptr<TProfile> & profile){
  int nbins = profile->GetNbinsX();
  float min = limit(3.); // zero suppression level .. something like S/Noise = 3
  for(int bin=1;bin<=nbins;++bin){
    if(profile->GetBinContent(bin)<min) { // set to zero
      profile->SetBinContent(bin,0.);
      profile->SetBinError(bin,0.);
      profile->SetBinEntries(bin,0);
    }
    else 
      profile->SetBinContent(bin,profile->GetBinEntries(bin)*correctMeasurement(profile->GetBinContent(bin),3.)); // correct the measurement 
  }
}

std::shared_ptr<TF1> fitProfile(const std::shared_ptr<TProfile> & prof, bool gaus = false){ 

  std::shared_ptr<TF1> pulse;
  if(gaus) { // gaussina fit
    pulse = std::shared_ptr<TF1>(new TF1(Form("Gaus_%s",prof->GetName()),"gaus(0)",-10,10));
    pulse->SetParameters(50,0,12);
  } else {// different fit for the profile
    pulse = std::shared_ptr<TF1>(TkPulseShape::GetDeconvFitter());
    pulse->SetName(Form("SignalFit_%s",prof->GetName()));
    pulse->SetParameters(0,0,3,50,15);
    pulse->FixParameter(0,0);
    pulse->FixParameter(3,50);
  }
  gROOT->SetBatch(1);
  prof->Fit(pulse.get(),"");
  std::cout << "Profile Name "<<prof->GetName()<<" "<<"Maximum at " << pulse->GetMaximumX(-10,10) << std::endl;
  return pulse;
}


void LayerPlots(const std::shared_ptr<TTree> & tree, 
		const std::shared_ptr<TTree> & map,
		const std::shared_ptr<TTree> & corrections,
		std::vector<std::shared_ptr<TProfile> > & TIBlayers,
		std::vector<std::shared_ptr<TProfile> > & TOBlayers,
		std::vector<std::shared_ptr<TProfile> > & TIDlayers,
		std::vector<std::shared_ptr<TProfile> > & TECPTlayers,
		std::vector<std::shared_ptr<TProfile> > & TECPtlayers,
		std::vector<std::shared_ptr<TProfile> > & TECMTlayers,
		std::vector<std::shared_ptr<TProfile> > & TECMtlayers) {

  std::cout<<"######################################################"<<std::endl;
  std::cout << "Preparing Layer plot for the fifferent subdetectors "<< std::endl; 
  std::cout<<"######################################################"<<std::endl;
  
  // TTree Reader appear not to be working with addFriend and EventList
  uint32_t detid;
  float    maxCharge, clCorrectedSignalOverNoise, clSignalOverNoise, clglobalX, clglobalY, clglobalZ, thickness;
  tree->SetBranchStatus("*",kFALSE);
  tree->SetBranchStatus("detid",kTRUE);
  tree->SetBranchStatus("runid",kTRUE);
  tree->SetBranchStatus("eventid",kTRUE);
  tree->SetBranchStatus("trackid0",kTRUE);
  tree->SetBranchStatus("maxCharge",kTRUE);
  tree->SetBranchStatus("clCorrectedSignalOverNoise",kTRUE);
  tree->SetBranchStatus("clSignalOverNoise",kTRUE);
  tree->SetBranchStatus("clglobalX",kTRUE);
  tree->SetBranchStatus("clglobalY",kTRUE);
  tree->SetBranchStatus("clglobalZ",kTRUE);
  tree->SetBranchStatus("thickness",kTRUE);
  tree->SetBranchAddress("detid",&detid);
  tree->SetBranchAddress("maxCharge",&maxCharge);
  tree->SetBranchAddress("clCorrectedSignalOverNoise",&clCorrectedSignalOverNoise);
  tree->SetBranchAddress("clSignalOverNoise",&clSignalOverNoise);
  tree->SetBranchAddress("clglobalX",&clglobalX);
  tree->SetBranchAddress("clglobalY",&clglobalY);
  tree->SetBranchAddress("clglobalZ",&clglobalZ);
  tree->SetBranchAddress("thickness",&thickness);

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
  
  // create vectors
  for(int ilayer = 0; ilayer < 4; ilayer++){
    std::shared_ptr<TProfile> temp (new TProfile(Form("TIB_layer%d",ilayer+1),"",20,-10,10,0,256));
    TIBlayers.push_back(temp);
  }
  for(int ilayer = 0; ilayer < 6; ilayer++){
    std::shared_ptr<TProfile> temp (new TProfile(Form("TOB_layer%d",ilayer+1),"",20,-10,10,0,256));
    TOBlayers.push_back(temp);
  }
  for(int ilayer = 0; ilayer < 3; ilayer++){
    std::shared_ptr<TProfile> temp (new TProfile(Form("TID_layer%d",ilayer+1),"",20,-10,10,0,256));
    TIDlayers.push_back(temp);
  }
  for(int ilayer = 0; ilayer < 9; ilayer++){
    std::shared_ptr<TProfile> temp (new TProfile(Form("TECP_layer%d_T",ilayer+1),"",20,-10,10,0,256));
    TECPTlayers.push_back(temp);
    temp = std::shared_ptr<TProfile>(new TProfile(Form("TECP_layer%d_t",ilayer+1),"",20,-10,10,0,256));
    TECPtlayers.push_back(temp);
  }
  for(int ilayer = 0; ilayer < 9; ilayer++){
    std::shared_ptr<TProfile> temp (new TProfile(Form("TECM_layer%d_T",ilayer+1),"",20,-10,10,0,256));
    TECMTlayers.push_back(temp);
    temp = std::shared_ptr<TProfile>(new TProfile(Form("TECM_layer%d_t",ilayer+1),"",20,-10,10,0,256));
    TECMtlayers.push_back(temp);
  }

  TEntryList* entryList = tree->GetEntryList();
  
  std::cout<<"Tree with nEntries "<<tree->GetEntries()<<", surviving selections "<<entryList->GetN()<<std::endl;
  if(tree->GetEntries() < entryList->GetN()){
    cout<<"[LayerPlots] Huston we have a problem with surviving events --> return "<<endl;
    return;
  }

  long int iEntry  = 0;
  long int iEvent  = 0;

  for( ; iEvent < entryList->GetN()/reductionFactor; iEvent++){    
    // take the selected event in the main tree and skip in case
    iEntry = entryList->GetEntry(iEvent);
    if(iEntry < 0) continue;
    // take the event
    tree->GetEntry(iEntry);
    // take the map delay from the detid
    map->GetEntryWithIndex(detid);
    // take the correction from the detid
    corrections->GetEntryWithIndex(detid);

    cout.flush();
    if(iEvent % 100000 == 0) cout<<"\r"<<"iEvent "<<100*double(iEvent)/(entryList->GetN()/reductionFactor)<<" % ";

    uint32_t subdetid    = int((detid-0x10000000)/0x2000000);
    uint32_t barrellayer = int((detid%33554432)/0x4000);
    uint32_t TIDlayer    = int((detid%33554432)/0x800)%4;
    uint32_t TECPlayer   = int((detid%33554432)/0x4000)-32;
    uint32_t TECMlayer   = int((detid%33554432)/0x4000)-16;
    float    R           = sqrt(clglobalX*clglobalX+clglobalY*clglobalY+clglobalZ*clglobalZ);

    if(subdetid == 3)      
      TIBlayers.at(barrellayer-1)->Fill(delay-correction,maxCharge*(clCorrectedSignalOverNoise)/(clSignalOverNoise));
    else if(subdetid == 5)
      TOBlayers.at(barrellayer-1)->Fill(delay-correction,maxCharge*(clCorrectedSignalOverNoise)/(clSignalOverNoise));
    
    else if(subdetid == 4)
      TIDlayers.at(TIDlayer-1)->Fill(delay-correction,maxCharge*(clCorrectedSignalOverNoise)/(clSignalOverNoise));
    
    else if(subdetid == 6  and clglobalZ > 0 and thickness > 400)
      TECPTlayers.at(TECPlayer-1)->Fill(delay-correction,maxCharge*(clCorrectedSignalOverNoise)/(clSignalOverNoise));
    else if(subdetid == 6  and clglobalZ > 0 and thickness < 400)
      TECPtlayers.at(TECPlayer-1)->Fill(delay-correction,maxCharge*(clCorrectedSignalOverNoise)/(clSignalOverNoise));
    else if(subdetid == 6  and clglobalZ < 0 and thickness > 400)
      TECMTlayers.at(TECMlayer-1)->Fill(delay-correction,maxCharge*(clCorrectedSignalOverNoise)/(clSignalOverNoise));
    else if(subdetid == 6  and clglobalZ < 0 and thickness < 400)
      TECMtlayers.at(TECMlayer-1)->Fill(delay-correction,maxCharge*(clCorrectedSignalOverNoise)/(clSignalOverNoise));    
  }

  std::cout<<std::endl;
  std::cout<<"Loop on events terminated"<<std::endl;
  
  // correct profiles and fit them with a gaussian
  std::cout<<"Analyze TIB profiles"<<std::endl;
  for(auto iprof : TIBlayers){
     correctProfile(iprof);
     fitProfile(iprof,isGaussian);
  }
  std::cout<<"Analyze TOB profiles"<<std::endl;
  for(auto iprof : TOBlayers){
     correctProfile(iprof);
     fitProfile(iprof,isGaussian);
  }     
  std::cout<<"Analyze TID profiles"<<std::endl;
  for(auto iprof : TIDlayers){
     correctProfile(iprof);
     fitProfile(iprof,isGaussian);
  }
  std::cout<<"Analyze TECPT profiles"<<std::endl;
  for(auto iprof : TECPTlayers){
     correctProfile(iprof);
     fitProfile(iprof,isGaussian);
  }
  std::cout<<"Analyze TECPt profiles"<<std::endl;
  for(auto iprof : TECPtlayers){
     correctProfile(iprof);
     fitProfile(iprof,isGaussian);
  }
  std::cout<<"Analyze TECMT profiles"<<std::endl;
  for(auto iprof : TECMTlayers){
     correctProfile(iprof);
     fitProfile(iprof,isGaussian);
  }
  std::cout<<"Analyze TECMt profiles"<<std::endl;
  for(auto iprof : TECMtlayers){
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
	    std::vector<std::shared_ptr<TProfile> > & TECMtrings){



  std::cout<<"######################################################"<<std::endl;
  std::cout << "Preparing Ring plot for the fifferent subdetectors "<< std::endl; 
  std::cout<<"######################################################"<<std::endl;

  cout<<"tree set branch status "<<endl;
  // TTree Reader appear not to be working with addFriend and EventList
  uint32_t detid;
  float    maxCharge, clCorrectedSignalOverNoise, clSignalOverNoise, clglobalX, clglobalY, clglobalZ, thickness;
  tree->SetBranchStatus("*",kFALSE);
  tree->SetBranchStatus("detid",kTRUE);
  tree->SetBranchStatus("runid",kTRUE);
  tree->SetBranchStatus("eventid",kTRUE);
  tree->SetBranchStatus("trackid0",kTRUE);
  tree->SetBranchStatus("maxCharge",kTRUE);
  tree->SetBranchStatus("clCorrectedSignalOverNoise",kTRUE);
  tree->SetBranchStatus("clSignalOverNoise",kTRUE);
  tree->SetBranchStatus("clglobalX",kTRUE);
  tree->SetBranchStatus("clglobalY",kTRUE);
  tree->SetBranchStatus("clglobalZ",kTRUE);
  tree->SetBranchStatus("thickness",kTRUE);
  tree->SetBranchAddress("detid",&detid);
  tree->SetBranchAddress("maxCharge",&maxCharge);
  tree->SetBranchAddress("clCorrectedSignalOverNoise",&clCorrectedSignalOverNoise);
  tree->SetBranchAddress("clSignalOverNoise",&clSignalOverNoise);
  tree->SetBranchAddress("clglobalX",&clglobalX);
  tree->SetBranchAddress("clglobalY",&clglobalY);
  tree->SetBranchAddress("clglobalZ",&clglobalZ);
  tree->SetBranchAddress("thickness",&thickness);

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

  
  cout<<"create profiles  "<<endl;
  // create vectors
  for(int iring = 0; iring < TIBRing.nDivision; iring++){
    std::shared_ptr<TProfile> temp (new TProfile(Form("TIB_ring%d",iring+1),"",20,-10,10,0,256));
    TIBrings.push_back(temp);
  }
  for(int iring = 0; iring < TOBRing.nDivision; iring++){
    std::shared_ptr<TProfile> temp (new TProfile(Form("TOB_ring%d",iring+1),"",20,-10,10,0,256));
    TOBrings.push_back(temp);
  }
  for(int iring = 0; iring < TIDRing.nDivision; iring++){
    std::shared_ptr<TProfile> temp (new TProfile(Form("TID_ring%d",iring+1),"",20,-10,10,0,256));
    TIDrings.push_back(temp);
  }
  for(int iring = 0; iring < TECRing.nDivision; iring++){
    std::shared_ptr<TProfile> temp (new TProfile(Form("TECP_ring%d_T",iring+1),"",20,-10,10,0,256));
    TECPTrings.push_back(temp);
    temp = std::shared_ptr<TProfile>(new TProfile(Form("TECP_ring%d_t",iring+1),"",20,-10,10,0,256));
    TECPtrings.push_back(temp);
  }
  for(int iring = 0; iring < TECRing.nDivision; iring++){
    std::shared_ptr<TProfile> temp (new TProfile(Form("TECM_ring%d_T",iring+1),"",20,-10,10,0,256));
    TECMTrings.push_back(temp);
    temp = std::shared_ptr<TProfile>(new TProfile(Form("TECM_ring%d_t",iring+1),"",20,-10,10,0,256));
    TECMtrings.push_back(temp);
  }

  cout<<"set entry list "<<endl;
  TEntryList* entryList = tree->GetEntryList();
  
  std::cout<<"Tree with nEntries "<<tree->GetEntries()<<", surviving selections "<<entryList->GetN()<<std::endl;
  if(tree->GetEntries() < entryList->GetN()){
    cout<<"[RPlots] Huston we have a problem with surviving events --> return "<<endl;
    return;
  }

  long int iEntry = 0;
  long int iEvent = 0;
  for( ; iEvent < entryList->GetN()/reductionFactor; iEvent++){    
    iEntry = entryList->GetEntry(iEvent);
    if(iEntry < 0) continue;
    tree->GetEntry(iEntry);
    // take the map delay from the detid
    map->GetEntryWithIndex(detid);
    // take the correction from the detid
    corrections->GetEntryWithIndex(detid);

    cout.flush();
    if(iEvent % 100000 == 0) cout<<"\r"<<"iEvent "<<100*double(iEvent)/(entryList->GetN()/reductionFactor)<<" % ";

    uint32_t subdetid    = int((detid-0x10000000)/0x2000000);
    float    R           = sqrt(clglobalX*clglobalX+clglobalY*clglobalY+clglobalZ*clglobalZ);

    if(subdetid == 3)
      TIBrings.at(int((R-TIBRing.rMin)/((TIBRing.rMax-TIBRing.rMin)/TIBRing.nDivision)))->Fill(delay-correction,maxCharge*(clCorrectedSignalOverNoise)/(clSignalOverNoise));    
    else if(subdetid == 5)
      TOBrings.at(int((R-TOBRing.rMin)/((TOBRing.rMax-TOBRing.rMin)/TOBRing.nDivision)))->Fill(delay-correction,maxCharge*(clCorrectedSignalOverNoise)/(clSignalOverNoise));  
    else if(subdetid == 4)
      TIDrings.at(int((R-TIDRing.rMin)/((TIDRing.rMax-TIDRing.rMin)/TIDRing.nDivision)))->Fill(delay-correction,maxCharge*(clCorrectedSignalOverNoise)/(clSignalOverNoise));    
    else if(subdetid == 6  and clglobalZ > 0 and thickness > 400)
      TECPTrings.at(int((R-TECRing.rMin)/((TECRing.rMax-TECRing.rMin)/TECRing.nDivision)))->Fill(delay-correction,maxCharge*(clCorrectedSignalOverNoise)/(clSignalOverNoise));
    else if(subdetid == 6  and clglobalZ > 0 and thickness < 400)
      TECPtrings.at(int((R-TECRing.rMin)/((TECRing.rMax-TECRing.rMin)/TECRing.nDivision)))->Fill(delay-correction,maxCharge*(clCorrectedSignalOverNoise)/(clSignalOverNoise));
    else if(subdetid == 6  and clglobalZ < 0 and thickness > 400)
      TECMTrings.at(int((R-TECRing.rMin)/((TECRing.rMax-TECRing.rMin)/TECRing.nDivision)))->Fill(delay-correction,maxCharge*(clCorrectedSignalOverNoise)/(clSignalOverNoise));
    else if(subdetid == 6  and clglobalZ < 0 and thickness < 400)
      TECMtrings.at(int((R-TECRing.rMin)/((TECRing.rMax-TECRing.rMin)/TECRing.nDivision)))->Fill(delay-correction,maxCharge*(clCorrectedSignalOverNoise)/(clSignalOverNoise));
  }

  std::cout<<std::endl;
  std::cout<<"Loop on events terminated"<<std::endl;
  
  // correct profiles and fit them with a gaussian
  std::cout<<"Analyze TIB profiles"<<std::endl;
  for(auto iprof : TIBrings){
     correctProfile(iprof);
     fitProfile(iprof,isGaussian);
  }
  std::cout<<"Analyze TOB profiles"<<std::endl;
  for(auto iprof : TOBrings){
     correctProfile(iprof);
     fitProfile(iprof,isGaussian);
  }     
  std::cout<<"Analyze TID profiles"<<std::endl;
  for(auto iprof : TIDrings){
     correctProfile(iprof);
     fitProfile(iprof,isGaussian);
  }
  std::cout<<"Analyze TECPT profiles"<<std::endl;
  for(auto iprof : TECPTrings){
     correctProfile(iprof);
     fitProfile(iprof,isGaussian);
  }
  std::cout<<"Analyze TECPt profiles"<<std::endl;
  for(auto iprof : TECPtrings){
     correctProfile(iprof);
     fitProfile(iprof,isGaussian);
  }
  std::cout<<"Analyze TECMT profiles"<<std::endl;
  for(auto iprof : TECMTrings){
     correctProfile(iprof);
     fitProfile(iprof,isGaussian);
  }
  std::cout<<"Analyze TECMt profiles"<<std::endl;
  for(auto iprof : TECMtrings){
    correctProfile(iprof);
    fitProfile(iprof,isGaussian);
  }  

  std::cout<<"#################################"<<std::endl;
  std::cout<<"##### End of Ring Analysis #####"<<std::endl;
  std::cout<<"#################################"<<std::endl;

}

std::shared_ptr<TCanvas> prepareCanvas(const char* name = ""){
  std::shared_ptr<TCanvas>  c (new TCanvas(name,name,600,625));
  c->cd();
  frame.GetYaxis()->SetTitle("corrected signal (ADC)");
  frame.GetXaxis()->SetTitle("delay (ns)");
  frame.GetXaxis()->SetTitleOffset(1.2);
  frame.GetYaxis()->SetTitleOffset(1.4);
  frame.SetMarkerSize(1.0);
  frame.SetMarkerStyle(20);       
  frame.Draw();
  CMS_lumi(c.get(),"");
  return c;
}

void plotAll(const std::shared_ptr<TCanvas> & canvas, const std::vector<std::shared_ptr<TProfile> > & curves){
  canvas->cd();
  float yMin = 10000.;
  float yMax = -1.;
  int   icolor = 1;
  for(std::vector<std::shared_ptr<TProfile> >::const_iterator it = curves.begin(); it != curves.end(); ++it) {
    if((*it)->Integral() == 0 or (*it)->GetFunction(Form("Gaus_%s",(*it)->GetName())) == 0) continue;
    (*it)->SetLineColor(icolor);
    (*it)->SetMarkerColor(icolor);
    (*it)->SetMarkerStyle(20);
    (*it)->SetMarkerSize(1);
    (*it)->GetFunction(Form("Gaus_%s",(*it)->GetName()))->SetLineColor(icolor);
    (*it)->GetFunction(Form("Gaus_%s",(*it)->GetName()))->SetLineWidth(2);
    (*it)->Draw("same");
    if((*it)->GetMaximum() > yMax)
      yMax = (*it)->GetMaximum();
    if((*it)->GetMinimum() < yMin)
      yMin = (*it)->GetMinimum();
    icolor++;
  }
  frame.GetYaxis()->SetRangeUser(yMin*0.75,yMax*1.25);
}

void plotMaxima(const std::shared_ptr<TCanvas> & canvas, const std::vector<std::shared_ptr<TProfile> > & curves, string outputDIR, string postfix){


  canvas->cd();
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->SetBottomMargin(0.21);
  
  std::shared_ptr<TH1F> graph (new TH1F(Form("graph_%s",postfix.c_str()),"",curves.size(),0,curves.size()+1));
  int i = 0;
  for(std::vector<std::shared_ptr<TProfile> >::const_iterator it = curves.begin(); it != curves.end(); ++it,++i) {
    if((*it)->GetListOfFunctions()->At(0)) {
      graph->SetBinContent(i+1,((TF1*)(*it)->GetListOfFunctions()->At(0))->GetMaximumX());
      graph->SetBinError(i+1,((TF1*)(*it)->GetListOfFunctions()->At(0))->GetParError(1) );
    }
  }

  
  i=0;
  for(std::vector<std::shared_ptr<TProfile> >::const_iterator it = curves.begin(); it < curves.end(); ++it,++i){
    graph->GetXaxis()->SetBinLabel(i+1,(*it)->GetName());
  }
  graph->GetYaxis()->SetTitle("delay (ns)");
  graph->GetXaxis()->LabelsOption("v");
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(1);
  graph->SetMarkerColor(kBlack);
  graph->SetLineColor(kBlack);
  graph->GetYaxis()->SetRangeUser(-10,10);
  graph->Draw();
  CMS_lumi(canvas.get(),"");
  graph->SetFillColor(kGreen+1);
  graph->Draw("E2same");
  std::auto_ptr<TF1> funz (new TF1("funz","0",0,graph->GetBinLowEdge(graph->GetNbinsX()+1)));
  funz->SetLineColor(kRed);
  funz->SetLineWidth(2);
  funz->Draw("same");
  graph->Draw("EPsame");
  canvas->RedrawAxis("sameaxis");
  canvas->Print(Form("%s/layersGraph_%s.root",outputDIR.c_str(),postfix.c_str()));
  
  return;
}

void delayValidation(string file0, 
		     string file1, 
		     string file2 = "nocorrection.root", 
		     bool isBOn          = true, 
		     bool plotPartitions = true, 
		     bool plotLayer      = true, 
		     bool plotSlices     = false, 
		     bool plotModules    = false,
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
  std::auto_ptr<TFile> _file2 (TFile::Open(file2.c_str()));
  std::shared_ptr<TTree> clusters   ((TTree*)_file0->Get("analysis/trackerDPG/clusters"));
  std::shared_ptr<TTree> readoutMap ((TTree*)_file1->Get("analysis/trackerDPG/readoutMap"));
  std::shared_ptr<TTree> delayCorrections ((TTree*)_file2->Get("delayCorrections"));
  // add tracks and event tree also as friends to apply the event selection by string
  clusters->AddFriend((TTree*)_file0->Get("analysis/trackerDPG/events"));
  clusters->AddFriend((TTree*)_file0->Get("analysis/trackerDPG/tracks0"));

  // apply common preselection cuts on events, track and cluster quantities
  std::cout<<"Select Events"<<std::endl;
  string eventSelection;
  string trackSelection;
  string clusterSelection;
  if(isBOn){
    eventSelection   = "L1TechnicalBits[0] && L1TechnicalBits[41] && !(L1TechnicalBits[36]||L1TechnicalBits[37]||L1TechnicalBits[38]||L1TechnicalBits[39])&&lowPixelProbabilityFraction[0]<0.4 && lowPixelProbabilityFraction[0]>-0.5 && nVertices>0";
    trackSelection   = "p>1 && p<10 && quality>2 && pterr/pt<0.1 && pt>1 && dedx1<10 && quality>2 && dedx1<3.5";
    clusterSelection = "onTrack && angle>0 && maxCharge<254";
  }
  else{
    eventSelection   = "nVertices>=0";
    trackSelection   = "pt>= 0";
    clusterSelection = "onTrack && angle>0 && maxCharge<254";     
  }
  
  clusters->SetEventList(0);  
  // select events -> make a single event list  
  TEventList* selectedEvents = new TEventList("selectedEvents","selectedEvents");
  std::cout << "Applying common event selection: "<<eventSelection+" && "+trackSelection+" && "+clusterSelection<<std::endl;
  clusters->Draw(">> selectedEvents",(eventSelection+" && "+trackSelection+" && "+clusterSelection).c_str());
  selectedEvents->Print();
  clusters->SetEventList(selectedEvents);  

  // create the plots per layer
  if(plotLayer){

    std::vector<std::shared_ptr<TProfile> > TIBlayers;
    std::vector<std::shared_ptr<TProfile> > TOBlayers;
    std::vector<std::shared_ptr<TProfile> > TIDlayers;
    std::vector<std::shared_ptr<TProfile> > TECPTlayers;
    std::vector<std::shared_ptr<TProfile> > TECPtlayers;
    std::vector<std::shared_ptr<TProfile> > TECMTlayers;
    std::vector<std::shared_ptr<TProfile> > TECMtlayers;
    LayerPlots(clusters,readoutMap,delayCorrections,TIBlayers,TOBlayers,TIDlayers,TECPTlayers,TECPtlayers,TECMTlayers,TECMtlayers);
  
    std::shared_ptr<TCanvas> c1 = prepareCanvas("TIB_layers");
    plotAll(c1,TIBlayers);
    c1->Print(Form("%s/TIB_layers.root",outputDIR.c_str()));

    std::shared_ptr<TCanvas> c2 = prepareCanvas("TID_layers");
    plotAll(c2,TIDlayers);
    c2->Print(Form("%s/TID_layers.root",outputDIR.c_str()));

    std::shared_ptr<TCanvas> c3 = prepareCanvas("TOB_layers");
    plotAll(c3,TOBlayers);
    c3->Print(Form("%s/TOB_layers.root",outputDIR.c_str()));

    std::shared_ptr<TCanvas> c4 = prepareCanvas("TECPt_layers");
    plotAll(c4,TECPtlayers);
    c4->Print(Form("%s/TECPt_layers.root",outputDIR.c_str()));

    std::shared_ptr<TCanvas> c5 = prepareCanvas("TECPT_layers");
    plotAll(c5,TECPTlayers);
    c5->Print(Form("%s/TECPT_layers.root",outputDIR.c_str()));

    std::shared_ptr<TCanvas> c6 = prepareCanvas("TECMt_layers");
    plotAll(c6,TECMtlayers);
    c6->Print(Form("%s/TECMt_layers.root",outputDIR.c_str()));

    std::shared_ptr<TCanvas> c7 = prepareCanvas("TECMT_layers");
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

    RPlots(clusters,readoutMap,delayCorrections,TIBrs,TOBrs,TIDrs,TECPTrs,TECPtrs,TECMTrs,TECMtrs);
    
    std::shared_ptr<TCanvas> c1b = prepareCanvas("TIB_distance");
    plotAll(c1b,TIBrs);
    c1b->Print(Form("%s/TIB_distance.root",outputDIR.c_str()));

    std::shared_ptr<TCanvas> c2b = prepareCanvas("TID_distance");
    plotAll(c2b,TIDrs);
    c2b->Print(Form("%s/TID_distance.root",outputDIR.c_str()));
    
    std::shared_ptr<TCanvas> c3b = prepareCanvas("TOB_distance");
    plotAll(c3b,TOBrs);
    c3b->Print(Form("%s/TOB_distance.root",outputDIR.c_str()));

    std::shared_ptr<TCanvas> c4b = prepareCanvas("TECPt_distance");
    plotAll(c4b,TECPtrs);
    c4b->Print(Form("%s/TECPt_distance.root",outputDIR.c_str()));

    std::shared_ptr<TCanvas> c5b = prepareCanvas("TECPT_distance");
    plotAll(c5b,TECPTrs);
    c5b->Print(Form("%s/TECPT_distance.root",outputDIR.c_str()));

    std::shared_ptr<TCanvas> c6b = prepareCanvas("TECMt_distance");
    plotAll(c6b,TECMtrs);
    c6b->Print(Form("%s/TECMt_distance.root",outputDIR.c_str()));

    std::shared_ptr<TCanvas> c7b = prepareCanvas("TECMT_distance");
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
    
    RPlots(clusters,readoutMap,delayCorrections,TIB,TOB,TID,TECPT,TECPt,TECMT,TECMt);
        
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

    std::shared_ptr<TCanvas> c1 = prepareCanvas("Partitions");
    plotAll(c1,allPartitions);
    c1->Print(Form("%s/Partitions.root",outputDIR.c_str()));
  }

  if(plotModules){ // useful for the new automatic random delay scan

    // top be developed --> one profile for each detid
  }  
}

