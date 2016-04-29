#include "TH1F.h"
#include "TF1.h"
#include "TProfile.h"
#include "TString.h"
#include "TkPulseShape.h"
#include "CMS_lumi.h"
#include "TFitResultPtr.h"

static int savePlotEvery = 10;

// basic profile for maxCharge                                                                                                                                                 
TProfile* frame;
TLegend*  legend;

/// define limit and binning for the different observables
void setLimitsAndBinning(const string & observable, float & xMin, float & xMax, int & nBin){
  if(observable == "maxCharge"){
    xMin = 10;
    xMax = 254;
    nBin = 40;
  }
  else if(observable == "clSignalOverNoise" or observable == "clCorrectedSignalOverNoise"){
    xMin = 0;
    xMax = 80;
    nBin = 40;
  }
  else if(observable == "delay"){
    xMin = -10*1.04;
    xMax = 10*1.04;
    nBin = 20;
  }
  else{
    xMin = 20;
    xMax = 250;
    nBin = 40;
  } 

  return;
}

void setLimitsAndBinning(const string & observable, vector<double> & limits){

  if(observable == "delay")
    limits = {-10.94,-9.9,-8.86,-7.82,-6.78,-5.74,-4.70,-3.66,-2.62,-1.58,-0.54,0.54,1.58,2.62,3.66,4.70,5.74,6.78,7.82,8.86,9.9,10.94};
  
  return;
}


// struct to handle ring definition for the different tracker partitions
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

/// Signal over Noise correction                                                                                                                                       
float limit(float SoNcut){
  return 3.814567e+00+8.336601e+00*SoNcut-1.511334e-01*pow(SoNcut,2);
}

// correct the measurement accordint to the SoN cut
float correctMeasurement(float mean, float SoNcut){
  if(mean>limit(SoNcut))
    return -8.124872e+00+9.860108e-01*mean-3.618158e-03*pow(mean,2)+2.037263e-05*pow(mean,3);
  else return 0.;
}

// Profile correction
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

  return;
}

//// fitting profiles
TFitResultPtr fitProfile(const std::shared_ptr<TProfile> & prof, bool gaus = false, string options = "", bool verbosity = false){
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
  TFitResultPtr result = prof->Fit(pulse.get(),(options+"S").c_str());
  if(verbosity)
    std::cout << "Profile Name "<<prof->GetName()<<" "<<"Maximum at " << pulse->GetMaximumX(-10,10) << std::endl;
  return result;
}

// prepare a canvas for the final plot
std::shared_ptr<TCanvas> prepareCanvas(const string & name = "",const string & observable = "maxCharge"){

  std::shared_ptr<TCanvas>  c (new TCanvas(name.c_str(),name.c_str(),600,625));
  c->cd();

  float yMin = 0.; 
  float yMax = 0.;
  int   nBinsY = 0.;

  vector<double> xBins;

  if(frame == 0 or frame == NULL){
    setLimitsAndBinning(observable,yMin,yMax,nBinsY);
    setLimitsAndBinning("delay",xBins);
    frame = new TProfile("frame","",xBins.size()-1,&xBins[0],yMin,yMax);
  }

  frame->GetYaxis()->SetTitle("corrected signal (ADC)");
  frame->GetXaxis()->SetTitle("delay (ns)");
  frame->GetXaxis()->SetTitleOffset(1.1);
  frame->GetYaxis()->SetTitleOffset(1.2);
  frame->SetMarkerSize(1.0);
  frame->SetMarkerStyle(20);
  frame->Draw();
  CMS_lumi(c.get(),"");
  return c;
}

// plot all the profiles on a canvas
void plotAll(const std::shared_ptr<TCanvas> & canvas, const std::vector<std::shared_ptr<TProfile> > & curves){

  canvas->cd();
  float yMin = 10000.;
  float yMax = -1.;
  int   icolor = 1;

  if(legend == 0 or legend == NULL)
    legend = new TLegend(0.55,0.7,0.85,0.92);
  else
    legend->Clear();

  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);  

  for(std::vector<std::shared_ptr<TProfile> >::const_iterator it = curves.begin(); it != curves.end(); ++it) {
    if((*it)->Integral() == 0 or (*it)->GetFunction(Form("Gaus_%s",(*it)->GetName())) == 0) continue;
    (*it)->SetLineColor(icolor);
    (*it)->SetMarkerColor(icolor);
    (*it)->SetMarkerStyle(20);
    (*it)->SetMarkerSize(1);
    legend->AddEntry((*it).get(),(*it)->GetName(),"EP");
    (*it)->GetFunction(Form("Gaus_%s",(*it)->GetName()))->SetLineColor(icolor);
    (*it)->GetFunction(Form("Gaus_%s",(*it)->GetName()))->SetLineWidth(2);
    (*it)->Draw("same");
    if((*it)->GetFunction(Form("Gaus_%s",(*it)->GetName()))->GetMaximum() > yMax)
      yMax = (*it)->GetFunction(Form("Gaus_%s",(*it)->GetName()))->GetMaximum();
    if((*it)->GetFunction(Form("Gaus_%s",(*it)->GetName()))->GetMinimum() < yMin)
      yMin = (*it)->GetFunction(Form("Gaus_%s",(*it)->GetName()))->GetMinimum();
    icolor++;
  }
  frame->GetYaxis()->SetRangeUser(yMin*0.75,yMax*1.50);
  legend->Draw("same");

  return;

}

// plot delay correspoding to the maximum of the profiles
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
  graph->GetYaxis()->SetRangeUser(-5,5);
  graph->Draw();
  CMS_lumi(canvas.get(),"");
  graph->SetFillColor(kGreen+1);
  std::auto_ptr<TF1> funz (new TF1("funz","0",0,graph->GetBinLowEdge(graph->GetNbinsX()+1)));
  funz->SetLineColor(kRed);
  funz->SetLineWidth(2);
  graph->Draw("HISTsame");
  funz->Draw("same");
  graph->Draw("EPsame");
  canvas->RedrawAxis("sameaxis");
  canvas->Print(Form("%s/layersGraph_%s.root",outputDIR.c_str(),postfix.c_str()));
  
  return;
}

void saveAll(const std::map<uint32_t,std::shared_ptr<TProfile> > channelMap, const string & outputDIR, const string & observable){

  std::shared_ptr<TFile> outputFile (new TFile((outputDIR+"/channelMap.root").c_str(),"RECREATE"));
  outputFile->cd();
  float yMin = 10000.;
  float yMax = -1.;

  cout<<"### saveAll channel distribution and fits "<<endl;
  long int imap = 0;
  for(auto itMap : channelMap){
    std::shared_ptr<TCanvas> c1 = prepareCanvas("channelMap",observable);
    c1->cd();
    if(imap % savePlotEvery == 0){ // save one every ten
      cout.flush();
      cout<<"\r"<<"iChannel "<<100*double(imap)/(channelMap.size())<<" % ";
      if(itMap.second->Integral() == 0 or itMap.second->GetFunction(Form("Gaus_%s",itMap.second->GetName())) == 0) continue;
      itMap.second->SetLineColor(kBlack);
      itMap.second->SetMarkerColor(kBlack);
      itMap.second->SetMarkerStyle(20);
      itMap.second->SetMarkerSize(1);
      itMap.second->GetFunction(Form("Gaus_%s",itMap.second->GetName()))->SetLineColor(kRed);
      itMap.second->GetFunction(Form("Gaus_%s",itMap.second->GetName()))->SetLineWidth(2);
      itMap.second->Draw("same");
      frame->GetYaxis()->SetRangeUser(itMap.second->GetFunction(Form("Gaus_%s",itMap.second->GetName()))->GetMinimum()*0.75,
				      itMap.second->GetFunction(Form("Gaus_%s",itMap.second->GetName()))->GetMaximum()*1.25);
      
      c1->Write(Form("profile_detid_%s",to_string(itMap.first).c_str()));
    }
    imap++;    
  }
  std::cout<<std::endl;
  outputFile->Close();  
  return;
  
  
} 
