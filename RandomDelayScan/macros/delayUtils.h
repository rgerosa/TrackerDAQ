#include "TH1F.h"
#include "TF1.h"
#include "TProfile.h"
#include "TString.h"
#include "TkPulseShape.h"
#include "CMS_lumi.h"

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

// prepare a canvas for the final plot
std::shared_ptr<TCanvas> prepareCanvas(const string & name = "",const string & observable = "maxCharge"){

  std::shared_ptr<TCanvas>  c (new TCanvas(name.c_str(),name.c_str(),600,625));
  c->cd();

  float yMin = 0.; 
  float yMax = 0.;
  int   nBinsY = 0.;

  float xMin = 0.; 
  float xMax = 0.;
  int   nBinsX = 0.;

  if(frame == 0 or frame == NULL){
    setLimitsAndBinning(observable,yMin,yMax,nBinsY);
    setLimitsAndBinning("delay",xMin,xMax,nBinsX);
    frame = new TProfile("frame","",nBinsX,xMin,xMax,yMin,yMax);
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
