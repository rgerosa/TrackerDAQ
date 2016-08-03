#include "TH1F.h"
#include "TF1.h"
#include "TProfile.h"
#include "TString.h"
#include "TkPulseShape.h"
#include "CMS_lumi.h"
#include "TFitResultPtr.h"
#include "TLegend.h"

using namespace std;

static int savePlotEvery = 10;

// basic profile for maxCharge                                                                                                                                                 
TH1F*     frame;
TLegend*  legend;

/// define limit and binning for the different observables
void setLimitsAndBinning(const string & observable, float & xMin, float & xMax, int & nBin){
  if(observable == "maxCharge"){
    xMin = 10;
    xMax = 254;
    nBin = 40;
  }
  else if(observable == "clSignalOverNoise" or observable == "clCorrectedSignalOverNoise"){
    xMin = 5;
    xMax = 70;
    nBin = 35;
  }
  else if(observable == "delay"){
    xMin = -10*1.04;
    xMax = 10*1.04;
    nBin = 20;
  }
  else{
    xMin = 20;
    xMax = 250;
    nBin = 150;
  } 

  return;
}

void setLimitsAndBinning(const string & observable, vector<double> & limits){

  if(observable == "delay")
    limits = {-10.94,-9.9,-8.86,-7.82,-6.78,-5.74,-4.70,-3.66,-2.62,-1.58,-0.54,0.54,1.58,2.62,3.66,4.70,5.74,6.78,7.82,8.86,9.9,10.94};
  
  return;
}

int getFilledBins(const std::shared_ptr<TProfile> & prof){
  int nfilled = 0;
  for(int iBin = 1; iBin <= prof->GetNbinsX(); iBin++){
    if(prof->GetBinContent(iBin) != 0) nfilled++;
  }
  return nfilled;
}

int getFilledBins(const std::shared_ptr<TH1F> & prof){
  int nfilled = 0;
  for(int iBin = 1; iBin <= prof->GetNbinsX(); iBin++){
    if(prof->GetBinContent(iBin) != 0) nfilled++;
  }
  return nfilled;
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

void correctHistogram(const std::shared_ptr<TH1F> & histo){
  int nbins = histo->GetNbinsX();
  float min = limit(3.); // zero suppression level .. something like S/Noise = 3                                                                                               
  for(int bin=1;bin<=nbins;++bin){
    if(histo->GetBinContent(bin)<min) { // set to zero                                                                                                                   
      histo->SetBinContent(bin,0.);
      histo->SetBinError(bin,0.);
    }
    else{
      histo->SetBinContent(bin,correctMeasurement(histo->GetBinContent(bin),3.)); // correct the measurement                                
    }
  }

  return;
}

//// fitting profiles
TFitResultPtr fitProfile(const std::shared_ptr<TProfile> & prof, bool gaus = false, string options = "", bool verbosity = false){
  std::shared_ptr<TF1> pulse;
  if(gaus) { // gaussina fit                                                                                                                                                   
    pulse = std::shared_ptr<TF1>(new TF1(Form("Gaus_%s",prof->GetName()),"gaus(0)",-10.5,10.5));
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
    std::cout << "Profile Name "<<prof->GetName()<<" "<<"Maximum at " << pulse->GetParameter(1) << std::endl;
  return result;
}

//// fitting profiles
TFitResultPtr fitHistogram(const std::shared_ptr<TH1F> & histo, bool gaus = false, string options = "", bool verbosity = false){
  std::shared_ptr<TF1> pulse;
  if(gaus) { // gaussina fit                                                                                                                                                   
    pulse = std::shared_ptr<TF1>(new TF1(Form("Gaus_%s",histo->GetName()),"gaus(0)",-10.5,10.5));
    pulse->SetParameters(50,0,12);
  } else {// different fit for the histoile                                                                                                                                     
    pulse = std::shared_ptr<TF1>(TkPulseShape::GetDeconvFitter());
    pulse->SetName(Form("SignalFit_%s",histo->GetName()));
    pulse->SetParameters(0,0,3,50,15);
    pulse->FixParameter(0,0);
    pulse->FixParameter(3,50);
  }
  gROOT->SetBatch(1);
  TFitResultPtr result = histo->Fit(pulse.get(),(options+"S").c_str());
  if(verbosity)
    std::cout << "Histogram Name "<<histo->GetName()<<" "<<"Maximum at " << pulse->GetParameter(1) << std::endl;
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
    frame = new TH1F("frame","",xBins.size()-1,&xBins[0]);
    frame->GetYaxis()->SetRangeUser(yMin,yMax);
  }

  if(observable == "clSignalOverNoise")
    frame->GetYaxis()->SetTitle("S/N");
  else if (observable == "clCorrectedSignalOverNoise")
    frame->GetYaxis()->SetTitle("corrected S/N");
  else
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
void plotAll(const std::shared_ptr<TCanvas> & canvas, const std::shared_ptr<TProfile>  & curves){

  canvas->cd();
  curves->SetLineColor(kBlack);
  curves->SetMarkerColor(kBlack);
  curves->SetMarkerStyle(20);
  curves->SetMarkerSize(1);
  curves->GetFunction(Form("Gaus_%s",curves->GetName()))->SetLineColor(kRed);
  curves->GetFunction(Form("Gaus_%s",curves->GetName()))->SetLineWidth(2);
  curves->Draw("same");
  frame->GetYaxis()->SetRangeUser(curves->GetFunction(Form("Gaus_%s",curves->GetName()))->GetMinimum()*0.75,curves->GetFunction(Form("Gaus_%s",curves->GetName()))->GetMaximum()*1.25);
  return;
}

void plotAll(const std::shared_ptr<TCanvas> & canvas, const std::shared_ptr<TH1F>  & curves){

  canvas->cd();
  curves->SetLineColor(kBlack);
  curves->SetMarkerColor(kBlack);
  curves->SetMarkerStyle(20);
  curves->SetMarkerSize(1);
  curves->GetFunction(Form("Gaus_%s",curves->GetName()))->SetLineColor(kRed);
  curves->GetFunction(Form("Gaus_%s",curves->GetName()))->SetLineWidth(2);
  curves->Draw("same");
  frame->GetYaxis()->SetRangeUser(curves->GetFunction(Form("Gaus_%s",curves->GetName()))->GetMinimum()*0.75,curves->GetFunction(Form("Gaus_%s",curves->GetName()))->GetMaximum()*1.25);
  return;
}


// plot all the profiles on a canvas
void plotAll(const std::shared_ptr<TCanvas> & canvas, const std::vector<std::shared_ptr<TProfile> > & curves, string  remove  = ""){

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
    TString legenEntry = Form("%s",(*it)->GetName());
    legenEntry.ReplaceAll("_"," ").ReplaceAll(remove.c_str(),"");
    legend->AddEntry((*it).get(),legenEntry,"EP");
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

// plot all the profiles on a canvas
void plotAll(const std::shared_ptr<TCanvas> & canvas, const std::vector<std::shared_ptr<TH1F> > & curves, string  remove  = ""){

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

  for(std::vector<std::shared_ptr<TH1F> >::const_iterator it = curves.begin(); it != curves.end(); ++it) {
    if((*it)->Integral() == 0 or (*it)->GetFunction(Form("Gaus_%s",(*it)->GetName())) == 0) continue;
    (*it)->SetLineColor(icolor);
    (*it)->SetMarkerColor(icolor);
    (*it)->SetMarkerStyle(20);
    (*it)->SetMarkerSize(1);
    TString legenEntry = Form("%s",(*it)->GetName());
    legenEntry.ReplaceAll("_"," ").ReplaceAll(remove.c_str(),"");
    legend->AddEntry((*it).get(),legenEntry,"EP");
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
      graph->SetBinContent(i+1,((TF1*)(*it)->GetListOfFunctions()->At(0))->GetParameter(1));
      graph->SetBinError(i+1,((TF1*)(*it)->GetListOfFunctions()->At(0))->GetParError(1) );
    }
  }

  i=0;
  for(std::vector<std::shared_ptr<TProfile> >::const_iterator it = curves.begin(); it < curves.end(); ++it,++i){
    TString label = Form("%s",(*it)->GetName());
    label.ReplaceAll("_"," ");
    graph->GetXaxis()->SetBinLabel(i+1,label);
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


void plotMaxima(const std::shared_ptr<TCanvas> & canvas, const std::vector<std::shared_ptr<TH1F> > & curves, string outputDIR, string postfix){


  canvas->cd();
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->SetBottomMargin(0.21);

  std::shared_ptr<TH1F> graph (new TH1F(Form("graph_%s",postfix.c_str()),"",curves.size(),0,curves.size()+1));
  int i = 0;
  for(std::vector<std::shared_ptr<TH1F> >::const_iterator it = curves.begin(); it != curves.end(); ++it,++i) {
    if((*it)->GetListOfFunctions()->At(0)) {
      graph->SetBinContent(i+1,((TF1*)(*it)->GetListOfFunctions()->At(0))->GetParameter(1));
      graph->SetBinError(i+1,((TF1*)(*it)->GetListOfFunctions()->At(0))->GetParError(1) );
    }
  }


  i=0;
  for(std::vector<std::shared_ptr<TH1F> >::const_iterator it = curves.begin(); it < curves.end(); ++it,++i){
    TString label = Form("%s",(*it)->GetName());
    label.ReplaceAll("_"," ");
    graph->GetXaxis()->SetBinLabel(i+1,label);
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

void saveAll(const std::map<uint32_t,std::shared_ptr<TProfile> > channelMap, const string & outputDIR, const string & observable, const string & postfix){

  std::shared_ptr<TFile> outputFile (new TFile((outputDIR+"/channelProfileMap_"+postfix+".root").c_str(),"RECREATE"));
  outputFile->cd();

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



void saveAll(const std::map<uint32_t,std::shared_ptr<TH1F> > channelHistoMap, const string & outputDIR, const string & observable, const string & postfix){

  std::shared_ptr<TFile> outputFile (new TFile((outputDIR+"/channelHistoMap_"+postfix+".root").c_str(),"RECREATE"));
  outputFile->cd();
  cout<<"### saveAll channel distribution for each delay "<<endl;
  long int idetid = 0;
  for(auto itMap : channelHistoMap){
    std::shared_ptr<TCanvas> c1 = prepareCanvas("channelMap",observable);
    c1->cd();
    if(idetid % savePlotEvery == 0){ 
      cout.flush();
      cout<<"\r"<<"idetid "<<100*double(idetid)/channelHistoMap.size()<<" % ";
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


      c1->Write(Form("histogram_detid_%s",to_string(itMap.first).c_str()));
    }
    idetid++;
  }
  std::cout<<std::endl;
  outputFile->Close();  
  return;
  
  
} 



Double_t langaufun(Double_t *x, Double_t *par) {

  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation),
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location

  // Control constants
  Double_t np = 200.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;


  // MP shift correction
  mpc = par[1] - mpshift * par[0];

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp-xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);

    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}
