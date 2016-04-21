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

//#include "tdrstyle.C"
#include "TkPulseShape.h"
 
static unsigned int id = 0;
static TProfile frame("frame","frame",40,-10,10,0,256);

float limit(float SoNcut)
{
  return 3.814567e+00+8.336601e+00*SoNcut-1.511334e-01*pow(SoNcut,2);
}

float correctMeasurement(float mean, float SoNcut)
{
  if(mean>limit(SoNcut))
      return -8.124872e+00+9.860108e-01*mean-3.618158e-03*pow(mean,2)+2.037263e-05*pow(mean,3);
  else return 0.;
}

void correctProfile(TProfile* profile)
{
  int nbins=profile->GetNbinsX();
  float min = limit(3.);
  for(int bin=1;bin<=nbins;++bin)
    if(profile->GetBinContent(bin)<min) {
      profile->SetBinContent(bin,0.);
      profile->SetBinError(bin,0.);
      profile->SetBinEntries(bin,0);
    }
    else {
      profile->SetBinContent(bin,
                             profile->GetBinEntries(bin)*
                             correctMeasurement(profile->GetBinContent(bin),3.));
    }
}

TF1* fitProfile(TProfile* prof, bool gaus=false)
{ 
   TF1* pulse = 0;
   if(gaus) {
     pulse = new TF1("fg","gaus(0)",-10,10);
     pulse->SetParameters(50,0,12);
   } else {
     pulse = TkPulseShape::GetDeconvFitter();
     pulse->SetName(Form("SignalFit_%s",prof->GetName()));
     pulse->SetParameters(0,0,3,50,15);
     pulse->FixParameter(0,0);
     pulse->FixParameter(3,50);
   }
   gROOT->SetBatch(1);
   prof->Fit(pulse,"");
   std::cout << "Maximum at " << pulse->GetMaximumX(-10,10) << std::endl;
   return pulse;
}

std::vector<TProfile*> LayerPlots(TTree* tree, unsigned int subdet, unsigned int nlayers, const char* extracut = "1") 
{
   std::cout << "Preparing plots for subdet " << subdet << " with extra cuts " << extracut << std::endl; 
   std::vector<TProfile*> output;
   for(unsigned int i=1;i<=nlayers;++i) {
     std::cout << "Processing layer " << i << std::endl;
     TProfile* prof = new TProfile(Form("Signal_s%d_l%d_id%d",subdet,i,id),Form("Signal_s%d_l%d_id%d",subdet,i,id),40,-10,10,0,256);
     output.push_back(prof);
     tree->Draw(Form("maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise:delay-correction>>Signal_s%d_l%d_id%d",subdet,i,id),
                Form("subdetid==%d && layer==%d && %s",subdet,i,extracut),
                "prof, goff",50000000);
     TFile file("tmp.root","CREATE");
     file.cd();
     prof->Write();
     correctProfile(prof);
     prof->Write("after_correction");
     fitProfile(prof,true);
     prof->Write("withFit");
     id++;
     file.Close();
   }
   return output;
}

std::vector<TProfile*> RPlots(TTree* tree, unsigned int subdet, unsigned int nslices = 10., float min = 0., float max = 300., const char* extracut = "1") 
{
   std::cout << "Preparing plots for subdet " << subdet << std::endl; 
   std::vector<TProfile*> output;
   float locmin = min;
   float step = (max-min)/nslices;
   float locmax = min+step;
   for(unsigned int i=1;i<=nslices;++i) {
     std::cout << "Processing slice [ " << locmin << " , " << locmax << " ]" << std::endl;
     TProfile* prof = new TProfile(Form("Signal_s%d_sl%d_id%d",subdet,i,id),Form("Signal_s%d_sl%d_id%d",subdet,i,id),40,-10,10,0,256);
     output.push_back(prof);
     tree->Draw(Form("maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise:delay-correction>>Signal_s%d_sl%d_id%d",subdet,i,id),
                Form("subdetid==%d && R>%f && R<%f && %s",subdet,locmin,locmax, extracut),
                "prof, goff",50000000);
     correctProfile(prof);
     fitProfile(prof,true);
     id++;
     locmin += step;
     locmax += step;
   }
   return output;
}

TCanvas* prepareCanvas(const char* name = "")
{
   TCanvas* c = new TCanvas(name,name);
   TProfile* p = (TProfile*)frame.DrawClone();
   c->SetGridx();
   p->GetYaxis()->SetRangeUser(0,100);
   p->GetYaxis()->SetTitle("corrected signal (ADC)");
   p->GetXaxis()->SetTitle("delay (ns)");
   p->GetXaxis()->SetTitleOffset(1.2);
   //TLatex* cmsPreliminary = new TLatex(-6., 25.,Form("CMS preliminary 2010 - 7TeV \\ %s",name));
   //cmsPreliminary->SetTextSize(20);
   //cmsPreliminary->Draw();
   return c;
}

void plotAll(TCanvas* canvas, std::vector<TProfile*> curves, unsigned int color)
{
   canvas->cd();
   for(std::vector<TProfile*>::const_iterator it = curves.begin(); it < curves.end(); ++it) {
     (*it)->SetLineColor(color);
     (*it)->SetMarkerColor(color);
     (*it)->Draw("same");
   }
}

TCanvas* plotMaxima(std::vector<TProfile*> curves)
{
   TGraphErrors* graph = new TGraphErrors;
   int i=0;
   int j=0;
   for(std::vector<TProfile*>::const_iterator it = curves.begin(); it < curves.end(); ++it,++i) {
     if((*it)->GetListOfFunctions()->At(0)) {
       graph->SetPoint(j,i+0.5,((TF1*)(*it)->GetListOfFunctions()->At(0))->GetMaximumX());
       graph->SetPointError(j, 0.,((TF1*)(*it)->GetListOfFunctions()->At(0))->GetParError(1) );
       ++j;
     }
   }
   TCanvas* c1 = new TCanvas;
   graph->Draw("AP");
   c1->SetGridy();
   new TCanvas;
   i=0;
   for(std::vector<TProfile*>::const_iterator it = curves.begin(); it < curves.end(); ++it,++i) {
     graph->GetXaxis()->SetBinLabel(graph->GetXaxis()->FindBin(i+0.5),(*it)->GetTitle());
   }
   graph->GetYaxis()->SetTitle("delay (ns)");
   graph->GetXaxis()->SetTitle("slice");

   return c1;
}

void delayValidation(const char* file0, const char* file1, const char* file2 = "nocorrection.root", string postfix = "prompt")
{
  // prepare style and load macros
  //setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gROOT->SetBatch(kTRUE);

  // open the file and prepare the cluster tree
  TFile *_file0 = TFile::Open(file0);
  TFile *_file1 = TFile::Open(file1);
  TFile *_file2 = TFile::Open(file2);
  TTree* clusters = (TTree*)_file0->Get("analysis/trackerDPG/clusters");
  clusters->AddFriend((TTree*)_file1->Get("analysis/trackerDPG/readoutMap"));
  clusters->AddFriend((TTree*)_file2->Get("delayCorrections"));
  clusters->AddFriend((TTree*)_file0->Get("analysis/trackerDPG/events"));
  clusters->AddFriend((TTree*)_file0->Get("analysis/trackerDPG/tracks0"));
  // setup the aliases
  // selection will be based on :
  //   - subdetector
  //   - layer
  //   - distance to IP
  clusters->SetAlias("subdetid","int((detid-0x10000000)/0x2000000)");
  clusters->SetAlias("barrellayer","int((detid%33554432)/0x4000)");
  clusters->SetAlias("TIDlayer","int((detid%33554432)/0x800)%4");
  clusters->SetAlias("TECPlayer","int((detid%33554432)/0x4000)-32");
  clusters->SetAlias("TECMlayer","int((detid%33554432)/0x4000)-16");
  clusters->SetAlias("R","sqrt(clglobalX**2+clglobalY**2+clglobalZ**2)");

  // apply common preselection cuts on events, track and cluster quantities
  //TCut eventSelection = "L1TechnicalBits[0]&&L1TechnicalBits[41]&&!(L1TechnicalBits[36]||L1TechnicalBits[37]||L1TechnicalBits[38]||L1TechnicalBits[39])&&lowPixelProbabilityFraction[0]<0.4 && lowPixelProbabilityFraction[0]>-0.5";//&& nVertices>0";
  TCut eventSelection = "";//&& nVertices>0";
  //TCut trackSelection = "p>1 && p<10 && quality>2 && pterr/pt<0.1 && pt>1 && dedx1<10 && quality>2 && dedx1<3.5";
  TCut trackSelection = "pt>=0";
  TCut clusterSelection = "onTrack && angle>0 && maxCharge<254";
  clusters->SetEventList(0);
  TEventList* selectedEventsA = new TEventList("selectedEventsA","selectedEventsA");
  TEventList* selectedEventsB = new TEventList("selectedEventsB","selectedEventsB");
  TEventList* selectedEventsC = new TEventList("selectedEventsC","selectedEventsC");
  std::cout << "Applying common event selection..." << std::endl;
  clusters->Draw(">>selectedEventsA",eventSelection,"");
  selectedEventsA->Print();
  clusters->SetEventList(selectedEventsA);
  clusters->Draw(">>selectedEventsB",trackSelection,"");
  selectedEventsB->Print();
  clusters->SetEventList(selectedEventsB);
  clusters->Draw(">>selectedEventsC", clusterSelection,"");
  selectedEventsC->Print();
  clusters->SetEventList(selectedEventsC);

  // create the plots per layer
  clusters->SetAlias("layer","barrellayer");
  std::vector<TProfile*> TIBlayers       = LayerPlots(clusters,3,4);
  std::vector<TProfile*> TOBlayers       = LayerPlots(clusters,5,6);
  clusters->SetAlias("layer","TIDlayer");
  std::vector<TProfile*> TIDlayers       = LayerPlots(clusters,4,3);
  clusters->SetAlias("layer","TECPlayer");
  std::vector<TProfile*> TECPthinlayers  = LayerPlots(clusters,6,9,"thickness<400");
  std::vector<TProfile*> TECPthicklayers = LayerPlots(clusters,6,9,"thickness>400");
  clusters->SetAlias("layer","TECMlayer");
  std::vector<TProfile*> TECMthinlayers  = LayerPlots(clusters,6,9,"thickness<400");
  std::vector<TProfile*> TECMthicklayers = LayerPlots(clusters,6,9,"thickness>400");
  TCanvas* c1 = prepareCanvas("TIBTID_layers");
  plotAll(c1,TIBlayers,1);
  plotAll(c1,TIDlayers,2);
  TCanvas* c2 = prepareCanvas("TOB_layers");
  plotAll(c2,TOBlayers,3);
  TCanvas* c3 = prepareCanvas("TEC_layers");
  plotAll(c3,TECPthinlayers,4);
  plotAll(c3,TECPthicklayers,5);
  plotAll(c3,TECMthinlayers,6);
  plotAll(c3,TECMthicklayers,7);
  std::vector<TProfile*> alllayers;
  alllayers.reserve(50);
  alllayers.insert(alllayers.end(),TIBlayers.begin(),TIBlayers.end());
  alllayers.insert(alllayers.end(),TIDlayers.begin(),TIDlayers.end());
  alllayers.insert(alllayers.end(),TOBlayers.begin(),TOBlayers.end());
  alllayers.insert(alllayers.end(),TECPthinlayers.begin(),TECPthinlayers.end());
  alllayers.insert(alllayers.end(),TECPthicklayers.begin(),TECPthicklayers.end());
  alllayers.insert(alllayers.end(),TECMthinlayers.begin(),TECMthinlayers.end());
  alllayers.insert(alllayers.end(),TECMthicklayers.begin(),TECMthicklayers.end());
  plotMaxima(alllayers)->Print(Form("layersGraph_%s.root",postfix.c_str()));
  c1->Print(Form("TIBTID_layers_%s.root",postfix.c_str()));
  c2->Print(Form("TOB_layers_%s.root",postfix.c_str()));
  c3->Print(Form("TEC_layers_%s.root",postfix.c_str()));

  // create the plots in R slices 
  std::vector<TProfile*> TIBrs   = RPlots(clusters,3,6,20,80);
  std::vector<TProfile*> TIDrs   = RPlots(clusters,4,4,80,120);
  std::vector<TProfile*> TOBrs   = RPlots(clusters,5,16,60,140);
  std::vector<TProfile*> TECpTrs  = RPlots(clusters,6,9,120,300,"thickness>400 && clglobalZ>0");
  std::vector<TProfile*> TECptrs  = RPlots(clusters,6,9,120,300,"thickness<400 && clglobalZ>0");
  std::vector<TProfile*> TECmTrs  = RPlots(clusters,6,9,120,300,"thickness>400 && clglobalZ<0");
  std::vector<TProfile*> TECmtrs  = RPlots(clusters,6,9,120,300,"thickness<400 && clglobalZ<0");
  TCanvas* c1b = prepareCanvas("TIBTID_distance");
  plotAll(c1b,TIBrs,1);
  plotAll(c1b,TIDrs,2);
  TCanvas* c2b = prepareCanvas("TOB_distance");
  plotAll(c2b,TOBrs,3);
  TCanvas* c3b = prepareCanvas("TEC_distance");
  plotAll(c3b,TECpTrs,4);
  plotAll(c3b,TECptrs,4);
  plotAll(c3b,TECmTrs,4);
  plotAll(c3b,TECmtrs,4);
  std::vector<TProfile*> allrs;
  allrs.reserve(TIBrs.size()+TIDrs.size()+TOBrs.size()+TECmTrs.size()+TECmtrs.size()+TECpTrs.size()+TECptrs.size());
  allrs.insert(allrs.end(),TIBrs.begin(),TIBrs.end());
  allrs.insert(allrs.end(),TIDrs.begin(),TIDrs.end());
  allrs.insert(allrs.end(),TOBrs.begin(),TOBrs.end());
  allrs.insert(allrs.end(),TECmTrs.begin(),TECmTrs.end());
  allrs.insert(allrs.end(),TECmtrs.begin(),TECmtrs.end());
  allrs.insert(allrs.end(),TECpTrs.begin(),TECpTrs.end());
  allrs.insert(allrs.end(),TECptrs.begin(),TECptrs.end());
  plotMaxima(allrs)->Print(Form("distanceGraph_%s.root",postfix.c_str()));
  c1b->Print(Form("TIBTID_distance_%s.root",postfix.c_str()));
  c2b->Print(Form("TOB_distance_%s.root",postfix.c_str()));
  c3b->Print(Form("TEC_distance_%s.root",postfix.c_str()));

  // create the plots per partition
  std::vector<TProfile*> TIBp   = RPlots(clusters,3,1,20,80);
  std::vector<TProfile*> TIDp   = RPlots(clusters,4,1,80,120);
  std::vector<TProfile*> TOBp   = RPlots(clusters,5,1,60,140);
  std::vector<TProfile*> TECTp  = RPlots(clusters,6,1,120,300,"thickness>400");
  std::vector<TProfile*> TECtp  = RPlots(clusters,6,1,120,300,"thickness<400");
  TCanvas* c4 = prepareCanvas("Partitions");
  plotAll(c4,TIBp,1);
  plotAll(c4,TIDp,2);
  plotAll(c4,TOBp,3);
  plotAll(c4,TECTp,4);
  plotAll(c4,TECtp,5);
  c4->Print(Form("Partitions_%s.root",postfix.c_str()));
}

