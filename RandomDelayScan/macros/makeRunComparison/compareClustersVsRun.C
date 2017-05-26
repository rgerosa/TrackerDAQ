#include "../delayUtils.h"

static int reductionFactor    = 50;
static string outputDirectory = "";
static TFile* outputFile_maxCharge = NULL;
static TFile* outputFile_SoN       = NULL;
static TCanvas* canvas = NULL;
static vector<int> colorList = {kBlack,kRed,kBlue,kGreen+1,kCyan+1,kViolet+1,kGray,kOrange+1,kYellow+1};

vector<TH1F*> histo_maxCharge_TIB_L1;
vector<TH1F*> histo_maxCharge_TIB_L2;
vector<TH1F*> histo_maxCharge_TIB_L3;
vector<TH1F*> histo_maxCharge_TIB_L4;

vector<TH1F*> histo_maxCharge_TOB_L1;
vector<TH1F*> histo_maxCharge_TOB_L2;
vector<TH1F*> histo_maxCharge_TOB_L3;
vector<TH1F*> histo_maxCharge_TOB_L4;
vector<TH1F*> histo_maxCharge_TOB_L5;
vector<TH1F*> histo_maxCharge_TOB_L6;

vector<TH1F*> histo_maxCharge_TID_D1;
vector<TH1F*> histo_maxCharge_TID_D2;
vector<TH1F*> histo_maxCharge_TID_D3;

vector<TH1F*> histo_maxCharge_TECP_D1;
vector<TH1F*> histo_maxCharge_TECP_D2;
vector<TH1F*> histo_maxCharge_TECP_D3;
vector<TH1F*> histo_maxCharge_TECP_D4;
vector<TH1F*> histo_maxCharge_TECP_D5;
vector<TH1F*> histo_maxCharge_TECP_D6;
vector<TH1F*> histo_maxCharge_TECP_D7;
vector<TH1F*> histo_maxCharge_TECP_D8;
vector<TH1F*> histo_maxCharge_TECP_D9;

vector<TH1F*> histo_maxCharge_TECM_D1;
vector<TH1F*> histo_maxCharge_TECM_D2;
vector<TH1F*> histo_maxCharge_TECM_D3;
vector<TH1F*> histo_maxCharge_TECM_D4;
vector<TH1F*> histo_maxCharge_TECM_D5;
vector<TH1F*> histo_maxCharge_TECM_D6;
vector<TH1F*> histo_maxCharge_TECM_D7;
vector<TH1F*> histo_maxCharge_TECM_D8;
vector<TH1F*> histo_maxCharge_TECM_D9;

vector<TH1F*> histo_SoN_TIB_L1;
vector<TH1F*> histo_SoN_TIB_L2;
vector<TH1F*> histo_SoN_TIB_L3;
vector<TH1F*> histo_SoN_TIB_L4;

vector<TH1F*> histo_SoN_TOB_L1;
vector<TH1F*> histo_SoN_TOB_L2;
vector<TH1F*> histo_SoN_TOB_L3;
vector<TH1F*> histo_SoN_TOB_L4;
vector<TH1F*> histo_SoN_TOB_L5;
vector<TH1F*> histo_SoN_TOB_L6;

vector<TH1F*> histo_SoN_TID_D1;
vector<TH1F*> histo_SoN_TID_D2;
vector<TH1F*> histo_SoN_TID_D3;

vector<TH1F*> histo_SoN_TECP_D1;
vector<TH1F*> histo_SoN_TECP_D2;
vector<TH1F*> histo_SoN_TECP_D3;
vector<TH1F*> histo_SoN_TECP_D4;
vector<TH1F*> histo_SoN_TECP_D5;
vector<TH1F*> histo_SoN_TECP_D6;
vector<TH1F*> histo_SoN_TECP_D7;
vector<TH1F*> histo_SoN_TECP_D8;
vector<TH1F*> histo_SoN_TECP_D9;

vector<TH1F*> histo_SoN_TECM_D1;
vector<TH1F*> histo_SoN_TECM_D2;
vector<TH1F*> histo_SoN_TECM_D3;
vector<TH1F*> histo_SoN_TECM_D4;
vector<TH1F*> histo_SoN_TECM_D5;
vector<TH1F*> histo_SoN_TECM_D6;
vector<TH1F*> histo_SoN_TECM_D7;
vector<TH1F*> histo_SoN_TECM_D8;
vector<TH1F*> histo_SoN_TECM_D9;

void createHistograms(const size_t & npos){

  float xMin_maxCharge, xMax_maxCharge, xMin_SoN, xMax_SoN;
  int nBin_maxCharge, nBin_SoN;
  
  setLimitsAndBinning("maxCharge",xMin_maxCharge,xMax_maxCharge,nBin_maxCharge);
  setLimitsAndBinning("clCorrectedSignalOverNoise",xMin_SoN,xMax_SoN,nBin_SoN);
    
  for(int ihist = 0; ihist < int(npos); ihist++){

    histo_maxCharge_TIB_L1.push_back(new TH1F(Form("histo_maxCharge_TIB_L1_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TIB_L1.back()->Sumw2();
    histo_maxCharge_TIB_L2.push_back(new TH1F(Form("histo_maxCharge_TIB_L2_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TIB_L2.back()->Sumw2();
    histo_maxCharge_TIB_L3.push_back(new TH1F(Form("histo_maxCharge_TIB_L3_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TIB_L3.back()->Sumw2();
    histo_maxCharge_TIB_L4.push_back(new TH1F(Form("histo_maxCharge_TIB_L4_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TIB_L4.back()->Sumw2();

    histo_maxCharge_TOB_L1.push_back(new TH1F(Form("histo_maxCharge_TOB_L1_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TOB_L1.back()->Sumw2();
    histo_maxCharge_TOB_L2.push_back(new TH1F(Form("histo_maxCharge_TOB_L2_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TOB_L2.back()->Sumw2();
    histo_maxCharge_TOB_L3.push_back(new TH1F(Form("histo_maxCharge_TOB_L3_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TOB_L3.back()->Sumw2();
    histo_maxCharge_TOB_L4.push_back(new TH1F(Form("histo_maxCharge_TOB_L4_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TOB_L4.back()->Sumw2();
    histo_maxCharge_TOB_L5.push_back(new TH1F(Form("histo_maxCharge_TOB_L5_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TOB_L5.back()->Sumw2();
    histo_maxCharge_TOB_L6.push_back(new TH1F(Form("histo_maxCharge_TOB_L6_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TOB_L6.back()->Sumw2();

    histo_maxCharge_TID_D1.push_back(new TH1F(Form("histo_maxCharge_TID_D1_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TID_D1.back()->Sumw2();
    histo_maxCharge_TID_D2.push_back(new TH1F(Form("histo_maxCharge_TID_D2_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TID_D2.back()->Sumw2();
    histo_maxCharge_TID_D3.push_back(new TH1F(Form("histo_maxCharge_TID_D3_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TID_D3.back()->Sumw2();

    histo_maxCharge_TECP_D1.push_back(new TH1F(Form("histo_maxCharge_TECP_D1_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TECP_D1.back()->Sumw2();
    histo_maxCharge_TECP_D2.push_back(new TH1F(Form("histo_maxCharge_TECP_D2_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TECP_D2.back()->Sumw2();
    histo_maxCharge_TECP_D3.push_back(new TH1F(Form("histo_maxCharge_TECP_D3_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TECP_D3.back()->Sumw2();
    histo_maxCharge_TECP_D4.push_back(new TH1F(Form("histo_maxCharge_TECP_D4_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TECP_D4.back()->Sumw2();
    histo_maxCharge_TECP_D5.push_back(new TH1F(Form("histo_maxCharge_TECP_D5_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TECP_D5.back()->Sumw2();
    histo_maxCharge_TECP_D6.push_back(new TH1F(Form("histo_maxCharge_TECP_D6_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TECP_D6.back()->Sumw2();
    histo_maxCharge_TECP_D7.push_back(new TH1F(Form("histo_maxCharge_TECP_D7_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TECP_D7.back()->Sumw2();
    histo_maxCharge_TECP_D8.push_back(new TH1F(Form("histo_maxCharge_TECP_D8_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TECP_D8.back()->Sumw2();
    histo_maxCharge_TECP_D9.push_back(new TH1F(Form("histo_maxCharge_TECP_D9_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TECP_D9.back()->Sumw2();

    histo_maxCharge_TECM_D1.push_back(new TH1F(Form("histo_maxCharge_TECM_D1_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TECM_D1.back()->Sumw2();
    histo_maxCharge_TECM_D2.push_back(new TH1F(Form("histo_maxCharge_TECM_D2_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TECM_D2.back()->Sumw2();
    histo_maxCharge_TECM_D3.push_back(new TH1F(Form("histo_maxCharge_TECM_D3_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TECM_D3.back()->Sumw2();
    histo_maxCharge_TECM_D4.push_back(new TH1F(Form("histo_maxCharge_TECM_D4_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TECM_D4.back()->Sumw2();
    histo_maxCharge_TECM_D5.push_back(new TH1F(Form("histo_maxCharge_TECM_D5_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TECM_D5.back()->Sumw2();
    histo_maxCharge_TECM_D6.push_back(new TH1F(Form("histo_maxCharge_TECM_D6_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TECM_D6.back()->Sumw2();
    histo_maxCharge_TECM_D7.push_back(new TH1F(Form("histo_maxCharge_TECM_D7_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TECM_D7.back()->Sumw2();
    histo_maxCharge_TECM_D8.push_back(new TH1F(Form("histo_maxCharge_TECM_D8_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TECM_D8.back()->Sumw2();
    histo_maxCharge_TECM_D9.push_back(new TH1F(Form("histo_maxCharge_TECM_D9_%d",ihist),"",nBin_maxCharge,xMin_maxCharge,xMax_maxCharge));
    histo_maxCharge_TECM_D9.back()->Sumw2();

    histo_SoN_TIB_L1.push_back(new TH1F(Form("histo_SoN_TIB_L1_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TIB_L1.back()->Sumw2();
    histo_SoN_TIB_L2.push_back(new TH1F(Form("histo_SoN_TIB_L2_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TIB_L2.back()->Sumw2();
    histo_SoN_TIB_L3.push_back(new TH1F(Form("histo_SoN_TIB_L3_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TIB_L3.back()->Sumw2();
    histo_SoN_TIB_L4.push_back(new TH1F(Form("histo_SoN_TIB_L4_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TIB_L4.back()->Sumw2();

    histo_SoN_TOB_L1.push_back(new TH1F(Form("histo_SoN_TOB_L1_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TOB_L1.back()->Sumw2();
    histo_SoN_TOB_L2.push_back(new TH1F(Form("histo_SoN_TOB_L2_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TOB_L2.back()->Sumw2();
    histo_SoN_TOB_L3.push_back(new TH1F(Form("histo_SoN_TOB_L3_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TOB_L3.back()->Sumw2();
    histo_SoN_TOB_L4.push_back(new TH1F(Form("histo_SoN_TOB_L4_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TOB_L4.back()->Sumw2();
    histo_SoN_TOB_L5.push_back(new TH1F(Form("histo_SoN_TOB_L5_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TOB_L5.back()->Sumw2();
    histo_SoN_TOB_L6.push_back(new TH1F(Form("histo_SoN_TOB_L6_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TOB_L6.back()->Sumw2();

    histo_SoN_TID_D1.push_back(new TH1F(Form("histo_SoN_TIDP_D1_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TID_D1.back()->Sumw2();
    histo_SoN_TID_D2.push_back(new TH1F(Form("histo_SoN_TIDP_D2_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TID_D2.back()->Sumw2();
    histo_SoN_TID_D3.push_back(new TH1F(Form("histo_SoN_TIDP_D3_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TID_D3.back()->Sumw2();

    histo_SoN_TECP_D1.push_back(new TH1F(Form("histo_SoN_TECP_D1_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TECP_D1.back()->Sumw2();
    histo_SoN_TECP_D2.push_back(new TH1F(Form("histo_SoN_TECP_D2_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TECP_D2.back()->Sumw2();
    histo_SoN_TECP_D3.push_back(new TH1F(Form("histo_SoN_TECP_D3_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TECP_D3.back()->Sumw2();
    histo_SoN_TECP_D4.push_back(new TH1F(Form("histo_SoN_TECP_D4_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TECP_D4.back()->Sumw2();
    histo_SoN_TECP_D5.push_back(new TH1F(Form("histo_SoN_TECP_D5_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TECP_D5.back()->Sumw2();
    histo_SoN_TECP_D6.push_back(new TH1F(Form("histo_SoN_TECP_D6_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TECP_D6.back()->Sumw2();
    histo_SoN_TECP_D7.push_back(new TH1F(Form("histo_SoN_TECP_D7_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TECP_D7.back()->Sumw2();
    histo_SoN_TECP_D8.push_back(new TH1F(Form("histo_SoN_TECP_D8_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TECP_D8.back()->Sumw2();
    histo_SoN_TECP_D9.push_back(new TH1F(Form("histo_SoN_TECP_D9_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TECP_D9.back()->Sumw2();

    histo_SoN_TECM_D1.push_back(new TH1F(Form("histo_SoN_TECM_D1_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TECM_D1.back()->Sumw2();
    histo_SoN_TECM_D2.push_back(new TH1F(Form("histo_SoN_TECM_D2_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TECM_D2.back()->Sumw2();
    histo_SoN_TECM_D3.push_back(new TH1F(Form("histo_SoN_TECM_D3_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TECM_D3.back()->Sumw2();
    histo_SoN_TECM_D4.push_back(new TH1F(Form("histo_SoN_TECM_D4_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TECM_D4.back()->Sumw2();
    histo_SoN_TECM_D5.push_back(new TH1F(Form("histo_SoN_TECM_D5_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TECM_D5.back()->Sumw2();
    histo_SoN_TECM_D6.push_back(new TH1F(Form("histo_SoN_TECM_D6_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TECM_D6.back()->Sumw2();
    histo_SoN_TECM_D7.push_back(new TH1F(Form("histo_SoN_TECM_D7_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TECM_D7.back()->Sumw2();
    histo_SoN_TECM_D8.push_back(new TH1F(Form("histo_SoN_TECM_D8_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TECM_D8.back()->Sumw2();
    histo_SoN_TECM_D9.push_back(new TH1F(Form("histo_SoN_TECM_D9_%d",ihist),"",nBin_SoN,xMin_SoN,xMax_SoN));
    histo_SoN_TECM_D9.back()->Sumw2();
  }

  return;
}

// laudau + gauss fit
void fitSignalShape(const vector<TH1F*> histo, const vector<float> & runList, TFile* outputFile, const string & observable){

  Double_t parameters[4];
  Double_t parametersHigh[4];
  Double_t parametersLow[4];
  const Double_t *fit_parameters;
  const Double_t *fit_parameters_error;
  Double_t chi2;
  Int_t    ndf;

  int ipos   = 0;

  canvas->cd();
  canvas->Clear();
  TH1* frame = (TH1*) histo.at(0)->Clone("frame");
  frame->Reset();
  frame->GetYaxis()->SetTitle("a.u.");
  frame->GetXaxis()->SetTitle(("Cluster "+observable).c_str());
  frame->GetYaxis()->SetTitleSize(0.045);
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetLabelSize(0.038);
  frame->GetXaxis()->SetLabelSize(0.038);
  frame->Draw();

  CMS_lumi(canvas,"",true,true);

  TLegend leg (0.5,0.5,0.95,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);

  float maxValue = 0.;

  for(auto ihist : histo){

    ihist->Scale(1./ihist->Integral()); // normalize to a.u.

   // Width of the landau distribution                                                                                                                                      
    parameters[0]     = 2.; 
    parametersHigh[0] = 50.; 
    parametersLow[0]  = 0.001;
    // MPV of landau peak                                                                                                                                                    
    parameters[1]     = ihist->GetMean();  
    parametersHigh[1] = ihist->GetBinLowEdge(1); 
    parametersLow[1]  = ihist->GetBinLowEdge(ihist->GetNbinsX()+1);
    // Total area                                                                                                                                                            
    parameters[2]     = ihist->Integral(); 
    parametersHigh[2] = ihist->Integral()*5; 
    parametersLow[2] = ihist->Integral()/5;
    // width of gaussian                                                                                                                                                     
    parameters[3] = 10; 
    parametersHigh[3] = 40; 
    parametersLow[3] = 1.;

   
    // create the function                                                                                                                                                   
    TF1 *    fitfunc = new TF1(Form("fit_%s",ihist->GetName()),langaufun,ihist->GetBinLowEdge(1),ihist->GetBinLowEdge(ihist->GetNbinsX()+1),4);
    if(observable == "maxCharge"){
      if(TString(ihist->GetName()).Contains("TIB"))
	fitfunc->SetRange(ihist->GetBinCenter(ihist->GetMaximumBin())-1.5*ihist->GetRMS(),ihist->GetBinCenter(ihist->GetMaximumBin())+ihist->GetRMS()*2.5);
      if(TString(ihist->GetName()).Contains("TID"))
	fitfunc->SetRange(ihist->GetBinCenter(ihist->GetMaximumBin())-1.5*ihist->GetRMS(),ihist->GetBinCenter(ihist->GetMaximumBin())+ihist->GetRMS()*2.5);
      else if(TString(ihist->GetName()).Contains("TOB"))
	fitfunc->SetRange(ihist->GetBinCenter(ihist->GetMaximumBin())-2.0*ihist->GetRMS(),ihist->GetBinCenter(ihist->GetMaximumBin())+ihist->GetRMS()*3.0);
      if(TString(ihist->GetName()).Contains("TEC"))
	fitfunc->SetRange(ihist->GetBinCenter(ihist->GetMaximumBin())-1.5*ihist->GetRMS(),ihist->GetBinCenter(ihist->GetMaximumBin())+ihist->GetRMS()*2.5);
    }    
    else
      fitfunc->SetRange(ihist->GetBinCenter(ihist->GetMaximumBin())-1.5*ihist->GetRMS(),ihist->GetBinCenter(ihist->GetMaximumBin())+ihist->GetRMS()*3.0);
    fitfunc->SetParameters(parameters);
    fitfunc->SetParNames("Width","MPV","Area","GSigma");
    fitfunc->SetParLimits(0,parametersLow[0],parametersHigh[0]);
    fitfunc->SetParLimits(0,parametersLow[1],parametersHigh[1]);
    fitfunc->SetParLimits(0,parametersLow[2],parametersHigh[2]);
    fitfunc->SetParLimits(0,parametersLow[3],parametersHigh[3]);
    // make fit and get parameters                                                                                                                                           
    TFitResultPtr fitResult = ihist->Fit(fitfunc,"RSQN");
    if(not fitResult.Get()) continue;
    fit_parameters = fitResult->GetParams();
    fit_parameters_error = fitResult->GetErrors();
    chi2   = fitResult->Chi2();
    ndf    = fitResult->Ndf();

    if(ipos < colorList.size()){
      ihist->SetMarkerColor(colorList.at(ipos));
      ihist->SetLineColor(colorList.at(ipos));
      fitfunc->SetLineColor(colorList.at(ipos));
    }
    else{
      ihist->SetMarkerColor(colorList.at(ipos-colorList.size()));
      ihist->SetLineColor(colorList.at(ipos-colorList.size()));
      fitfunc->SetLineColor(colorList.at(ipos-colorList.size()));
    }
    ihist->SetMarkerStyle(20);
    ihist->SetMarkerSize(0.8);
    fitfunc->SetLineWidth(2);
    outputFile->cd();
    ihist->Write(ihist->GetName());

    if(ihist->GetMaximum() > maxValue)
      maxValue = ihist->GetMaximum();
    
    canvas->cd();
    ihist->Draw("same");
    fitfunc->Draw("Lsame");
    
    TLegendEntry* l1 = leg.AddEntry((TObject*)0,Form("Run = %d",int(runList.at(ipos))),"");
    if(ipos < colorList.size())
      l1->SetTextColor(colorList.at(ipos));
    else
      l1->SetTextColor(colorList.at(ipos-colorList.size()));
    l1 = leg.AddEntry((TObject*)0,Form("Mean       = %.2f #pm %.2f",ihist->GetMean(),ihist->GetMeanError()),"");
    if(ipos < colorList.size())
      l1->SetTextColor(colorList.at(ipos));
    else
      l1->SetTextColor(colorList.at(ipos-colorList.size()));
    l1 = leg.AddEntry((TObject*)0,Form("MPV        = %.2f #pm %.2f",fitfunc->GetMaximumX(),fit_parameters_error[1]),"");      
    if(ipos < colorList.size())
      l1->SetTextColor(colorList.at(ipos));
    else
      l1->SetTextColor(colorList.at(ipos-colorList.size()));

    ipos++;

  }

  frame->GetYaxis()->SetRangeUser(0.,maxValue*1.2);
  leg.Draw("same");
  
  canvas->SaveAs((outputDirectory+"/"+histo.at(0)->GetName()+"_"+observable+".pdf").c_str(),"pdf");
  canvas->SaveAs((outputDirectory+"/"+histo.at(0)->GetName()+"_"+observable+".png").c_str(),"png");

  return;
}


void compareClustersVsRun(string inputDIR, vector<float> runList, string outputDIR, string postfix = "tree"){


  outputDirectory = outputDIR;

  // prepare style and load macros                                                                                                                                            
  setTDRStyle();
  gROOT->SetBatch(kTRUE);

  system(("mkdir -p "+outputDIR).c_str());

  // open the file and prepare the cluster tree, adding the other trees as frined --> memory consuming                                                                       
  std::cout<<"#########################################"<<std::endl;
  std::cout<<"###### compareCluster proprierties ######"<<std::endl;
  std::cout<<"#########################################"<<std::endl;

  std::cout<<"### Create basic histograms "<<endl;
  // to allocate all the histograms
  size_t size = runList.size();
  createHistograms(size);


  std::cout<<"### Make input file list"<<std::endl;
  system(("find "+inputDIR+" -name \"*"+postfix+"*.root\" > file.temp").c_str());
  std::ifstream infile;
  string line;
  vector<string> fileList;
  infile.open("file.temp",ifstream::in);
  if(infile.is_open()){
    while(!infile.eof()){
      getline(infile,line);
      if(line != "" and TString(line).Contains(".root"))
        fileList.push_back(line);
    }
  }
  system("rm file.temp");
  std::sort(fileList.begin(),fileList.end());

  std::cout<<"### Build the input TTree Vector"<<std::endl;
  vector<TFile*> files;
  vector<TTree*> clusters;
  for(auto ifile : fileList){
    files.push_back(TFile::Open(ifile.c_str()));
    clusters.push_back((TTree*) (files.back()->FindObjectAny("clusters")));
  }

  // create the histogram list
  cout<<"### Start looping on events"<<std::endl;
  for(int iTree = 0; iTree < clusters.size(); iTree++){

    std::cout<<"Analyzing tree "<<iTree<<" of "<<clusters.size()<<std::endl;

    // set branches for the cluster, readoutmap and no corrections trees                                                                                                        
    uint32_t detid, runid;
    float    clCorrectedSignalOverNoise, clSignalOverNoise, maxCharge, clglobalZ;
    clusters.at(iTree)->SetBranchStatus("*",kFALSE);
    clusters.at(iTree)->SetBranchStatus("detid",kTRUE);
    clusters.at(iTree)->SetBranchStatus("clCorrectedSignalOverNoise",kTRUE);
    clusters.at(iTree)->SetBranchStatus("clSignalOverNoise",kTRUE);
    clusters.at(iTree)->SetBranchStatus("maxCharge",kTRUE);
    clusters.at(iTree)->SetBranchStatus("runid",kTRUE);
    clusters.at(iTree)->SetBranchStatus("clglobalZ",kTRUE);
    clusters.at(iTree)->SetBranchAddress("detid",&detid);
    clusters.at(iTree)->SetBranchAddress("runid",&runid);
    clusters.at(iTree)->SetBranchAddress("clCorrectedSignalOverNoise",&clCorrectedSignalOverNoise);
    clusters.at(iTree)->SetBranchAddress("clSignalOverNoise",&clSignalOverNoise);
    clusters.at(iTree)->SetBranchAddress("maxCharge",&maxCharge);
    clusters.at(iTree)->SetBranchAddress("clglobalZ",&clglobalZ);

    for(long int iEvent = 0; iEvent < clusters.at(iTree)->GetEntries()/reductionFactor; iEvent++){
      
      cout.flush();
      if(iEvent % 10000 == 0) cout<<"\r"<<"iEvent "<<100*double(iEvent)/(clusters.at(iTree)->GetEntries()/reductionFactor)<<" % ";      
      // take the event                                                                                                                                                         
      clusters.at(iTree)->GetEntry(iEvent);
      
      uint32_t subdetid    = int((detid-0x10000000)/0x2000000);
      uint32_t barrellayer = int((detid%33554432)/0x4000);
      uint32_t TIDlayer    = int((detid%33554432)/0x800)%4;
      uint32_t TECPlayer   = int((detid%33554432)/0x4000)-32;
      uint32_t TECMlayer   = int((detid%33554432)/0x4000)-16;
      
      size_t pos = std::find(runList.begin(),runList.end(),float(runid))-runList.begin();
      if(pos == runList.size()) continue;
      
      if(subdetid == 3){ // TIB
	if(barrellayer == 1){
	  histo_maxCharge_TIB_L1.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TIB_L1.at(pos)->Fill(clCorrectedSignalOverNoise);
	}
	else if(barrellayer == 2){
	  histo_maxCharge_TIB_L2.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TIB_L2.at(pos)->Fill(clCorrectedSignalOverNoise);
	}
	else if(barrellayer == 3){
	  histo_maxCharge_TIB_L3.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TIB_L3.at(pos)->Fill(clCorrectedSignalOverNoise);
	}
	else if(barrellayer == 4){
	  histo_maxCharge_TIB_L4.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TIB_L4.at(pos)->Fill(clCorrectedSignalOverNoise);
	}
      }
      else if(subdetid == 5){ // TOB
	if(barrellayer == 1){
	  histo_maxCharge_TOB_L1.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TOB_L1.at(pos)->Fill(clCorrectedSignalOverNoise);
	}
	else if(barrellayer == 2){
	  histo_maxCharge_TOB_L2.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TOB_L2.at(pos)->Fill(clCorrectedSignalOverNoise);
	}
	else if(barrellayer == 3){
	  histo_maxCharge_TOB_L3.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TOB_L3.at(pos)->Fill(clCorrectedSignalOverNoise);
	}
	else if(barrellayer == 4){
	  histo_maxCharge_TOB_L4.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TOB_L4.at(pos)->Fill(clCorrectedSignalOverNoise);
	}
	else if(barrellayer == 5){
	  histo_maxCharge_TOB_L5.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TOB_L5.at(pos)->Fill(clCorrectedSignalOverNoise);
	}
	else if(barrellayer == 6){
	  histo_maxCharge_TOB_L6.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TOB_L6.at(pos)->Fill(clCorrectedSignalOverNoise);
	}
      }    
      else if(subdetid == 4){ // TID
	if(TIDlayer == 1){
	histo_maxCharge_TID_D1.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	histo_SoN_TID_D1.at(pos)->Fill(clCorrectedSignalOverNoise);	  
	}
	else if(TIDlayer == 2){
	  histo_maxCharge_TID_D2.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TID_D2.at(pos)->Fill(clCorrectedSignalOverNoise);	  	  
	}
	else if(TIDlayer == 3){
	  histo_maxCharge_TID_D3.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TID_D3.at(pos)->Fill(clCorrectedSignalOverNoise);	  	  
	}	
      }
      else if(subdetid == 6){ // TEC
	// TECP
	if( TECPlayer == 1 and clglobalZ > 0){
	  histo_maxCharge_TECP_D1.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TECP_D1.at(pos)->Fill(clCorrectedSignalOverNoise);	  	  
	}
	else if( TECPlayer == 2 and clglobalZ > 0){
	  histo_maxCharge_TECP_D2.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TECP_D2.at(pos)->Fill(clCorrectedSignalOverNoise);	  	  
	}
	else if( TECPlayer == 3 and clglobalZ > 0){
	  histo_maxCharge_TECP_D3.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TECP_D3.at(pos)->Fill(clCorrectedSignalOverNoise);	  	  
	}
	else if( TECPlayer == 4 and clglobalZ > 0){
	  histo_maxCharge_TECP_D4.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TECP_D4.at(pos)->Fill(clCorrectedSignalOverNoise);	  	  
	}
	else if( TECPlayer == 5 and clglobalZ > 0){
	  histo_maxCharge_TECP_D5.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TECP_D5.at(pos)->Fill(clCorrectedSignalOverNoise);	  	  
	}
	else if( TECPlayer == 6 and clglobalZ > 0){
	  histo_maxCharge_TECP_D6.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TECP_D6.at(pos)->Fill(clCorrectedSignalOverNoise);	  	  
	}
	else if( TECPlayer == 7 and clglobalZ > 0){
	  histo_maxCharge_TECP_D7.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TECP_D7.at(pos)->Fill(clCorrectedSignalOverNoise);	  	  
	}
	else if( TECPlayer == 8 and clglobalZ > 0){
	  histo_maxCharge_TECP_D8.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TECP_D8.at(pos)->Fill(clCorrectedSignalOverNoise);	  	  
	}
	else if( TECPlayer == 9 and clglobalZ > 0){
	  histo_maxCharge_TECP_D9.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TECP_D9.at(pos)->Fill(clCorrectedSignalOverNoise);	  	  
	}
	// TECM 
	if( TECMlayer == 1 and clglobalZ < 0){
	  histo_maxCharge_TECM_D1.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TECM_D1.at(pos)->Fill(clCorrectedSignalOverNoise);	  	  
	}
	else if( TECMlayer == 2 and clglobalZ < 0){
	  histo_maxCharge_TECM_D2.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TECM_D2.at(pos)->Fill(clCorrectedSignalOverNoise);	  	  
	}
	else if( TECMlayer == 3 and clglobalZ < 0){
	  histo_maxCharge_TECM_D3.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TECM_D3.at(pos)->Fill(clCorrectedSignalOverNoise);	  	  
	}
	else if( TECMlayer == 4 and clglobalZ < 0){
	  histo_maxCharge_TECM_D4.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TECM_D4.at(pos)->Fill(clCorrectedSignalOverNoise);	  	  
	}
	else if( TECMlayer == 5 and clglobalZ < 0){
	  histo_maxCharge_TECM_D5.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TECM_D5.at(pos)->Fill(clCorrectedSignalOverNoise);	  	  
	}
	else if( TECMlayer == 6 and clglobalZ < 0){
	  histo_maxCharge_TECM_D6.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TECM_D6.at(pos)->Fill(clCorrectedSignalOverNoise);	  	  
	}
	else if( TECMlayer == 7 and clglobalZ < 0){
	  histo_maxCharge_TECM_D7.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TECM_D7.at(pos)->Fill(clCorrectedSignalOverNoise);	  	  
	}
	else if( TECMlayer == 8 and clglobalZ < 0){
	  histo_maxCharge_TECM_D8.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
	  histo_SoN_TECM_D8.at(pos)->Fill(clCorrectedSignalOverNoise);	  	  
	}
	else if( TECMlayer == 9 and clglobalZ < 0){
	  histo_maxCharge_TECM_D9.at(pos)->Fill(maxCharge*clCorrectedSignalOverNoise/clSignalOverNoise);
          histo_SoN_TECM_D9.at(pos)->Fill(clCorrectedSignalOverNoise);	  	  
	}	
      }
    }
    files.at(iTree)->Close(); // prevent memory leaks
    std::cout<<std::endl;
  }
  /////////////////////////  
  outputFile_maxCharge = new TFile((outputDIR+"/distributionPerLayer_maxCharge.root").c_str(),"RECREATE");
  outputFile_SoN = new TFile((outputDIR+"/distributionPerLayer_SoN.root").c_str(),"RECREATE");

  canvas = new TCanvas("canvas","",600,625);
  
  ///// fitting part
  fitSignalShape(histo_maxCharge_TIB_L1,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TIB_L2,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TIB_L3,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TIB_L4,runList,outputFile_maxCharge,"maxCharge");
    
  fitSignalShape(histo_maxCharge_TOB_L1,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TOB_L2,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TOB_L3,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TOB_L4,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TOB_L5,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TOB_L6,runList,outputFile_maxCharge,"maxCharge");
  
  fitSignalShape(histo_maxCharge_TID_D1,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TID_D2,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TID_D3,runList,outputFile_maxCharge,"maxCharge");
  
  fitSignalShape(histo_maxCharge_TECP_D1,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TECP_D2,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TECP_D3,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TECP_D4,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TECP_D5,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TECP_D6,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TECP_D7,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TECP_D8,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TECP_D9,runList,outputFile_maxCharge,"maxCharge");
  
  fitSignalShape(histo_maxCharge_TECM_D1,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TECM_D2,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TECM_D3,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TECM_D4,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TECM_D5,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TECM_D6,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TECM_D7,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TECM_D8,runList,outputFile_maxCharge,"maxCharge");
  fitSignalShape(histo_maxCharge_TECM_D9,runList,outputFile_maxCharge,"maxCharge");
  
  ///
  fitSignalShape(histo_SoN_TIB_L1,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TIB_L2,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TIB_L3,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TIB_L4,runList,outputFile_SoN,"SoN");
  
  fitSignalShape(histo_SoN_TOB_L1,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TOB_L2,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TOB_L3,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TOB_L4,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TOB_L5,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TOB_L6,runList,outputFile_SoN,"SoN");
  
  fitSignalShape(histo_SoN_TID_D1,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TID_D2,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TID_D3,runList,outputFile_SoN,"SoN");
  
  fitSignalShape(histo_SoN_TECP_D1,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TECP_D2,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TECP_D3,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TECP_D4,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TECP_D5,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TECP_D6,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TECP_D7,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TECP_D8,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TECP_D9,runList,outputFile_SoN,"SoN");

  fitSignalShape(histo_SoN_TECM_D1,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TECM_D2,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TECM_D3,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TECM_D4,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TECM_D5,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TECM_D6,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TECM_D7,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TECM_D8,runList,outputFile_SoN,"SoN");
  fitSignalShape(histo_SoN_TECM_D9,runList,outputFile_SoN,"SoN");
  
  return;
  
}
