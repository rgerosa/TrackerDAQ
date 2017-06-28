#include "../CMS_lumi.h"

void plotStripCanvas(string inputFileName, string outputDIR){

  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());
  gROOT->SetBatch(kTRUE);

  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  TList* listKeys = inputFile->GetListOfKeys();
  TIter next(listKeys);
  while(TKey *obj = (TKey*) next()){
    if(string(obj->GetClassName()) == "TCanvas"){
      TCanvas* canvas = (TCanvas*) obj->ReadObj();
      canvas->Draw(); 
      canvas->Modified(); 
      canvas->Update();
      canvas->SaveAs((outputDIR+"/"+string(obj->GetName())+".png").c_str(),"png");
      
    }
  }

  
}
