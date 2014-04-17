#ifndef PtDistrRadCase_h
#define PtDistrRadCase_h

// system include files
#include <iostream>
#include <string>

#include <TFile.h>
#include <TH1.h>
#include <TProfile.h>
#include <TH3.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TBranch.h>
#include <TTree.h>
#include <sstream>
#include <vector>
#include "TKey.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TopTreeProducer/interface/TRootJet.h" //not really needed

using namespace std;

class PtDistrRadCase{

 public:
  ~PtDistrRadCase();
  PtDistrRadCase(TString);
  void FillPlots(bool ContainsR, double pt);
  void Write();
  void WritePS(TString);

 private:
  TString *trueName_; 
  TString *falseName_; 
  TH1D *Pt_true_;
  TH1D *Pt_false_;
  
};

#endif
