#ifndef WorkingPointBin_h
#define WorkingPointBin_h

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

using namespace std;

class WorkingPointBin{ 

 public:

  ~WorkingPointBin();//destructor should delete pointers
  WorkingPointBin(int);
  void DefineHistos(int,double,double,int,double,double,int,double,double,int,double,double,int,double,double,int,double,double,int,double,double,int,double,double);
  void FillHistos(double*);
  void WriteHistos(TFile*);
 private:

  TString *genericName_; 
  TH1D *TH1D_Eff_;
  TString *title_TH1D_Eff_;
  TH1D *TH1D_EffMC_;
  TString *title_TH1D_EffMC_;
  TH1D *TH1D_EffErr_;
  TString *title_TH1D_EffErr_;

  TH1D *TH1D_EffDiff_;
  TString *title_TH1D_EffDiff_;
  TH1D *TH1D_EffPull_;
  TString *title_TH1D_EffPull_;

  TH1D *TH1D_FitParam_[2];
  TString *title_TH1D_FitParam_[2];

  TH1D *TH1D_FRatio_;
  TString *title_TH1D_FRatio_;
  TH1D *TH1D_FRatioMC_;
  TString *title_TH1D_FRatioMC_;

  TH1D *TH1D_FRatioErr_;
  TString *title_TH1D_FRatioErr_;
  TH1D *TH1D_FRatioMCErr_;
  TString *title_TH1D_FRatioMCErr_;

  TH1D *TH1D_FRatioDiff_;
  TString *title_TH1D_FRatioDiff_;
  TH1D *TH1D_FRatioPull_;
  TString *title_TH1D_FRatioPull_;

};

#endif
