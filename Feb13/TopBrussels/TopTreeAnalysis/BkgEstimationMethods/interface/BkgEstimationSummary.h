#ifndef BkgEstimationSummary_h
#define BkgEstimationSummary_h

#include <iostream>
#include <vector>
#include <iomanip>
#include <Math/VectorUtil.h>

#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TString.h"

#include "TopTreeAnalysis/Tools/interface/PlottingTools.h"

using namespace std;

	/**
	//	Aim: Creating summary plots (one now)  with MC-expectation, 
	//		data-driven estimation and data themself.
	//		It's actually done for QCD, WJets and Ttbar
	//
	//	Example of usage:
	//	- call constructor
	//	- SetColor with default colors are not nice
	//	- Draw()
	//	- Write()
	//
	//	Missing pieces:
	//	- update it for SingleTop, Zjets and VV
	*/

class BkgEstimationSummary {
  
 public:
  
  BkgEstimationSummary();
  BkgEstimationSummary(TString Name, TString XaxisTitle, TString YaxisTitle, TH1F* hData, TH1F* hQCDEstim, TH1F* hWJetEstim, TH1F* hTtJetEstim=0);
  BkgEstimationSummary(TString Name, TString XaxisTitle, TString YaxisTitle, TH1F* hData, TH1F* hQCDEstim, TH1F* hQCDMC, TH1F* hWJetEstim, TH1F* hWJetMC, TH1F* hTtJetEstim=0, TH1F* hTtJetMC=0);
  ~BkgEstimationSummary();

  void SetColor(int DataColor, int QCDColor, int WJetColor, int TtJetColor) {DataColor_=DataColor; QCDColor_=QCDColor; WJetColor_=WJetColor; TtJetColor_=TtJetColor;}

  void Draw();
  void Write(TFile* fout, string label = string(""));

  ////////////////////////////
  //  Access to histo
  ////////////////////////////
  // Data
  TH1F* GetData() const { return hData_; };
  // Expectation from MC
  TH1F* GetQCDMC() const {return hQCDMC_;}; 
  TH1F* GetWJetMC() const {return hWJetMC_;}; 
  TH1F* GetTtJetMC() const {return hTtJetMC_;}; 
  TH1F* GetCompleteMC() const {return hCompleteMC_;}
  // Estimation
  TH1F* GetQCDEstimation() const {return hQCDEstim_;}; 
  TH1F* GetWJetEstimation() const {return hWJetEstim_;}; 
  TH1F* GetTtJetEstimation() const {return hTtJetEstim_;}; 
  TH1F* GetCompleteEstimation() const {return hCompleteEstim_;}

 private:
  

  TString Name_;
  TString XaxisTitle_;
  TString YaxisTitle_;

  int DataColor_;
  int QCDColor_;
  int WJetColor_;
  int TtJetColor_;
  int QCDLColor_;
  int WJetLColor_;
  int TtJetLColor_;
  int LWidth_;

  // data
  TH1F* hData_;
  
  // estimation
  TH1F* hQCDEstim_;
  TH1F* hWJetEstim_;
  TH1F* hTtJetEstim_;
  TH1F* hCompleteEstim_;

  // expectation from MC
  TH1F* hQCDMC_;
  TH1F* hWJetMC_;
  TH1F* hTtJetMC_;
  TH1F* hCompleteMC_;
  
  // Used in Draw() method
  TLegend* legend;
  TCanvas* SumCanvas;
  THStack* MCStackHisto;
};

#endif
