
#ifndef PtEtaBin_h
#define PtEtaBin_h

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
#include "TObjArray.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TRandom3.h"
#include "TFractionFitter.h"
#include "TStyle.h"

#include "TopTreeProducer/interface/TRootJet.h"

//from Greg's http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/NewPhysicsBrussels/TopTreeAnalysis/BkgEstimationMethods/interface/VJetEstimation.h?revision=1.1.2.10&view=markup&pathrev=CMSSW_36X
// RooFit librairies
#include "RooArgSet.h"
#include "RooAddition.h"
#include "RooCategory.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooHistPdf.h"
//#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooGenericPdf.h"
//#include "RooSimultaneous.h"

#include "RooMCStudy.h"
#include "RooMinuit.h"
//#include "RooMinimizer.h"
//#include "RooNLLVar.h"
//#include "RooProfileLL.h"
#include "RooPlot.h"
#include "RooChi2Var.h"

#include "RooFitResult.h"

#include <sys/stat.h>

class PtEtaBin{ 
  
 public:
  
  //inline PtEtaBin(){};
  //~PtEtaBin(){};
  ~PtEtaBin();
  PtEtaBin(int, int, int, int, double, double, double, double, bool, bool);

  void DefineSignalSamplePlots(int, int, double, double, int, double, double, int, double, double, int, double, double, int, double, double, int*, double*, double*);
  void DefineControlSamplePlots(int,int, double, double, int, double, double, int, double, double, int*, double*, double*);
  //void FillSignalSamplePlots(double, int, double, double, double, double, int, double, double, double, double, double, double, double);

  void FillSignalSamplePlots(double, double, int, bool, bool, double, double, double*, double, double, double, double, double, double, double, double,double);
  void FillControlSamplePlots(double, int, bool, bool, double, double, double, double, double, double, double, double, double);

	void FillXStemplates(double weight, string datasetname, int partonflavour, double btag, double *btagCuts, double controlVar0, double m3,double lowCutVar0, double centralLowCutVar0, double centralUpCutVar0, double upCutVar0);

	void SetVarBins(vector<float> rangesbTag);

  void WriteHistosToRootFile(bool,bool,bool,bool,bool,bool,bool,bool,bool,bool,bool,bool);
  void WriteHistoToPSFile(TString*,bool);
  void SetErrorsSignalSamples();
  void SetErrorsControlSamples();

  void MakeSoverSBPlots();
  void MakeMCEffPlots();
  void Make1DXplots();
  void Make1DYplots();
  void Make1DYVar2plots();
  void SetError1DXplots();
  void SetError1DYplots();
  void SetError1DYVar2plots();
  void MakeProfileXplots();
  void MakeSCprojectionPlots();

  void MakeXRatioPlot(bool);
  
void GetLeftRightPars(double*, double*, double*);
  void SetLeftRightPars(double*, double*, double*);
  
  void Make2DRatioPlot();
  void MakeSCVar12RatioPlot();
  void ReweighLeft();
  void ReweighRight();
  //obsolete void ReweighRightChangeFitParams(double, double);
  //void FillReweighRight(bool, bool, double, int, double, double, double, double, int, double, double, double, double, double, double, double);//they have the same input (but not alle variables ared used in each function)
  void FillReweighRight(bool, bool, double, int, double, double, double, double, double, double, double, double, double, double, double, double);//they have the same input (but not alle variables ared used in each function)
  //void FillReweighControl(bool, bool, double, int, double, double, double, double, int, double, double, double, double, double, double, double);//they have the same input (but not alle variables ared used in each function)
  void FillReweighControl(double,double*,bool, double, int, bool, bool, double, double, double, double, double, double, double, double);//they have the same input (but not alle variables ared used in each function)
  void MakeReweighRatio();
  void GetLRratio(bool, bool, double, bool, double);  
  void GetPercentiles(double,double,double,int*);
  void MeasureEff(bool);
  void MeasureEffLR(bool);
    void MeasureEffRR(bool);

    void MeasureMistagEffRR(bool);

  void GetWPEff(bool, bool, double, double*, bool, bool, double);
  void CoutWPEff(bool, bool, double, double*, bool, bool, int, int, double);

  vector<double> doMLJTemplateFit(string chi2cut,int mode, string data_postfi,int nSystematic);
	
	std::map<std::string,float> GetEffCalcDetails() const { return EffCalcDetails_; }
	
	std::vector<string> GetFitPlotPaths() const { return FitPlotPaths; }
	
  double getmljMeanVal();
  double getmljMeannoRWVal();
  double getmljMeanMCVal();
  double getmljSigmaVal();
  double getmljSigmanoRWVal();
  double getmljSigmaMCVal();

  double getmlj_W_MeanVal();
  double getmlj_W_MeannoRWVal();
  double getmlj_W_MeanMCVal();
  double getmlj_W_SigmaVal();
  double getmlj_W_SigmanoRWVal();
  double getmlj_W_SigmaMCVal();

  double getmlj_R_MeanVal();
  double getmlj_R_MeannoRWVal();
  double getmlj_R_MeanMCVal();
  double getmlj_R_SigmaVal();
  double getmlj_R_SigmanoRWVal();
  double getmlj_R_SigmaMCVal();
	
	float nTTbar () const { return nTTbar_;}
    
    void setLumi(double l) { lumi_ = l; }
    
    bool findTemplates(string chi2cut,int mode, string data_postfix);
    void loadTemplates(std::map<string,TH1D*> &h, double &lumi, string chi2cut,int mode, string data_postfix);
	vector<float> doTemplateFit (TH1D* ttbar, TH1D* vvmc, TH1D* vvdata,TH1D* data, TString PrefixPlot);

private:
    
    void writeTemplates(string chi2cut,int mode, string data_postfix);
    
    void FillTemplates(string process,double weight, int partonflavour, double btag, double *btagCuts, double controlVar0, double m3);

    int nBinsVar0_;
    float lowRangeVar0_;
    float upRangeVar0_;
    
    double lumi_;
    double templateLumi_;

  void EffCalculation(bool, TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*);
  void CalculateReweigh(TH2D*,TH1D*,TH1D*);
  void DiffCalculation(TH1D*,TH1D*,TH1D*);
  void WriteSignalSamplePlots();
  void WriteControlSamplePlots();
  void WriteSoverSBPlot();
  void WriteMCEffPlots(bool);
  void Write1DXplots(bool);
  void Write1DYplots();
  void Write1DYVar2plots();
  void WriteEffPlots();
  void WriteXRatioPlot();
  void WriteSCPlot();
  void Write2DRatioPlot();
  void WriteEffPlotsLR();
  void WriteEffPlotsRR();
  void SetError(TH1D*);
  void WriteProfileXplots();

  void DrawbtagNR(TString*);
  void DrawEffNR(TString*);
  void DrawEffLR(TString*);
  void DrawEffRR(TString*);
  void DrawEffDiffNR(TString*);
  void DrawEffDiffLR(TString*);
  void DrawEffDiffRR(TString*);
  void DrawBasicFPlot(TString*);
  void DrawControlFPlot(TString*);
  void DrawControlRWFPlot(TString*);
  void DrawSignalFPlot(TString*);
  void DrawControlOnlyRPlot(TString*);
  void DrawControlOnlyWPlot(TString*);
  void DrawPtDistro(TString*);
  void DrawSoverSB(TString*);
  void DrawControlRWXPlot(TString*);
  void DrawControlRWYPlot(TString*);
  void DrawControlWSignalControlXPlot(TString*);
  void DrawControlWSignalControlYPlot(TString*);
  void DrawControlRSignalControlXPlot(TString*);
  void DrawControlRSignalControlYPlot(TString*);

  //TString GiveName(TString);
  void GiveName(TString*);

  int debug_;
  int nVar0_;
  int nVar1_;
  int nBdisc_;
  double ptbinlow_;
  double etabinlow_; 
  double ptbinup_;
  double etabinup_; 
  bool varBinSize_;
  bool doShift_;
  TString *genericName; 

  double rangesbTag_[101];//number of bins and ranges are copied by hand.
  double rangesVar1_[21];
	
	std::map<TString,TH1D*> histo1D;
	std::map<TString,TH2D*> histo2D;

    TH1D * TH1Data_PartonFlavor;
    TString title_PartonFlavor_;

    
  TH2D * TH2Sng;
  TH2D * TH2Bkg;
  TH2D * TH2Data;
  TString titleSng_;
  TString titleBkg_;
  TString titleData_;
 
  TH1D * TH2Sng_ProfileX;
  TH1D * TH2Bkg_ProfileX;
  TH1D * TH2Data_ProfileX;
  TString titleSng_ProfileX_;
  TString titleBkg_ProfileX_;
  TString titleData_ProfileX_;
  
    TH2D * TH2Sng_Left;
    TH2D * TH2Bkg_Left;
    TH2D * TH2Data_Left;
    TString titleSng_Left_;
    TString titleBkg_Left_;
    TString titleData_Left_;
    
    TH2D * TH2Sng_Right;
    TH2D * TH2Bkg_Right;
    TH2D * TH2Data_Right;
    TString titleSng_Right_;
    TString titleBkg_Right_;
    TString titleData_Right_;
    
    TH2D * TH2Sng_LeftControl;
    TH2D * TH2Bkg_LeftControl;
    TH2D * TH2Data_LeftControl;
    TString titleSng_LeftControl_;
    TString titleBkg_LeftControl_;
    TString titleData_LeftControl_;
    
    TH2D * TH2Sng_RightControl;
    TH2D * TH2Bkg_RightControl;
    TH2D * TH2Data_RightControl;
    TString titleSng_RightControl_;
    TString titleBkg_RightControl_;
    TString titleData_RightControl_;

  TH2D * TH2SngVar12_Left; //var12 stands for var 1 versus var 2
  TH2D * TH2BkgVar12_Left;
  TH2D * TH2DataVar12_Left;
  TString titleSngVar12_Left_;
  TString titleBkgVar12_Left_;
  TString titleDataVar12_Left_;
	
  TH2D * TH2SngVar12_Right;
  TH2D * TH2BkgVar12_Right;
  TH2D * TH2DataVar12_Right;
  TString titleSngVar12_Right_;
  TString titleBkgVar12_Right_;
  TString titleDataVar12_Right_;

  TH1D * TH1Sng_Var0;
  TH1D * TH1Bkg_Var0;
  TH1D * TH1Bkg_W_Var0;
  TH1D * TH1Bkg_R_Var0;
    TH1D * TH1Data_Var0;
    TH1D * TH1Data_Var0_XS;
  TString titleSng_Var0_;
  TString titleBkg_Var0_;
  TString titleBkg_W_Var0_;
  TString titleBkg_R_Var0_;
    TString titleData_Var0_;
    TString titleData_Var0_XS_;
	
	TH1D * TH1Data_Pt_Left;
	TH1D * TH1Sng_Pt_Left;
	TH1D * TH1Data_Pt_Right;
	TH1D * TH1Data_Pt_RightReweigh;
	TH1D * TH1Data_Pt_Control;
	TH1D * TH1Data_Pt_ControlReweigh;
    TH1D * TH1Data_Eta_Control;
	TH1D * TH1Data_Eta_ControlReweigh;
    
	TString titleData_Pt_Left;
	TString titleSng_Pt_Left;
	TString titleData_Pt_Right;
	TString titleData_Pt_RightReweigh;
	TString titleData_Pt_Control;
	TString titleData_Pt_ControlReweigh;
	TString titleData_Eta_Control;
	TString titleData_Eta_ControlReweigh;


  TString titleSng_chisq_;
  TString titleBkg_chisq_;
  TString titleData_chisq_;
  TString titleSng_chisqControl_;
  TString titleBkg_chisqControl_;
  TString titleData_chisqControl_;
  TH1D *TH1Sng_chisq;
  TH1D *TH1Bkg_chisq;
  TH1D *TH1Data_chisq;
  TH1D *TH1Sng_chisqControl;
  TH1D *TH1Bkg_chisqControl;
  TH1D *TH1Data_chisqControl;

  TH1D * TH1SoverSB_Var0;
  TString titleSoverSB_Var0_;
  TH1D * TH1SoverSB_ControlVar;
  TString titleSoverSB_ControlVar_;

  TH1D * TH1Sng_BtagAll;
  TH1D * TH1Bkg_BtagAll;
  TH1D * TH1Data_BtagAll;
  TString titleSng_BtagAll_;
  TString titleBkg_BtagAll_;
  TString titleData_BtagAll_;
  TH1D * TH1Sng_BtagEffAll;
  TH1D * TH1Bkg_BtagEffAll;
  TH1D * TH1Data_BtagEffAll;
  TString titleSng_BtagEffAll_;
  TString titleBkg_BtagEffAll_;
  TString titleData_BtagEffAll_;

  TH1D * TH1Sng_ControlVar;
  TH1D * TH1Bkg_ControlVar;
  TH1D * TH1Bkg_W_ControlVar;
  TH1D * TH1Bkg_R_ControlVar;
  TH1D * TH1Data_ControlVar;
  TString titleSng_ControlVar_;
  TString titleBkg_ControlVar_;
  TString titleBkg_W_ControlVar_;
  TString titleBkg_R_ControlVar_;
  TString titleData_ControlVar_;

  //plots for the SC reweiging
  TH2D * TH2Sng_Var12;
  TH2D * TH2Bkg_Var12;
  TH2D * TH2Bkg_W_Var12;
  TH2D * TH2Bkg_R_Var12;
  TH2D * TH2Data_Var12;
  TString titleSng_Var12_;
  TString titleBkg_Var12_;
  TString titleBkg_W_Var12_;
  TString titleBkg_R_Var12_;
  TString titleData_Var12_;
  TH1D * TH2Sng_Var12_1DX;
  TH1D * TH2Bkg_Var12_1DX;
  TH1D * TH2Bkg_W_Var12_1DX;
  TH1D * TH2Bkg_R_Var12_1DX;
  TH1D * TH2Data_Var12_1DX;
  TString titleSng_Var12_1DX_;
  TString titleBkg_Var12_1DX_;
  TString titleBkg_W_Var12_1DX_;
  TString titleBkg_R_Var12_1DX_;
  TString titleData_Var12_1DX_;
  TH1D * TH2Sng_Var12_1DY;
  TH1D * TH2Bkg_Var12_1DY;
  TH1D * TH2Bkg_W_Var12_1DY;
  TH1D * TH2Bkg_R_Var12_1DY;
  TH1D * TH2Data_Var12_1DY;
  TString titleSng_Var12_1DY_;
  TString titleBkg_Var12_1DY_;
  TString titleBkg_W_Var12_1DY_;
  TString titleBkg_R_Var12_1DY_;
  TString titleData_Var12_1DY_;

  TH2D * TH2Sng_ControlVar12;
  TH2D * TH2Bkg_ControlVar12;
  TH2D * TH2Bkg_W_ControlVar12;
  TH2D * TH2Bkg_R_ControlVar12;
  TH2D * TH2Data_ControlVar12;
  TString titleSng_ControlVar12_;
  TString titleBkg_ControlVar12_;
  TString titleBkg_W_ControlVar12_;
  TString titleBkg_R_ControlVar12_;
  TString titleData_ControlVar12_;
  TH1D * TH2Sng_ControlVar12_1DX;
  TH1D * TH2Bkg_ControlVar12_1DX;
  TH1D * TH2Bkg_W_ControlVar12_1DX;
  TH1D * TH2Bkg_R_ControlVar12_1DX;
  TH1D * TH2Data_ControlVar12_1DX;
  TString titleSng_ControlVar12_1DX_;
  TString titleBkg_ControlVar12_1DX_;
  TString titleBkg_W_ControlVar12_1DX_;
  TString titleBkg_R_ControlVar12_1DX_;
  TString titleData_ControlVar12_1DX_;
  TH1D * TH2Sng_ControlVar12_1DY;
  TH1D * TH2Bkg_ControlVar12_1DY;
  TH1D * TH2Bkg_W_ControlVar12_1DY;
  TH1D * TH2Bkg_R_ControlVar12_1DY;
  TH1D * TH2Data_ControlVar12_1DY;
  TString titleSng_ControlVar12_1DY_;
  TString titleBkg_ControlVar12_1DY_;
  TString titleBkg_W_ControlVar12_1DY_;
  TString titleBkg_R_ControlVar12_1DY_;
  TString titleData_ControlVar12_1DY_;

  TH2D * TH2Data_SignalControlVar12;
  TString titleData_SignalControlVar12_; 
  TH1D * TH2Data_SignalControlVar12_1DX;
  TString titleData_SignalControlVar12_1DX_; 
  TH1D * TH2Data_SignalControlVar12_1DY;
  TString titleData_SignalControlVar12_1DY_; 

  //plots for btag efficiency
  TH1D * TH1Sng_BtagEffMC;
  TH1D * TH1Bkg_BtagEffMC;
  TH1D * TH1Data_BtagEffMC;
  TString titleSng_BtagEffMC_;
  TString titleBkg_BtagEffMC_;
  TString titleData_BtagEffMC_;
  
  //plots for 1D of pt/btag plots
    TH1D * TH2Sng_Left1DX;
    TH1D * TH2Bkg_Left1DX;
    TH1D * TH2Data_Left1DX;
    TString titleSng_Left1DX_;
    TString titleBkg_Left1DX_;
    TString titleData_Left1DX_;
    
    TH1D * TH2Sng_Right1DX;
    TH1D * TH2Bkg_Right1DX;
    TH1D * TH2Data_Right1DX;
    TString titleSng_Right1DX_;
    TString titleBkg_Right1DX_;
    TString titleData_Right1DX_;
    
    TH1D * TH2Sng_Left1DXControl;
    TH1D * TH2Bkg_Left1DXControl;
    TH1D * TH2Data_Left1DXControl;
    TString titleSng_Left1DXControl_;
    TString titleBkg_Left1DXControl_;
    TString titleData_Left1DXControl_;
    
    TH1D * TH2Sng_Right1DXControl;
    TH1D * TH2Bkg_Right1DXControl;
    TH1D * TH2Data_Right1DXControl;
    TString titleSng_Right1DXControl_;
    TString titleBkg_Right1DXControl_;
    TString titleData_Right1DXControl_;
 
  TH1D * TH2SngVar12_Left1DY;
  TH1D * TH2BkgVar12_Left1DY;
  TH1D * TH2DataVar12_Left1DY;
  TString titleSngVar12_Left1DY_;
  TString titleBkgVar12_Left1DY_;
  TString titleDataVar12_Left1DY_;

  TH1D * TH2SngVar12_Right1DY;
  TH1D * TH2BkgVar12_Right1DY;
  TH1D * TH2DataVar12_Right1DY;
  TString titleSngVar12_Right1DY_;
  TString titleBkgVar12_Right1DY_;
  TString titleDataVar12_Right1DY_;

  TH1D * TH2SngVar12_Right1DYReweigh;
  TH1D * TH2BkgVar12_Right1DYReweigh;
  TH1D * TH2DataVar12_Right1DYReweigh;
  TString titleSngVar12_Right1DYReweigh_;
  TString titleBkgVar12_Right1DYReweigh_;
  TString titleDataVar12_Right1DYReweigh_;

  TH1D * TH2SngVar12_Right1DYReweighRatio;
  TH1D * TH2BkgVar12_Right1DYReweighRatio;
  TH1D * TH2DataVar12_Right1DYReweighRatio;
  TString titleSngVar12_Right1DYReweighRatio_;
  TString titleBkgVar12_Right1DYReweighRatio_;
  TString titleDataVar12_Right1DYReweighRatio_;

  TH1D * TH2Data_RightLeft1DX;
    TH1D * TH2Data_LeftRight1DX;
    TString titleData_RightLeft1DX_;
    TString titleData_LeftRight1DX_;
    TF1 * TFData_LeftRight1DXFit;
    TString titleData_LeftRight1DXFit_;
    
    TH1D * TH2Data_LeftRight1DXControl;
    TString titleData_LeftRight1DXControl_;
    TF1 * TFData_LeftRight1DXControlFit;
    TString titleData_LeftRight1DXControlFit_;

  TH2D * TH2DataVar12_LeftRight;
  TString titleDataVar12_LeftRight_; 

  TH1D * TH1Sng_ControlVarReweigh;
  TH1D * TH1Bkg_ControlVarReweigh;
  TH1D * TH1Bkg_W_ControlVarReweigh;
  TH1D * TH1Bkg_R_ControlVarReweigh;
  TH1D * TH1Data_ControlVarReweigh;
  TString titleSng_ControlVarReweigh_;
  TString titleBkg_ControlVarReweigh_;
  TString titleBkg_W_ControlVarReweigh_;
  TString titleBkg_R_ControlVarReweigh_;
  TString titleData_ControlVarReweigh_;

  TH1D * TH2Data_Left1DY;
  TH1D * TH2Sng_Left1DY;
  TH1D * TH2Bkg_Left1DY;
  TString titleData_Left1DY_;
  TString titleSng_Left1DY_;
  TString titleBkg_Left1DY_;
  
  TH1D * TH2Data_Right1DY;
  TH1D * TH2Sng_Right1DY;
  TH1D * TH2Bkg_Right1DY;
  TString titleData_Right1DY_;
  TString titleSng_Right1DY_;
  TString titleBkg_Right1DY_;

  //  TH1D * TH2Data_Left1DYReweigh;
  //TString titleData_Left1DYReweigh_;

  TH1D * TH2Bkg_Right1DYReweigh;
  TH1D * TH2Sng_Right1DYReweigh;
  TH1D * TH2Data_Right1DYReweigh;
  TString titleBkg_Right1DYReweigh_;
  TString titleSng_Right1DYReweigh_;
  TString titleData_Right1DYReweigh_;
    
  TH1D * TH2Bkg_Right1DYReweighRatio;
  TH1D * TH2Sng_Right1DYReweighRatio;
  TH1D * TH2Data_Right1DYReweighRatio;
  TString titleBkg_Right1DYReweighRatio_;
  TString titleSng_Right1DYReweighRatio_;
  TString titleData_Right1DYReweighRatio_;
 
  TH1D * TH2Bkg_Right1DY_LeftRightRatio;
  TH1D * TH2Sng_Right1DY_LeftRightRatio;
  TH1D * TH2Data_Right1DY_LeftRightRatio;
  TString titleBkg_Right1DY_LeftRightRatio_;
  TString titleSng_Right1DY_LeftRightRatio_;
  TString titleData_Right1DY_LeftRightRatio_;
  
  TH1D * TH2Bkg_Right1DYReweigh_LeftRightRatio;
  TH1D * TH2Sng_Right1DYReweigh_LeftRightRatio;
  TH1D * TH2Data_Right1DYReweigh_LeftRightRatio;
  TString titleBkg_Right1DYReweigh_LeftRightRatio_;
  TString titleSng_Right1DYReweigh_LeftRightRatio_;
  TString titleData_Right1DYReweigh_LeftRightRatio_;

  TH1D * TH2Bkg_Left1DYReweigh;
  TH1D * TH2Sng_Left1DYReweigh;
  TH1D * TH2Data_Left1DYReweigh;
  TString titleBkg_Left1DYReweigh_;
  TString titleSng_Left1DYReweigh_;
  TString titleData_Left1DYReweigh_;
    
  TH1D * TH2Bkg_Left1DYReweighRatio;
  TH1D * TH2Sng_Left1DYReweighRatio;
  TH1D * TH2Data_Left1DYReweighRatio;
  TString titleBkg_Left1DYReweighRatio_;
  TString titleSng_Left1DYReweighRatio_;
  TString titleData_Left1DYReweighRatio_;
 
  TH1D * TH2Bkg_Left1DYReweigh_LeftRightRatio;
  TH1D * TH2Sng_Left1DYReweigh_LeftRightRatio;
  TH1D * TH2Data_Left1DYReweigh_LeftRightRatio;
  TString titleBkg_Left1DYReweigh_LeftRightRatio_;
  TString titleSng_Left1DYReweigh_LeftRightRatio_;
  TString titleData_Left1DYReweigh_LeftRightRatio_;
    
    TH1D* TH2Bkg_Left1DYControl;
    TH1D* TH2Sng_Left1DYControl;
    TH1D* TH2Data_Left1DYControl;
    TString titleBkg_Left1DYControl_;
    TString titleSng_Left1DYControl_;
    TString titleData_Left1DYControl_;

    TH1D* TH2Bkg_Left1DYControlReweigh;
    TH1D* TH2Sng_Left1DYControlReweigh;
    TH1D* TH2Data_Left1DYControlReweigh;
    TString titleBkg_Left1DYControlReweigh_;
    TString titleSng_Left1DYControlReweigh_;
    TString titleData_Left1DYControlReweigh_;

    TH1D* TH2Bkg_Right1DYControl;
    TH1D* TH2Sng_Right1DYControl;
    TH1D* TH2Data_Right1DYControl;
    TString titleBkg_Right1DYControl_;
    TString titleSng_Right1DYControl_;
    TString titleData_Right1DYControl_;
    
    TH1D* TH2Bkg_Right1DYControlReweigh;
    TH1D* TH2Sng_Right1DYControlReweigh;
    TH1D* TH2Data_Right1DYControlReweigh;
    TString titleBkg_Right1DYControlReweigh_;
    TString titleSng_Right1DYControlReweigh_;
    TString titleData_Right1DYControlReweigh_;


  int lowerBin_;
  int upperBin_;
  int centralLowBin_;
  int centralUpBin_;
  double leftNonB_;
  double rightNonB_;
  double leftNonBMC_;
  double rightNonBMC_;
  double leftNonBReweigh_;
  double rightNonBReweigh_;
  double leftNonBe_;
  double rightNonBe_;
  double leftNonBMCe_;
  double rightNonBMCe_;
  double leftNonBReweighe_;
  double rightNonBReweighe_;
  double F_;
  double FMC_;
  double FReweigh_;
  double Fe_;
  double FMCe_;
  double FReweighe_;

  TH1D * TH1Data_BtagMeasured;
  TString titleData_BtagMeasured_;
  TH1D * TH1Data_BtagMCMeasured;
  TString titleData_BtagMCMeasured_;
  
  TH1D * TH1Data_BtagEffMeasured;
  TString titleData_BtagEffMeasured_;
  TH1D * TH1Data_BtagEffMCMeasured;
  TString titleData_BtagEffMCMeasured_;
 
  TH1D * TH1Data_BtagEffMeasuredDiff;
  TString titleData_BtagEffMeasuredDiff_;
  TH1D * TH1Data_BtagEffMCMeasuredDiff;
  TString titleData_BtagEffMCMeasuredDiff_;
 
  TH1D * TH1Data_BtagMeasuredLR;//LR aka Left region b-tag plot Reweighted
  TString titleData_BtagMeasuredLR_;
  TH1D * TH1Data_BtagMCMeasuredLR;
  TString titleData_BtagMCMeasuredLR_;
  
  TH1D * TH1Data_BtagEffMeasuredLR;
  TString titleData_BtagEffMeasuredLR_;
  TH1D * TH1Data_BtagEffMCMeasuredLR;
  TString titleData_BtagEffMCMeasuredLR_;  

  TH1D * TH1Data_BtagEffMeasuredLRDiff;
  TString titleData_BtagEffMeasuredLRDiff_;
  TH1D * TH1Data_BtagEffMCMeasuredLRDiff;
  TString titleData_BtagEffMCMeasuredLRDiff_;

    TH1D * TH1Data_BtagMeasuredRR;//RR aka Right region b-tag plot Reweighted
    TString titleData_BtagMeasuredRR_;
    TH1D * TH1Data_BtagShapeMC_MeasuredRR;
    TString titleData_BtagShapeMC_MeasuredRR_;
    
    TH1D * TH1Data_BtagMCMeasuredRR;
    TString titleData_BtagMCMeasuredRR_;
    TH1D * TH1Data_BtagMC_ShapeMC_MeasuredRR;
    TString titleData_BtagMC_ShapeMC_MeasuredRR_;
    
    TH1D * TH1Data_BtagEffMeasuredRR;
    TString titleData_BtagEffMeasuredRR_;
    TH1D * TH1Data_BtagEffShapeMC_MeasuredRR;
    TString titleData_BtagEffShapeMC_MeasuredRR_;
    
    TH1D * TH1Data_BtagEffMCMeasuredRR;
    TString titleData_BtagEffMCMeasuredRR_;
    TH1D * TH1Data_BtagEffMC_ShapeMC_MeasuredRR;
    TString titleData_BtagEffMC_ShapeMC_MeasuredRR_;
    
    TH1D * TH1Data_BtagEffMeasuredRRDiff;
    TString titleData_BtagEffMeasuredRRDiff_;
    TH1D * TH1Data_BtagEffShapeMC_MeasuredRRDiff;
    TString titleData_BtagEffShapeMC_MeasuredRRDiff_;
    
    TH1D * TH1Data_BtagEffMCMeasuredRRDiff;
    TString titleData_BtagEffMCMeasuredRRDiff_;
    TH1D * TH1Data_BtagEffMC_ShapeMC_MeasuredRRDiff;
    TString titleData_BtagEffMC_ShapeMC_MeasuredRRDiff_;
	
    // mistag
    TH1D * TH1Data_MistagMeasuredRR;//RR aka Right region b-tag plot Reweighted
    TString titleData_MistagMeasuredRR_;
    TH1D * TH1Data_MistagShapeMC_MeasuredRR;
    TString titleData_MistagShapeMC_MeasuredRR_;
    
    TH1D * TH1Data_MistagEffMeasuredRR;
    TString titleData_MistagEffMeasuredRR_;
    TH1D * TH1Data_MistagEffShapeMC_MeasuredRR;
    TString titleData_MistagEffShapeMC_MeasuredRR_;
        
    TH1D * TH1Data_MistagEffMeasuredRRDiff;
    TString titleData_MistagEffMeasuredRRDiff_;
    TH1D * TH1Data_MistagEffShapeMC_MeasuredRRDiff;
    TString titleData_MistagEffShapeMC_MeasuredRRDiff_;
    

	std::map<std::string, float> EffCalcDetails_;
	
	vector<string> FitPlotPaths;
	
	float nTTbar_;

    int fitMode; // template fit variable 0: mlj 1: M3 2: 2D (mlj,m3)

	string data_postfix_;
    
    int nSystematic_;
};


#endif
