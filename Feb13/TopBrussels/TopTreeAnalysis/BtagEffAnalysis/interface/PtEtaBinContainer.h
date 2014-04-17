#ifndef PtEtaBinContainer_h
#define PtEtaBinContainer_h

// system include files
#include <iostream>
#include <string>

#include <TFile.h>
#include <TH1.h>
#include <TProfile.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TBranch.h>
#include <TTree.h>
#include <sstream>
#include <vector>
#include <TCanvas.h>

#include "TopTreeProducer/interface/TRootJet.h"
#include "../interface/PtEtaBin.h"

class PtEtaBinContainer{

 public:
  PtEtaBinContainer(){};
  PtEtaBinContainer(int,int,int,int,bool,bool,bool,int,int,int); 
  //~PtEtaBinContainer(){}; 
  ~PtEtaBinContainer(); 

  void DefineSignalSamplePlots(int, double, double, int, double, double, int, double, double, int, double, double, int, double, double, int[], double[], double[]);
  void DefineControlSamplePlots(int, double, double,int, double, double,int, double, double,int[], double[], double[]);
  //void FillSignalSamplePlots(double, int, double*, double, double, double, int, double, double, double, double, double, double, double);
  void FillSignalSamplePlots(double, double, int, bool, bool, double, double*, std::map<int,vector<double> > WPMap, double, double, double, double, double, double, double,double);
  void FillControlSamplePlots(double, int, bool, bool, double, double*, double, double, double, double, double, double, double);

  void FillXStemplates(double weight, string datasetname, int partonFlavour, double* btag,std::map<int,vector<double> > WPMap, double controlVar0, double m3, double lowCutVar0, double centralLowCutVar0, double centralUpCutVar0, double upCutVar0);

	void SetVarBins(std::map<int,vector<float> > rangesbTag);
    
    void setLumi(double lumi);

  void WriteContainerToRootFile(TFile *,bool,bool,bool,bool,bool,bool,bool,bool,bool,bool,bool,bool);
  void WriteHistoToPSFile(TString *);
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

  void MakeXRatioPlot(bool);
  void ChangeLeftRightPars();
  void Make2DRatioPlot();
  void MakeSCprojectionPlots();
  void MakeSCVar12RatioPlot();
  void ReweighLeft();
  void ReweighRight();   
  //void ReweighRightChangeFitParams(double, double); //obsolete
  //void FillReweighRight(bool, bool, double, int, double*, double, double, double, int, double, double, double, double, double, double, double);
  void FillReweighRight(bool, bool, double, int, double*, double, double, double, double, double, double, double);
  //void FillReweighControl(bool, bool, double, int, double*, double, double, double, int, double, double, double, double, double, double, double);
  void FillReweighControl(double*,double*,bool, double, int, bool, bool, double, double, double, double, double, double, double, double);
  void MakeReweighRatio();
  void GetLRratio(double, bool, double*, double*);
  void GetPercentiles(double,double,double,int*);
  void MeasureEff(bool);
  void MeasureEffLR(bool);
  void MeasureEffRR(bool);
    void MeasureMistagEffRR(bool);

  void GetWPEff(bool, bool, double, int, double*, bool, bool, int, int, bool, double);
  void CoutWPEff(bool, bool, double, int, double*, bool, bool, int, int, bool, double);
	
	std::map<std::string,float> GetEffCalcDetails(float wp, int tagger, int ptbin, int etabin, bool ptetabin);
	
	std::vector<string> GetFitPlotPaths(int tagger);

	std::map<int,vector<double> > doMLJTemplateFit(string chi2cut,int mode, string data_postfix="",int nSystematic=0);

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

 private:
  int debug_;
  int nVar1_;
  int nVar0_;
  int nBdiscrAlgos_;
  
  
  //double ptBins_[9];
  //double etaBins_[9];
  double ptBins_[6];
  double etaBins_[6];
  int nPtBins_;
  int nEtaBins_;

  //double ptBins_[21];
  //double etaBins_[2];

  //double ptBins_[1];
  //double etaBins_[21];

 
  double ptUpCut_;
  double ptLowCut_;

  double ptUpCutComplement_;
  double ptLowCutComplement_;

  
  PtEtaBin ** binGlobal_;
  //PtEtaBin ** binVector_[10]; //pt5 + eta5
  //PtEtaBin ** binVector_[16]; //pt8 + eta8
  PtEtaBin ** binVector_[21]; //20+1 or 1+20
  
  /* void FillContainer(float, TRootJet *, float, float);   
  void WriteContainerToRootFile(TFile *, int, int, bool);

  void SetSumw2();
  void MakeSons();
  void ReBinTH1(int, int, bool);

 private:  
  int nBdiscrAlgos_;
  int nVar_;
  int nBinsX_; 
  float lowRangeX_;
  float upRangeX_;
  int nBinsY_; 
  float lowRangeY_;
  float upRangeY_;
  int nDaughter_;

  PtEtaBin binVector_[5][2][1]; //now the number of bins and the number of btaggers is hardcoded, this could change
                                 //dimension 1: pt bin, dim 2: eta bin, dim 3: btagger
  PtEtaBin binGlobal_[1];
  float ptBins_[4];//This can be done with a vector of pairs (lower limit, upperlimit)
  int nPtBins_;
  float etaBins_[3];
  int nEtaBins_;*/


};



#endif
