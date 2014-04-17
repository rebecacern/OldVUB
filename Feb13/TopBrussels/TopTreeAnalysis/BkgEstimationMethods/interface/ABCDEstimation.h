#ifndef ABCDEstimation_h
#define ABCDEstimation_h

#include <iostream>
#include <vector>
#include <iomanip>
#include <Math/VectorUtil.h>

#include "TFile.h"
  //#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TString.h"

#include "TopTreeAnalysis/Tools/interface/PlottingTools.h"

using namespace std;

/**
	//	Aim: Estimation of the number of QCD bkg events in our Signal Region with ABCD method.
	//
	//	It takes 2 variables are compute the usual ABCD method.\n
	//	This was extended to work as function of a 3rd variable like Mttbar or MET in order to obtain a estimation of the shape of this variable for QCD.\n
	//	The methods/variables "*3D" or "*Z" are refering to this extension.\n
	//
	//	Example of usage:
	//	- call constructor
	//	- Fill() in a loop over events
	//	- ComputeEstimate()
	//	- Print2DEstimate() (can be called later or not called at all)
	//	- Draw() 
	//	- Write()
	//
	//	Different modes:
	//	- Fill() or Fill3D() in loop over events
	//	- Use Contructor(TH2F*) or Constructor(TH3F*) to take directly the histo
	// 		without looping over events
	//
	//	All methods containing QCD refer to method where only the data flagged as QCD
	//	are considered. It's a MC approach to test the performance of the method without
	//	being spoiled by other processus.

*/

class ABCDEstimation {
  
 public:
  
  ABCDEstimation();
  /** default constructor */
  ABCDEstimation(TString hName, TString hTitle, TString XaxisTitle, TString YaxisTitle, int NXbins, float Xmin, float Xmax, int NYbins, float Ymin, float Ymax);
  /** initialize 1D-2D histograms */
  ABCDEstimation(TString hName, TString hTitle, TString XaxisTitle, TString YaxisTitle, int NXbins, float Xmin, float Xmax, int NYbins, float Ymin, float Ymax, int NZbins, float Zmin, float Zmax);
  /**initialize 1D-2D-3D histograms */
  ABCDEstimation(TString hName, TString hTitle, TString XaxisTitle, TString YaxisTitle, int NXbins, float *binsX, int NYbins, float *binsY, int NZbins, float *binsZ);
  /**initialize 1D-2D-3D histograms with varying binning for Z axis */
  ABCDEstimation(TH2F *histo, TH2F *histo_QCD = 0);
  /**initialize directly with a 2D histo (no need to use Fill() method )*/
  ABCDEstimation(TH3F *histo, TH3F *histo_QCD = 0);
  /**initialize directly with a 3D histo (no need to use Fill3D() method )*/
  ~ABCDEstimation();
  /** destructor*/

  ABCDEstimation(const ABCDEstimation&);
  /**copy constructor*/

  void Fill(float X, float Y, float weight = 1., bool isQCD = false);
  /* X and Y are the 2 uncorrelated variables filled event by event. The event can have a weight. The boolean isQCD is for MC-matching. */
  void Fill3D(float X, float Y, float Z, float weight = 1., bool isQCD = false);
  /* X and Y are the 2 uncorrealted variables filled event by event, Z being the 3rd variable. The event can have a weight. The boolean isQCD is for MC-matching. */
	
  void Draw(bool dofit = false);
  /** if dofit==true,   will be fitted */
  void Write(TFile* fout, string label = string(""));
  /** string could be the name of the dataset, the name of the 3rd variable or wathever. */

  ////////////////////////////
  //  Access to histo
  ////////////////////////////
  //all (data)
  TH1F* GetShapePrediction()    const { return hShapePred_; }; /** Access to the predicted shape with computed errors*/
  TH2F* Get2DHistogram()        const { return h2D_; }; /** Access to the 2D histogram (X,Y)*/
  TH3F* Get3DHistogram()        const { return h3D_; }; /** Access to the 3D histogram (X,Y,Z) */
  TGraphAsymmErrors* GetXProf() const { return hXProf_; }; /** Stability of the prediction while varying the cut on X axis */
  TGraphAsymmErrors* GetYProf() const { return hYProf_; }; /** Stability of the prediction while varying the cut on Y axis */

  //QCD only (MC)
  TH1F* GetShapePredictionQCD()    const { return hShapePredQCD_; }; /** Access to the predicted shape with computed errors for QCD events only (MC)*/
  TH2F* Get2DHistogramQCD()        const { return h2DQCD_; }; /** Access to the 2D histogram (X,Y) for QCD events only (MC)*/
  TH3F* Get3DHistogramQCD()        const { return h3DQCD_; }; /** Access to the 3D histogram (X,Y,Z) for QCD events only (MC)*/
  TGraphAsymmErrors* GetXProfQCD() const { return hXProfQCD_; }; /** Stability of the prediction while varying the cut on X axis for QCD events only (MC)*/
  TGraphAsymmErrors* GetYProfQCD() const { return hYProfQCD_; }; /** Stability of the prediction while varying the cut on Y axis for QCD events only (MC)*/

  ////////////////////////////
  //  Set histo
  ////////////////////////////
  void Set2DHistogram(TH2F* histo) { h2D_ = histo; };
  void Set3DHistogram(TH3F* histo) { h3D_ = histo; };
  void Set1DHistogramQCD(TH1F* histo) { h1DQCD_ = histo; };
  void Set2DHistogramQCD(TH2F* histo) { h2DQCD_ = histo; };
  void Set3DHistogramQCD(TH3F* histo) { h3DQCD_ = histo; };

  void SetHistogramName(TString name) { Name_ = name; };
  TString GetHistogramName(TString name) { return (Name_.IsNull() ? "DefaultName" : Name_); };

  void SetLumi(float lumi) {Lumi_ = lumi;};

  float GetCorrelation() const  { return h2D_->GetCorrelationFactor(); }; /** Correlation for X,Y variables*/ 
  float GetCorrelationQCD() const  { return h2DQCD_->GetCorrelationFactor(); }; /** Correlation for X,Y variables for QCD events only (MC)*/
  
  /////////////////////////////////////////////////////////////////////////
  //  "Main" method which compute the ABCD estimation for all the regions
  /////////////////////////////////////////////////////////////////////////
  void ComputeEstimate(float Xmin, float X0, float X1, float Xmax, float Ymin, float Y0, float Y1, float Ymax, int region);
  /** Range for the 4 regions are:\n 
  //1:[Xmin-X0,Ymin-Y0] \n
  //2:[X1-Xmax,Ymin-Y0] \n
  //3:[Xmin-X0,Y1-Ymax] \n
  //4:[X1-Xmax,Y1-Ymax] \n 
  //region is the number of one of this 4 regions*/

  void Print2DEstimate();/** Print results of the estimation */
  void Print3DEstimate();/** Print results of the estimation in case of 3rd variable */

  float GetQCDData(int region) const { return (region>0 && region<5 ? estimateQCD_[region-1][0] : -1);};/** Nof of QCD events (MC) in a given region*/
  float GetQCDDataError(int region) const { return (region>0 && region<5 ? estimateQCD_[region-1][1] : -1);};/** statistical error associated*/
  
  float GetData(int region) const { return (region>0 && region<5 ? estimate_[region-1][0] : -1);};/** Nof of events from data in a given region*/
  float GetDataError(int region) const { return (region>0 && region<5 ? estimate_[region-1][1] : -1);};/** statistical error associated*/

  //only taking QCD into account (MC)
  float GetQCDEstimate(int region)  const {return (region>0 && region<5 ? estimateQCD_[region-1][2] : -1);};
  float GetQCDEstimateError(int region)  const {return (region>0 && region<5 ? estimateQCD_[region-1][3] : -1);};
  
  float GetQCDEstimateZ(int bin) const {return (bin<h3DQCD_->GetNbinsZ() ? estimateQCDZ_[bin][2] : -1);};
  float GetQCDEstimateZError(int bin) const {return (bin<h3DQCD_->GetNbinsZ() ? estimateQCDZ_[bin][3] : -1);};
  
  //taking all data into account
  float GetEstimate(int region)  const {return (region>0 && region<5 ? estimate_[region-1][2] : -1);};
  float GetEstimateError(int region)  const {return (region>0 && region<5 ? estimate_[region-1][3] : -1);};

  float GetEstimateZ(int bin) const {return (bin<h3D_->GetNbinsZ() ? estimateZ_[bin][2] : -1);};
  float GetEstimateZError(int bin) const {return (bin<h3D_->GetNbinsZ() ? estimateZ_[bin][3] : -1);};

  ABCDEstimation& operator +=(const ABCDEstimation& );
  /** Will add 1D & 2D & 3D plots, but take care while using it !
  // ComputeEstimate(), Draw().. have to be called after summing */
  
  void ReScale(float factor);/** Rescale the events in the 4 regions by a certain factor, it's equivalent to changing the luminosity*/


		//private methods

  /////////////////////////////////////////////////
  //  "Hidden methods" called by ComputeEstimate()
  /////////////////////////////////////////////////

  void  ABCDEstimator(float a, float b, float c, float d, float** tab);//tab can be estimate_ or estimateQCD_ 
  void  ABCDEstimator(float a, float b, float c, float d, int region, int iBinZ, float** tabZ);
  
  float ABCDErrorCalculator(float a, float b, float c, float d, int region);
  float ABCDRelErrorCalculator(float a, float b, float c, float d, int region);

  //void  ABCDErrorCalc(float a, float b, float c, float d, float** tab);//tab can be estimate_ or estimateQCD_ 
  //void  ABCDErrorCalc(float a, float b, float c, float d, int region, int iBinZ, float** tabZ); 
  // for a bin in Z variable - tabZ can be estimateZ_ or estimateQCDZ_

  private: 

  float    CoVar_XY_Z(float X, float Y, float Z, float Ntot);
  float	     CoVar_XY(float X, float Y, float Ntot);
  float      Sigma_XY(float X, float Y, float Ntot);
  float       Sigma_X(float X, float Ntot);
  float    Sigma_EffX(float X, float Ntot);

  float WilsonScoreIntervalHigh(float Non, float Ntot);
  float WilsonScoreIntervalLow( float Non, float Ntot);
		//private data
 private:
  
  float Xmin_;  
  float X0_;
  float X1_;
  float Xmax_;

  float Ymin_;
  float Y0_;
  float Y1_;
  float Ymax_;

  TString Name_;
  TString XaxisTitle_;
  TString YaxisTitle_;

  //for QCD (MC based)
  TH1F* h1DQCD_; /**shape expected - correspond to variable on Z axis of the 3D plot*/
  TH2F* h2DQCD_; /**plot used to perform ABCD method*/
  TH3F* h3DQCD_; /**plot 2DxVariable*/
  TH1F* hShapePredQCD_; /**shape predicted*/
  TGraphAsymmErrors* hXProfQCD_;
  TGraphAsymmErrors* hYProfQCD_;
  /* TCanvas* cXProfQCD_; */
  /* TCanvas* cYProfQCD_; */

  // all (data based)
  TH2F* h2D_; /**plot used to perform ABCD method*/
  TH3F* h3D_; /**plot 2DxVariable*/
  TH1F* hShapePred_; /**shape predicted*/
  TGraphAsymmErrors* hXProf_;
  TGraphAsymmErrors* hYProf_;
  /* TCanvas* cXProf_; */
  /* TCanvas* cYProf_; */

  TH1F* hEstim_vs_Xcut;
  TH1F* hEstim_vs_Ycut;

  int region_; /**bin of the Z axis*/
  float** estimate_;/**contains the ABCD estimation values & errors associated for each region*/
  float** estimateZ_;
  //idem for QCD only (MC - not contamination)
  float** estimateQCD_;
  float** estimateQCDZ_;

  // Used in Draw() method
  TF1* func;
  TLine* line;
  /* TCanvas* canvasRatioEstim_Exp; */
  /* TCanvas* canvasRatioQCDEstim_Exp; */ 
  TH1F* plotRatioEstim_Exp;
  TH1F* plotRatioQCDEstim_Exp;
  TGraph* hErrorVsLuminosity_;
  TGraph* hErrorVsLuminosityQCD_;
  /* TCanvas* cErrorVsLuminosity_; */
  /* TCanvas* cErrorVsLuminosityQCD_; */
  float Lumi_;/** luminosity */


};

#endif
