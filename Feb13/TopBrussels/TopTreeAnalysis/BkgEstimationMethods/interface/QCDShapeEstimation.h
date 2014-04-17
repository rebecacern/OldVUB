#ifndef QCDShapeEstimation_h
#define QCDShapeEstimation_h


#include <iostream>
#include <iomanip>
#include <cmath>
#include "TH1F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"

//user
#include "TopTreeAnalysis/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysis/Content/interface/Dataset.h"

using namespace std;

        /**
	//	Aim: Estimation of the shape for a given observable for QCD background using a control region 
	//	and normalized to the number of events estimated with ABCD method.
	//
	//
	//	Example of usage:
	//	- Call the constructor
	//	- SetControlShape
	//	- Fill() in the loop over events
	//	- Normalized() to #events expected
	//	- Draw()
	//	- PrintResults() //can be called later or not called at all
	//	- Write()
	//
	//	The estimation is performed over all data and/or for QCD (MC) only.\n
	//	A comparison between both methods can be done.
	*/

class QCDShapeEstimation{

	
	public:
	
		QCDShapeEstimation();
		QCDShapeEstimation(int Nbins, float Xmin, float Xmax, string xlabel, string varLabel = string(""));
		QCDShapeEstimation(int Nbins, float* binsX, string xlabel, string varLabel = string(""));
		~QCDShapeEstimation();
		
		void SetControlShape(vector<float> Xbins); //bins of relIso
		void SetControlShape(vector<float> Xbins, string); //bins of the chosen variable
		void Fill(float value, float weight, float controlValue, bool isQCD, bool isInSR, bool isInCR); //controlValue = relIso
		void Fill(float value, float weight, bool isQCD, bool isInSR, bool isInCR);
		void Draw(bool dofit = false, string label = string(""));
		void PrintResults();
		void Write(TFile* fout, string label = string(""));
		
		float NofQCDEventsExpected() const { return NofQCDEventsExpected_;};
		void Normalize(float nofQCDExp);
		
		void AddNormalisationError(float factor);//this factor is a fraction %
		void ComputeSystematicError(bool add = true);
		void AddSystematicError(TH1F* SystError); //coming from MC ...
		TH1F* GetSystematicError() {return h_SysError_;};

		//access to histograms
		TH1F* GetQCDMCShapeSR() { return h_QCDMCShape_SR_;};
		TH1F* GetQCDMCShapeCR() { return h_QCDMCShape_CR_;};
		TH1F* GetShapeCR() { return h_Shape_CR_;};
		TH1F* GetQCDShapeEstimated() { return h_QCDShapeEstimated_;};
		

	private:
		float NofQCDEventsExpected_;
		
		TH1F* h_QCDMCShape_SR_;
		TH1F* h_QCDMCShape_CR_;
		TH1F* h_Shape_CR_;
		TH1F* h_QCDShapeEstimated_;
		TH1F* h_SysError_;
		TProfile* h_ControlShapeProfile_;
		vector<pair<float,TH1F*> > h_ControlShape_;
		TLegend* legControlShape_;
		TH1F** h_Projection_; 
		
	
		//for drawing 
		TLine* line;
		TF1 *func;
	        TH1F* plotRatioCR_SR;
		TH1F* plotRatioEstim_Exp;
		TCanvas* canvasRatioCR_SR;
		TCanvas* canvasRatioEstim_Exp;
		TCanvas* canvasControlShape;
	  		
};

#endif
