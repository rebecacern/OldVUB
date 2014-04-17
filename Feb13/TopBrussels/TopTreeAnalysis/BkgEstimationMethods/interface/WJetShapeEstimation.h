#ifndef WJetShapeEstimation_h
#define WJetShapeEstimation_h


#include <iostream>
#include <iomanip>
#include <cmath>
#include "TH1F.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"

//user
#include "TopTreeAnalysis/Tools/interface/PlottingTools.h"

using namespace std;


        /**
        //
	//	Aim: Estimation of the shape of a variable for W+jets background
	//	using the distribution of this variable for 0 and 1 b-tagged jets bins 
	//	and the number of W & ttbar events estimated with WJetEstimation in both bins.
	//
	//	It's a iteratif processus. 20 iterations should be enough in principle.
	//
	//	Example of usage:
	//	- call the default constructor
	//	- Initialisation(...)
	//	- IterativeHistoSubstraction()
	//	- Scale(..) // take as arguments the #events estimated from WJetEstimation (W,ttbar) 
	//	- PrintResults() // can be called later or not at all
	//	- Draw()
	//	- Write()
	//
	//	Missing pieces:
	//	- histo for QCD expectation is not taken into account and substracted
	*/
	
class WJetShapeEstimation{

	public:
	
		WJetShapeEstimation(TAxis* axis, int NofIterations, string VarLabel = string(""));
		WJetShapeEstimation(int nbins, float xmin, float xmax, int NofIterations, string VarLabel = string(""));
		~WJetShapeEstimation();
		
		//Fill methods
		void Fill(int nofbjets, float variable, float weight = 1.);
		void FillMC(bool isWJets, bool isTtJets, float variable, float weight = 1.);
		void AddMCHitos(TH1F* histo_wjets, TH1F* histo_ttjets);
		
		//Method to call to run the algorithm itself
		void IterativeHistoSubstraction(float RatioNbEvts_Wlike_1bjet, float RatioNbEvts_TTlike_0bjet, float Nevts0bj, float Nevts1bj);

		//Scale method as to be called before Draw, PrintResults and Write
		void Scale(float Ntt, float Nw);
		void Draw(string AxisLabel = string(""));
		void PrintResults();
		void Write(TFile* fout, string label = string(""));
		
		void AddNormalisationError(float factor);//this factor is a fraction %
		void ComputeSystematicError(bool add = true);
		void AddSystematicError(TH1F* SystError); //coming from MC ...
		TH1F* GetSystematicError() {return h_SysError_;};
				

		//access to histograms
		TH1F* GetHisto0bjet() { return histo_0bj_;};
		TH1F* GetHisto1bjet() { return histo_1bj_;};
		TH1F* GetHistoWLikeEstimated() { return histo_Wlike_;};
		TH1F* GetHistoTTLikeEstimated() { return histo_TTlike_;};
		TH1F* GetHistoInclusif() { return histo_inclusif_;};
	
			//private methods
	private:
		void SubstractHisto(TH1F*& histo, TH1F* histoToSubstract);

			//privat data
	private:
	
		int NofIterations_;
		float RatioNbEvts_Wlike_1bjet_;
		float RatioNbEvts_TTlike_0bjet_;
		float Nevts0bj_; 
		float Nevts1bj_;
		//histo for b-jets bins
		TH1F* histo_0bj_;
		TH1F* histo_1bj_;
		//histo estimators
 		TH1F* histo_Wlike_;
		TH1F* histo_TTlike_;
		//histo for estimate eff
 		TH1F* histo_Wlike_EstEff_;
		TH1F* histo_TTlike_EstEff_;
		//reference histo: MC
		TH1F* histo_Wjets_;
		TH1F* histo_ttjets_;
		//histo of reference: inclusive distribution
		TH1F* histo_inclusif_;
		//histo convergence
		TH1F*  histo_Wlike_cvg_;
		TH1F*  histo_TTlike_cvg_;
		//histo for systematic error
		TH1F* h_SysError_;
		// Canvas to hold the histos
		// - comparison MC and estimation
		TCanvas* c_Wjets_;
		TCanvas* c_ttjets_;
		// - comparison and ratio between MC and estimation
		TCanvas* c_Wjets_Split;
		TCanvas* c_ttjets_Split;
		// - comparison MC and N0bjet, N1bjet
		TCanvas* c_Wjets_Nbjet;
		TCanvas* c_ttjets_Nbjet;
		// - comparison and ratio between MC and N0bjet, N1bjet
		TCanvas* c_Wjets_Split_Nbjet;
		TCanvas* c_ttjets_Split_Nbjet;
        // Legend
		TLegend* l_Wjets_;
		TLegend* l_ttjets_;
		TLegend* l_Wjets_Nbjet;
		TLegend* l_ttjets_Nbjet;
};

#endif
