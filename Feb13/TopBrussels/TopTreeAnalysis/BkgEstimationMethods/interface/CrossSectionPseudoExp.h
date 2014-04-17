#ifndef CrossSectionPseudoExp_h
#define CrossSectionPseudoExp_h

#include <iostream>
#include <iomanip>
#include <vector>

//ROOT
#include "TH1F.h"
#include "TFile.h"

//Our Package
#include "TopTreeAnalysis/Selection/interface/Selection.h"
#include "TopTreeAnalysis/Content/interface/MCExpectation.h"

        /**
	//	Aim: Saving information for each MC Pseudo Experiment performed for the Cross-Section measurement\n
	//	Plots are obtained from this values and write in a root file.\n
	//	
	//	Usage:\n
	//	- Call the constructor
	//	- Fill(...) at the end of the loop over the datasets
	//	- Compute()  at the end of the macro
	//	- Draw()     at the end of the macro
	//	- Write()    at the end of the macro
	*/

	//missing
	// - using ABCDRegions
	// - adding more plots
	// - implementing Compute and Draw methods

struct ABCDRegions{
	float a;
	float b;
	float c;
	float d;
};

class CrossSectionPseudoExp{

  public:
	CrossSectionPseudoExp();
	CrossSectionPseudoExp(int NofBins, float QCDMean, float VJetMean, float TTJetMean);/** Mean are used for histograms[0-Mean] */
	~CrossSectionPseudoExp();

	void Fill(MCExpectation& mcexp, pair<float,float> QCDEstimation, pair<float,float> VJetsEstimation, pair<float,float> TTJetsEstimation, pair<float,float> XsectionEstimation);
	void Compute();
	void Draw();
	void Write(TFile* fout, string label = string(""));

  private:
  	//A entry corresponds to one event in the vectors
	//info about events passing selection cuts
	vector<MCExpectation> MCExp;
	//info about estimation: first = estimation - second = error
	vector<pair<float,float> > QCDEstimation;
	vector<pair<float,float> > VJetEstimation;
	//ABCD
	vector<ABCDRegions> ABCDInfos;/** the numbers of events in the 4 regions*/
	//V+Jets estimation methods
	vector<vector<TH2F*> > BtagInfos;/** second dimension is the number of working point - TH2 contain jets vs bjets multiplicities*/

	///////////
	//Plots  //
	///////////
	//ABCD
	TH1F* h_ABCD_Estimation;/** Nof QCD events estimated*/
	TH1F* h_ABCD_Diff;/** Estimation - Expectation(MC) */
	TH1F* h_ABCD_Pull;/** Pull for ABCD estimation */
	//V+jets
	TH1F* h_VJets_Estimation;/** Nof V+Jets events estimated*/
	TH1F* h_VJets_Diff;/** Estimation - Expectation(MC) */
	TH1F* h_VJets_Pull;/**Pull for V+jets estimation */
	TH1F* h_TTJets_Estimation;/** Nof TT+Jets events estimated*/
	TH1F* h_TTJets_Diff;/** Estimation - Expectation(MC) */
	TH1F* h_TTJets_Pull;/**Pull for TT+jets estimation */
	TH1F* h_Xsection_Estimation;/** Xsection estimated*/
	TH1F* h_Xsection_Diff;/** Estimation - Expectation(MC) */
	TH1F* h_Xsection_Pull;/**Pull for X-section estimation */
	TH1F* h_Xsection_RelError;
};

#endif
