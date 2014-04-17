#ifndef TtJetEstPseudoExp_h
#define TtJetEstPseudoExp_h

#include <iostream>
#include <iomanip>
#include <vector>

#include "TRandom3.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TProfile.h"
//#include "TMathBase.h"

#include "TtJetEstimation.h"

using namespace std;


class TtJetEstPseudoExp {

  public:
	TtJetEstPseudoExp(int NbOfPseudoExp, TtJetEstimation * ttj);
	~TtJetEstPseudoExp();

	void RollTheDice(bool verbose = false,  bool addEstimError = false, vector<float> EstimError = vector<float>());
	
	
	double GetNTtJetError(bool doFit=false) const;
	double GetNTTBias(bool doFit=false) const;

	void Write(TFile* fout, string label);
	void PrintResults();
	
  private:
	int NbOfPseudoExp_;
	
	TtJetEstimation *ttjEst_;
	TRandom3* rand;
        
	//normalised expectation for all datasets
	vector<TH1F*> ExpPlotCRBkg;
	vector<TH1F*> ExpPlotSRBkg;
	vector<TH1F*> ExpPlotCRTtJet;
	vector<TH1F*> ExpPlotSRTtJet;
	vector<TH1F*> ExpPlotCRNP;
	vector<TH1F*> ExpPlotSRNP;
	
	//For a given pseudo-experiment (PE)
	vector<TH1F*> PEPlotCRBkg;
	vector<TH1F*> PEPlotSRBkg;
	vector<TH1F*> PEPlotBkgEstim;
	vector<TH1F*> PEPlotCRTtJet;
	vector<TH1F*> PEPlotSRTtJet;
	vector<TH1F*> PEPlotCRNP;
	vector<TH1F*> PEPlotSRNP;
 	TH1F* hPEPlotEstim;
 	TH1F* hPEPlotExpect;
	TH1F* hPEData;
	
	TH1F*  hNttEst;/** histo of Nttbar estimated for each pseudo-experiment*/
	TH1F*  hNttEstMinusNttExp;/** histo of (Nttbar estimated - Nttbar expected) for each pseudo-experiment*/
	TH1F*  hPull;/** histo of (Nttbar estimated - Nttbar expected)/Sigma for each pseudo-experiment*/
	TProfile* hEstMinusExp;/** TProfile with Estimation - Expectation */
	TH1F* hSystError;/** Mean value of the TProfile hEstMinusExp*/
};

#endif
