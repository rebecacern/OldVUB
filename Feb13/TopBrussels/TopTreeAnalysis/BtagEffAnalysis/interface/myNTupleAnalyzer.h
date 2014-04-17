//file to define the proto-type of my analyzer which is containing the b-tag measurement
#ifndef myNTupleAnalyzer_h
#define myNTupleAnalyzer_h

#include "TString.h"
#include <stdio.h>
#include <time.h>
#include <ctime>
#include <iomanip>
#include "TopTreeProducer/interface/TRootMuon.h"
#include "TopTreeProducer/interface/TRootElectron.h"
#include "TopTreeProducer/interface/TRootJet.h"
#include "TopTreeProducer/interface/TRootMET.h"
#include "TopTreeProducer/interface/TRootGenEvent.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootParticle.h"
#include "TopTreeProducer/interface/TRootMCParticle.h"

#include "../interface/TRootNTuple.h"
#include "../interface/PtEtaBinContainer.h"
#include "../interface/PtEtaBin.h"
#include "../interface/WorkingPointBin.h"
#include "../interface/MakerJetOrigin.h"

#include "TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"

#include <TFile.h>
#include <TObject.h>
#include <TH1.h>
#include <TProfile.h>
#include <TH2.h>
#include "TCanvas.h"
#include <TBranch.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TLine.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "TClassTable.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TColor.h"
#include "TPaveStats.h"

#include "../../macros/Style.C"

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

using namespace std;

class myNTupleAnalyzer{ 
	
public:
	
	~myNTupleAnalyzer();
	myNTupleAnalyzer(TString *inDir, TString *outDir, int *runSamples_, int nRunSamples_, double chisqCut, int inBin_, int outBin_, double, int, bool);
	
	void run(int verbosity, 
			 int leftlimit, 
			 int centerleftlimit, 
			 int centerrightlimit, 
			 int rightlimit, 
			 bool doSCreweigh, 
			 bool doTwoLights,
			 bool useFit, 
			 bool do2D, 
			 bool do2Dcontrol,
			 bool doPtEtaBin,
			 bool doJESchange, 
			 double JESfactor,
			 bool doNewF, double, double, double,
			 double runNb,
			 bool doFfromMC,
			 int nSystematic,
			 int decay,
			 int fitMode,
			 int nBtag_);
	void getEff(double*);
	void getEffVal(double*);
	void getEffMCVal(double*);
	void getEffRRVal(double*);
	void getEffRRMCVal(double*);
	//void getFMCVal(double*);
	double getFMCVal();
	double getFVal();
	double getFnoRWVal();
	double getFMCerrVal();
	double getFerrVal();
	double getFnoRWerrVal();
	
	double getmljMeannoRWVal();
	double getmljMeanVal();
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
	
	double get_lepb_b_Counter();
	double get_lepb_nonb_Counter();
	double get_lepb_hadqq_Counter();
	double get_lepb_radq_Counter();
	double get_lepbtag_b_Counter();
	double get_lepbtag_nonb_Counter();
	double get_lepbtag_hadqq_Counter();
	double get_lepbtag_radq_Counter();
	double get_q1q2_b_Counter();
	double get_q1q2_nonb_Counter();
	double get_q1q2_hadqq_Counter();
	double get_q1q2_radq_Counter();
	
	void getPercentiles(int*);
	
	void setFMCBias(double);
	void setBackgroundFraction(double);
	
private: 
	
	string datasetName;
	float origLumi;
	int nTaggers;
	
	int nSystematic;
	
	TString inDir_;
	TString outDir_;
	int *runSamples_;
	int nRunSamples_;
	double matchChiSquareCut_;
	
	int nWP;
	double eff[18];
	double effVal[18];
	double effMCVal[18];
	double effRRVal[18];
	double effRRMCVal[18];
	double FMCVal;
	double FVal;
	double FMCerrVal;
	double FerrVal;
	
	double mljMeannoRWVal;
	double mljMeanVal;
	double mljMeanMCVal;
	double mljSigmaVal;
	double mljSigmanoRWVal;
	double mljSigmaMCVal;
	double mlj_W_MeanVal;
	double mlj_W_MeannoRWVal;
	double mlj_W_MeanMCVal;
	double mlj_W_SigmaVal;
	double mlj_W_SigmanoRWVal;
	double mlj_W_SigmaMCVal;
	double mlj_R_MeanVal;
	double mlj_R_MeannoRWVal;
	double mlj_R_MeanMCVal;
	double mlj_R_SigmaVal;
	double mlj_R_SigmanoRWVal;
	double mlj_R_SigmaMCVal;
	
	double lepb_b_Counter;
	double lepb_nonb_Counter;
	double lepb_hadqq_Counter;
	double lepb_radq_Counter;
	double lepbtag_b_Counter;
	double lepbtag_nonb_Counter;
	double lepbtag_hadqq_Counter;
	double lepbtag_radq_Counter;
	double q1q2_b_Counter;
	double q1q2_nonb_Counter;
	double q1q2_hadqq_Counter;
	double q1q2_radq_Counter;
	
	
	double FMCBias;
	bool changeFbias;
	double backgroundFraction;
	bool changeBackgroundFraction; 
	
	int inBin_;
	int outBin_;
	
	int nIncreaseSample_;
	int *increaseSample_;
	double increaseWeightFraction_;
	
	double desiredIntLum_;
	double dataLum_;
	int nPseudoExp_;
	bool doPseudoExp_; //to do the pseudo-exps

	bool useTTJetsExcl_;	

	double* newFpt;
	double* newFeta;
	int* percentiles_;
	
	std::map<string,TH1D*> chi2CutHistos;
	
	std::map<string,TH1D*> histos1D;
	std::map<string,TH2D*> histos2D;
	std::map<string,TGraphErrors*> graphs;
	std::map<string,TCanvas*> canvas_2Dmeas;
	
	std::map<string,TH1D*> pullHistos;
	std::map<string,TH2D*> pullHistos2D;
	std::map<string,TH3D*> pullHistos3D;
	
	std::map<string,TH2D*> XSHistos;
	
	std::map<float,vector<std::pair<float,float> > > bTagPseudoExpResults;
	std::map<float,vector<std::pair<float,float> > > misTagPseudoExpResults;
	
	std::map<string,vector<float> > bTagPseudoExpResultsForF;
	std::map<string,vector<float> > bTagPseudoExpResultsForXS;
	std::vector<float> bTagEffBiasResults;
	
	double nTTbarBeforeChiSq;
	double nTTbarAfterChiSq;
	
	double nTTbarBeforeMLBCUT;
	double nTTbarAfterMLBCUT;
	
	double nTTbarBeforeRefSel;
	double nTTbarAfterRefSel;
    int nTTbarAfterRefSel_nPV[100];
	double nTTbarAfterRefSel2;
	
	std::map<std::string, std::vector<float> > nSelected_;
	
	vector<Dataset*> datasets;
	vector<int> datasets_color;
	vector<string> datasets_title;
	
	std::map<double,TH1D*> MLj_shift_templates;
	
	
};

#endif
