#ifndef InclFourthGenSearchTools_h
#define InclFourthGenSearchTools_h

#include "TopTreeProducer/interface/TRootMuon.h"
#include "TopTreeProducer/interface/TRootJet.h"
#include "TopTreeProducer/interface/TRootElectron.h"
#include "TopTreeProducer/interface/TRootMET.h"
#include "TopTreeProducer/interface/TRootGenEvent.h"
#include "TopTreeProducer/interface/TRootMCParticle.h"

#include "TopTreeAnalysis/Content/interface/Dataset.h"
#include "TopTreeAnalysis/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysis/MCInformation/interface/JetPartonMatching.h"
//for Kinematic Fit
#include "TopTreeAnalysis/MCInformation/interface/ResolutionFit.h"
#include "TopTreeAnalysis/KinFitter/interface/TKinFitter.h"
#include "TopTreeAnalysis/KinFitter/interface/TFitConstraintM.h"
#include "TopTreeAnalysis/KinFitter/interface/TFitParticleEtThetaPhi.h"


#include <iostream>
#include <iomanip>
#include <string.h>
#include <sstream>
#include <map>
#include "TH1F.h"
#include "TFile.h"
#include "TMatrixD.h"


using namespace std;
using namespace TopTree;

std::string IntToStr( int n )
{
	std::ostringstream result;
	result << n;
	return result.str();
}

class InclFourthGenSearchTools
{
 public:
   InclFourthGenSearchTools(bool semiMuon, bool semiElectron, vector<Dataset* > datasets, float Luminosity, bool doKinematicFit);
   InclFourthGenSearchTools(const InclFourthGenSearchTools & i);
   ~InclFourthGenSearchTools();
   void SetResolutionFit(ResolutionFit* resFitLightJets);
   void FillPlots(int dataset, int nbOfBtags, int nbOfWs, float HT, vector<TRootMuon*> selectedMuons, vector<TRootElectron*> selectedElectrons, vector<TRootMET*> mets, vector<TRootJet*> selectedJets, float scaleFactor);
   void FillPlots(int dataset, int nbOfBtags, int nbOfWs, float HT, vector<TLorentzVector> selectedMuons, vector<TLorentzVector> selectedElectrons, float met, vector<TLorentzVector> selectedJets, float scaleFactor); 
	 void CalculateTopMass(TRootJet* WJet1, TRootJet* WJet2, TRootJet* HadBJet);
   void CalculateTopMass(TLorentzVector WJet1, TLorentzVector WJet2, TLorentzVector HadBJet);
	 float GetHT() const {return HT_;};
   float GetMtop() const {return Mtop_;};
   float GetMtop_kinfit() const {return Mtop_kinfit_;};
   void FillMassPlots(int d, int nbOfBtags, int nbOfWs, float scaleFactor);
   void WritePlots(TFile* fout, TDirectory* th1dir, bool savePNG, string pathPNG);
   void TestPurityGoodCombinations(int d, int nbOfBtags, int nbOfWs, TRootGenEvent* genEvt, vector<TRootMCParticle*> mcParticles, vector<TRootJet*> selectedJets_MVAinput, bool TprimeEvaluation, float scaleFactor);
   void PrintPurityGoodCombinations();
   
 private:
   bool semiMuon_;
   bool semiElectron_;
   vector<Dataset* > datasets_;
   float Luminosity_;
   bool doKinematicFit_;
   
   ResolutionFit * resFitLightJets_;
   
   map<string,TH1F*> histo1D;
   map<string,TH2F*> histo2D;   
   map<string,MultiSamplePlot*> MSPlot;
   
   float HT_;
   float Mtop_;
   float Mtop_kinfit_;
   
   int counter_all4JetsMatched_MCdef;
   int counter_hadronictopJetsMatched_MCdef;
   int counter_hadronictopJetsMatched_largestMass;
   int counter_hadronictopJetsMatched_smallestMass;
   int counter_hadronictopJetsMatched_largestPt;
   int counter_hadronictopJetsMatched_smallestPt;
   int counter_hadronictopJetsMatched_largestMass_WJetsFixed;
   int counter_hadronictopJetsMatched_smallestMass_WJetsFixed;
   int counter_hadronictopJetsMatched_largestPt_WJetsFixed;
   int counter_hadronictopJetsMatched_smallestPt_WJetsFixed;
   
};

#endif
