#ifndef TopFCNC_KinFit_h
#define TopFCNC_KinFit_h

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cstdio>
#include <algorithm>
#include <sys/stat.h>

#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "TH2F.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "TopTreeAnalysis/MCInformation/interface/ResolutionFit.h"
#include "TopTreeAnalysis/KinFitter/interface/TKinFitter.h"
#include "TopTreeAnalysis/KinFitter/interface/TFitConstraintM.h"
#include "TopTreeAnalysis/KinFitter/interface/TFitConstraintEp.h"
//#include "TopTreeAnalysis/KinFitter/interface/TFitParticleEtThetaPhiEMomFix.h"
#include "TopTreeAnalysis/KinFitter/interface/TFitParticleEtThetaPhi.h"
#include "TopTreeAnalysis/Reconstruction/interface/FactorizedJetCorrector.h"
#include "TopTreeAnalysis/Reconstruction/interface/JetCorrectorParameters.h"
#include "TopTreeAnalysis/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysis/TopFCNC/interface/TopFCNC_GenEvt.h"
#include "TopTreeAnalysis/TopFCNC/interface/TopFCNC_Evt.h"

#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeProducer/interface/TRootJet.h"
#include "TopTreeProducer/interface/TRootMuon.h"

using namespace std;
using namespace TopTree;

class TopFCNC_KinFit {

  public:
    TopFCNC_KinFit(Dataset* dataset, ResolutionFit *resFitLeptons, ResolutionFit *resFitBJets, ResolutionFit *resFitQJets, ResolutionFit *resFitLightJets, float WMass = 80.4, float Zmass = 91.2, float topMass = 172.5);
    ~TopFCNC_KinFit();
    void FitEvent(TopFCNC_GenEvt *topFCNC_GenEvt);
    void FitEvent(TopFCNC_Evt *topFCNC_Evt);
    void FitEvent(TLorentzVector &lepton1, TLorentzVector &lepton2, TLorentzVector &qJet, TLorentzVector &lightJet1, TLorentzVector &lightJet2, TLorentzVector &bJet);
  
    void Write(TFile* fout, bool savePNG = false, string pathPNG = string(""));

    Double_t GetProb() {return prob_;};
    Double_t GetChi2() {return chi2_;};
    Double_t GetNdof() {return ndof_;};
    
    void SetMaxNbIter(Int_t maxNbIter)    {maxNbIter_ = maxNbIter;}
    void SetMaxDeltaS(Double_t maxDeltaS) {maxDeltaS_ = maxDeltaS;}
    void SetMaxF(Double_t maxF)           {maxF_ = maxF;}

    void SetConstraints(vector<string> &constraints);
    void SetVerbosity(Bool_t verbose) {verbose_ = verbose;};
    void SetFitVerbosity(Bool_t verbosity_fit) {verbosity_fit_ = verbosity_fit;};
  
  private:
    map<string, TH1F*> histo1D_;
    Dataset       *dataset_;

    TKinFitter *kinfit_;

    TFitParticleEtThetaPhi *lepton1_;
    TFitParticleEtThetaPhi *lepton2_;

    TFitParticleEtThetaPhi *lightJet1_;
    TFitParticleEtThetaPhi *lightJet2_;
    TFitParticleEtThetaPhi *bJet_;
    TFitParticleEtThetaPhi *qJet_;

    TFitConstraintM  *consHadW_;
    TFitConstraintM  *consLepZ_;
    TFitConstraintM  *consHadTop_;
    TFitConstraintM  *consFcncTop_;
    TFitConstraintM  *consEqualTop_;
    TFitConstraintEp *consSumPx_;
    TFitConstraintEp *consSumPy_;

    ResolutionFit *resFitLeptons_;
    ResolutionFit *resFitBJets_;
    ResolutionFit *resFitQJets_;
    ResolutionFit *resFitLightJets_;

    Double_t       prob_;
    Double_t       chi2_;
    Int_t          ndof_;
    Int_t          maxNbIter_;  // Maximum number of iterations
    Double_t       maxDeltaS_;  // Convergence criterium for deltaS
    Double_t       maxF_;       // Convergence criterium for F
  
    Bool_t         constrainSumPt_;
    Bool_t         verbose_;
    Int_t          verbosity_fit_;
};

#endif
