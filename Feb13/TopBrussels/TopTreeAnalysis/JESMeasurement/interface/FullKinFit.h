#ifndef FullKinFit_h
#define FullKinFit_h

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

#include "TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"
#include "TopTreeAnalysisBase/KinFitter/interface/TKinFitter.h"
#include "TopTreeAnalysisBase/KinFitter/interface/TFitConstraintM.h"
#include "TopTreeAnalysisBase/KinFitter/interface/TFitParticleEtThetaPhiEMomFix.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"

#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeProducer/interface/TRootJet.h"
#include "TopTreeProducer/interface/TRootMuon.h"

using namespace std;
using namespace TopTree;

class FullKinFit {

  public:
    FullKinFit(Dataset* d, ResolutionFit *resFitLightJets, ResolutionFit *resFitBJets, bool measureTopMass = false, bool measureTopMassDiff = false, bool dummy = false);
    ~FullKinFit();
    void SetResFitL7(ResolutionFit *resFitLightJetsL7, ResolutionFit *resFitBJetsL7);
    void SetJets(vector<TRootJet*> jets);
    void SetMVAStuff(pair<float, vector<unsigned int> > MVAvals);
    TH2F* FitEvent(TRootEvent* event, float WMass = 80.4, float topMass = 172.5, bool writePNG = false, int jetCombi=0);
    TH2F* DummyMonster(int jetCombi=0);
    float* EstimateTopMass(TRootEvent* event, float WMass = 80.4, bool writePNG = false, int jetCombi = 0, bool correctCombi = false);
    void Write(TFile* fout, bool savePNG = false, string pathPNG = string(""));
  
  private:
    map<string, TH1F*> histo1D_;
    vector<TRootJet*> jets_;
    ResolutionFit *resFitLightJets_;
    ResolutionFit *resFitBJets_;
    ResolutionFit *resFitLightJetsL7_; // reso's after L7 corrections
    ResolutionFit *resFitBJetsL7_; // reso's after L7 corrections
    string MVAMethod_;
    Dataset* dataset_;
    
    bool measureTopMass_;
    bool measureTopMassDiff_;
    
    pair<float, vector<unsigned int> > MVAvals_;
    
    unsigned int nBinsCorrB_;
    float stepCorrB_;
    float startCorrB_;
    unsigned int nBinsCorrLight_;
    float stepCorrLight_;
    float startCorrLight_;
    unsigned int nBinsTopMass_;
    float stepTopMass_;
    float startTopMass_;
    
    // some control stuff
    int nTotalJetCombis_;
    int nParabolaFitJetCombis_;
    int nFinalKinFitJetCombis_;
    int nFinalJetCombis_;
};

#endif
