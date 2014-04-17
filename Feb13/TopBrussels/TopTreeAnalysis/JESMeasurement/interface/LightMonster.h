#ifndef LightMonster_h
#define LightMonster_h

#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "Rtypes.h"
#include "TObject.h"
#include "TVector3.h"
#include "TRef.h"
#include "TH2F.h"
#include "TLorentzVector.h"

using namespace std;

class LightMonster : public TObject
{
 public:
	
  LightMonster() :
    TObject()
    ,eventID_(0)
    ,runID_(0)
    ,lumiBlockID_(0)
    ,idParton1_(-9999)
    ,xParton1_(-9999.)
    ,idParton2_(-9999)
    ,xParton2_(-9999.)
    ,factorizationScale_(-9999.)
    ,nPV_(0)
    ,nPUBXm1_(0)
    ,nPU_(0)
    ,nPUBXp1_(0)
    ,nTruePU_(0)
    ,topMass_(0)
    ,antiTopMass_(0)
    ,selectedSemiMu_(false)
    ,semiMuDecay_(false)
    ,semiElDecay_(false)
    ,topDecayedLept_(false)
    ,all4JetsMCMatched_(false)
    ,allHadronicJetsMCMatched_(false)
    ,mvaVals_()
    ,mvaResults_()
    ,eventWeight_(0)
    ,mTopFit_()
    ,sigmaMTopFit_()
    ,chi2MTopFit_()
    ,hadrBJet_(0)
    ,hadrLJet1_(0)
    ,hadrLJet2_(0)
    ,leptBJet_(0)
    ,MET_()
    ,selectedJets_()
    ,bTagCSV_()
    ,lepton_()
    ,leptonCharge_(0)
    ,leptonPFRelIso_(-9999.)
    ,hadrBQuark_()
    ,hadrLQuark1_()
    ,hadrLQuark2_()
    ,leptBQuark_()
   {;}
  
  ~LightMonster() {;}
  
  unsigned int eventID() const { return eventID_; }
  unsigned int runID() const { return runID_; }
  unsigned int lumiBlockID() const { return lumiBlockID_; }
  int idParton1() const { return idParton1_; }
  float xParton1() const { return xParton1_; }
  int idParton2() const { return idParton2_; }
  float xParton2() const { return xParton2_; }
  float factorizationScale() const { return factorizationScale_; }
  unsigned int nPV() const { return nPV_; }
  unsigned int nPUBXm1() const { return nPUBXm1_; }
  unsigned int nPU() const { return nPU_; }
  unsigned int nPUBXp1() const { return nPUBXp1_; }
  unsigned int nTruePU() const { return nTruePU_; }
  float topMass() const { return topMass_; }
  float antiTopMass() const { return antiTopMass_; }
  bool selectedSemiMu() const { return selectedSemiMu_; }
  bool semiMuDecay() const { return semiMuDecay_; }
  bool semiElDecay() const { return semiElDecay_; }
  bool topDecayedLept() const { return topDecayedLept_; }
  bool all4JetsMCMatched() const { return all4JetsMCMatched_; }
  bool allHadronicJetsMCMatched() const { return allHadronicJetsMCMatched_; }
  float* mvaVals() { return mvaVals_; }
  float mvaVal(int i) const { return mvaVals_[i]; }
  unsigned int* mvaResult(int i) { return mvaResults_[i]; }
  float eventWeight() const { return eventWeight_; }
  float mTopFit(int iCombi) { return mTopFit_[iCombi]; } 
  float sigmaMTopFit(int iCombi) { return sigmaMTopFit_[iCombi]; }
  float chi2MTopFit(int iCombi) { return chi2MTopFit_[iCombi]; }
  int hadrBJet() const { return hadrBJet_; } // index, according to MC
  int hadrLJet1() const { return hadrLJet1_; } // index, according to MC
  int hadrLJet2() const { return hadrLJet2_; } // index, according to MC
  int leptBJet() const { return leptBJet_; } // index, according to MC
  TLorentzVector met() const { return MET_; }
  vector<TLorentzVector> selectedJets() const { return selectedJets_; }
  TLorentzVector selectedJet(int i) const { return selectedJets_[i]; }
  vector<float> bTagCSV() const { return bTagCSV_; }
  TLorentzVector lepton() const { return lepton_; }
  int leptonCharge() const { return leptonCharge_; }
  float leptonPFRelIso() const { return leptonPFRelIso_; }
  TLorentzVector hadrBQuark() const { return hadrBQuark_; }
  TLorentzVector hadrLQuark1() const { return hadrLQuark1_; }
  TLorentzVector hadrLQuark2() const { return hadrLQuark2_; }
  TLorentzVector leptBQuark() const { return leptBQuark_; }
  
  void setEventID(unsigned int eventID) { eventID_ = eventID; }
  void setRunID(unsigned int runID) { runID_ = runID; }
  void setLumiBlockID(unsigned int lumiBlockID) { lumiBlockID_ = lumiBlockID; }
  void setIdParton1(int Id1) { idParton1_ = Id1; }
  void setXParton1(float x1) { xParton1_ = x1; }
  void setIdParton2(int Id2) { idParton2_ = Id2; }
  void setXParton2(float x2) { xParton2_ = x2; }
  void setFactorizationScale(float scale) { factorizationScale_ = scale; }
  void setNPV(unsigned int nPV) { nPV_ = nPV; }
  void setNPUBXm1(unsigned int nPUBXm1) { nPUBXm1_ = nPUBXm1; }
  void setNPU(unsigned int nPU) { nPU_ = nPU; }
  void setNPUBXp1(unsigned int nPUBXp1) { nPUBXp1_ = nPUBXp1; }
  void setNTruePU(unsigned int nTruePU) { nTruePU_ = nTruePU; }
  void setTopMass(float topMass) { topMass_ = topMass; }
  void setAntiTopMass(float antiTopMass) { antiTopMass_ = antiTopMass; }
  void setSelectedSemiMu(bool selectedSemiMu) { selectedSemiMu_ = selectedSemiMu; }
  void setSemiMuDecay(bool semiMuDecay) { semiMuDecay_ = semiMuDecay; }
  void setSemiElDecay(bool semiElDecay) { semiElDecay_ = semiElDecay; }
  void setTopDecayedLept(bool topDecayedLept) { topDecayedLept_ = topDecayedLept; }
  void setAll4JetsMCMatched(bool all4JetsMCMatched) { all4JetsMCMatched_ = all4JetsMCMatched; }
  void setAllHadronicJetsMCMatched(bool allHadronicJetsMCMatched) { allHadronicJetsMCMatched_ = allHadronicJetsMCMatched; }
  void setMvaVals(vector<float> mvaVals)
  {
    for(unsigned int i=0; i<12; i++)
      mvaVals_[i] = mvaVals[i];
  }
  void setMvaResults(vector< vector<unsigned int> > mvaResults)
  {
    for(unsigned int iCombi=0; iCombi<12; iCombi++)
      for(unsigned int iJet=0; iJet<4; iJet++)
        mvaResults_[iCombi][iJet] = mvaResults[iCombi][iJet];
  }
  void setEventWeight(float eventWeight) { eventWeight_ = eventWeight; }
  void setMTopFitResults(vector< vector<float> > mTopFitResults)
  {
    for(unsigned int i=0; i<12; i++)
    {
      mTopFit_[i] = mTopFitResults[i][0];
      sigmaMTopFit_[i] = mTopFitResults[i][1];
      chi2MTopFit_[i] = mTopFitResults[i][2];
    }
  }
  void setHadrBJet(int hadrBJet) { hadrBJet_ = hadrBJet; }
  void setHadrLJet1(int hadrLJet1) { hadrLJet1_ = hadrLJet1; }
  void setHadrLJet2(int hadrLJet2) { hadrLJet2_ = hadrLJet2; }
  void setLeptBJet(int leptBJet) { leptBJet_ = leptBJet; }
  void setMET(TLorentzVector MET) { MET_ = MET; }
  void setSelectedJets(vector<TLorentzVector> selectedJets) { selectedJets_ = selectedJets; }
  void setBTagCSV(vector<float> bTagCSV) { bTagCSV_ = bTagCSV; }
  void setLepton(TLorentzVector lepton) { lepton_ = lepton; }
  void setLeptonCharge(int leptonCharge) { leptonCharge_ = leptonCharge; }
  void setLeptonPFRelIso(float leptonPFRelIso) { leptonPFRelIso_ = leptonPFRelIso; }
  void setHadrBQuark(TLorentzVector hadrBQuark) { hadrBQuark_ = hadrBQuark; }
  void setHadrLQuark1(TLorentzVector hadrLQuark1) { hadrLQuark1_ = hadrLQuark1; }
  void setHadrLQuark2(TLorentzVector hadrLQuark2) { hadrLQuark2_ = hadrLQuark2; }
  void setLeptBQuark(TLorentzVector leptBQuark) { leptBQuark_ = leptBQuark; }
  
 protected:
  
  unsigned int eventID_;
  unsigned int runID_;
  unsigned int lumiBlockID_;
  int idParton1_;
  float xParton1_;
  int idParton2_;
  float xParton2_;
  float factorizationScale_;
  unsigned int nPV_;
  unsigned int nPUBXm1_;
  unsigned int nPU_;
  unsigned int nPUBXp1_;
  unsigned int nTruePU_;
  float topMass_;
  float antiTopMass_;
  bool selectedSemiMu_;
  bool semiMuDecay_;
  bool semiElDecay_;
  bool topDecayedLept_;
  bool all4JetsMCMatched_;
  bool allHadronicJetsMCMatched_;
  float mvaVals_[12];
  unsigned int mvaResults_[12][4]; // jet indices
  float eventWeight_;
  float mTopFit_[12];
  float sigmaMTopFit_[12];
  float chi2MTopFit_[12];
  int hadrBJet_; //index according to MC
  int hadrLJet1_;
  int hadrLJet2_;
  int leptBJet_;
  TLorentzVector MET_;
  vector<TLorentzVector> selectedJets_; // all selected jets
  vector<float> bTagCSV_; // indices like selectedJets indices
  TLorentzVector lepton_;
  int leptonCharge_;
  float leptonPFRelIso_;
  TLorentzVector hadrBQuark_;
  TLorentzVector hadrLQuark1_;
  TLorentzVector hadrLQuark2_;
  TLorentzVector leptBQuark_;
  
  ClassDef (LightMonster,2);
};

#endif

