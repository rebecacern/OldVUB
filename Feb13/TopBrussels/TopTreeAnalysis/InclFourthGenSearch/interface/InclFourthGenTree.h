#ifndef InclFourthGenTree_h
#define InclFourthGenTree_h

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
#include "TLorentzVector.h"

using namespace std;

class InclFourthGenTree : public TObject
{
 public:
	
  InclFourthGenTree() :
    TObject()
    ,eventID_(0)
    ,runID_(0)
    ,lumiBlockID_(0)
    ,nPV_(0)
    ,nPUBXm1_(0)
    ,nPU_(0)
    ,nPUBXp1_(0)
		,flavorHistoryPath_(0)
    ,SelectedSingleLepton_(false)
		,SelectedSingleMu_(false)
		,SelectedSingleEl_(false)
		,SelectedSSLepton_(false)
		,SelectedSSMu_(false)
		,SelectedSSEl_(false)
		,SelectedSSMuEl_(false)
		,SelectedTriLepton_(false)
		,SelectedMuMuMu_(false)
		,SelectedMuMuEl_(false)
		,SelectedMuElEl_(false)
		,SelectedElElEl_(false)
		,EE_(false)
		,EB_(false)
		,EEEE_(false)
		,EBEE_(false)
		,EBEB_(false)
    ,semiMuDecay_(false)
    ,semiElDecay_(false)
		,Wbosonpartonsmatched_(false)
		,WMassmatched_(0)
		,eventWeight_(0)
		,MET_()
    ,selectedJets_()
    ,bTagTCHE_()
    ,bTagTCHP_()
		,partonFlavourJet_()
		,selectedForwardJets_()
		,InitJets_()
		,InitJetsbTagTCHE_()
    ,InitJetsbTagTCHP_()
    ,selectedMuons_()
		,selectedElectrons_()
		,selectedMuonsRelIso_()
		,selectedElectronsRelIso_()
		,chargeMisIdRateBarrel_(0)
		,chargeMisIdRateEndcap_(0)
//    ,all4JetsMCMatched_(false)
//    ,allHadronicJetsMCMatched_(false)
//    ,hadrLJet1_(9999)
//    ,hadrLJet2_(9999)
//		,hadrBJet_(9999)
//		,leptBJet_(9999)
//		,MatchedJetsIndices_()
		,mcQuarksForMatching_()
//    ,hadrBQuark_()
//		,leptBQuark_()
//    ,hadrLQuark1_()
//    ,hadrLQuark2_()
//    ,topDecayedLept_(false)
//    ,mvaVals_()
//    ,mvaResults_()
		,quarksFromW_()
   {;}
  
  ~InclFourthGenTree() {;}
  
  unsigned int eventID() const { return eventID_; }
  unsigned int runID() const { return runID_; }
  unsigned int lumiBlockID() const { return lumiBlockID_; }
  unsigned int nPV() const { return nPV_; }
  unsigned int nPUBXm1() const { return nPUBXm1_; }
  unsigned int nPU() const { return nPU_; }
  unsigned int nPUBXp1() const { return nPUBXp1_; }
	int flavorHistoryPath() const { return flavorHistoryPath_; }
  bool SelectedSingleLepton() const { return SelectedSingleLepton_; }
	bool SelectedSingleMu() const { return SelectedSingleMu_; }
	bool SelectedSingleEl() const { return SelectedSingleEl_; }
	bool SelectedSSLepton() const { return SelectedSSLepton_; }
	bool SelectedSSMu() const { return SelectedSSMu_; }
	bool SelectedSSEl() const { return SelectedSSEl_; }
	bool SelectedSSMuEl() const { return SelectedSSMuEl_; }
	bool SelectedTriLepton() const { return SelectedTriLepton_; }
	bool SelectedMuMuMu() const { return SelectedMuMuMu_; }
	bool SelectedMuMuEl() const { return SelectedMuMuEl_; }
	bool SelectedMuElEl() const { return SelectedMuElEl_; }
	bool SelectedElElEl() const { return SelectedElElEl_; }	
	bool isEE() const { return EE_; }	
	bool isEB() const { return EB_; }	
	bool isEEEE() const { return EEEE_; }	
	bool isEBEE() const { return EBEE_; }	
	bool isEBEB() const { return EBEB_; }	
  bool semiMuDecay() const { return semiMuDecay_; }
  bool semiElDecay() const { return semiElDecay_; }
	bool Wbosonpartonsmatched() const { return Wbosonpartonsmatched_; }
	float WMassmatched() const { return WMassmatched_; }
  float eventWeight() { return eventWeight_; }
	TLorentzVector met() const { return MET_; }
  vector<TLorentzVector> selectedJets() const { return selectedJets_; }
  TLorentzVector selectedJet(int i) const { return selectedJets_[i]; }
  vector<float> bTagTCHE() const { return bTagTCHE_; }
  vector<float> bTagTCHP() const { return bTagTCHP_; }
	vector<int> partonFlavourJet() const { return partonFlavourJet_; }
	vector<TLorentzVector> selectedForwardJets() const { return selectedForwardJets_; }
	vector<TLorentzVector> InitJets() const { return InitJets_; }
  vector<float> InitJetsbTagTCHE() const { return InitJetsbTagTCHE_; }
  vector<float> InitJetsbTagTCHP() const { return InitJetsbTagTCHP_; }
	vector<TLorentzVector> selectedMuons() const { return selectedMuons_; }
  vector<TLorentzVector> selectedElectrons() const { return selectedElectrons_; }
	vector<float> selectedMuonsRelIso() const { return selectedMuonsRelIso_; }
  vector<float> selectedElectronsRelIso() const { return selectedElectronsRelIso_; }
//	bool all4JetsMCMatched() const { return all4JetsMCMatched_; }
//  bool allHadronicJetsMCMatched() const { return allHadronicJetsMCMatched_; }
//  unsigned int hadrLJet1() const { return hadrLJet1_; } // index, according to MC
//  unsigned int hadrLJet2() const { return hadrLJet2_; } // index, according to MC
//	unsigned int hadrBJet() const { return hadrBJet_; } // index, according to MC
//	unsigned int leptBJet() const { return leptBJet_; } // index, according to MC
//	vector<unsigned int> MatchedJetsIndices() const { return MatchedJetsIndices_; } // indices, according to MC
	vector<TLorentzVector> mcQuarksForMatching() const { return mcQuarksForMatching_; }
//  TLorentzVector hadrBQuark() const { return hadrBQuark_; }
//	TLorentzVector leptBQuark() const { return leptBQuark_; }
//  TLorentzVector hadrLQuark1() const { return hadrLQuark1_; }
//  TLorentzVector hadrLQuark2() const { return hadrLQuark2_; }	
  float chargeMisIdRateBarrel() { return chargeMisIdRateBarrel_; }
  float chargeMisIdRateEndcap() { return chargeMisIdRateEndcap_; }
	vector<TLorentzVector> quarksFromW() const { return quarksFromW_; }

	void setEventID(unsigned int eventID) { eventID_ = eventID; }
  void setRunID(unsigned int runID) { runID_ = runID; }
  void setLumiBlockID(unsigned int lumiBlockID) { lumiBlockID_ = lumiBlockID; }  
  void setNPV(unsigned int nPV) { nPV_ = nPV; }
  void setNPUBXm1(unsigned int nPUBXm1) { nPUBXm1_ = nPUBXm1; }
  void setNPU(unsigned int nPU) { nPU_ = nPU; }
  void setNPUBXp1(unsigned int nPUBXp1) { nPUBXp1_ = nPUBXp1; }
	void setFlavorHistoryPath(int flavorHistoryPath) { flavorHistoryPath_ = flavorHistoryPath; }
	void setSelectedSingleLepton(bool SelectedSingleLepton) { SelectedSingleLepton_ = SelectedSingleLepton; }
	void setSelectedSingleMu(bool SelectedSingleMu) { SelectedSingleMu_ = SelectedSingleMu; }
	void setSelectedSingleEl(bool SelectedSingleEl) { SelectedSingleEl_ = SelectedSingleEl;}
	void setSelectedSSLepton(bool SelectedSSLepton) { SelectedSSLepton_ = SelectedSSLepton;}
	void setSelectedSSMu(bool SelectedSSMu) { SelectedSSMu_ = SelectedSSMu; }
	void setSelectedSSEl(bool SelectedSSEl) { SelectedSSEl_ = SelectedSSEl; }
	void setSelectedSSMuEl(bool SelectedSSMuEl) { SelectedSSMuEl_ = SelectedSSMuEl; }
	void setSelectedTriLepton(bool SelectedTriLepton) { SelectedTriLepton_ = SelectedTriLepton; }
	void setSelectedMuMuMu(bool SelectedMuMuMu) { SelectedMuMuMu_ = SelectedMuMuMu; }
	void setSelectedMuMuEl(bool SelectedMuMuEl) { SelectedMuMuEl_ = SelectedMuMuEl; }
	void setSelectedMuElEl(bool SelectedMuElEl) { SelectedMuElEl_ = SelectedMuElEl; }
	void setSelectedElElEl(bool SelectedElElEl) { SelectedElElEl_ = SelectedElElEl; }	
	void setSelectedEE(bool isEE) { EE_ = isEE; }	
	void setSelectedEB(bool isEB) { EB_ = isEB; }	
	void setSelectedEEEE(bool isEEEE) { EEEE_ = isEEEE; }	
	void setSelectedEBEE(bool isEBEE) { EBEE_ = isEBEE; }	
	void setSelectedEBEB(bool isEBEB) { EBEB_ = isEBEB; }	
  void setSemiMuDecay(bool semiMuDecay) { semiMuDecay_ = semiMuDecay; }
  void setSemiElDecay(bool semiElDecay) { semiElDecay_ = semiElDecay; }
	void setWbosonpartonsmatched(bool Wbosonpartonsmatched) { Wbosonpartonsmatched_ = Wbosonpartonsmatched; }
	void setWMassmatched(float WMassmatched) { WMassmatched_ = WMassmatched; }
	void setEventWeight(float eventWeight) { eventWeight_ = eventWeight; }
	void setMET(TLorentzVector MET) { MET_ = MET; }
	void setSelectedJets(vector<TLorentzVector> selectedJets) { selectedJets_ = selectedJets; }
  void setBTagTCHE(vector<float> bTagTCHE) { bTagTCHE_ = bTagTCHE; }
  void setBTagTCHP(vector<float> bTagTCHP) { bTagTCHP_ = bTagTCHP; }
	void setpartonFlavourJet(vector<int> partonFlavourJet) { partonFlavourJet_ = partonFlavourJet; }
	void setSelectedForwardJets(vector<TLorentzVector> selectedForwardJets) { selectedForwardJets_ = selectedForwardJets; }
	void setInitJets(vector<TLorentzVector> InitJets) { InitJets_ = InitJets; }
	void setInitJetsBTagTCHE(vector<float> InitJetsbTagTCHE) { InitJetsbTagTCHE_ = InitJetsbTagTCHE; }
  void setInitJetsBTagTCHP(vector<float> InitJetsbTagTCHP) { InitJetsbTagTCHP_ = InitJetsbTagTCHP; }
	void setMuons(vector<TLorentzVector> selectedMuons) { selectedMuons_ = selectedMuons; }
	void setElectrons(vector<TLorentzVector> selectedElectrons) { selectedElectrons_ = selectedElectrons; }
	void setMuonsRelIso(vector<float> selectedMuonsRelIso) { selectedMuonsRelIso_ = selectedMuonsRelIso; }
	void setElectronsRelIso(vector<float> selectedElectronsRelIso) { selectedElectronsRelIso_ = selectedElectronsRelIso; }
//	void setAll4JetsMCMatched(bool all4JetsMCMatched) { all4JetsMCMatched_ = all4JetsMCMatched; }
//  void setAllHadronicJetsMCMatched(bool allHadronicJetsMCMatched) { allHadronicJetsMCMatched_ = allHadronicJetsMCMatched; }
//  void setHadrLJet1(int hadrLJet1) { hadrLJet1_ = hadrLJet1; }
//  void setHadrLJet2(int hadrLJet2) { hadrLJet2_ = hadrLJet2; }
//  void setHadrBJet(int hadrBJet) { hadrBJet_ = hadrBJet; }
//	void setLeptBJet(int leptBJet) { leptBJet_ = leptBJet; }
//	void setMatchedJetsIndices(vector<unsigned int> MatchedJetsIndices) { MatchedJetsIndices_ = MatchedJetsIndices; }
  void setmcQuarksForMatching(vector<TLorentzVector> mcQuarksForMatching) { mcQuarksForMatching_ = mcQuarksForMatching; }
//  void setHadrBQuark(TLorentzVector hadrBQuark) { hadrBQuark_ = hadrBQuark; }
//	void setLeptBQuark(TLorentzVector leptBQuark) { leptBQuark_ = leptBQuark; }
//  void setHadrLQuark1(TLorentzVector hadrLQuark1) { hadrLQuark1_ = hadrLQuark1; }
//  void setHadrLQuark2(TLorentzVector hadrLQuark2) { hadrLQuark2_ = hadrLQuark2; }
	void setChargeMisIdRateBarrel(float chargeMisIdRateBarrel) { chargeMisIdRateBarrel_ = chargeMisIdRateBarrel; }	
	void setChargeMisIdRateEndcap(float chargeMisIdRateEndcap) { chargeMisIdRateEndcap_ = chargeMisIdRateEndcap; }	
	void setQuarksFromW(vector<TLorentzVector> quarksFromW) { quarksFromW_ = quarksFromW; }
	
	//bool topDecayedLept() const { return topDecayedLept_; }
  //float* mvaVals() { return mvaVals_; }
  //float mvaVal(int i) const { return mvaVals_[i]; }
  //unsigned int* mvaResult(int i) { return mvaResults_[i]; }
	  
  //void setSemiMuDecay(bool semiMuDecay) { semiMuDecay_ = semiMuDecay; }
  //void setSemiElDecay(bool semiElDecay) { semiElDecay_ = semiElDecay; }
  //void setTopDecayedLept(bool topDecayedLept) { topDecayedLept_ = topDecayedLept; }
  //void setMvaVals(vector<float> mvaVals)
  //{
  //  for(unsigned int i=0; i<12; i++)
  //    mvaVals_[i] = mvaVals[i];
  //}
  //void setMvaResults(vector< vector<unsigned int> > mvaResults)
  //{
  //  for(unsigned int iCombi=0; iCombi<12; iCombi++)
  //    for(unsigned int iJet=0; iJet<4; iJet++)
  //      mvaResults_[iCombi][iJet] = mvaResults[iCombi][iJet];
  //}
  
 protected:
  
  unsigned int eventID_;
  unsigned int runID_;
  unsigned int lumiBlockID_;
  unsigned int nPV_;
  unsigned int nPUBXm1_;
  unsigned int nPU_;
  unsigned int nPUBXp1_;
	int flavorHistoryPath_;
  bool SelectedSingleLepton_;
	bool SelectedSingleMu_;
	bool SelectedSingleEl_;
	bool SelectedSSLepton_;
	bool SelectedSSMu_;
	bool SelectedSSEl_;
	bool SelectedSSMuEl_;
	bool SelectedTriLepton_;
	bool SelectedMuMuMu_;
	bool SelectedMuMuEl_;
	bool SelectedMuElEl_;
	bool SelectedElElEl_;	
	bool EE_;	
	bool EB_;	
	bool EEEE_;	
	bool EBEE_;	
	bool EBEB_;	
  bool semiMuDecay_;
  bool semiElDecay_;
	bool Wbosonpartonsmatched_;
	float WMassmatched_;
	float eventWeight_;
	TLorentzVector MET_;
	vector<TLorentzVector> selectedJets_; // all selected jet
	vector<float> bTagTCHE_; // indices like selectedJets indices
  vector<float> bTagTCHP_;
	vector<int> partonFlavourJet_;
	vector<TLorentzVector> selectedForwardJets_; // all selected forward jet
	vector<TLorentzVector> InitJets_; // all initial jets, needed for MC trigger efficiency reweighting
  vector<float> InitJetsbTagTCHE_; // indices like InitJets indices
  vector<float> InitJetsbTagTCHP_;
	vector<TLorentzVector> selectedMuons_;
	vector<TLorentzVector> selectedElectrons_;
	vector<float> selectedMuonsRelIso_;
	vector<float> selectedElectronsRelIso_;	
//	bool all4JetsMCMatched_;
//  bool allHadronicJetsMCMatched_;
//  int hadrLJet1_;
//  int hadrLJet2_;
//	int hadrBJet_; //index according to MC
//  int leptBJet_;
//	vector<unsigned int> MatchedJetsIndices_;
  float chargeMisIdRateBarrel_;
  float chargeMisIdRateEndcap_;
	  vector<TLorentzVector> mcQuarksForMatching_;
//  TLorentzVector hadrBQuark_;
//	TLorentzVector leptBQuark_;
//  TLorentzVector hadrLQuark1_;
//  TLorentzVector hadrLQuark2_;
	
  //bool topDecayedLept_;
  //float mvaVals_[12];
  //unsigned int mvaResults_[12][4]; // jet indices
	vector<TLorentzVector> quarksFromW_;
  
	
	ClassDef (InclFourthGenTree,1);
};

#endif
