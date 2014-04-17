#ifndef TRootNTuple_h
#define TRootNTuple_h

#include "Rtypes.h"
#include "TObject.h"
#include <string>

#include <vector>
#include <iostream>
#include <iomanip>
#include <map>
#include <cmath>

using namespace std;

class TRootNTuple : public TObject 
{
	
public:
	TRootNTuple() :
	  TObject()
	  ,eventId_(-999)
	  ,runId_(-999)
	  ,lumiBlockId_(-999)
	  ,selSemiMu_(0)
	  ,trueSemiMu_(0)
	  ,selSemiEl_(0)
	  ,trueSemiEl_(0)
	  ,nPV_(-999)
	  ,mlj_(-999.) //property of signal sample jet
	  ,m3_(-999.) //property of signal sample jet
	  ,pt_(-999.) //property of signal sample jet
	  ,E_(-999.) //property of signal sample jet
	  ,eta_(-999) //property of signal sample jet
	  ,weight_(-999.) //general event property
	  ,scalefactor_(-999.) //general event property
	  ,scalefactor3D_(-999.) //general event property
	  ,dataSetNumber_(-1) //general event property
	  ,partonFlavour_(-999) //property of signal sample jet
	  ,nConstituents_(-999) //property of signal sample jet
	  ,n90_(-999) //property of signal sample jet
	  ,n60_(-999) //property of signal sample jet
	  ,jetArea_(-999.) //property of signal sample jet
	  ,pileupEnergy_(-999.) //property of signal sample jet
	  ,maxDistance_(-999.) //property of signal sample jet
	  ,ecalEnergyFraction_(-999.) //property of signal sample jet
	  ,hcalEnergyFraction_(-999.) //property of signal sample jet
	  ,maxEInEmTowers_(-999.) //property of signal sample jet
	  ,maxEInHadTowers_(-999.) //property of signal sample jet
	  ,towersArea_(-999.)  //property of signal sample jet
	  ,mljControl_(-999.) //property of control sample jet
	  ,mljControl2_(-999.) //property of control sample jet
	  ,mljControl3_(-999.) //property of control sample jet
	  ,chiSq_(-999.) //event property
	  ,ptControl_(-999.) //property of control sample jet
	  ,ptControl2_(-999.) //property of control sample jet
	  ,EControl_(-999.) //property of control sample jet
	  ,EControl2_(-999.) //property of control sample jet
	  ,ptControl3_(-999.) //property of control sample jet
	  ,pthadb_(-999.) //property of control sample jet
	  ,ptMuon_(-999.)
	  ,EMuon_(-999.)
	  ,etaMuon_(-999.)
	  ,chargeMuon_(-999.)
	  ,mvaTrigID_(-999.)
	  ,delRlj_(-999.)
	  ,partonFlavourControl_(-999)
	  ,partonFlavourControl2_(-999)
	  ,partonFlavourControl3_(-999)
	  ,bTagControl_(-999.)
	  ,bTagControl2_(-999.)
	  ,bTagControl3_(-999.)
	  ,etaControl_(-999.)
	  ,etaControl2_(-999.)
	  ,etaControl3_(-999.)  
	  ,delRljControl_(-999.)
	  ,delRljControl2_(-999.)
	  ,btag_trackCountingHighEffBJetTags_hadb_(-999.)
	  ,nJets_(-999)
	  ,nMuons_(-999)
	  ,jetsOverlap_(-999.)
	  ,jetMuonOverlap_(-999.)
	  ,jet_is_b_(false)
	  ,jet_is_nonb_(false)
	  ,jet_is_radq_(false)
	  ,jet_is_hadqq_(false)
	  ,jetControl_is_b_(false)
	  ,jetControl_is_nonb_(false)
	  ,jetControl_is_radq_(false)
	  ,jetControl_is_hadqq_(false)
	  ,jetControl2_is_b_(false)
	  ,jetControl2_is_nonb_(false)
	  ,jetControl2_is_radq_(false)
	  ,jetControl2_is_hadqq_(false)
	  ,jethadb_is_b_(false)
	  ,jethadb_is_nonb_(false)
	  ,jethadb_is_radq_(false)
	  ,jethadb_is_hadqq_(false)
	  ,HadWCandMass_(-999.)
	  ,HadTopCandMass_(-999.)
	  ,LepBlCandMass_(-999.)
	  ,R_inAll_(false)
	  ,R_inHad_(false)
	  ,R_or_lepb_inHad_(false)
	  ,allMatched_(false)
	  ,delOmegalj_(-999.)
	  ,delOmegaljhad_(-999.)
	  ,delOmegaljControl_(-999.)
	  ,delOmegaljControl2_(-999.)
	  ,pthadtop_(-999.)
	  ,delPhitt_(-999.)
	  ,etafifth_(-999.)
	  ,ptfifth_(-999.)
	  ,MET_(-999.)
	  {;}
	  
	  
	  ~TRootNTuple(){;}
	
	  //output methods

	  Int_t eventId() const { return eventId_; }
	  Int_t runId() const { return runId_; }
	  Int_t lumiBlockId() const { return lumiBlockId_; }

	  Bool_t semiMu() const {return selSemiMu_;}
	  Bool_t truesemiMu() const {return trueSemiMu_;}
	  Bool_t semiEl() const {return selSemiEl_;}
	  Bool_t truesemiEl() const {return trueSemiEl_;}

	  Int_t nPV ( ) const { return nPV_; }

	Double_t mlj() const { return mlj_;}
	Double_t m3() const { return m3_;}
	
	Double_t Btag(unsigned int index) const { return Btag_[index]; }
	Double_t BtagCS1(unsigned int index) const { return BtagC1_[index]; }
	Double_t BtagCS2(unsigned int index) const { return BtagC2_[index]; }
	
	Double_t pt() const { return pt_;}
	Double_t E() const { return E_;}
	Double_t eta() const { return eta_;}
	Double_t weight() const { return weight_;}
	Double_t scalefactor() const { return scalefactor_;}
	Double_t scalefactor3D() const { return scalefactor3D_;}
	Int_t dataSetNumber() const { return dataSetNumber_;}
	string dataSetName() const { return dataSetName_;}
	Int_t partonFlavour() const { return partonFlavour_;}
	Int_t nConstituents() const { return nConstituents_; }
	Int_t n90() const { return n90_; }
	Int_t n60() const { return n60_; }
	Float_t jetArea() const { return jetArea_; }
	Float_t pileupEnergy() const { return pileupEnergy_; }
	Float_t maxDistance() const { return maxDistance_; }
	Float_t ecalEnergyFraction() const { return ecalEnergyFraction_; }
	Float_t hcalEnergyFraction() const { return hcalEnergyFraction_; }
	Float_t maxEInEmTowers() const { return maxEInEmTowers_; }
	Float_t maxEInHadTowers() const { return maxEInHadTowers_;}
	Float_t towersArea() const { return towersArea_;} 
	Double_t mljControl() const { return mljControl_;}
	Double_t mljControl2() const { return mljControl2_;}
	Double_t mljControl3() const { return mljControl3_;}
	Double_t chiSq() const { return chiSq_;}
	Double_t ptControl() const { return ptControl_;}
	Double_t ptControl2() const { return ptControl2_;}
	Double_t EControl() const { return EControl_;}
	Double_t EControl2() const { return EControl2_;}
	Double_t ptControl3() const { return ptControl3_;}
	Double_t pthadb() const { return pthadb_;}
	Double_t ptMuon() const { return ptMuon_;}
	Double_t EMuon() const { return EMuon_;}
	Double_t etaMuon() const { return etaMuon_;}
	Double_t chargeMuon() const { return chargeMuon_;}
	Double_t mvaTrigID() const {return mvaTrigID_;}
	Double_t delRlj() const { return delRlj_;}  
	Int_t partonFlavourControl() const { return partonFlavourControl_;}
	Int_t partonFlavourControl2() const { return partonFlavourControl2_;}
	Int_t partonFlavourControl3() const { return partonFlavourControl3_;}
	Double_t bTagControl() const { return bTagControl_;}
	Double_t bTagControl2() const { return bTagControl2_;}
	Double_t bTagControl3() const { return bTagControl3_;}
	Double_t etaControl() const { return etaControl_;}
	Double_t etaControl2() const { return etaControl2_;}
	Double_t etaControl3() const { return etaControl3_;}
	Double_t delRljControl() const { return delRljControl_;}  
	Double_t delRljControl2() const { return delRljControl2_;}  
	Bool_t jet_is_b() {return jet_is_b_;}
	Bool_t jet_is_nonb() {return jet_is_nonb_;}
	Bool_t jet_is_radq() {return jet_is_radq_;}
	Bool_t jet_is_hadqq() {return jet_is_hadqq_;}
	Bool_t jetControl_is_b() {return jetControl_is_b_;}
	Bool_t jetControl_is_nonb() {return jetControl_is_nonb_;}
	Bool_t jetControl_is_radq() {return jetControl_is_radq_;}
	Bool_t jetControl_is_hadqq() {return jetControl_is_hadqq_;}
	Bool_t jetControl2_is_b() {return jetControl2_is_b_;}
	Bool_t jetControl2_is_nonb() {return jetControl2_is_nonb_;}
	Bool_t jetControl2_is_radq() {return jetControl2_is_radq_;}
	Bool_t jetControl2_is_hadqq() {return jetControl2_is_hadqq_;}
	Bool_t jethadb_is_b() {return jethadb_is_b_;}
	Bool_t jethadb_is_nonb() {return jethadb_is_nonb_;}
	Bool_t jethadb_is_radq() {return jethadb_is_radq_;}
	Bool_t jethadb_is_hadqq() {return jethadb_is_hadqq_;}
	Double_t HadWCandMass() {return HadWCandMass_;}
	Double_t HadTopCandMass() {return HadTopCandMass_;}
	Double_t LepBlCandMass() {return LepBlCandMass_;}
	Bool_t R_inAll() {return R_inAll_;}
	Bool_t R_inHad() {return R_inHad_;}
	Bool_t R_or_lepb_inHad() {return R_or_lepb_inHad_;}
	Bool_t allMatched() {return allMatched_;}
	Double_t MET() {return MET_;}
	
	//some extra info for further event selection
	Double_t btag_trackCountingHighEffBJetTags_hadb() const { return btag_trackCountingHighEffBJetTags_hadb_;}
	Int_t nJets() const { return nJets_;}
	Int_t nMuons() const { return nMuons_;}
	Double_t jetsOverlap() const { return jetsOverlap_;}
	Double_t jetMuonOverlap() const { return jetMuonOverlap_;}
	Double_t delOmegalj() const { return delOmegalj_;}
	Double_t delOmegaljhad() const { return delOmegaljhad_;}
	Double_t delOmegaljControl() const { return delOmegaljControl_;}
	Double_t delOmegaljControl2() const { return delOmegaljControl2_;}
	Double_t pthadtop() const { return pthadtop_;}
	Double_t delPhitt() const { return delPhitt_;}
	Double_t etafifth() const { return etafifth_;}
	Double_t ptfifth() const { return ptfifth_;}

	// pdf unc

	Float_t getx1() const { return x1_;}
	Float_t getx2() const { return x2_;}
	Float_t getq2() const { return q2_;}
	Int_t getid1() const { return id1_;}
	Int_t getid2() const { return id2_;}

	
	//input methods

	void setEventId( int i ) { eventId_ = i ; }
	void setRunId( int i ) { runId_ = i; }
	void setLumiBlockId( int i ) { lumiBlockId_ = i; }

	void setDecay(bool selSemiMu, bool trueSemiMu, bool selSemiEl, bool trueSemiEl) {
	  selSemiMu_=selSemiMu;
	  trueSemiMu_=trueSemiMu;
	  selSemiEl_=selSemiEl;
	  trueSemiEl_=trueSemiEl;
	}

	void setnPV ( int i ) { nPV_ = i; }
	
	void setMlj(Double_t mlj) { mlj_=mlj;}
	void setM3(Double_t m3) { m3_=m3;}
	
	void setBtag(unsigned int index, Double_t btag) { Btag_[index] = btag; }
	void setBtagCS1(unsigned int index, Double_t btag) { BtagC1_[index] = btag; }
	void setBtagCS2(unsigned int index, Double_t btag) { BtagC2_[index] = btag; }
	
	void setPt(Double_t pt) { pt_=pt;}
	void setE(Double_t E) { E_=E;}
	void setEta(Double_t eta) { eta_=eta;}
	void setWeight(Double_t weight) { weight_=weight;}
	void setScaleFactor(Double_t f) { scalefactor_=f;}
	void setScaleFactor3D(Double_t f) { scalefactor3D_=f;}
	void setDataSetNumber(Int_t dataSetNumber) { dataSetNumber_=dataSetNumber;} 
	void setDataSetName(string dataSetName) { dataSetName_=dataSetName;} 
	void setPartonFlavour(Int_t partonFlavour) { partonFlavour_=partonFlavour;}
	void setNConstituents(Int_t nConstituents) { nConstituents_ = nConstituents; }
	void setN90(Int_t n90) { n90_ = n90; }
	void setN60(Int_t n60) { n60_ = n60; }
	void setJetArea(Float_t jetArea) { jetArea_ = jetArea; }
	void setPileupEnergy(Float_t pileupEnergy) { pileupEnergy_ = pileupEnergy; }
	void setMaxDistance(Float_t maxDistance) { maxDistance_ = maxDistance; }
	void setMaxEInEmTowers(Float_t maxEInEmTowers) { maxEInEmTowers_ = maxEInEmTowers; }
	void setMaxEInHadTowers(Float_t maxEInHadTowers) { maxEInHadTowers_ = maxEInHadTowers; }
	void setTowersArea(Float_t towersArea) {towersArea_ = towersArea; }
	void setEcalEnergyFraction(Float_t ecalEnergyFraction) { ecalEnergyFraction_ = ecalEnergyFraction; }
	void setHcalEnergyFraction(Float_t hcalEnergyFraction) { hcalEnergyFraction_ = hcalEnergyFraction; }
	void setMljControl(Double_t mljControl) { mljControl_=mljControl;}
	void setMljControl2(Double_t mljControl2) { mljControl2_=mljControl2;}
	void setMljControl3(Double_t mljControl3) { mljControl3_=mljControl3;}
	void setChiSq(Double_t chiSq) { chiSq_=chiSq;}
	void setPtControl(Double_t ptControl) { ptControl_=ptControl;}
	void setPtControl2(Double_t ptControl2) { ptControl2_=ptControl2;}
	void setEControl(Double_t EControl) { EControl_=EControl;}
	void setEControl2(Double_t EControl2) { EControl2_=EControl2;}
	void setPtControl3(Double_t ptControl3) { ptControl3_=ptControl3;}
	void setPthadb(Double_t pthadb) { pthadb_=pthadb;}
	void setPtMuon(Double_t ptMuon) {ptMuon_=ptMuon;}
	void setEMuon(Double_t EMuon) {EMuon_=EMuon;}
	void setEtaMuon(Double_t etaMuon) {etaMuon_=etaMuon;}
	void setChargeMuon(Double_t chargeMuon) {chargeMuon_=chargeMuon;}
	void setmvaTrigID(Double_t i) {mvaTrigID_=i;}
	void setDelRlj(Double_t delRlj) {delRlj_=delRlj;}  
	void setPartonFlavourControl(Int_t partonFlavourControl) { partonFlavourControl_=partonFlavourControl;}
	void setPartonFlavourControl2(Int_t partonFlavourControl2) { partonFlavourControl2_=partonFlavourControl2;}
	void setPartonFlavourControl3(Int_t partonFlavourControl3) { partonFlavourControl3_=partonFlavourControl3;}
	void setbTagControl(Double_t bTagControl) { bTagControl_=bTagControl;}
	void setbTagControl2(Double_t bTagControl2) { bTagControl2_=bTagControl2;}
	void setbTagControl3(Double_t bTagControl3) { bTagControl3_=bTagControl3;}
	void setEtaControl(Double_t etaControl) { etaControl_=etaControl;}
	void setEtaControl2(Double_t etaControl2) { etaControl2_=etaControl2;}
	void setEtaControl3(Double_t etaControl3) { etaControl3_=etaControl3;}
	void setDelRljControl(Double_t delRljControl) {delRljControl_=delRljControl;}  
	void setDelRljControl2(Double_t delRljControl2) {delRljControl2_=delRljControl2;}  
	void setBtag_trackCountingHighEffBJetTags_hadb(Double_t btag_trackCountingHighEffBJetTags_hadb) { btag_trackCountingHighEffBJetTags_hadb_ = btag_trackCountingHighEffBJetTags_hadb;}
	void setnJets(Int_t nJets) { nJets_=nJets;}
	void setnMuons(Int_t nMuons) { nMuons_=nMuons;}
	void setJetsOverlap(Double_t jetsOverlap) { jetsOverlap_=jetsOverlap;}
	void setJetMuonOverlap(Double_t jetMuonOverlap) { jetMuonOverlap_=jetMuonOverlap;}
	void setJet_is_b(Bool_t jet_is_b) {jet_is_b_ = jet_is_b;}
	void setJet_is_nonb(Bool_t jet_is_nonb) {jet_is_nonb_ = jet_is_nonb;}
	void setJet_is_radq(Bool_t jet_is_radq) {jet_is_radq_ = jet_is_radq;}
	void setJet_is_hadqq(Bool_t jet_is_hadqq) {jet_is_hadqq_ = jet_is_hadqq;}
	void setJetControl_is_b(Bool_t jetControl_is_b) {jetControl_is_b_ = jetControl_is_b;}
	void setJetControl_is_nonb(Bool_t jetControl_is_nonb) {jetControl_is_nonb_ = jetControl_is_nonb;}
	void setJetControl_is_radq(Bool_t jetControl_is_radq) {jetControl_is_radq_ = jetControl_is_radq;}
	void setJetControl_is_hadqq(Bool_t jetControl_is_hadqq) {jetControl_is_hadqq_ = jetControl_is_hadqq;}
	void setJetControl2_is_b(Bool_t jetControl2_is_b) {jetControl2_is_b_ = jetControl2_is_b;}
	void setJetControl2_is_nonb(Bool_t jetControl2_is_nonb) {jetControl2_is_nonb_ = jetControl2_is_nonb;}
	void setJetControl2_is_radq(Bool_t jetControl2_is_radq) {jetControl2_is_radq_ = jetControl2_is_radq;}
	void setJetControl2_is_hadqq(Bool_t jetControl2_is_hadqq) {jetControl2_is_hadqq_ = jetControl2_is_hadqq;}
	void setJethadb_is_b(Bool_t jethadb_is_b) {jethadb_is_b_ = jethadb_is_b;}
	void setJethadb_is_nonb(Bool_t jethadb_is_nonb) {jethadb_is_nonb_ = jethadb_is_nonb;}
	void setJethadb_is_radq(Bool_t jethadb_is_radq) {jethadb_is_radq_ = jethadb_is_radq;}
	void setJethadb_is_hadqq(Bool_t jethadb_is_hadqq) {jethadb_is_hadqq_ = jethadb_is_hadqq;}
	//virtual TString typeName() const { return "TRootNTuples"; }
	void setHadWCandMass(Double_t HadWCandMass) {HadWCandMass_ = HadWCandMass;}
	void setHadTopCandMass(Double_t HadTopCandMass) {HadTopCandMass_ = HadTopCandMass;}
	void setLepBlCandMass(Double_t LepBlCandMass) {LepBlCandMass_ = LepBlCandMass;}
	void setR_inAll(Bool_t R_inAll) {R_inAll_ = R_inAll;}
	void setR_inHad(Bool_t R_inHad) {R_inHad_ = R_inHad;}
	void setR_or_lepb_inHad(Bool_t R_or_lepb_inHad) {R_or_lepb_inHad_ = R_or_lepb_inHad;}
	void setAllMatched(Bool_t allMatched) {allMatched_ = allMatched;}
	void setDelOmegalj(Double_t delOmegalj) { delOmegalj_ = delOmegalj;}
	void setDelOmegaljhad(Double_t delOmegaljhad) { delOmegaljhad_ = delOmegaljhad;}
	void setDelOmegaljControl(Double_t delOmegaljControl) { delOmegaljControl_ = delOmegaljControl;}
	void setDelOmegaljControl2(Double_t delOmegaljControl2) { delOmegaljControl2_ = delOmegaljControl2;}
	void setPthadtop(Double_t pthadtop) { pthadtop_ = pthadtop;}
	void setDelPhitt(Double_t delPhitt) { delPhitt_ = delPhitt;}
	void setEtafifth(Double_t etafifth) { etafifth_ = etafifth;}
	void setPtfifth(Double_t ptfifth) { ptfifth_ = ptfifth;}
	
	void setx1(Float_t x1) { x1_ = x1;}
	void setx2(Float_t x2) { x2_ = x2;}
	void setq2(Float_t q2) { q2_ = q2;}
	void setid1(Int_t id1) { id1_ = id1;}
	void setid2(Int_t id2) { id2_ = id2;}
	
	void setMET(Double_t MET) { MET_ = MET; } 
	
private:

	Int_t eventId_;
	Int_t runId_;
	Int_t lumiBlockId_;

	bool selSemiMu_;
	bool trueSemiMu_;
	bool selSemiEl_;
	bool trueSemiEl_;

	Int_t nPV_;

	Double_t mlj_;
	Double_t m3_;
	
	Double_t Btag_[20];
	Double_t BtagC1_[20];
	Double_t BtagC2_[20];

 
	Double_t pt_;
	Double_t E_;
	Double_t eta_;
	Double_t weight_;
	Double_t scalefactor_;
	Double_t scalefactor3D_;
	Int_t dataSetNumber_;
	string dataSetName_;
	Int_t partonFlavour_;
	Int_t nConstituents_;               // Number of constituents of the jet (calotowers for CaloJet / PFParticles for PFJet)
	Int_t n90_;                         // Number of constituents of the jet carrying 90% of tje jet energy
	Int_t n60_;                         // Number of constituents of the jet carrying 60% of tje jet energy
	Float_t jetArea_;                   // Jet area
	Float_t pileupEnergy_;               // Pile-up Energy
	Float_t maxDistance_;               // Maximum distance from jet to constituent
	Float_t ecalEnergyFraction_;        // ECAL Energy Fraction
	Float_t hcalEnergyFraction_;        // HCAL Energy Fraction
	Float_t maxEInEmTowers_;
	Float_t maxEInHadTowers_;
	Float_t towersArea_; 
	Double_t mljControl_;
	Double_t mljControl2_;
	Double_t mljControl3_;
	Double_t chiSq_;
	Double_t ptControl_;
	Double_t ptControl2_;
	Double_t EControl_;
	Double_t EControl2_;
	Double_t ptControl3_;
	Double_t pthadb_;
	Double_t ptMuon_;
	Double_t EMuon_;
	Double_t etaMuon_;
	Double_t chargeMuon_;
	Double_t mvaTrigID_;
	Double_t delRlj_;
	Int_t partonFlavourControl_;
	Int_t partonFlavourControl2_;
	Int_t partonFlavourControl3_;
	Double_t bTagControl_;
	Double_t bTagControl2_;
	Double_t bTagControl3_;
	Double_t etaControl_;
	Double_t etaControl2_;
	Double_t etaControl3_;
	Double_t delRljControl_;
	Double_t delRljControl2_;
	Double_t btag_trackCountingHighEffBJetTags_hadb_;
	Int_t nJets_;
	Int_t nMuons_;
	Double_t jetsOverlap_;
	Double_t jetMuonOverlap_;
	Bool_t jet_is_b_;
	Bool_t jet_is_nonb_;
	Bool_t jet_is_radq_;
	Bool_t jet_is_hadqq_;
	Bool_t jetControl_is_b_;
	Bool_t jetControl_is_nonb_;
	Bool_t jetControl_is_radq_;
	Bool_t jetControl_is_hadqq_;
	Bool_t jetControl2_is_b_;
	Bool_t jetControl2_is_nonb_;
	Bool_t jetControl2_is_radq_;
	Bool_t jetControl2_is_hadqq_;
	Bool_t jethadb_is_b_;
	Bool_t jethadb_is_nonb_;
	Bool_t jethadb_is_radq_;
	Bool_t jethadb_is_hadqq_;
	Double_t HadWCandMass_;
	Double_t HadTopCandMass_; 
	Double_t LepBlCandMass_; 
	Bool_t R_inAll_;
	Bool_t R_inHad_; 
	Bool_t R_or_lepb_inHad_; 
	Bool_t allMatched_;
	Double_t delOmegalj_;
	Double_t delOmegaljhad_;
	Double_t delOmegaljControl_;
	Double_t delOmegaljControl2_;
	Double_t pthadtop_;
	Double_t delPhitt_;
	Double_t etafifth_;
	Double_t ptfifth_;
	Double_t MET_;
	// for PDF unc

	Float_t x1_;
	Float_t x2_;
	Int_t id1_;
	Int_t id2_;
	Float_t q2_;

	ClassDef (TRootNTuple,1);
};

#endif
