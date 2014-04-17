#ifndef TopFCNC_Evt_h
#define TopFCNC_Evt_h

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include "TLorentzVector.h"
#include "TopTreeProducer/interface/TRootJet.h"
#include "TopTreeProducer/interface/TRootMET.h"
#include "TopTreeAnalysis/Reconstruction/interface/Combination.h"
#include "TopTreeAnalysis/Reconstruction/interface/MEzCalculator.h"
#include "TopTreeAnalysis/Selection/interface/Selection.h"

using namespace std;
using namespace stdcomb;

class TopFCNC_Evt : public TObject
{
	public: 
		// W/Z boson leptonic decay channel
		enum LeptonType {kNone, kMuon, kElec, kTau};

		TopFCNC_Evt() :
			TObject(),
      eventID_(0),
      runID_(0),
      lumiBlockID_(0),
      idParton1_(-1),
      idParton2_(-1),
      xParton1_(-1.),
      xParton2_(-1.),
      factorizationScale_(-1.),
      nPV_(0),
      nTruePU_(-1.),
      eventWeight_(1.),
			zLeptonicChannel_(kNone),
			wLeptonicChannel_(kNone),
			isDiLeptonic_(false),
			isTriLeptonic_(false),
			verbose_(false)
	  {;}

		TopFCNC_Evt(LeptonType type) :
			TObject(),
      eventID_(0),
      runID_(0),
      lumiBlockID_(0),
      idParton1_(-1),
      idParton2_(-1),
      xParton1_(-1.),
      xParton2_(-1.),
      factorizationScale_(-1.),
      nPV_(0),
      nTruePU_(-1.),
      eventWeight_(1.),
			zLeptonicChannel_(type),
			wLeptonicChannel_(kNone),
			isDiLeptonic_(true),
			isTriLeptonic_(false),
			verbose_(false)
	  {;}

		TopFCNC_Evt(LeptonType type1, LeptonType type2) :
			TObject(),
      eventID_(0),
      runID_(0),
      lumiBlockID_(0),
      idParton1_(-1),
      idParton2_(-1),
      xParton1_(-1.),
      xParton2_(-1.),
      factorizationScale_(-1.),
      nPV_(0),
      nTruePU_(-1.),
      eventWeight_(1.),
			zLeptonicChannel_(type1),
			wLeptonicChannel_(type2),
			isDiLeptonic_(false),
			isTriLeptonic_(true),
			verbose_(false)
		{;}

		TopFCNC_Evt(const TopFCNC_Evt& evt) :
			TObject(evt),
      eventID_(evt.eventID_),
      runID_(evt.runID_),
      lumiBlockID_(evt.lumiBlockID_),
      idParton1_(evt.idParton1_),
      idParton2_(evt.idParton2_),
      xParton1_(evt.xParton1_),
      xParton2_(evt.xParton2_),
      factorizationScale_(evt.factorizationScale_),
      nPV_(evt.nPV_),
      nTruePU_(evt.nTruePU_),
      eventWeight_(evt.eventWeight_),
			zLeptonicChannel_(evt.zLeptonicChannel_),
			wLeptonicChannel_(evt.wLeptonicChannel_),
			isDiLeptonic_(evt.isDiLeptonic_),
			isTriLeptonic_(evt.isTriLeptonic_),	
			smDecayTop_(evt.smDecayTop_),
			W_(evt.W_),
			B_(evt.Z_),
			fcncDecayTop_(evt.fcncDecayTop_),
			Z_(evt.Z_),
			Q_(evt.Q_),
			leptonFromW_(evt.leptonFromW_),
			neutrino_(evt.neutrino_),
			quark1FromW_(evt.quark1FromW_),
			quark2FromW_(evt.quark2FromW_),
			lepton1FromZ_(evt.lepton1FromZ_),
			lepton2FromZ_(evt.lepton2FromZ_),
			selectedJets_(evt.selectedJets_),
			ISR_(evt.ISR_),
			smDecayTopRadiation_(evt.smDecayTopRadiation_),
			fcncDecayTopRadiation_(evt.fcncDecayTopRadiation_),
			verbose_(evt.verbose_)
		{;}
	
		~TopFCNC_Evt(){;}

    const UInt_t eventID()     const { return eventID_;}
    const UInt_t runID()       const { return runID_;}
    const UInt_t lumiBlockID() const { return lumiBlockID_;}

    const Int_t idParton1() const { return idParton1_;}
    const Int_t idParton2() const { return idParton2_;}

    const Float_t xParton1()           const { return xParton1_;}
    const Float_t xParton2()           const { return xParton2_;}
    
    const Float_t factorizationScale() const { return factorizationScale_;}

    const UInt_t  nPV()     const { return nPV_;}
    const Float_t nTruePU() const { return nTruePU_;}

    const Float_t eventWeight() const { return eventWeight_;}

		LeptonType zLeptonicChannel() const { return zLeptonicChannel_;}
		LeptonType wLeptonicChannel() const { return wLeptonicChannel_;}
	  
		Bool_t isDiLeptonic()  const {return isDiLeptonic_;}
		Bool_t isTriLeptonic() const {return isTriLeptonic_;}
		Bool_t isDiLeptonic(LeptonType t)                 const { return (zLeptonicChannel_==t ? true : false); }
		Bool_t isTriLeptonic(LeptonType t1,LeptonType t2) const { return((wLeptonicChannel_==t1 && zLeptonicChannel_==t2) ? true : false);}

		const TRootParticle leptonFromW()   const { return leptonFromW_;}
		const TRootParticle neutrino()      const { return neutrino_;}
		const TRootJet quark1FromW()        const { return quark1FromW_;}
		const TRootJet quark2FromW()        const { return quark2FromW_;}

		const TRootParticle lepton1FromZ()  const { return lepton1FromZ_;}
		const TRootParticle lepton2FromZ()  const { return lepton2FromZ_;}

		const TRootParticle smDecayTop()    const { return smDecayTop_;}
		const TRootParticle W()             const { return W_;}
		const TRootJet B()                  const { return B_;}
		const TRootParticle fcncDecayTop()  const { return fcncDecayTop_;}
		const TRootParticle Z()             const { return Z_;}
		const TRootJet Q()                  const { return Q_;}

    const std::vector<TRootJet> selectedJets()             const { return selectedJets_;}
		const std::vector<TRootJet> ISR()                      const { return ISR_;}
		const std::vector<TRootJet> smDecayTopRadiation()      const { return smDecayTopRadiation_;}
		const std::vector<TRootJet> fcncDecayTopRadiation()    const { return fcncDecayTopRadiation_;}
		
		const TRootMET met() const { return met_;}

    void SetEventID(UInt_t eventID)         { eventID_ = eventID; }
    void SetRunID(UInt_t runID)             { runID_ = runID; }
    void SetLumiBlockID(UInt_t lumiBlockID) { lumiBlockID_ = lumiBlockID; }

    void SetIdParton1(Int_t idParton1) { idParton1_ = idParton1; }
    void SetIdParton2(Int_t idParton2) { idParton2_ = idParton2; }

    void SetxParton1(Float_t xParton1)           { xParton1_ = xParton1; }
    void SetxParton2(Float_t xParton2)           { xParton2_ = xParton2; }
    void SetFactorizationScale(Float_t factorizationScale) { factorizationScale_ = factorizationScale; }

    void SetnPV(UInt_t nPV)         { nPV_ = nPV; }
    void SetnTruePU(Float_t nTruePU){ nTruePU_ = nTruePU; }

    void SetEventWeight(Float_t eventWeight) { eventWeight_ = eventWeight;}

		void SetDiLeptonicChannel(LeptonType type)
		{
			isDiLeptonic_ = true;
			isTriLeptonic_= false;
			zLeptonicChannel_ = type;
		}
		void SetTriLeptonicChannel(LeptonType type1, LeptonType type2)
		{
			isDiLeptonic_  = false;
			isTriLeptonic_ = true;
			zLeptonicChannel_ = type1;
			wLeptonicChannel_ = type2;
		}

		void SetTLorentzVector(TRootParticle &smDecayTop, TRootParticle &W, TRootJet &B, TRootParticle &fcncDecayTop, TRootParticle &Z, TRootJet &Q, TRootParticle &leptonFromW, TRootParticle &neutrino, TRootJet &quark1FromW, TRootJet &quark2FromW, TRootParticle &lepton1FromZ, TRootParticle &lepton2FromZ)
		{
			smDecayTop_    = smDecayTop;
			W_             = W;
			B_             = B;
			fcncDecayTop_  = fcncDecayTop;
			Z_             = Z;
			Q_             = Q;

			leptonFromW_   = leptonFromW;
			neutrino_      = neutrino;
			quark1FromW_   = quark1FromW;
			quark2FromW_   = quark2FromW;

			lepton1FromZ_  = lepton1FromZ;
			lepton2FromZ_  = lepton2FromZ;

		}

		void SetSmDecayTop(TRootParticle smDecayTop){ smDecayTop_ = smDecayTop; }
		void SetW(TRootParticle W){ W_ = W; }
		void SetB(TRootJet B){ B_ = B; }
		void SetFcncDecayTop(TRootParticle fcncDecayTop){ fcncDecayTop_ = fcncDecayTop; }
		void SetZ(TRootParticle Z){ Z_ = Z; }
		void SetQ(TRootJet Q){ Q_ = Q; }

		void SetLeptonFromW(TRootParticle leptonFromW){ leptonFromW_ = leptonFromW; }
		void SetNeutrino(TRootParticle neutrino){ neutrino_ = neutrino; }
		void SetQuark1FromW(TRootJet quark1FromW){ quark1FromW_ = quark1FromW; }
		void SetQuark2FromW(TRootJet quark2FromW){ quark2FromW_ = quark2FromW; }

		void SetLepton1FromZ(TRootParticle lepton1FromZ){ lepton1FromZ_ = lepton1FromZ; }
		void SetLepton2FromZ(TRootParticle lepton2FromZ){ lepton2FromZ_ = lepton2FromZ; }

		void SetRadiation(std::vector<TRootJet> smDecayTopRadiation, std::vector<TRootJet> fcncDecayTopRadiation, std::vector<TRootJet> ISR)
		{
			smDecayTopRadiation_   = smDecayTopRadiation; 
			fcncDecayTopRadiation_ = fcncDecayTopRadiation;
			ISR_ = ISR; 
		}
    void SetSelectedJets(std::vector<TRootJet*> selectedJets){
      selectedJets_.clear();
			for(UInt_t i=0;i<selectedJets.size();++i) selectedJets_.push_back(*selectedJets[i]);
    }
    void SetMET(TRootMET met){ met_ = met;}
    
		void SetISR(std::vector<TRootJet> ISR){ISR_ = ISR;}
		void SetSmDecayTopRadiation(std::vector<TRootJet> smDecayTopRadiation){ smDecayTopRadiation_ = smDecayTopRadiation; }
		void SetFcncDecayTopRadiation(std::vector<TRootJet> fcncDecayTopRadiation){ fcncDecayTopRadiation_ = fcncDecayTopRadiation; }

		void SetVerbosity(Bool_t verbose)                  { verbose_ = verbose; }

		void ReconstructEvt()
		{
			if(isDiLeptonic_ == isTriLeptonic_)
			{
				cerr << "Specify if top FCNC event candidate is di- or tri-leptonic before event reconstruction" << endl;
				exit(1);
			}
			else if(isDiLeptonic_)  ReconstructDiLeptEvt();
			else if(isTriLeptonic_) ReconstructTriLeptEvt();
		}

		void ReconstructDiLeptEvt()
		{
			if(!isDiLeptonic_)
			{
				cerr << "Top FCNC event candidate is not di-leptonic. Cannot be reconstructed as such" << endl;
				exit(1);
			}
			W_ = quark1FromW_ + quark2FromW_ ;
			smDecayTop_ = B_ + W_ ;
			Z_ = lepton1FromZ_ + lepton2FromZ_ ;
			fcncDecayTop_ = Q_ + Z_ ;
		}

		void ReconstructTriLeptEvt()
		{
			if(!isTriLeptonic_)
			{
				cerr << "Top FCNC event candidate is not tri-leptonic. Cannot be reconstructed as such" << endl;
				exit(1);
			}
			W_ = leptonFromW_+ neutrino_ ;
			smDecayTop_ = B_ + W_ ;
			Z_ = lepton1FromZ_ + lepton2FromZ_ ;
			fcncDecayTop_ = Q_ + Z_ ;
		}

		virtual TString typeName() const { return "TopFCNC_Evt"; }

		friend std::ostream& operator<< (std::ostream& stream, const TopFCNC_Evt& fcncEvt)
		{
			stream << "TopFCNC_Evt - is ";
			if(fcncEvt.isDiLeptonic()) {
				stream <<" DiLeptonic (";
				if(fcncEvt.zLeptonicChannel()==0) stream <<" no channel )";
				else if(fcncEvt.zLeptonicChannel()==1) stream <<" muonic )";
				else if(fcncEvt.zLeptonicChannel()==2) stream <<" electronic )";
				else if(fcncEvt.zLeptonicChannel()==3) stream <<" tauonic )";
			}
			else if(fcncEvt.isTriLeptonic()){
				stream <<" TriLeptonic (";
				if(fcncEvt.wLeptonicChannel()==0) stream << "W : no channel and ";
				else if(fcncEvt.wLeptonicChannel()==1) stream << "W : muonic and ";
				else if(fcncEvt.wLeptonicChannel()==2) stream << "W : electronic and ";
				else if(fcncEvt.wLeptonicChannel()==3) stream << "W : tauonic and ";

				if(fcncEvt.zLeptonicChannel()==0) stream << "Z : no channel )";
				else if(fcncEvt.zLeptonicChannel()==1) stream << "Z : muonic )";
				else if(fcncEvt.zLeptonicChannel()==2) stream << "Z : electronic )";
				else if(fcncEvt.zLeptonicChannel()==3) stream << "Z : tauonic )";
			}
			stream << std::endl;
			stream << "Nof ISR: "<< fcncEvt.ISR().size() << std::endl;
			stream << "Nof Top radiations: "<< fcncEvt.smDecayTopRadiation().size() + fcncEvt.fcncDecayTopRadiation().size() << std::endl;
			return stream;
		};


//	private:
	protected:

    UInt_t eventID_;
    UInt_t runID_;
    UInt_t lumiBlockID_;

    Int_t idParton1_;
    Int_t idParton2_;

    Float_t xParton1_;
    Float_t xParton2_;
    Float_t factorizationScale_;

    UInt_t  nPV_;
/*
    UInt_t  nPUBXm1_;
    UInt_t  nPU_;
    UInt_t  nPUBXp1_;
*/
    Float_t nTruePU_;

    Float_t eventWeight_;

		LeptonType zLeptonicChannel_;
		LeptonType wLeptonicChannel_;

		Bool_t isDiLeptonic_;
		Bool_t isTriLeptonic_;

		TRootParticle  smDecayTop_;
		TRootParticle  W_;
		TRootJet       B_;
		TRootParticle  fcncDecayTop_;
		TRootParticle  Z_;
		TRootJet       Q_;

		TRootParticle leptonFromW_;
		TRootParticle neutrino_;
		TRootJet      quark1FromW_;
		TRootJet      quark2FromW_;

		TRootParticle lepton1FromZ_;
		TRootParticle lepton2FromZ_;

		std::vector<TRootJet> selectedJets_;
		std::vector<TRootJet> ISR_;
		std::vector<TRootJet> smDecayTopRadiation_;
		std::vector<TRootJet> fcncDecayTopRadiation_;
		
		TRootMET met_;

		Bool_t verbose_;

		ClassDef (TopFCNC_Evt,2);
	};

#endif
