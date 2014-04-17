#ifndef StopPair_Evt_h
#define StopPair_Evt_h

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include "TLorentzVector.h"
#include "TopTreeProducer/interface/TRootJet.h"
#include "TopTreeProducer/interface/TRootMET.h"
#include "TopTreeAnalysis/Reconstruction/interface/MEzCalculator.h"
#include "TopTreeAnalysis/Selection/interface/Selection.h"

using namespace std;

class StopPair_Evt : public TObject
{
	public: 
		// W boson leptonic decay channel
		enum LeptonType {kNone, kMuon, kElec, kTau};

		StopPair_Evt() :
			TObject(),
			LeptonicChannel_(kNone)
			{;}

		StopPair_Evt(LeptonType type) :
			TObject(),
			LeptonicChannel_(type)
			{;}

		StopPair_Evt(const StopPair_Evt& evt) :
			TObject(evt),
			LeptonicChannel_(evt.LeptonicChannel_),
			hadronicStopCand_(evt.hadronicStopCand_),
			hadronicWCand_(evt.hadronicWCand_),
			hadronicCharginoCand_(evt.hadronicCharginoCand_),
			hadronicQCand_(evt.hadronicQCand_),
			hadronicQbarCand_(evt.hadronicQbarCand_),
			hadronicNeutralinoCand_(evt.hadronicNeutralinoCand_),
			leptonicStopCand_(evt.leptonicStopCand_),
			leptonicWCand_(evt.leptonicWCand_),
			leptonicCharginoCand_(evt.leptonicCharginoCand_),
			leptonicMuonCand_(evt.leptonicMuonCand_),
			leptonicElectronCand_(evt.leptonicElectronCand_),
			leptonicNeutrinoCand_(evt.leptonicNeutrinoCand_),
			leptonicNeutralinoCand_(evt.leptonicNeutralinoCand_),
			totMET_(evt.totMET_),
			ISR_(evt.ISR_),
			hadDecayTopRadiation_(evt.hadDecayTopRadiation_),
			lepDecayTopRadiation_(evt.lepDecayTopRadiation_)
			{;}
	
		~StopPair_Evt(){;}

		LeptonType LeptonicChannel() const { return LeptonicChannel_;}

		const TRootParticle hadronicStopCand()       const { return hadronicStopCand_;}
		const TRootParticle hadronicWCand()          const { return hadronicWCand_;}
		const TRootParticle hadronicCharginoCand()   const { return hadronicCharginoCand_;}
		const TRootJet hadronicQCand()               const { return hadronicQCand_;}
		const TRootJet hadronicQbarCand()            const { return hadronicQbarCand_;}

		const TRootParticle leptonicStopCand()       const { return leptonicStopCand_;}
		const TRootParticle leptonicWCand()          const { return leptonicWCand_;}
		const TRootParticle leptonicCharginoCand()   const { return leptonicCharginoCand_;}
		const TRootMuon     leptonicMuonCand()       const { return leptonicMuonCand_;}
		const TRootElectron leptonicElectronCand()   const { return leptonicElectronCand_;}
		const TRootParticle leptonicNeutrinoCand()   const { return leptonicNeutrinoCand_;}
		const TRootParticle leptonicNeutralinoCand() const { return leptonicNeutralinoCand_;}

		void SetLeptonicChannel(LeptonType type) {LeptonicChannel_ = type;}

		void SetHadronicStopCand(TRootParticle hadronicStopCand){ hadronicStopCand_ = hadronicStopCand; }
		void SetHadronicWCand(TRootParticle hadronicWCand){ hadronicWCand_ = hadronicWCand; }
		void SetHadronicCharginoCand(TRootParticle hadronicCharginoCand){ hadronicCharginoCand_ = hadronicCharginoCand; }
		void SetHadronicQCand(TRootJet hadronicQCand){ hadronicQCand_ = hadronicQCand;}
		void SetHadronicQbarCand(TRootJet hadronicQbarCand){ hadronicQbarCand_ = hadronicQbarCand;}

		void SetLeptonicStopCand(TRootParticle leptonicStopCand){ leptonicStopCand_ = leptonicStopCand;}
		void SetLeptonicWCand(TRootParticle leptonicWCand){ leptonicWCand_ = leptonicWCand;}
		void SetLeptonicCharginoCand(TRootParticle leptonicCharginoCand){ leptonicCharginoCand_ = leptonicCharginoCand;}
		void SetLeptonicMuonCand(TRootMuon leptonicMuonCand){ leptonicMuonCand_ = leptonicMuonCand;}
		void SetLeptonicElectronCand(TRootElectron leptonicElectronCand){ leptonicElectronCand_ = leptonicElectronCand;}
		void SetLeptonicNeutrinoCand(TRootParticle leptonicNeutrinoCand){ leptonicNeutrinoCand_ = leptonicNeutrinoCand;}
		void SetLeptonicNeutralinoCand(TRootParticle leptonicNeutralinoCand){ leptonicNeutralinoCand_ = leptonicNeutralinoCand;}

	protected:
		LeptonType LeptonicChannel_;

		TRootParticle hadronicStopCand_;
		TRootParticle hadronicWCand_;
		TRootParticle hadronicCharginoCand_;
		TRootJet      hadronicQCand_;
		TRootJet      hadronicQbarCand_;
		TRootParticle hadronicNeutralinoCand_;

		TRootParticle leptonicStopCand_;
		TRootParticle leptonicWCand_;
		TRootParticle leptonicCharginoCand_;
		TRootMuon     leptonicMuonCand_;
		TRootElectron leptonicElectronCand_;
		TRootParticle leptonicNeutrinoCand_;
		TRootParticle leptonicNeutralinoCand_;

		TRootMET      totMET_;

		std::vector<TRootJet> ISR_;
		std::vector<TRootJet> hadDecayTopRadiation_;
		std::vector<TRootJet> lepDecayTopRadiation_;

		ClassDef (StopPair_Evt,1);
	};

#endif
