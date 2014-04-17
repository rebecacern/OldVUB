#include "../interface/TopFCNC_GenEvt.h"

//ClassImp(TopFCNC_GenEvt);

void TopFCNC_GenEvt::ReconstructEvt(const std::vector<TRootMCParticle*> &mcParticles, bool debug){

  wLeptonicChannel_ = kNone;
  zLeptonicChannel_ = kNone;

  for(unsigned int i=0; i<mcParticles.size(); i++){
	  if( mcParticles[i]->status() != 3) continue;
	  if( fabs(mcParticles[i]->type()) == 6 )
	  {
		  if( fabs(mcParticles[i]->dauOneId()) == 24 || fabs(mcParticles[i]->dauTwoId()) == 24 || fabs(mcParticles[i]->dauThreeId()) == 24 || fabs(mcParticles[i]->dauFourId()) == 24) {   smDecayTop_ = *mcParticles[i]; continue; }
		  else if( fabs(mcParticles[i]->dauOneId()) == 23 || fabs(mcParticles[i]->dauTwoId()) == 23 || fabs(mcParticles[i]->dauThreeId()) == 23 || fabs(mcParticles[i]->dauFourId()) == 23) { fcncDecayTop_ = *mcParticles[i]; continue; }
		  else {
		  	cerr << "Cannot assign top quark according to its decay..." << endl;
		  	exit(1);
		  }
  	}

  	if( fabs(mcParticles[i]->motherType()) == 6 )
  	{
	  	if( fabs(mcParticles[i]->type()) == 24 )
	  	{
	  		W_ = *mcParticles[i];
		  	continue;
		  }
		  else if( fabs(mcParticles[i]->type()) ==  5 )
		  {
			  B_ = *mcParticles[i];
			  continue;
		  }
		  else if( fabs(mcParticles[i]->type()) == 23 )
		  {
			  Z_ = *mcParticles[i];
			  continue;
		  }
		  else if( fabs(mcParticles[i]->type()) >= 1  && fabs(mcParticles[i]->type()) <= 4 )
		  {
		  	Q_ = *mcParticles[i];
		  	continue;
		  }
/*
		else {
			cerr << "Cannot assign top quark decay product. Type =" << mcParticles[i]->type() << endl;
			exit(1);
		}
*/
  	}
  	else if( fabs(mcParticles[i]->motherType()) == 24 && fabs(mcParticles[i]->grannyType()) == 6 )
  	{
  		if( fabs(mcParticles[i]->type()) >= 1  && fabs(mcParticles[i]->type()) <= 4 )
  		{
  			if(quark1FromW_.P() == 0){
  				quark1FromW_ = *mcParticles[i];
  				continue;
  			}
  			else
  			{
  				quark2FromW_ = *mcParticles[i];
  				continue;
  			}
  		}
  		else if( fabs(mcParticles[i]->type()) == 11 || fabs(mcParticles[i]->type()) == 13 || fabs(mcParticles[i]->type()) == 15 )
  		{
  			if(      fabs(mcParticles[i]->type()) == 11 ) wLeptonicChannel_ = kElec;
  			else if( fabs(mcParticles[i]->type()) == 13 ) wLeptonicChannel_ = kMuon;
  			else if( fabs(mcParticles[i]->type()) == 15 ) wLeptonicChannel_ = kTau;
  
  			leptonFromW_ = *mcParticles[i];
  			continue;
  		}
  		else if( fabs(mcParticles[i]->type()) == 12 || fabs(mcParticles[i]->type()) == 14 || fabs(mcParticles[i]->type()) == 16 )
  		{
  			neutrino_ = *mcParticles[i];
  			continue;
  		}
  /*
  		else {
  			cerr << "Cannot assign W boson decay. Type =" << mcParticles[i]->type() << endl;
  			exit(1);
  		}
  */
  	}
  	else if( fabs(mcParticles[i]->motherType()) == 23 && fabs(mcParticles[i]->grannyType()) == 6 )
  	{
  		if( fabs(mcParticles[i]->type()) == 11 || fabs(mcParticles[i]->type()) == 13 || fabs(mcParticles[i]->type()) == 15 )
  		{
  			if(      fabs(mcParticles[i]->type()) == 11 ) zLeptonicChannel_ = kElec;
  			else if( fabs(mcParticles[i]->type()) == 13 ) zLeptonicChannel_ = kMuon;
  			else if( fabs(mcParticles[i]->type()) == 15 ) zLeptonicChannel_ = kTau;
  
  			if(lepton1FromZ_.E() == 0){
  				lepton1FromZ_ = *mcParticles[i];
  				continue;
  			}
  			else{
  				lepton2FromZ_ = *mcParticles[i];
  				continue;
  			}
  		}
  /*
  		else {
  			cerr << "Cannot assign Z boson decay. Type =" << mcParticles[i]->type() << endl;
  			exit(1);
  		}
  */
  	}
  }
  //if(debug) cout<<(*this)<<endl;
}

void TopFCNC_GenEvt::MatchJetsToPartons(const std::vector<TRootJet*> &jets, const int algorithm, const bool useMaxDist, const bool useDeltaR, const double maxDist){
	if(zLeptonicChannel_ == kNone)
	{
		cerr << "Cannot find the Z leptonic decay channel. Use the method ReconstructEvt first." << endl;
		exit(1);
	}

	std::vector<TLorentzVector> partons;
	std::vector<TLorentzVector> tljets;

	for(unsigned int i=0;i<jets.size();i++)
	{
		tljets.push_back((TLorentzVector)*jets[i]);
	}

	if(wLeptonicChannel_ == kNone)// dileptonic event tt->Wb+Zq->qqb+llq
	{
		partons.push_back(quark1FromW_);
		partons.push_back(quark2FromW_);
		partons.push_back(B_);
		partons.push_back(Q_);
	}
	else // trileptonic events tt->Wb+Zq->lvb+llq
	{
		partons.push_back(B_);
		partons.push_back(Q_);
	}
	
	JetPartonMatching myJetPartonMatcher = JetPartonMatching(partons,tljets,algorithm,useMaxDist,useDeltaR,maxDist);
	//myJetPartonMatcher.print();

	if(wLeptonicChannel_ == kNone){ // dileptonic event tt->Wb+Zq->qqb+llq
	  HadQ1Idx_ = myJetPartonMatcher.getMatchForParton(0,0);
	  HadQ2Idx_ = myJetPartonMatcher.getMatchForParton(1,0);
	  HadBIdx_  = myJetPartonMatcher.getMatchForParton(2,0);
	  HadQIdx_  = myJetPartonMatcher.getMatchForParton(3,0);

		if(HadQ1Idx_ != -1 ) matchedQuark1FromW_ = *jets[HadQ1Idx_];
		if(HadQ2Idx_ != -1 ) matchedQuark2FromW_ = *jets[HadQ2Idx_];
		if(HadBIdx_  != -1 ) matchedB_           = *jets[HadBIdx_];
		if(HadQIdx_  != -1 ) matchedQ_           = *jets[HadQIdx_];
	}
	else{
	  HadBIdx_  = myJetPartonMatcher.getMatchForParton(0,0);
	  HadQIdx_  = myJetPartonMatcher.getMatchForParton(1,0);
	  
		if(HadBIdx_  != -1 ) matchedB_           = *jets[HadBIdx_];
		if(HadQIdx_  != -1 ) matchedQ_           = *jets[HadQIdx_];
	}
}

void TopFCNC_GenEvt::MatchLeptonsToZ(const std::vector<TRootMuon*> &leptons, const int algorithm, const bool useMaxDist, const bool useDeltaR, const double maxDist){
 
	std::vector<TLorentzVector> genlept;
	std::vector<TLorentzVector> tleptons;

	for(unsigned int i=0;i<leptons.size();i++)
	{
		tleptons.push_back((TLorentzVector)*leptons[i]);
	}

  genlept.push_back(lepton1FromZ_);
  genlept.push_back(lepton2FromZ_);

	JetPartonMatching myLeptonMatcher = JetPartonMatching(genlept,tleptons,algorithm,useMaxDist,useDeltaR,maxDist);

  Lep1Idx_ = myLeptonMatcher.getMatchForParton(0,0);
  Lep2Idx_ = myLeptonMatcher.getMatchForParton(1,0);
  
  if(Lep1Idx_ != -1) matchedLepton1FromZ_ = *leptons[Lep1Idx_];
  if(Lep2Idx_ != -1) matchedLepton2FromZ_ = *leptons[Lep2Idx_];

  if(Lep1Idx_ != -1 && Lep2Idx_ != -1)
  	matchedZ_ = matchedLepton1FromZ_+matchedLepton2FromZ_;
}
void TopFCNC_GenEvt::MatchLeptonsToZ(const std::vector<TRootElectron*> &leptons, const int algorithm, const bool useMaxDist, const bool useDeltaR, const double maxDist){
 
	std::vector<TLorentzVector> genlept;
	std::vector<TLorentzVector> tleptons;

	for(unsigned int i=0;i<leptons.size();i++){
		tleptons.push_back((TLorentzVector)*leptons[i]);
	}

  genlept.push_back(lepton1FromZ_);
  genlept.push_back(lepton2FromZ_);

	JetPartonMatching myLeptonMatcher = JetPartonMatching(genlept,tleptons,algorithm,useMaxDist,useDeltaR,maxDist);

  Lep1Idx_ = myLeptonMatcher.getMatchForParton(0,0);
  Lep2Idx_ = myLeptonMatcher.getMatchForParton(1,0);

  if(Lep1Idx_ != -1) matchedLepton1FromZ_ = *leptons[Lep1Idx_];
  if(Lep2Idx_ != -1) matchedLepton2FromZ_ = *leptons[Lep2Idx_];

  if(Lep1Idx_ != -1 && Lep2Idx_ != -1)
  	matchedZ_ = matchedLepton1FromZ_+matchedLepton2FromZ_;
}

void TopFCNC_GenEvt::FillResolutions(ResolutionFit* resFitMuons, ResolutionFit* resFitElectrons, ResolutionFit* resFitBJets, ResolutionFit* resFitQJets, ResolutionFit* resFitLightJets){

  if(zLeptonicChannel_ == kMuon)
  {
    if(Lep1Idx_ != -1) resFitMuons->Fill(&matchedLepton1FromZ_,&lepton1FromZ_);
    if(Lep2Idx_ != -1) resFitMuons->Fill(&matchedLepton2FromZ_,&lepton2FromZ_);
  }
  else if(zLeptonicChannel_ == kElec)
  {
    if(Lep1Idx_ != -1) resFitElectrons->Fill(&matchedLepton1FromZ_,&lepton1FromZ_);
    if(Lep2Idx_ != -1) resFitElectrons->Fill(&matchedLepton2FromZ_,&lepton2FromZ_);
  }

  if(HadBIdx_ != -1) resFitBJets->Fill(&matchedB_,&B_);
  if(HadQIdx_ != -1) resFitQJets->Fill(&matchedQ_,&Q_);

	if(wLeptonicChannel_ == kNone) // dileptonic event tt->Wb+Zq->qqb+llq
	{
		if(HadQ1Idx_ != -1) resFitLightJets->Fill(&matchedQuark1FromW_,&quark1FromW_);
		if(HadQ2Idx_ != -1) resFitLightJets->Fill(&matchedQuark2FromW_,&quark2FromW_);
  }
}
