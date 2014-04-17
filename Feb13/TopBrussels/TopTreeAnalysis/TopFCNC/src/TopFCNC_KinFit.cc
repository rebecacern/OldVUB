#include "../interface/TopFCNC_KinFit.h"

TopFCNC_KinFit::TopFCNC_KinFit(Dataset* dataset, ResolutionFit *resFitLeptons, ResolutionFit *resFitBJets, ResolutionFit *resFitQJets, ResolutionFit *resFitLightJets, float wMass, float zMass, float topMass){

  // Set up datasets
  dataset_ = dataset;
  
  TMatrixD empty3x3(3,3);
  // Set up leptons ___________________________________________________________________________________________________________
  lepton1_   = new TFitParticleEtThetaPhi("Lepton1",   "Lepton1",   0, &empty3x3); // Lepton 1 from Z boson decay
  lepton2_   = new TFitParticleEtThetaPhi("Lepton2",   "Lepton2",   0, &empty3x3); // Lepton 2 from Z boson decay
  
  // Set up jets ______________________________________________________________________________________________________________
  lightJet1_ = new TFitParticleEtThetaPhi("LightJet1", "LightJet1", 0, &empty3x3); // Light Jet 1 from hadronic W boson decay
  lightJet2_ = new TFitParticleEtThetaPhi("LightJet2", "LightJet2", 0, &empty3x3); // Light Jet 2 from hadronic W boson decay

  bJet_      = new TFitParticleEtThetaPhi("BJet",      "BJet",      0, &empty3x3); // B-Jet from hadronic top quark decay
  qJet_      = new TFitParticleEtThetaPhi("QJet",      "QJet",      0, &empty3x3); // (u/c)-Jet from FCNC top quark decay
  
  // Set up constraints ________________________________________________________________________________________________________
  consHadW_     = new TFitConstraintM("HadWMass",      "MassConstraint", 0, 0, wMass);
  consLepZ_     = new TFitConstraintM("LepZMass",      "MassConstraint", 0, 0, zMass);
  consHadTop_   = new TFitConstraintM("Had_TopMass",   "MassConstraint", 0, 0, topMass);
  consFcncTop_  = new TFitConstraintM("FCNC_TopMass",  "MassConstraint", 0, 0, topMass);
  consEqualTop_ = new TFitConstraintM("EqualTopMasses","EqualTopMasses", 0, 0, 0.);
  consSumPx_    = new TFitConstraintEp("SumPx",        "SumPx", 0, TFitConstraintEp::pX, 0.);
  consSumPy_    = new TFitConstraintEp("SumPy",        "SumPy", 0, TFitConstraintEp::pY, 0.);

  // - W boson mass constraint
  consHadW_->addParticles1(lightJet1_,lightJet2_);
  // - Z boson mass constraint
  consLepZ_->addParticles1(lepton1_,lepton2_);
  
  // - top quark mass constraints
  consHadTop_->addParticles1(bJet_,lightJet1_,lightJet2_);
  consFcncTop_->addParticles1(qJet_,lepton1_,lepton2_);
  
  // - equal top quark mass constraints
  consEqualTop_->addParticles1(bJet_,lightJet1_,lightJet2_);
  consEqualTop_->addParticles2(qJet_,lepton1_,lepton2_);
  
  // - additional kinematic constraints
  consSumPx_->addParticles(bJet_,lightJet1_,lightJet2_,qJet_,lepton1_,lepton2_);
  consSumPy_->addParticles(bJet_,lightJet1_,lightJet2_,qJet_,lepton1_,lepton2_);

  // Set up resolutions
  resFitLeptons_   = resFitLeptons;
  resFitBJets_     = resFitBJets;
  resFitQJets_     = resFitQJets;
  resFitLightJets_ = resFitLightJets;

  // Set up fitter ____________________________________________________________________________________________________________
  kinfit_ = new TKinFitter("topFCNC_Fit", "topFCNC_Fit");
  // - add measured particles
  kinfit_->addMeasParticles(bJet_,lightJet1_,lightJet2_,qJet_,lepton1_,lepton2_);

/*
  // - add constraints
  kinfit_->addConstraint(consHadW_);
  kinfit_->addConstraint(consLepZ_);
  kinfit_->addConstraint(consHadTop_);
  kinfit_->addConstraint(consFcncTop_);
  kinfit_->addConstraint(consEqualTop_);
  kinfit_->addConstraint(consSumPx_);
  kinfit_->addConstraint(consSumPy_);
*/
  
  // Set up fit parameters _____________________________________________________________________________________________________
  prob_ = -1.;
  chi2_ =  0.;
  ndof_ =  0;
  maxNbIter_ = 200;
  maxDeltaS_ = 0.00005;
  maxF_      = 0.0001;

  constrainSumPt_= false;
  verbose_       = false;
  verbosity_fit_ = 0;

  // Set up plot parameters ____________________________________________________________________________________________________
  if(dataset_)
  {
    // Make some plots
/*
    string mWUnCorrtitle = "mWUnCorr_"+dataset_->Name();
    string mTopUnCorrtitle = "mTopUnCorr_"+dataset_->Name();
    string mWCorrtitle = "mWCorr_"+dataset_->Name();
    string mTopCorrtitle = "mTopCorr_"+dataset_->Name();
    histo1D_[mWUnCorrtitle] = new TH1F(mWUnCorrtitle.c_str(),mWUnCorrtitle.c_str(),50,0,150);
    histo1D_[mTopUnCorrtitle] = new TH1F(mTopUnCorrtitle.c_str(),mTopUnCorrtitle.c_str(),50,100,250);
    histo1D_[mWCorrtitle] = new TH1F(mWCorrtitle.c_str(),mWCorrtitle.c_str(),50,0,150);
    histo1D_[mTopCorrtitle] = new TH1F(mTopCorrtitle.c_str(),mTopCorrtitle.c_str(),50,100,250);
*/
  }
}

TopFCNC_KinFit::~TopFCNC_KinFit(){
  // Delete fitter
  delete kinfit_;
  
  // Delete leptons
  delete lepton1_;
  delete lepton2_;
  
  // Delete jets
  delete lightJet1_;
  delete lightJet2_;
  
  delete bJet_;
  delete qJet_;
  
  // Delete constraints
  delete consHadW_;
  delete consLepZ_;
  delete consHadTop_;
  delete consFcncTop_;
  delete consEqualTop_;
  delete consSumPx_;
  delete consSumPy_;
}

void TopFCNC_KinFit::SetConstraints(vector<string> &constraints){
  // - add constraints
  for(unsigned int idx=0;idx<constraints.size();idx++){
    if(     constraints[idx].find("kHadWMass")       != string::npos) kinfit_->addConstraint(consHadW_);
    else if(constraints[idx].find("kLepZMass")       != string::npos) kinfit_->addConstraint(consLepZ_);
    else if(constraints[idx].find("kHadTopMass")     != string::npos) kinfit_->addConstraint(consHadTop_);
    else if(constraints[idx].find("kFcncTopMass")    != string::npos) kinfit_->addConstraint(consFcncTop_);
    else if(constraints[idx].find("kEqualTopMasses") != string::npos) kinfit_->addConstraint(consEqualTop_);
    else if(constraints[idx].find("kSumPx")          != string::npos) { kinfit_->addConstraint(consSumPx_); constrainSumPt_ = true; }
    else if(constraints[idx].find("kSumPy")          != string::npos) { kinfit_->addConstraint(consSumPy_); constrainSumPt_ = true; }
    else {
      cerr << "Unkown constraint : " << constraints[idx] << endl;
      cerr << "Available constraints are: " << endl;
      cerr << " - kHadWMass " << endl;
      cerr << " - kLepZMass " << endl;
      cerr << " - kHadTopMass " << endl;
      cerr << " - kFcncTopMass " << endl;
      cerr << " - kEqualTopMasses " << endl;
      cerr << " - kSumPx " << endl;
      cerr << " - kSumPy " << endl;
      exit(1);
    }
  }
}

void TopFCNC_KinFit::FitEvent(TopFCNC_Evt *topFCNC_Evt)
{

  Double_t prob_tmp = 0.;
  prob_ = -1.;
  chi2_ =  0.;
  ndof_ =  0;

  // Set up kinfit parameters
  kinfit_->setVerbosity(verbosity_fit_);
  kinfit_->setMaxNbIter(maxNbIter_);
  kinfit_->setMaxDeltaS(maxDeltaS_);
  kinfit_->setMaxF(maxF_);

  // Set up kinfit mass constraints
/*
  consHadW_->setMassConstraint(wMass);
  consLepZ_->setMassConstraint(zMass);
  consHadTop_->setMassConstraint(topMass);
  consFcncTop_->setMassConstraint(topMass);
*/  
  // Declare TLorentzVector for fitted particles
  TLorentzVector lepton1FromZ_fitted, lepton2FromZ_fitted, qJet_fitted, lightJet1_fitted, lightJet2_fitted, bJet_fitted;
  
  // Leptons ____________________________________________________________
  
  // Set up lepton kinematics
  TLorentzVector lepton1FromZ = topFCNC_Evt->lepton1FromZ();
  lepton1_->setIni4Vec(&lepton1FromZ);
  TLorentzVector lepton2FromZ = topFCNC_Evt->lepton2FromZ();
  lepton2_->setIni4Vec(&lepton2FromZ);

  // Set up lepton resolutions
  TMatrixD Ml1(3,3), Ml2(3,3);

  Ml1.Zero();
  Ml1(0,0) = pow(resFitLeptons_->EtResolution(&lepton1FromZ), 2);
  Ml1(1,1) = pow(resFitLeptons_->ThetaResolution(&lepton1FromZ), 2);
  Ml1(2,2) = pow(resFitLeptons_->PhiResolution(&lepton1FromZ), 2);
  lepton1_->setCovMatrix(&Ml1);
  
  Ml2.Zero();
  Ml2(0,0) = pow(resFitLeptons_->EtResolution(&lepton2FromZ), 2);
  Ml2(1,1) = pow(resFitLeptons_->ThetaResolution(&lepton2FromZ), 2);
  Ml2(2,2) = pow(resFitLeptons_->PhiResolution(&lepton2FromZ), 2);
  lepton2_->setCovMatrix(&Ml2);

  if(verbose_){
    Ml1.Print();
    Ml2.Print();
  }
  
  // Jets ____________________________________________________________

  // Set up seleted jet
  vector<TRootJet> jets  = topFCNC_Evt->selectedJets();
  if(verbose_)
	  cout<<"- Number of selected jets : "<<jets.size()<<endl;
  // Set up resolution matrices
  TMatrixD Mb(3,3),  Mq(3,3), Mj1(3,3), Mj2(3,3);

  // Set up permutations
  UInt_t *numbers = new UInt_t[jets.size()];
  for(UInt_t i=0;i<jets.size();i++) numbers[i]=i;
  UInt_t *Comb = 0;
  UInt_t *Permutation = 0;

  if(topFCNC_Evt->isDiLeptonic()){
    // Topology to reconstruct : tt-> bW + qZ -> bqq + qll

    UInt_t NofJets = 4 ;

    Comb = new UInt_t[NofJets];
    for(UInt_t i=0;i <NofJets;i++) Comb[i]=i;
	  Permutation = new UInt_t[NofJets];
    for(UInt_t i=0;i<NofJets;i++) Permutation[i]=i;

    do{ // - loop over all combinations
      if(verbose_)
		    cout<<"-- Jet Combination considered : "<<Comb[0]<<"/"<<Comb[1]<<"/"<<Comb[2]<<"/"<<Comb[3]<<endl;

		  prob_tmp = 0.;

	    do{ // - loop over all permutations for a given combination
		    if(verbose_)
		  	  cout<<"--- Permutations : "<<Comb[0]<<"/"<<Comb[1]<<"/"<<Comb[2]<<"/"<<Comb[3]<<endl;
    
        // - set up current jet configuration
        TLorentzVector bJet      = jets[Comb[0]];
        TLorentzVector qJet      = jets[Comb[1]];
        TLorentzVector lightJet1 = jets[Comb[2]];
        TLorentzVector lightJet2 = jets[Comb[3]];
        // - set up jet kinematics
        bJet_->setIni4Vec(&bJet);
        qJet_->setIni4Vec(&qJet);
        lightJet1_->setIni4Vec(&lightJet1);
        lightJet2_->setIni4Vec(&lightJet2);
        
        // - reset resolution matrices
        Mb.Zero();
        Mq.Zero();
        Mj1.Zero();
        Mj2.Zero();

        // - fill resolution matrices
        Mb(0,0)  = pow(resFitBJets_->EtResolution(&bJet), 2);
        Mb(1,1)  = pow(resFitBJets_->ThetaResolution(&bJet), 2);
        Mb(2,2)  = pow(resFitBJets_->PhiResolution(&bJet), 2);
		    if(verbose_)
          Mb.Print();

        Mq(0,0)  = pow(resFitQJets_->EtResolution(&qJet), 2);
        Mq(1,1)  = pow(resFitQJets_->ThetaResolution(&qJet), 2);
        Mq(2,2)  = pow(resFitQJets_->PhiResolution(&qJet), 2);
		    if(verbose_)
          Mq.Print();

        Mj1(0,0) = pow(resFitLightJets_->EtResolution(&lightJet1), 2);
        Mj1(1,1) = pow(resFitLightJets_->ThetaResolution(&lightJet1), 2);
        Mj1(2,2) = pow(resFitLightJets_->PhiResolution(&lightJet1), 2);
		    if(verbose_)
          Mj1.Print();
    
        Mj2(0,0) = pow(resFitLightJets_->EtResolution(&lightJet2), 2);
        Mj2(1,1) = pow(resFitLightJets_->ThetaResolution(&lightJet2), 2);
        Mj2(2,2) = pow(resFitLightJets_->PhiResolution(&lightJet2), 2);
		    if(verbose_)
          Mj2.Print();

        // - set up jet resolutions
        bJet_->setCovMatrix(&Mb);
        qJet_->setCovMatrix(&Mq);
        lightJet1_->setCovMatrix(&Mj1);
        lightJet2_->setCovMatrix(&Mj2);

        // - set up Px and Py constraint for curent event configuration so that sum Pt will be conserved
        if(constrainSumPt_){
          consSumPx_->setConstraint( lepton1FromZ.Px() + lepton2FromZ.Px() + qJet.Px() + lightJet1.Px() + lightJet2.Px() + bJet.Px() );
          consSumPy_->setConstraint( lepton1FromZ.Py() + lepton2FromZ.Py() + qJet.Py() + lightJet1.Py() + lightJet2.Py() + bJet.Py() );
        }
		    if(verbose_)
		  	  cout<<"---- Jet covariance matrices instantiated "<<endl;

        // - perform the fit!
        kinfit_->fit();
        if(kinfit_->getStatus() == 0) // if the fitter converged
          prob_tmp = TMath::Prob(kinfit_->getS(), kinfit_->getNDF());
        else
          prob_tmp = -1.;

        if(prob_tmp>prob_){
          for(UInt_t i=0;i<NofJets;i++) Permutation[i] = Comb[i];
          prob_ = prob_tmp;
          chi2_ = kinfit_->getS();
          ndof_ = kinfit_->getNDF();
          lepton1FromZ_fitted = *lepton1_  ->getCurr4Vec();
          lepton2FromZ_fitted = *lepton2_  ->getCurr4Vec();
          qJet_fitted         = *qJet_     ->getCurr4Vec();
          lightJet1_fitted    = *lightJet1_->getCurr4Vec();
          lightJet2_fitted    = *lightJet2_->getCurr4Vec();
          bJet_fitted         = *bJet_     ->getCurr4Vec();
        }
			}
		  while(next_permutation(Comb,Comb+NofJets));
   	}
    while(next_combination(numbers,numbers+jets.size(),Comb,Comb+NofJets));

    topFCNC_Evt->SetB(bJet_fitted);
    topFCNC_Evt->SetQ(qJet_fitted);
    topFCNC_Evt->SetQuark1FromW(lightJet1_fitted);
    topFCNC_Evt->SetQuark2FromW(lightJet2_fitted);

    delete Permutation;
    delete Comb;
  }
  else{
    // Topology to reconstruct : tt-> bW + qZ -> blv + qll
    prob_ = 0;
    chi2_ = 0;
    ndof_ = 0;
  }
}

void TopFCNC_KinFit::FitEvent(TopFCNC_GenEvt *topFCNC_GenEvt)
{
  
  prob_ = -1.;
  chi2_ =  0.;
  ndof_ =  0;
  
  // Set up kinfit parameters
  kinfit_->setVerbosity(verbosity_fit_);
  kinfit_->setMaxNbIter(maxNbIter_);
  kinfit_->setMaxDeltaS(maxDeltaS_);
  kinfit_->setMaxF(maxF_);
  
  
  // Declare TLorentzVector for fitted particles
  TLorentzVector lepton1FromZ_fitted, lepton2FromZ_fitted, qJet_fitted, lightJet1_fitted, lightJet2_fitted, bJet_fitted;
  
  // Leptons ____________________________________________________________
  
  // Set up lepton kinematics
  TLorentzVector lepton1FromZ = topFCNC_GenEvt->matchedLepton1FromZ();
  if(!lepton1FromZ.E()>0){
    cout << "Couldn't find matched B jet" << endl;
    return;
  }
  lepton1_->setIni4Vec(&lepton1FromZ);
  TLorentzVector lepton2FromZ = topFCNC_GenEvt->matchedLepton2FromZ();
  if(!lepton2FromZ.E()>0){
    cout << "Couldn't find matched B jet" << endl;
    return;
  }
  lepton2_->setIni4Vec(&lepton2FromZ);
  
  // Set up lepton resolutions
  TMatrixD Ml1(3,3), Ml2(3,3);
  
  Ml1.Zero();
  Ml1(0,0) = pow(resFitLeptons_->EtResolution(&lepton1FromZ), 2);
  Ml1(1,1) = pow(resFitLeptons_->ThetaResolution(&lepton1FromZ), 2);
  Ml1(2,2) = pow(resFitLeptons_->PhiResolution(&lepton1FromZ), 2);
  lepton1_->setCovMatrix(&Ml1);
  
  Ml2.Zero();
  Ml2(0,0) = pow(resFitLeptons_->EtResolution(&lepton2FromZ), 2);
  Ml2(1,1) = pow(resFitLeptons_->ThetaResolution(&lepton2FromZ), 2);
  Ml2(2,2) = pow(resFitLeptons_->PhiResolution(&lepton2FromZ), 2);
  lepton2_->setCovMatrix(&Ml2);
  
  if(verbose_){
    Ml1.Print();
    Ml2.Print();
  }
  
  // Jets ____________________________________________________________
  
  // Set up (b/q)jet kinematics
  TLorentzVector bJet = topFCNC_GenEvt->matchedB();
  if(!bJet.E()>0){
    cout << "Couldn't find matched B jet" << endl;
    return;
  }
  bJet_->setIni4Vec(&bJet);
  TLorentzVector qJet = topFCNC_GenEvt->matchedQ();
  if(!qJet.E()>0){
    cout << "Couldn't find matched Q jet" << endl;
    return;
  }
  qJet_->setIni4Vec(&qJet);
  
  // Set up resolution matrices
  TMatrixD Mb(3,3),  Mq(3,3), Mj1(3,3), Mj2(3,3);
  // - reset resolution matrices
  Mb.Zero();
  Mq.Zero();
  Mj1.Zero();
  Mj2.Zero();
  
  if(topFCNC_GenEvt->isDiLeptonic()){
    // Topology to reconstruct : tt-> bW + qZ -> bqq + qll
    
    // Set up light-jet kinematics
    TLorentzVector lightJet1 = topFCNC_GenEvt->matchedQuark1FromW();
    if(!lightJet1.E()>0.){
      cout << "Couldn't find matched light jet 1" << endl;
      return;
    }
    lightJet1_->setIni4Vec(&lightJet1);
    
    TLorentzVector lightJet2 = topFCNC_GenEvt->matchedQuark2FromW();
    if(!lightJet2.E()>0.){
      cout << "Couldn't find matched light jet 2" << endl;
      return;
    }
    lightJet2_->setIni4Vec(&lightJet2);
    
    // - fill resolution matrices
    Mb(0,0)  = pow(resFitBJets_->EtResolution(&bJet), 2);
    Mb(1,1)  = pow(resFitBJets_->ThetaResolution(&bJet), 2);
    Mb(2,2)  = pow(resFitBJets_->PhiResolution(&bJet), 2);
		if(verbose_)
      Mb.Print();
    
    Mq(0,0)  = pow(resFitQJets_->EtResolution(&qJet), 2);
    Mq(1,1)  = pow(resFitQJets_->ThetaResolution(&qJet), 2);
    Mq(2,2)  = pow(resFitQJets_->PhiResolution(&qJet), 2);
		if(verbose_)
      Mq.Print();
    
    Mj1(0,0) = pow(resFitLightJets_->EtResolution(&lightJet1), 2);
    Mj1(1,1) = pow(resFitLightJets_->ThetaResolution(&lightJet1), 2);
    Mj1(2,2) = pow(resFitLightJets_->PhiResolution(&lightJet1), 2);
    if(verbose_)
      Mj1.Print();
    
    Mj2(0,0) = pow(resFitLightJets_->EtResolution(&lightJet2), 2);
    Mj2(1,1) = pow(resFitLightJets_->ThetaResolution(&lightJet2), 2);
    Mj2(2,2) = pow(resFitLightJets_->PhiResolution(&lightJet2), 2);
    if(verbose_)
      Mj2.Print();
    
    // - set up jet resolutions
    bJet_->setCovMatrix(&Mb);
    qJet_->setCovMatrix(&Mq);
    lightJet1_->setCovMatrix(&Mj1);
    lightJet2_->setCovMatrix(&Mj2);
    
    // - set up Px and Py constraint for curent event configuration so that sum Pt will be conserved
    if(constrainSumPt_){
      consSumPx_->setConstraint( lepton1FromZ.Px() + lepton2FromZ.Px() + qJet.Px() + lightJet1.Px() + lightJet2.Px() + bJet.Px() );
      consSumPy_->setConstraint( lepton1FromZ.Py() + lepton2FromZ.Py() + qJet.Py() + lightJet1.Py() + lightJet2.Py() + bJet.Py() );
    }
    if(verbose_)
      cout<<"---- Jet covariance matrices instantiated "<<endl;
    
    // - perform the fit!
    kinfit_->fit();
    if(kinfit_->getStatus() == 0){ // if the fitter converged
      chi2_ = kinfit_->getS();
      ndof_ = kinfit_->getNDF();
      prob_ = TMath::Prob(chi2_, ndof_);
      lepton1FromZ_fitted = *lepton1_  ->getCurr4Vec();
      lepton2FromZ_fitted = *lepton2_  ->getCurr4Vec();
      qJet_fitted         = *qJet_     ->getCurr4Vec();
      lightJet1_fitted    = *lightJet1_->getCurr4Vec();
      lightJet2_fitted    = *lightJet2_->getCurr4Vec();
      bJet_fitted         = *bJet_     ->getCurr4Vec();
      
    }
    else
      prob_ = -1.;
    
  }
  else{
    // Topology to reconstruct : tt-> bW + qZ -> blv + qll
    prob_ = -1;
    chi2_ = 0;
    ndof_ = 0;
  }
}
void TopFCNC_KinFit::FitEvent(TLorentzVector &lepton1, TLorentzVector &lepton2, TLorentzVector &qJet, TLorentzVector &lightJet1, TLorentzVector &lightJet2, TLorentzVector &bJet)
{
  // Topology to reconstruct : tt-> bW + qZ -> bqq + qll
 
  prob_ = -1.;
  chi2_ =  0.;
  ndof_ =  0;
  
  // Set up kinfit parameters
  kinfit_->setVerbosity(verbosity_fit_);
  kinfit_->setMaxNbIter(maxNbIter_);
  kinfit_->setMaxDeltaS(maxDeltaS_);
  kinfit_->setMaxF(maxF_);
  
  
  // Declare TLorentzVector for fitted particles
  TLorentzVector lepton1FromZ_fitted, lepton2FromZ_fitted, qJet_fitted, lightJet1_fitted, lightJet2_fitted, bJet_fitted;
  
  // Leptons ____________________________________________________________
  
  // Set up lepton kinematics
  if(!lepton1.E()>0){
    cout << "Couldn't find matched B jet" << endl;
    return;
  }
  lepton1_->setIni4Vec(&lepton1);
  if(!lepton2.E()>0){
    cout << "Couldn't find matched B jet" << endl;
    return;
  }
  lepton2_->setIni4Vec(&lepton2);
  
  // Set up lepton resolutions
  TMatrixD Ml1(3,3), Ml2(3,3);
  
  Ml1.Zero();
  Ml1(0,0) = pow(resFitLeptons_->EtResolution(&lepton1), 2);
  Ml1(1,1) = pow(resFitLeptons_->ThetaResolution(&lepton1), 2);
  Ml1(2,2) = pow(resFitLeptons_->PhiResolution(&lepton1), 2);
  lepton1_->setCovMatrix(&Ml1);
  
  Ml2.Zero();
  Ml2(0,0) = pow(resFitLeptons_->EtResolution(&lepton2), 2);
  Ml2(1,1) = pow(resFitLeptons_->ThetaResolution(&lepton2), 2);
  Ml2(2,2) = pow(resFitLeptons_->PhiResolution(&lepton2), 2);
  lepton2_->setCovMatrix(&Ml2);
  
  if(verbose_){
    Ml1.Print();
    Ml2.Print();
  }
  
  // Jets ____________________________________________________________
  
  // Set up (b/q)jet kinematics
  if(!bJet.E()>0){
    cout << "Couldn't find matched B jet" << endl;
    return;
  }
  bJet_->setIni4Vec(&bJet);

  if(!qJet.E()>0){
    cout << "Couldn't find matched Q jet" << endl;
    return;
  }
  qJet_->setIni4Vec(&qJet);

  // Set up light-jet kinematics
  if(!lightJet1.E()>0.){
    cout << "Couldn't find matched light jet 1" << endl;
    return;
  }
  lightJet1_->setIni4Vec(&lightJet1);
  
  if(!lightJet2.E()>0.){
    cout << "Couldn't find matched light jet 2" << endl;
    return;
  }
  lightJet2_->setIni4Vec(&lightJet2);
  
  // Set up resolution matrices
  TMatrixD Mb(3,3),  Mq(3,3), Mj1(3,3), Mj2(3,3);
  // - reset resolution matrices
  Mb.Zero();
  Mq.Zero();
  Mj1.Zero();
  Mj2.Zero();
  
  // - fill resolution matrices
  Mb(0,0)  = pow(resFitBJets_->EtResolution(&bJet), 2);
  Mb(1,1)  = pow(resFitBJets_->ThetaResolution(&bJet), 2);
  Mb(2,2)  = pow(resFitBJets_->PhiResolution(&bJet), 2);
  if(verbose_)
    Mb.Print();
  
  Mq(0,0)  = pow(resFitQJets_->EtResolution(&qJet), 2);
  Mq(1,1)  = pow(resFitQJets_->ThetaResolution(&qJet), 2);
  Mq(2,2)  = pow(resFitQJets_->PhiResolution(&qJet), 2);
  if(verbose_)
    Mq.Print();
  
  Mj1(0,0) = pow(resFitLightJets_->EtResolution(&lightJet1), 2);
  Mj1(1,1) = pow(resFitLightJets_->ThetaResolution(&lightJet1), 2);
  Mj1(2,2) = pow(resFitLightJets_->PhiResolution(&lightJet1), 2);
  if(verbose_)
    Mj1.Print();
  
  Mj2(0,0) = pow(resFitLightJets_->EtResolution(&lightJet2), 2);
  Mj2(1,1) = pow(resFitLightJets_->ThetaResolution(&lightJet2), 2);
  Mj2(2,2) = pow(resFitLightJets_->PhiResolution(&lightJet2), 2);
  if(verbose_)
    Mj2.Print();
  
  // - set up jet resolutions
  bJet_->setCovMatrix(&Mb);
  qJet_->setCovMatrix(&Mq);
  lightJet1_->setCovMatrix(&Mj1);
  lightJet2_->setCovMatrix(&Mj2);
  
  // - set up Px and Py constraint for curent event configuration so that sum Pt will be conserved
  if(constrainSumPt_){
    consSumPx_->setConstraint( lepton1.Px() + lepton2.Px() + qJet.Px() + lightJet1.Px() + lightJet2.Px() + bJet.Px() );
    consSumPy_->setConstraint( lepton1.Py() + lepton2.Py() + qJet.Py() + lightJet1.Py() + lightJet2.Py() + bJet.Py() );
  }
  if(verbose_)
    cout<<"---- Jet covariance matrices instantiated "<<endl;
  
  // - perform the fit!
  kinfit_->fit();
  if(kinfit_->getStatus() == 0){ // if the fitter converged
    chi2_ = kinfit_->getS();
    ndof_ = kinfit_->getNDF();
    prob_ = TMath::Prob(chi2_, ndof_);
    lepton1FromZ_fitted = *lepton1_  ->getCurr4Vec();
    lepton2FromZ_fitted = *lepton2_  ->getCurr4Vec();
    qJet_fitted         = *qJet_     ->getCurr4Vec();
    lightJet1_fitted    = *lightJet1_->getCurr4Vec();
    lightJet2_fitted    = *lightJet2_->getCurr4Vec();
    bJet_fitted         = *bJet_     ->getCurr4Vec();
    
  }
  else
    prob_ = -1.;
  
}

void TopFCNC_KinFit::Write(TFile* fout, bool savePNG, string pathPNG){

}
