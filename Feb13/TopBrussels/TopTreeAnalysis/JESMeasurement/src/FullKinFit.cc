#include "../interface/FullKinFit.h"

#include <fstream>
#include <time.h>

FullKinFit::FullKinFit(Dataset* d, ResolutionFit *resFitLightJets, ResolutionFit *resFitBJets, bool measureTopMass, bool measureTopMassDiff, bool dummy)
{
  resFitLightJets_ = resFitLightJets;
  resFitBJets_ = resFitBJets;
  resFitLightJetsL7_ = 0;
  resFitBJetsL7_ = 0;
  nBinsCorrB_ = 41;
  stepCorrB_ = 0.02;
  startCorrB_ = 0.6;
  nBinsCorrLight_ = 41;
  stepCorrLight_ = 0.02;
  startCorrLight_ = 0.6;
  measureTopMass_ = measureTopMass;
  measureTopMassDiff_ = measureTopMassDiff;
  stepTopMass_ = 2;
  if(measureTopMass_)
  {
    nBinsTopMass_ = 26; //original: 26 bins, start at 150, width = 2
    startTopMass_ = 150.;
  }
  else if(measureTopMassDiff_)
  {
    nBinsCorrLight_ = 1;
    startCorrLight_ = 1.;
    nBinsTopMass_ = 51;
    startTopMass_ = 125.;
  }
  nTotalJetCombis_ = 0;
  nParabolaFitJetCombis_ = 0;
  nFinalKinFitJetCombis_ = 0;
  nFinalJetCombis_ = 0;
  
  dataset_ = d;
  
  if(dataset_)
  {
    if(measureTopMassDiff_)
    {
      string wJetsSFTitle = "wJetsSF_"+dataset_->Name();
      histo1D_[wJetsSFTitle] = new TH1F(wJetsSFTitle.c_str(),wJetsSFTitle.c_str(),100,0,2);
      string startTopMassTitle = "startTopMass_"+dataset_->Name();
      histo1D_[startTopMassTitle] = new TH1F(startTopMassTitle.c_str(),startTopMassTitle.c_str(),120,0,600);
      string massUncTitle = "massUnc_"+dataset_->Name();
      histo1D_[massUncTitle] = new TH1F(massUncTitle.c_str(),massUncTitle.c_str(),100,0,100);
    }
    
    // Make some plots
    string mWUnCorrtitle = "mWUnCorr_"+dataset_->Name();
    string mTopUnCorrtitle = "mTopUnCorr_"+dataset_->Name();
    string mWCorrtitle = "mWCorr_"+dataset_->Name();
    string mTopCorrtitle = "mTopCorr_"+dataset_->Name();
    histo1D_[mWUnCorrtitle] = new TH1F(mWUnCorrtitle.c_str(),mWUnCorrtitle.c_str(),50,0,150);
    histo1D_[mTopUnCorrtitle] = new TH1F(mTopUnCorrtitle.c_str(),mTopUnCorrtitle.c_str(),50,100,250);
    histo1D_[mWCorrtitle] = new TH1F(mWCorrtitle.c_str(),mWCorrtitle.c_str(),50,0,150);
    histo1D_[mTopCorrtitle] = new TH1F(mTopCorrtitle.c_str(),mTopCorrtitle.c_str(),50,100,250);
  }
  
  string MonsterDir = "Monsters/";
  mkdir(MonsterDir.c_str(),0777);

  mkdir("PlotsJES",0777);
  string path = "PlotsJES/Monsters/";
  mkdir(path.c_str(),0777);
}
  
FullKinFit::~FullKinFit()
{
}

void FullKinFit::SetResFitL7(ResolutionFit *resFitLightJetsL7, ResolutionFit *resFitBJetsL7)
{
  resFitLightJetsL7_ = resFitLightJetsL7;
  resFitBJetsL7_ = resFitBJetsL7;
}

void FullKinFit::SetJets(vector<TRootJet*> jets)
{
  jets_.clear();
  jets_ = jets;
}

void FullKinFit::SetMVAStuff(pair<float, vector<unsigned int> > MVAvals)
{
  MVAvals_ = MVAvals;
}

TH2F* FullKinFit::FitEvent(TRootEvent* event, float WMass, float topMass, bool writePNG, int jetCombi)
{
  int hadrBJetIndex = MVAvals_.second[2], lightJet1Index = MVAvals_.second[0], lightJet2Index = MVAvals_.second[1];
  TLorentzVector bJet = (*jets_[hadrBJetIndex]);
  TLorentzVector lightJet1 = (*jets_[lightJet1Index]);
	TLorentzVector lightJet2 = (*jets_[lightJet2Index]);
	
  if(measureTopMass_ || measureTopMassDiff_)
  {
    float mWUnCorr = ( lightJet1 + lightJet2 ).M();
    float mTopUnCorr = ( lightJet1 + lightJet2 + bJet ).M();

    float lightCorr1 = 1 - resFitLightJets_->EtCorrection(&lightJet1);
    float lightCorr2 = 1 - resFitLightJets_->EtCorrection(&lightJet2);
    float bCorr = 1 - resFitBJets_->EtCorrection(&bJet);
    
//    cout << lightCorr1 << " " << lightCorr2 << " " << bCorr << endl;
//    cout << "test: " << resFitLightJetsL7_->EtResolution(&lightJet1) << endl;
    
    lightJet1 = lightCorr1 * lightJet1;
    lightJet2 = lightCorr2 * lightJet2;
    bJet = bCorr * bJet;
    
    float mWCorr = ( lightJet1 + lightJet2 ).M();
    float mTopCorr = ( lightJet1 + lightJet2 + bJet ).M();
    
    string mWUnCorrtitle = "mWUnCorr_"+dataset_->Name();
    string mTopUnCorrtitle = "mTopUnCorr_"+dataset_->Name();
    string mWCorrtitle = "mWCorr_"+dataset_->Name();
    string mTopCorrtitle = "mTopCorr_"+dataset_->Name();
    histo1D_[mWUnCorrtitle]->Fill(mWUnCorr);
    histo1D_[mTopUnCorrtitle]->Fill(mTopUnCorr);
    histo1D_[mWCorrtitle]->Fill(mWCorr);
    histo1D_[mTopCorrtitle]->Fill(mTopCorr);
  }
  
  TH2F *histo = DummyMonster(jetCombi);

  // prepare everything for the Kinematic Fit
  TMatrixD Ml1(3,3), Ml2(3,3), Mb(3,3);
  Ml1.Zero(); Ml2.Zero(); Mb.Zero();
  if(measureTopMass_ || measureTopMassDiff_)
  {
    Ml1(0,0) = pow(resFitLightJetsL7_->EtResolution(&lightJet1, true), 2);
    Ml1(1,1) = pow(resFitLightJetsL7_->ThetaResolution(&lightJet1, true), 2);
    Ml1(2,2) = pow(resFitLightJetsL7_->PhiResolution(&lightJet1, true), 2);
    Ml2(0,0) = pow(resFitLightJetsL7_->EtResolution(&lightJet2, true), 2);
    Ml2(1,1) = pow(resFitLightJetsL7_->ThetaResolution(&lightJet2, true), 2);
    Ml2(2,2) = pow(resFitLightJetsL7_->PhiResolution(&lightJet2, true), 2);
    Mb(0,0) = pow(resFitBJetsL7_->EtResolution(&bJet, true), 2);
    Mb(1,1) = pow(resFitBJetsL7_->ThetaResolution(&bJet, true), 2);
    Mb(2,2) = pow(resFitBJetsL7_->PhiResolution(&bJet, true), 2);
  }
  else
  {
    Ml1(0,0) = pow(resFitLightJets_->EtResolution(&lightJet1, true), 2);
    Ml1(1,1) = pow(resFitLightJets_->ThetaResolution(&lightJet1, true), 2);
    Ml1(2,2) = pow(resFitLightJets_->PhiResolution(&lightJet1, true), 2);
    Ml2(0,0) = pow(resFitLightJets_->EtResolution(&lightJet2, true), 2);
    Ml2(1,1) = pow(resFitLightJets_->ThetaResolution(&lightJet2, true), 2);
    Ml2(2,2) = pow(resFitLightJets_->PhiResolution(&lightJet2, true), 2);
    Mb(0,0) = pow(resFitBJets_->EtResolution(&bJet, true), 2);
    Mb(1,1) = pow(resFitBJets_->ThetaResolution(&bJet, true), 2);
    Mb(2,2) = pow(resFitBJets_->PhiResolution(&bJet, true), 2);
  }
  
  if(measureTopMass_ || measureTopMassDiff_)
  {
    for(unsigned int binCorrLight=0; binCorrLight<nBinsCorrLight_; binCorrLight++)
    {
      for(unsigned int binTopMass=0; binTopMass<nBinsTopMass_; binTopMass++)
      {
        float DE = 1;
        if(measureTopMass_)
          DE = DE * ( startCorrLight_ + binCorrLight * stepCorrLight_ );
        float DEb = 1;
        float DEl = DE;
        
        TLorentzVector bJetKinFit = bJet * DEb;
        TLorentzVector lightJet1KinFit = lightJet1 * DEl;
        TLorentzVector lightJet2KinFit = lightJet2 * DEl;
        
        float shiftedTopMass = startTopMass_ + binTopMass * stepTopMass_;

//        cout << "DEb: " << DEb << "  DEl: " << DEl << "  shiftedTopMass: " << shiftedTopMass << endl;
        
        TKinFitter *theFitter = new TKinFitter("hadtopFit", "hadtopFit");
        theFitter->setVerbosity(0);
        
        TFitParticleEtThetaPhiEMomFix *fitLight1 = new TFitParticleEtThetaPhiEMomFix("lightJet1", "lightJet1", &lightJet1KinFit, &Ml1);
        TFitParticleEtThetaPhiEMomFix *fitLight2 = new TFitParticleEtThetaPhiEMomFix("lightJet2", "lightJet2", &lightJet2KinFit, &Ml2);
        TFitParticleEtThetaPhiEMomFix *fitB = new TFitParticleEtThetaPhiEMomFix("bJet", "bJet", &bJetKinFit, &Mb);
        theFitter->addMeasParticles(fitLight1,fitLight2,fitB);
        
        TFitConstraintM *consW = new TFitConstraintM(WMass, "MassConstraint", 0, 0, WMass );
        TFitConstraintM *consTop = new TFitConstraintM(topMass, "MassConstraint", 0, 0, shiftedTopMass );
        consW->addParticles1(fitLight1,fitLight2);
        consTop->addParticles1(fitB,fitLight1,fitLight2);

        theFitter->addConstraint(consW);
        theFitter->addConstraint(consTop);
        theFitter->setMaxNbIter(30);
        theFitter->setMaxDeltaS(5e-5);
        theFitter->setMaxF(1e-4);
        
        //do the fit!
        theFitter->fit();
        
        if(measureTopMassDiff_)
        {
          if (theFitter->getStatus() == 0) // if the fitter converged
            histo->Fill(DE, shiftedTopMass, theFitter->getS());
//          else
//            histo->Fill(DE, shiftedTopMass, 9999);
        }
        else if (theFitter->getStatus() == 0) // if the fitter converged
          histo->Fill(DE, shiftedTopMass, TMath::Prob(theFitter->getS(), theFitter->getNDF()) );
        
        delete theFitter;
        delete fitLight1;
        delete fitLight2;
        delete fitB;
        delete consW;
        delete consTop;
      }
    }
  }
  else // Measure b and light JEC
  {
    for(unsigned int binCorrB=0; binCorrB<nBinsCorrB_; binCorrB++)
    {
      for(unsigned int binCorrLight=0; binCorrLight<nBinsCorrLight_; binCorrLight++)
      {
        float DEb = startCorrB_ + binCorrB * stepCorrB_;
        float DElight = startCorrLight_ + binCorrLight * stepCorrLight_;
        TLorentzVector bJetFit = bJet * DEb;
        TLorentzVector lightJet1Fit = lightJet1 * DElight;
        TLorentzVector lightJet2Fit = lightJet2 * DElight;
        
        TKinFitter *theFitter = new TKinFitter("hadtopFit", "hadtopFit");
        theFitter->setVerbosity(0);
        
        TFitParticleEtThetaPhiEMomFix *fitLight1 = new TFitParticleEtThetaPhiEMomFix("lightJet1", "lightJet1", &lightJet1Fit, &Ml1);
        TFitParticleEtThetaPhiEMomFix *fitLight2 = new TFitParticleEtThetaPhiEMomFix("lightJet2", "lightJet2", &lightJet2Fit, &Ml2);
        TFitParticleEtThetaPhiEMomFix *fitB = new TFitParticleEtThetaPhiEMomFix("bJet", "bJet", &bJetFit, &Mb);
        theFitter->addMeasParticles(fitLight1,fitLight2,fitB);
   	    
        TFitConstraintM *consW = new TFitConstraintM(WMass, "MassConstraint", 0, 0, WMass );
        TFitConstraintM *consTop = new TFitConstraintM(topMass, "MassConstraint", 0, 0, topMass );
        consW->addParticles1(fitLight1,fitLight2);
        consTop->addParticles1(fitB,fitLight1,fitLight2);
        
        theFitter->addConstraint(consW);
        theFitter->addConstraint(consTop);
        theFitter->setMaxNbIter(30);
        theFitter->setMaxDeltaS(5e-5);
        theFitter->setMaxF(1e-4);
        
        //do the fit!
        theFitter->fit();
        if (theFitter->getStatus() == 0) // if the fitter converged
          histo->Fill(DElight, DEb, TMath::Prob(theFitter->getS(), theFitter->getNDF()) );
        //else
        //  cout << "FIT NOT CONVERGED" << endl;
        
        delete theFitter;
        delete fitLight1;
        delete fitLight2;
        delete fitB;
        delete consW;
        delete consTop;
        
      } // loop over light corrections
    } // loop over b corrections
  }
  
  if(writePNG)
  {
    stringstream s1; s1 << event->runId();
    stringstream s2; s2 << event->lumiBlockId();
    stringstream s3; s3 << event->eventId();
    stringstream s4; s4 << jetCombi;
    string monsterName = "Monster_Data_" + s1.str() + "_" + s2.str() + "_" + s3.str() + "_" + s4.str();
    
    mkdir("PlotsJES",0777);
    string path = "PlotsJES/Monsters/";
    mkdir(path.c_str(),0777);
    
    if(measureTopMassDiff_)
    {
      TH1F* histo1D = (TH1F*) histo->ProjectionY(monsterName.c_str());
      TCanvas* tempCanvas = TCanvasCreator(histo1D, monsterName);
//      tempCanvas->SaveAs( (path + monsterName + ".root").c_str() );
      tempCanvas->SaveAs( (path + monsterName + ".png").c_str() );
      delete histo1D;
    }
    else
    {
      TCanvas* tempCanvas = TCanvasCreator(histo, monsterName);
//      tempCanvas->SaveAs( (path + monsterName + ".root").c_str() );
      tempCanvas->SaveAs( (path + monsterName + ".png").c_str() );
    }
  }
  
  return histo;
}

TH2F* FullKinFit::DummyMonster(int jetCombi)
{
  TH2F* histo=0;
	stringstream s1;
	s1 << jetCombi;
	string histoName = "KinFitResults_" + s1.str();
  if(measureTopMass_)
    histo = new TH2F(histoName.c_str(),histoName.c_str(), nBinsCorrLight_, startCorrLight_ - stepCorrLight_/2, startCorrLight_ + stepCorrLight_*nBinsCorrLight_ - stepCorrLight_/2, nBinsTopMass_, startTopMass_ - stepTopMass_/2, startTopMass_ + stepTopMass_*nBinsTopMass_ - stepTopMass_/2);
  else if(measureTopMassDiff_)
    histo = new TH2F(histoName.c_str(),histoName.c_str(), 1, 0.99, 1.01, nBinsTopMass_, startTopMass_ - stepTopMass_/2, startTopMass_ + stepTopMass_*nBinsTopMass_ - stepTopMass_/2);
  else
    histo = new TH2F(histoName.c_str(),histoName.c_str(), nBinsCorrLight_, startCorrLight_ - stepCorrLight_/2, startCorrLight_ + stepCorrLight_*nBinsCorrLight_ - stepCorrLight_/2, nBinsCorrB_,startCorrB_ - stepCorrB_/2, startCorrB_ + stepCorrB_*nBinsCorrB_ - stepCorrB_/2);
  return histo;
}

float* FullKinFit::EstimateTopMass(TRootEvent* event, float WMass, bool writePNG, int jetCombi, bool correctCombi)
{
  nTotalJetCombis_++;
  
//  cout << "kinFit!!!" << endl;
  
  float* result = new float[3];
  int hadrBJetIndex = MVAvals_.second[2], lightJet1Index = MVAvals_.second[0], lightJet2Index = MVAvals_.second[1];
  TLorentzVector bJet = (*jets_[hadrBJetIndex]);
  TLorentzVector lightJet1 = (*jets_[lightJet1Index]);
	TLorentzVector lightJet2 = (*jets_[lightJet2Index]);
	
	float mTopRaw = (bJet+lightJet1+lightJet2).M();
	float mWRaw = (lightJet1+lightJet2).M();
	
  float lightCorr1 = 1 - resFitLightJets_->EtCorrection(&lightJet1);
  float lightCorr2 = 1 - resFitLightJets_->EtCorrection(&lightJet2);
  float bCorr = 1 - resFitBJets_->EtCorrection(&bJet);
  
  lightJet1 = lightCorr1 * lightJet1;
  lightJet2 = lightCorr2 * lightJet2;
  bJet = bCorr * bJet;
  
  float mTopCorr = (bJet+lightJet1+lightJet2).M();
  float mWCorr = (lightJet1+lightJet2).M();
  
  // prepare the resolutions for the Kinematic Fit
  TMatrixD Ml1(3,3), Ml2(3,3), Mb(3,3);
  Ml1.Zero(); Ml2.Zero(); Mb.Zero();
  Ml1(0,0) = pow(resFitLightJetsL7_->EtResolution(&lightJet1, true), 2);
  Ml1(1,1) = pow(resFitLightJetsL7_->ThetaResolution(&lightJet1, true), 2);
  Ml1(2,2) = pow(resFitLightJetsL7_->PhiResolution(&lightJet1, true), 2);
  Ml2(0,0) = pow(resFitLightJetsL7_->EtResolution(&lightJet2, true), 2);
  Ml2(1,1) = pow(resFitLightJetsL7_->ThetaResolution(&lightJet2, true), 2);
  Ml2(2,2) = pow(resFitLightJetsL7_->PhiResolution(&lightJet2, true), 2);
  Mb(0,0) = pow(resFitBJetsL7_->EtResolution(&bJet, true), 2);
  Mb(1,1) = pow(resFitBJetsL7_->ThetaResolution(&bJet, true), 2);
  Mb(2,2) = pow(resFitBJetsL7_->PhiResolution(&bJet, true), 2);
  
//  cout << "resoLight1 :  " << lightJet1.Pt() << " " << resFitLightJetsL7_->EtResolution(&lightJet1, true) << "  | resoLight2 :  " << lightJet2.Pt() << " " << resFitLightJetsL7_->EtResolution(&lightJet2, true) << "  | resoB :  " << bJet.Pt() << " " << resFitBJetsL7_->EtResolution(&bJet, true) << endl;
  
  // mass to start from
  float wJetsSF = WMass/( (lightJet1+lightJet2).M() ); // factor to scale up/down the light jets
  float startTopMass = ( ( (lightJet1+lightJet2)*wJetsSF ) + bJet ).M();
  histo1D_["wJetsSF_"+dataset_->Name()]->Fill(wJetsSF);
  histo1D_["startTopMass_"+dataset_->Name()]->Fill(startTopMass);
  
  if(correctCombi)
  {
    histo1D_["mWUnCorr_"+dataset_->Name()]->Fill(mWRaw);
    histo1D_["mTopUnCorr_"+dataset_->Name()]->Fill(mTopRaw);
    histo1D_["mWCorr_"+dataset_->Name()]->Fill(mWCorr);
    histo1D_["mTopCorr_"+dataset_->Name()]->Fill(mTopCorr);
  }
  
  Float_t mTopConstraint[51], chi2[51];
  float massMinChi2 = -1., minChi2 = 9999999., maxChi2 = 0;
  int nGoodFits = 0;
  for(int i=0; i<51; i++)
  {
    float shiftedTopMass = startTopMass - 50 + i * stepTopMass_;
    if( shiftedTopMass < 0 ) continue;
    
    TKinFitter *theFitter = new TKinFitter("hadtopFit", "hadtopFit");
    theFitter->setVerbosity(0);
    
    TFitParticleEtThetaPhiEMomFix *fitLight1 = new TFitParticleEtThetaPhiEMomFix("lightJet1", "lightJet1", &lightJet1, &Ml1);
    TFitParticleEtThetaPhiEMomFix *fitLight2 = new TFitParticleEtThetaPhiEMomFix("lightJet2", "lightJet2", &lightJet2, &Ml2);
    TFitParticleEtThetaPhiEMomFix *fitB = new TFitParticleEtThetaPhiEMomFix("bJet", "bJet", &bJet, &Mb);
    theFitter->addMeasParticles(fitLight1,fitLight2,fitB);
    
    TFitConstraintM *consW = new TFitConstraintM(WMass, "MassConstraint", 0, 0, WMass );
    TFitConstraintM *consTop = new TFitConstraintM(shiftedTopMass, "MassConstraint", 0, 0, shiftedTopMass );
    consW->addParticles1(fitLight1,fitLight2);
    consTop->addParticles1(fitB,fitLight1,fitLight2);
    
    theFitter->addConstraint(consW);
    theFitter->addConstraint(consTop);
    theFitter->setMaxNbIter(30);
    theFitter->setMaxDeltaS(5e-5);
    theFitter->setMaxF(1e-4);
    
    //do the fit!
    theFitter->fit();
    
    if (theFitter->getStatus() == 0) // if the fitter converged
    {
      chi2[i] = theFitter->getS();
      mTopConstraint[i] = shiftedTopMass;
      nGoodFits++;
      if(chi2[i] < minChi2)
      {
        massMinChi2 = shiftedTopMass;
        minChi2 = chi2[i];
      }
      if(chi2[i] > maxChi2) maxChi2 = chi2[i];
    }
    else
    {
      chi2[i] = 99999;
      mTopConstraint[i] = 99999;
    }
    
    delete theFitter;
    delete fitLight1;
    delete fitLight2;
    delete fitB;
    delete consW;
    delete consTop;
  }
  
  int nGoodFitsInFitWindow = 0;
  for(int i=0; i<51; i++)
    if( mTopConstraint[i] > massMinChi2-15 && mTopConstraint[i] < massMinChi2+15 ) nGoodFitsInFitWindow++;
  
  result[0] = 99999;
  result[1] = 99999;
  result[2] = 99999;
  if( nGoodFitsInFitWindow > 5 )
  {
    nParabolaFitJetCombis_++;
    
    // now extract the fitted top mass
    TGraph* graph = new TGraph(51,mTopConstraint,chi2);
    TF1* parabola = new TF1("parabola","pol2",massMinChi2-15,massMinChi2+15);
//    TF1* parabola = new TF1("parabola","[2]*(x-[1])^2+[0]",massMinChi2-15,massMinChi2+15);
    parabola->SetParLimits(2,0,9999);
    parabola->SetParameter(2,100);
    parabola->SetParName(2,"a");
  //  parabola->SetParLimits(1,-9999,0);
  //  parabola->SetParameter(1,-10);
    parabola->SetParName(1,"b");
  //  parabola->SetParLimits(0,0,9999);
  //  parabola->SetParameter(0,100);
    parabola->SetParName(0,"c");
    parabola->SetLineWidth(2);
    parabola->SetLineColor(kRed);
    graph->Fit("parabola","QR");
    
    float a = parabola->GetParameter(2);
    float b = parabola->GetParameter(1);
    float c = parabola->GetParameter(0);

    float minim = - b/(2*a); // best top mass
    float minimchi2 = a*minim*minim + b*minim +c;
    float minplus = (-b + sqrt(b*b-4*a*(c-minimchi2-1)))/(2*a); //upper bound right
    float minmin  = (-b - sqrt(b*b-4*a*(c-minimchi2-1)))/(2*a); //lower bound left
    float uncert1 = minplus-minim;
    float uncert2 = minim-minmin;
    float uncert = uncert2;
    if(uncert1 > uncert2) uncert = uncert1;
    result[0] = minim;
    result[1] = uncert;
    
    if( minimchi2 < 5 )
      histo1D_["massUnc_"+dataset_->Name()]->Fill(uncert);
    
    if( minim > 0 && minim > startTopMass-51 && minim < startTopMass+51 && a > 0 ) //&& minChi2 < 40 )
    {
      nFinalKinFitJetCombis_++;
      
      if( a < 0 || !(uncert<2000) || minim <= 50 )
      {
        cout << event->runId() << " " << event->lumiBlockId() << " " << event->eventId() << " " << jetCombi << " |  wJetsSF: " << wJetsSF
          << "  startTopMass: " << startTopMass << "  massMinChi2: " << massMinChi2 <<  "  finalTopMass: " << minim << " " << nGoodFits << " " << minChi2 << " " << nGoodFitsInFitWindow << endl;
        writePNG = true;
      }
      
      // now redo kinFit with this top mass as constraint
      TKinFitter *theFitter = new TKinFitter("hadtopFit", "hadtopFit");
      theFitter->setVerbosity(0);
      
      TFitParticleEtThetaPhiEMomFix *fitLight1 = new TFitParticleEtThetaPhiEMomFix("lightJet1", "lightJet1", &lightJet1, &Ml1);
      TFitParticleEtThetaPhiEMomFix *fitLight2 = new TFitParticleEtThetaPhiEMomFix("lightJet2", "lightJet2", &lightJet2, &Ml2);
      TFitParticleEtThetaPhiEMomFix *fitB = new TFitParticleEtThetaPhiEMomFix("bJet", "bJet", &bJet, &Mb);
      theFitter->addMeasParticles(fitLight1,fitLight2,fitB);
      
      TFitConstraintM *consW = new TFitConstraintM(WMass, "MassConstraint", 0, 0, WMass );
      TFitConstraintM *consTop = new TFitConstraintM(result[0], "MassConstraint", 0, 0, result[0] );
      consW->addParticles1(fitLight1,fitLight2);
      consTop->addParticles1(fitB,fitLight1,fitLight2);
      
      theFitter->addConstraint(consW);
      theFitter->addConstraint(consTop);
      theFitter->setMaxNbIter(30);
      theFitter->setMaxDeltaS(5e-5);
      theFitter->setMaxF(1e-4);
      theFitter->fit();      //do the fit!
      
      if (theFitter->getStatus() == 0) // if the fitter converged
      {
        result[2] = theFitter->getS();
        nFinalJetCombis_++;
      }
      else
        result[2] = 99999;
      delete theFitter;
      delete fitLight1;
      delete fitLight2;
      delete fitB;
      delete consW;
      delete consTop;
    }
    if( writePNG )
    {
      stringstream s1; s1 << event->runId();
      stringstream s2; s2 << event->lumiBlockId();
      stringstream s3; s3 << event->eventId();
      stringstream s4; s4 << jetCombi;
      string monsterName = "Monster_" + dataset_->Name() + "_" + s1.str() + "_" + s2.str() + "_" + s3.str() + "_" + s4.str();
      string path = "PlotsJES/Monsters/";
      graph->GetXaxis()->SetLimits(startTopMass-55,startTopMass+55);
      graph->GetYaxis()->SetRangeUser(minChi2-2,maxChi2+5);
      TCanvas* tempCanvas = TCanvasCreator(graph, monsterName);
      tempCanvas->SaveAs( (path + monsterName + ".root").c_str() );
      tempCanvas->SaveAs( (path + monsterName + ".png").c_str() );
      delete tempCanvas;
    }
    delete parabola;
    delete graph;
  }
  return result;
}

void FullKinFit::Write(TFile* fout, bool savePNG, string pathPNG)
{
  cout << "FullKinFit: " << dataset_->Name() << endl;
  cout << "nTotalJetCombis: " << nTotalJetCombis_ << " nParabolaFitJetCombis: " << nParabolaFitJetCombis_ << " nFinalKinFitJetCombis: " << nFinalKinFitJetCombis_ << " nFinalJetCombis: " << nFinalJetCombis_ << endl;
  
  string mWCorrtitle = "mWCorr_"+dataset_->Name();
  string mTopCorrtitle = "mTopCorr_"+dataset_->Name();
  
  if(histo1D_[mWCorrtitle]->GetEntries() > 5)
  {
    TF1* fit_func_W = new TF1("W mass","gaus");
    double left = histo1D_[mWCorrtitle]->GetMean() - histo1D_[mWCorrtitle]->GetRMS();
    double right = histo1D_[mWCorrtitle]->GetMean() + histo1D_[mWCorrtitle]->GetRMS();
    fit_func_W->SetRange(left,right);
    histo1D_[mWCorrtitle]->Fit(fit_func_W,"RQ");

    TF1* fit_func_Top = new TF1("Top mass","gaus");
    left = histo1D_[mTopCorrtitle]->GetMean() - histo1D_[mTopCorrtitle]->GetRMS();
    right = histo1D_[mTopCorrtitle]->GetMean() + histo1D_[mTopCorrtitle]->GetRMS();
    fit_func_Top->SetRange(left,right);
    histo1D_[mTopCorrtitle]->Fit(fit_func_Top,"RQ");
  }
  
  mkdir(pathPNG.c_str(),0777);
  string newPathPNG = pathPNG + dataset_->Name() + "/";
  mkdir(newPathPNG.c_str(),0777);
  
  fout->cd();
  string dirname = "FullKinFit_" + dataset_->Name();
	if(fout->Get(dirname.c_str())==0)
		fout->mkdir(dirname.c_str());
	fout->cd(dirname.c_str());
	
	for(std::map<std::string,TH1F*>::const_iterator it = histo1D_.begin(); it != histo1D_.end(); it++)
	{
		TH1F *temp = it->second;
		temp->Write();
		TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
		tempCanvas->SaveAs( (newPathPNG+it->first+".png").c_str() );
	}
}
	 
