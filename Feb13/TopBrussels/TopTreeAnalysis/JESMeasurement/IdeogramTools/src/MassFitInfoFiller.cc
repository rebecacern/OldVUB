#include "TopTreeAnalysis/JESMeasurement/IdeogramTools/interface/MassFitInfoFiller.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

void fillMassFitInfoHeader(Dataset* dataSet, ideogram::MassFitInfoHeader& header, bool isSemiMu, bool containsPDFsyst)
{
  header.channel_ = isSemiMu;
  string dataSetName = dataSet->Name();
  if( dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 || dataSetName.find("DAtaDriven") != string::npos )
    header.isMC_ = false;
  else header.isMC_ = true;
  header.datasetName_ = dataSetName;
  header.datasetTitle_ = dataSet->Title();
  header.datasetCrossSection_ = dataSet->Xsection();
  header.datasetIntegratedLuminosity_ = dataSet->EquivalentLumi();
  header.pdf_ = containsPDFsyst;
}

void fillMassFitInfoEvent(LightMonster* monster, const ideogram::MassFitInfoHeader& header, ideogram::MassFitInfoEvent& event, float PUweight, float PUup, float PUdown, bool onlyHighestWeight, int bestCombi)
{
  event.runNumber = monster->runID();
  event.luminositySection = monster->lumiBlockID();
  event.eventNumber = monster->eventID();
  vector<TLorentzVector> selJets = monster->selectedJets();
  vector<float> bTag = monster->bTagCSV();
  event.nJet = selJets.size();
  if(header.channel_ == false)
  {
    event.nElectron = 1;
    event.nMuon = 0;
  }
  else
  {
    event.nElectron = 0;
    event.nMuon = 1;
  }
  event.nVertex = monster->nPV();
  event.nPU = monster->nPU();
  event.leptonCharge = monster->leptonCharge();
  TLorentzVector lepton = monster->lepton();
  event.leptonPt = lepton.Pt();
  event.leptonEta = lepton.Eta();
  event.leptonPhi = lepton.Phi();
  TLorentzVector met =  monster->met();
  event.METx = met.Px();
  event.METy = met.Py();
  double Ht = 0;
  for (size_t j = 0 ; j != selJets.size() ; ++j)
  {
    Ht += selJets[j].Pt();
    if(j < ideogram::MassFitInfoEvent::kNJET)
    {
      event.jetPt[j] = selJets[j].Pt();
      event.jetEta[j] = selJets[j].Eta();
      event.jetPhi[j] = selJets[j].Phi();
      event.jetBTagDiscriminant[j] = bTag[j];
    }
  }
  event.HT = Ht;
  event.eventWeight = monster->eventWeight();
  event.lumiWeight = PUweight;
  event.PUdown = PUdown;
  event.PUup = PUup;
  // calculate pz neutrino (for mttbar)
  double M_W  = 80.4;
  double M_mu =  0.10566;
  double pznu = 0.;
  
  double a = M_W*M_W - M_mu*M_mu + 2.0*lepton.Px()*met.Px() + 2.0*lepton.Py()*met.Py();
  double A = 4.0*(pow(lepton.E(),2)- pow(lepton.Pz(),2));
  double B = -4.0*a*lepton.Pz();
  double C = 4.0*pow(lepton.E(),2)*(pow(met.Px(),2) + pow(met.Py(),2)) - a*a;
  double tmproot = B*B - 4.0*A*C;
  
  if(tmproot < 0) pznu = - B/(2*A); // take real part for complex roots
  else
  {
    double tmpsol1 = (-B + TMath::Sqrt(tmproot))/(2.0*A);
    double tmpsol2 = (-B - TMath::Sqrt(tmproot))/(2.0*A);
    pznu = TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ? tmpsol1 : tmpsol2;
  }
  met.SetPxPyPzE(met.Px(), met.Py(), pznu, sqrt(met.Px()*met.Px() + met.Py()*met.Py() + pznu*pznu ) );
  event.Mttbar = (selJets[0]+selJets[1]+selJets[2]+selJets[3]+lepton+met).M();
  event.PtTTbar = (selJets[0]+selJets[1]+selJets[2]+selJets[3]+lepton+met).Pt();
  
  fillMassFitResult(monster, event, onlyHighestWeight, bestCombi);
}

void fillMassFitResult(LightMonster* monster, ideogram::MassFitInfoEvent& event, bool onlyHighestWeight, int bestCombi)
{
  // remove what is already stored in the event
  event.fitInfo.clear();
  event.fitSortedChi2Index.clear();
  size_t nFit = 0;
  for(int iCombi=0; iCombi<12; iCombi++)
  {
    if(!onlyHighestWeight || iCombi==bestCombi)
    {
      unsigned int* combi = monster->mvaResult(iCombi);
      vector<Int_t> jetType;
      for (size_t j = 0 ; j!= ideogram::MassFitInfoEvent::kNJET ; ++j) jetType.push_back(combi[j]);
      
      event.fitInfo.push_back( ideogram::MassFitInfo(jetType, monster->chi2MTopFit(iCombi), exp(-0.5*monster->chi2MTopFit(iCombi)), 1., true, monster->mTopFit(iCombi), monster->sigmaMTopFit(iCombi)) );
      nFit++;
    }
  }
  event.nFit = nFit;
  event.nFitConverge = nFit; // Implement ad-hoc calculation of fit convergence criteria here
  ideogram::normalizeFitWeights(event);
  
  // fill the sorted chi-squared index
  std::vector< std::pair<size_t, double> > Chi2Index;
  for (size_t f = 0 ; f != nFit ; ++f)
    Chi2Index.push_back(std::pair<size_t,double>(f,event.fitInfo[f].fitChi2));
  std::sort(Chi2Index.begin(),Chi2Index.end(),ideogram::IndexedQuantityAbsLessThan<double>);
  for (size_t f = 0 ; f != nFit ; ++f)
    event.fitSortedChi2Index.push_back(Chi2Index[f].first);
}

int highestWeightCombi(LightMonster* monster, double chi2Cut, double bTagCut, double bTagEff, double misTagRate)
{
  double maxWeight = -1;
  int maxWeightCombi = -1;
  vector<float> bTag = monster->bTagCSV();
  for(unsigned int iCombi=0; iCombi<12; iCombi++)
  {
    if(monster->chi2MTopFit(iCombi) >= chi2Cut) continue;
    double weight = exp(-0.5*monster->chi2MTopFit(iCombi));
    unsigned int* combi = monster->mvaResult(iCombi); // 0 = light1, 1 = light2, 2 = hadrB, 3 = leptB
    
    for (UInt_t j = 0 ; j != ideogram::MassFitInfoEvent::kNJET ; ++j)
    {
      if(bTag[j] > bTagCut)
      {
        if( combi[0] == j || combi[1] == j ) weight *= misTagRate;
        else if( combi[2] == j || combi[3] ==  j ) weight *= bTagEff;
      }
      else
      {
        if( combi[0] == j || combi[1] == j ) weight *= (1.-misTagRate);
        else if( combi[2] == j || combi[3] ==  j ) weight *= (1.-bTagEff);
      }
    }
    if(weight > maxWeight)
    {
      maxWeight = weight;
      maxWeightCombi = iCombi;
    }
  }
  return maxWeightCombi;
}

