#include "freyaLooper.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "inputs.h"
//#include <iostream>
#include "Riostream.h"
#include <vector>
#include <string>
#include "TFile.h"
#include "TChain.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"

void freyaLooper::Loop(){
  // running default loop:
  myLoop(0,0,0);
}
void freyaLooper::myLoop(int nsel, int mode, bool silent)
{



  
  double SFval = 0.95;
  bool SFplus = false;
  bool SFminus = false;
  // samples used
  //  double x_sec = 0.;
  char plotName[300];
  sprintf(plotName,"test");
  
  if (nsel == 0)                	{sprintf(plotName,"tt");}
  else if (nsel == 1)   		{sprintf(plotName,"twdr");}
  else if (nsel == -1)   		{sprintf(plotName,"twds");}
  else if (nsel == 2)			{sprintf(plotName, "st");}
  else if (nsel == 3)   		{sprintf(plotName,"wjets");}
  else if (nsel == 4)   		{sprintf(plotName,"zjets");}
  else if (nsel == 5)   		{sprintf(plotName,"dy");}
  else if (nsel == 6)   		{sprintf(plotName,"qcd_mu");}
  else if (nsel == 7)   		{sprintf(plotName,"di");}
  else if (nsel == 77)                	{sprintf(plotName,"non_tt");}
  else if (nsel == 66)                	{sprintf(plotName,"others");}
  else if (nsel == 666)                	{sprintf(plotName,"data");}
  else if (nsel == 555)                	{sprintf(plotName,"mc");}
  else if (nsel == 777)                 {sprintf(plotName,"all");}
  
  // if (mode != 0 &&  mode !=1 && mode !=2) mode = 0;
  if (!silent){
    if      (mode == 0) 	std::cout << " Electron-Muon Mixed channel " << std::endl;
    else if (mode == 1) 	std::cout << " Di-Muon channel " << std::endl;
    else if (mode == 2) 	std::cout << " Di-Electron channel " << std::endl;
    std::cout << "*************************************************" << std::endl;
  }
  
  bool nosf = true;
  
  if(!silent){
    std::cout << std::endl;
    std::cout << "******************************************" << std::endl;
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "Starting the analysis: " << plotName <<  std::endl;
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "Number of Raw events: " <<  fChain->GetEntries() << std::endl;
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "******************************************" << std::endl;
  }


  
  char newRootFile[300];
  double lumi = luminosity;
  if (mode == 0 )        lumi = 4626.297;
  else if ( mode == 1)   lumi = 4534.871;
  else if ( mode == 2)   lumi = 4593.348;
  else if ( mode == 3)   lumi = 4.5;
  sprintf(newRootFile,"results/an_%dpb_%d.root", (int)lumi, mode);
  
  TFile f_var(newRootFile, "UPDATE");
  std::cout << "******************************************" << std::endl;
  std::cout << "creating new root file " << newRootFile << std::endl;
  std::cout << "******************************************" << std::endl;
  
  //////////
  char title[300];
  sprintf(title,"cuts_%s",plotName);
  TH1F* histo = new TH1F( title, " ", 11,  -0.5, 10.5);
  histo->Sumw2();
  
  sprintf(title,"met_%s",plotName);
  TH1F* histo_met = new TH1F( title, " ", 100,  0, 200 );
  histo_met->Sumw2();
  
  sprintf(title,"met_high_%s",plotName);
  TH1F* histo_met_high = new TH1F( title, " ", 100,  0, 200 );
  histo_met_high->Sumw2();
  
  sprintf(title,"met_low_%s",plotName);
  TH1F* histo_met_low = new TH1F( title, " ", 100,  0, 200 );
  histo_met_low->Sumw2();
  
  sprintf(title,"promet_%s",plotName);
  TH1F* histo_promet = new TH1F( title, " ", 100,  0, 200 );
  histo_promet->Sumw2();
  
  sprintf(title,"met_cut_%s",plotName);
  TH1F* histo_met_cut = new TH1F( title, " ", 100,  0, 200 );
  histo_met_cut->Sumw2();
  
  sprintf(title,"met_bt_%s",plotName);
  TH1F* histo_met_bt = new TH1F( title, " ", 100,  0, 200 );
  histo_met_bt->Sumw2();
  
  sprintf(title,"mll_%s",plotName);
  TH1F* histo_mll = new TH1F( title, " ", 100,  0, 200 );
  histo_mll->Sumw2();
  
  sprintf(title,"mll_after_%s",plotName);
  TH1F* histo_mll_after = new TH1F( title, " ", 100,  0, 200 );
  histo_mll_after->Sumw2();
  
  sprintf(title,"njets_%s",plotName);
  TH1F* histo_njets = new TH1F( title, " ", 10,  0, 10 );
  histo_njets->Sumw2();
  
  sprintf(title,"njets_cut_%s",plotName);
  TH1F* histo_njets_cut = new TH1F( title, " ", 10,  -0.5, 9.5 );
  histo_njets_cut->Sumw2();
  
  sprintf(title,"njetsbt_%s",plotName);
  TH1F* histo_njetsbt = new TH1F( title, " ", 10,  -0.5, 9.5 );
  histo_njetsbt->Sumw2();
  
  sprintf(title,"njetsbt_cut_%s",plotName);
  TH1F* histo_njetsbt_cut = new TH1F( title, " ", 10,   -0.5, 9.5 );
  histo_njetsbt_cut->Sumw2();
  
  sprintf(title,"njetsbt_high_%s",plotName);
  TH1F* histo_njetsbt_high = new TH1F( title, " ", 10,   -0.5, 9.5 );
  histo_njetsbt_high->Sumw2();

  sprintf(title,"njetsbt_low_%s",plotName);
  TH1F* histo_njetsbt_low = new TH1F( title, " ", 10,   -0.5, 9.5 );
  histo_njetsbt_low->Sumw2();
  
  sprintf(title,"njets_high_%s",plotName);
  TH1F* histo_njets_high = new TH1F( title, " ", 10,   -0.5, 9.5 );
  histo_njets_high->Sumw2();

  sprintf(title,"njets_low_%s",plotName);
  TH1F* histo_njets_low = new TH1F( title, " ", 10,   -0.5, 9.5 );
  histo_njets_low->Sumw2();
  
  sprintf(title,"ptsys_%s",plotName);
  TH1F* histo_ptsys = new TH1F( title, " ", 100,  0, 200 );
  histo_ptsys->Sumw2();
  
  sprintf(title,"ptsys_high_%s",plotName);
  TH1F* histo_ptsys_high = new TH1F( title, " ", 100,  0, 200 );
  histo_ptsys_high->Sumw2();
  
  sprintf(title,"ptsys_low_%s",plotName);
  TH1F* histo_ptsys_low = new TH1F( title, " ", 100,  0, 200 );
  histo_ptsys_low->Sumw2();
  
  sprintf(title,"ht_%s",plotName);
  TH1F* histo_ht = new TH1F( title, " ", 300,  0, 600 );
  histo_ht->Sumw2();
  
  sprintf(title,"ht_high_%s",plotName);
  TH1F* histo_ht_high = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_high->Sumw2();
  
  sprintf(title,"ht_low_%s",plotName);
  TH1F* histo_ht_low = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_low->Sumw2();
  
  sprintf(title,"ht_cut_%s",plotName);
  TH1F* histo_ht_cut = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_cut->Sumw2();
  
  sprintf(title,"pt_max_%s",plotName);
  TH1F* histo_pt_max = new TH1F( title, " ", 100,  0, 200 );
  histo_pt_max->Sumw2();
  
  sprintf(title,"pt_min_%s",plotName);
  TH1F* histo_pt_min = new TH1F( title, " ", 100,  0, 200 );
  histo_pt_min->Sumw2();
  
  sprintf(title,"pt_leading_%s",plotName);
  TH1F* histo_pt_leading = new TH1F( title, " ", 100,  0, 200 );
  histo_pt_leading->Sumw2();
  
  sprintf(title,"btagHE_%s",plotName);
  TH2F* histo_btagHE = new TH2F( title, " ", 300,  -200, 100, 100, -2, 7);
  histo_btagHE->Sumw2();
  
  sprintf(title,"btagHP_%s",plotName);
  TH2F* histo_btagHP = new TH2F( title, " ", 300,  -200, 100, 100, -2, 7);
  histo_btagHP->Sumw2();
  
  sprintf(title,"etalepton_%s",plotName);
  TH1F* histo_etalepton = new TH1F( title, " ", 101,  -3, 3);
  histo_etalepton->Sumw2();
  
  sprintf(title,"ht_bf_%s",plotName);
  TH1F* histo_ht_bf = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_bf->Sumw2();
  
  sprintf(title,"ptsys_bf_%s",plotName);
  TH1F* histo_ptsys_bf = new TH1F( title, " ", 100,  0, 200 );
  histo_ptsys_bf->Sumw2();
  
  sprintf(title,"met_sqrtHT_%s",plotName);
  TH1F* histo_met_sqrtHT = new TH1F( title, " ", 150,  -10, 20 );
  histo_met_sqrtHT->Sumw2();
  
  
  sprintf(title,"met_sqrtHT_control_%s",plotName);
  TH1F* histo_met_sqrtHT_control = new TH1F( title, " ", 150,  -10, 20 );
  histo_met_sqrtHT_control->Sumw2();
  
  
  sprintf(title,"nvertex_%s",plotName);
  TH1F* histo_nvertex = new TH1F( title, " ", 30,   -0.5, 29.5 );
  histo_nvertex->Sumw2();\


  sprintf(title,"nvertex_control_%s",plotName);
  TH1F* histo_nvertex_control = new TH1F( title, " ", 30,   -0.5, 29.5 );
  histo_nvertex_control->Sumw2();
  
  sprintf(title,"npu_%s",plotName);
  TH1F* histo_npu = new TH1F( title, " ", 30,   -0.5, 29.5 );
  histo_npu->Sumw2();


  sprintf(title,"npixhitsmissed_%s",plotName);
  TH1F* histo_npixhitsmissed = new TH1F( title, " ", 5,   -0.5, 4.5 );
  histo_npixhitsmissed->Sumw2();


  sprintf(title,"mht_%s",plotName);
  TH1F* histo_mht = new TH1F( title, " ", 100,  0 , 400 );
  histo_mht->Sumw2();


  sprintf(title,"mht_met_%s",plotName);
  TH1F* histo_mht_met = new TH1F( title, " ", 100, -200 , 200 );
  histo_mht_met->Sumw2();

  sprintf(title,"mht_div_met_%s",plotName);
  TH1F* histo_mht_div_met = new TH1F( title, " ", 100, 0 , 20 );
  histo_mht_div_met->Sumw2();

  sprintf(title,"sigmaIetaIeta_%s",plotName);
  TH1F* histo_sigmaIetaIeta = new TH1F( title, " ", 100,   0, 0.05 );
  histo_sigmaIetaIeta->Sumw2();

  sprintf(title,"flavLepton_%s",plotName);
  TH1F* histo_flavLepton = new TH1F( title, " ", 30,   -15, 15 );
  histo_flavLepton->Sumw2();
  
  sprintf(title,"ecalRelIso_%s",plotName);
  TH1F* histo_ecalRelIso = new TH1F( title, " ", 300,   0, 0.15 );
  histo_ecalRelIso->Sumw2();


  sprintf(title,"hcalRelIso_%s",plotName);
  TH1F* histo_hcalRelIso = new TH1F( title, " ", 300,   0, 0.15 );
  histo_hcalRelIso->Sumw2();
  sprintf(title,"totalRelIso_%s",plotName);
  TH1F* histo_totalRelIso = new TH1F( title, " ", 300,   0, 0.15 );
  histo_totalRelIso->Sumw2();

  sprintf(title,"trkRelIso_%s",plotName);
  TH1F* histo_trkRelIso = new TH1F( title, " ", 300,   0, 0.3 );
  histo_trkRelIso->Sumw2();

  sprintf(title,"d0_%s",plotName);
  TH1F* histo_d0 = new TH1F( title, " ", 100,   -0.02, 0.02 );
  histo_d0->Sumw2();
  
  sprintf(title,"etaEle_%s",plotName);
  TH1F* histo_etaEle = new TH1F( title, " ", 30,   -3, 3 );
  histo_etaEle->Sumw2();
  
  sprintf(title,"phiEle_%s",plotName);
  TH1F* histo_phiEle = new TH1F( title, " ", 100,   -TMath::Pi(), TMath::Pi() );
  histo_phiEle->Sumw2();

  sprintf(title,"eOverpEle_%s",plotName);
  TH1F* histo_eOverpEle = new TH1F( title, " ", 100,   0, 5 );
  histo_eOverpEle->Sumw2();



  sprintf(title,"deltaEtaSCEle_%s",plotName);
  TH1F* histo_deltaEtaSCEle = new TH1F( title, " ", 100,   -0.25, 0.25);
  histo_deltaEtaSCEle->Sumw2();

  sprintf(title,"deltaPhiSCEle_%s",plotName);
  TH1F* histo_deltaPhiSCEle = new TH1F( title, " ", 100,   -0.25, 0.25 );
  histo_deltaPhiSCEle->Sumw2();

  sprintf(title,"convDcotEle_%s",plotName);
  TH1F* histo_convDcotEle = new TH1F( title, " ", 100,   -0.25, 0.25 );
  histo_convDcotEle->Sumw2();
  
  sprintf(title,"convDistEle_%s",plotName);
  TH1F* histo_convDistEle = new TH1F( title, " ", 100,   -0.25, 0.25 );
  histo_convDistEle->Sumw2();

  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //    if(jentry%100000==0 || (jentry<100000 && jentry%1000==0))
    //      std::cout << jentry << "/" << nentries << std::endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;


    if (lumi != lum && nsel != 666 && mode !=3){
      if (jentry == 0)std::cout << "Warning: This tree was made with a different luminosity (" << lum << ") than " << lumi << std::endl;
      xlWeight*=(lumi/lum);
    }

    if(ptLepton->size() != 2){
      std::cout << "Something is wrong, your Tree is not correctly filled" << std::endl;
      break;

    }
    

    histo->Fill(0.,xlWeight);
    // reject bad electrons on-the-fly:
    bool foundconv=false;
    for(int ilep=0; ilep<ptLepton->size() && !foundconv; ilep++){
      if(abs(flavLepton->at(ilep))!=11)
	continue;
      
      if(fabs(convDcotLepton->at(ilep))<0.02 && fabs(convDistLepton->at(ilep))<0.20){
	//	std::cout << convDcotLepton->at(ilep) << " " << convDistLepton->at(ilep) << std::endl;
	foundconv=true;
      }
      if(numPixHitsLepton->at(ilep)!=0){
	foundconv=true;
	//	std::cout << numPixHitsLepton->at(ilep) << std::endl;
      }
    }
    if(foundconv){
      //      std::cout << "rejecting electron!" << std::endl;
      continue;
    }


    histo->Fill(1,xlWeight);
    int ntightele=0;
    for(int ilep=0; ilep<ptLepton->size() && !foundconv; ilep++){
      if(abs(flavLepton->at(ilep))!=11){
	ntightele++;
	continue;
      }

      if(hcalRelIsoLepton->at(ilep)>0.03)
	continue;
      if(ecalRelIsoLepton->at(ilep)>0.04)
	continue;
      if(trkRelIsoLepton->at(ilep)>0.04)
	continue;
      ntightele++;
    }    
    if(ntightele<2)
      continue;
    
    histo->Fill(2.,xlWeight);

    Double_t mht_alljets_x=0;
    Double_t mht_alljets_y=0;
    Double_t mht_alljets=0;
    Double_t ht_alljets=0;
    
    TLorentzVector lepton0(pxLepton->at(0),pyLepton->at(0), pzLepton->at(0), eLepton->at(0));
    TLorentzVector lepton1(pxLepton->at(1),pyLepton->at(1), pzLepton->at(1), eLepton->at(1)); 
    TLorentzVector pair = lepton0+lepton1;
    
    for(size_t ijet=0; ijet<pxJet->size(); ijet++){
      if(ptJet->at(ijet)<30.)
	continue;
      mht_alljets_x-=pxJet->at(ijet);
      ht_alljets+=ptJet->at(ijet);
      mht_alljets_y-=pyJet->at(ijet);
    }
    for(size_t iilep=0; iilep<2; iilep++){
      mht_alljets_x-=pxLepton->at(iilep);
      mht_alljets_y-=pyLepton->at(iilep);
      mht_alljets+=ptLepton->at(iilep);
    }
    
    mht_alljets = sqrt(mht_alljets_x*mht_alljets_x+mht_alljets_y*mht_alljets_x);

    double phipairmet_t = 0;
    double pi_m = 3.1416/2;
    phipairmet_t = pi_m;
    
    TVector3 vmet(metPx, metPy, 0);
    
    double promet = metPt*sin(phipairmet_t);
    TVector3 m0(pxLepton->at(0),pyLepton->at(0), pzLepton->at(0));
    TVector3 m1(pxLepton->at(1),pyLepton->at(1), pzLepton->at(1)); 
    if (fabs(m0.DeltaPhi(vmet)) < phipairmet_t) phipairmet_t = fabs(m0.DeltaPhi(vmet));
    if (fabs(m1.DeltaPhi(vmet)) < phipairmet_t) phipairmet_t = fabs(m1.DeltaPhi(vmet));
    
    if (phipairmet_t == pi_m) promet = metPt;
    
    if (pair.M() <= 20)
      continue;


    histo->Fill(3, xlWeight);

    int nJetsBT = 0;
    int nJets = 0;
    bool bTagged = false;
    int iJet = -5;
    int iSF;
    double tempSF = SFval;
    if (SFminus)  tempSF = SFval - SFval*10/100;
    if (SFplus)   tempSF = SFval + SFval*10/100;
    int SFvalue = tempSF*100 +1;
    //
      

    if ( nsel == 666 || !nosf){
      for (size_t i =0; i < ptJet->size(); i ++){ 
	TLorentzVector tempJet(pxJet->at(i),pyJet->at(i), pzJet->at(i), eJet->at(i));
	if (ptJet->at(i) >= 30 && TMath::Min(fabs(lepton0.DeltaR(tempJet)), fabs(lepton1.DeltaR(tempJet))) > 0.3) {
	  nJets++;
	  if (iJet == -5) iJet = i;
	  if (btSSVHEJet->at(i) > 1.74){
	    bTagged = true;
	    nJetsBT++;
	  } 
	} else if (btSSVHEJet->at(i) > 1.74) nJetsBT++;
      }
    } 

    else {
      //// Regular SF
      if (SFvalue < 101){
	for (size_t i =0; i < ptJet->size(); i ++){ 
	  TLorentzVector tempJet(pxJet->at(i),pyJet->at(i), pzJet->at(i), eJet->at(i));
	  if (ptJet->at(i) >= 30 && TMath::Min(fabs(lepton0.DeltaR(tempJet)), fabs(lepton1.DeltaR(tempJet))) > 0.3) {
	    nJets++;
	    if (iJet == -5) iJet = i;
	    if (btSSVHEJet->at(i) > 1.74){
	      iSF = rand() % 101;
	      if (iSF < SFvalue ){
		bTagged = true;
		nJetsBT++;
	      } 
	    } 
	  } else if (btSSVHEJet->at(i) > 1.74){
	    iSF = rand() % 101;
	    if (iSF < SFvalue ) nJetsBT++;
	  }
	}
      }
      else {
	//// Large SF
	     
	for (size_t i =0; i < ptJet->size(); i ++){ 
	  TLorentzVector tempJet(pxJet->at(i),pyJet->at(i), pzJet->at(i), eJet->at(i));
	  if (ptJet->at(i) >= 30 && TMath::Min(fabs(lepton0.DeltaR(tempJet)), fabs(lepton1.DeltaR(tempJet))) > 0.3) {
	    nJets++;
	    if (iJet == -5) iJet = i;
	    if (btSSVHEJet->at(i) > 1.74 ){
	      bTagged = true;
	      nJetsBT++;
	    } 
	    else {
	      iSF = rand() % 101;
	      if (iSF < abs(100 - SFvalue)){
		nJetsBT++;
		bTagged = true;
	      }
	    }
	  }
	  else if (btSSVHEJet->at(i) > 1.74){
	    nJetsBT++;
	  } 
	  else {
	    iSF = rand() % 101;
	    if (iSF < abs(100 - SFvalue)) nJetsBT++;
	  }
	}
      }
    }
    
	 
    histo_pt_max->Fill(TMath::Max(lepton0.Pt(), lepton1.Pt()), xlWeight);
    histo_pt_min->Fill(TMath::Min(lepton0.Pt(), lepton1.Pt()), xlWeight);
    histo_njets->Fill(nJets,  xlWeight);
    histo_njetsbt->Fill(nJetsBT,  xlWeight);
    histo_mll->Fill(pair.M(),  xlWeight);
    histo_met->Fill(metPt,  xlWeight);
    histo_promet->Fill(promet, xlWeight);
	
    if (nvertex > 5){
      histo_met_high->Fill(metPt,  xlWeight);
      histo_njets_high->Fill(nJets,  xlWeight);
      histo_njetsbt_high->Fill(nJetsBT,  xlWeight);
	  
    } else {
      histo_met_low->Fill(metPt,  xlWeight);
      histo_njets_low->Fill(nJets,  xlWeight);
      histo_njetsbt_low->Fill(nJetsBT,  xlWeight);
	  
    }

    if (nJets) histo_pt_leading->Fill(ptJet->at(0), xlWeight);
	
    if (nJets == 1){
      histo_etalepton->Fill(lepton0.Eta(), xlWeight);
      TLorentzVector jet(pxJet->at(iJet),pyJet->at(iJet), pzJet->at(iJet), eJet->at(iJet));
		
      double ptSysPx1 = lepton0.Px() + lepton1.Px() + jet.Px() + metPx;
      double ptSysPy1 = lepton0.Py() + lepton1.Py() + jet.Py() + metPy;
      double ptSystem1 = sqrt(ptSysPx1*ptSysPx1 + ptSysPy1*ptSysPy1);
      double ht1 = lepton0.Pt() + lepton1.Pt() + jet.Pt() + metPt; 
      histo_ptsys_bf->Fill(ptSystem1, xlWeight);
      histo_ht_bf->Fill(ht1, xlWeight);
    }

    bool invMass = false;
    if      (mode == 0) invMass = true;
    else if (mode == 1  && (pair.M() > invMax || pair.M() < invMin)) invMass = true;
    else if (mode == 2 && (pair.M() > invMax || pair.M() < invMin)) invMass = true;
    
    if (invMass){
      histo->Fill(4, xlWeight);
      histo_mll_after->Fill(pair.M(),  xlWeight);
      histo_met_cut->Fill(metPt,  xlWeight);
      
      if(metPt >= metCut || mode ==0 ){
     //if (promet >= metCut || mode ==0){
	histo->Fill(5, xlWeight);
	histo_njets_cut->Fill(nJets, xlWeight);
	if(nJets==2 && nJetsBT>0){// control region

	  histo_nvertex_control->Fill(nvertex, xlWeight);
	  histo_met_sqrtHT_control->Fill((metPt-mht_alljets)/sqrt(ht_alljets),xlWeight);



		for(int ilep=0; ilep<2; ilep++){
	      
		  histo_flavLepton->Fill((*flavLepton)[ilep],xlWeight);
		  if(abs((*flavLepton)[ilep])!=11){
		    continue;
		  }
		  histo_sigmaIetaIeta->Fill((*sigmaIetaIetaLepton)[ilep],xlWeight);
		  histo_ecalRelIso->Fill((*ecalRelIsoLepton)[ilep],xlWeight);
		  histo_hcalRelIso->Fill((*hcalRelIsoLepton)[ilep],xlWeight);
		  histo_trkRelIso->Fill((*trkRelIsoLepton)[ilep],xlWeight);
		  histo_totalRelIso->Fill((*totalRelIsoLepton)[ilep],xlWeight);
		  histo_d0->Fill((*d0Lepton)[ilep],xlWeight);
		  histo_etaEle->Fill((*etaLepton)[ilep],xlWeight);
		  histo_phiEle->Fill((*phiLepton)[ilep],xlWeight);
		  
		  histo_eOverpEle->Fill((*eOverpLepton)[ilep],xlWeight);
		  histo_npixhitsmissed->Fill((*numPixHitsLepton)[ilep],xlWeight);
		  
		  histo_deltaPhiSCEle->Fill((*deltaPhiSCTrkAtVtxLepton)[ilep],xlWeight);
		  histo_deltaEtaSCEle->Fill((*deltaEtaSCTrkAtVtxLepton)[ilep],xlWeight);
		  histo_convDcotEle->Fill((*convDcotLepton)[ilep],xlWeight);
		  histo_convDistEle->Fill((*convDistLepton)[ilep],xlWeight);
		  
		}
	}
	if (nJets == 1){
	  histo->Fill(6, xlWeight);
	  histo_njetsbt_cut->Fill(nJetsBT, xlWeight);
	  

	  if (nJetsBT == 1 && bTagged){
	    histo->Fill(7, xlWeight);

	    TLorentzVector jet(pxJet->at(iJet),pyJet->at(iJet), pzJet->at(iJet), eJet->at(iJet));
	    
	    double ptSysPx = lepton0.Px() + lepton1.Px() + jet.Px() + metPx;
	    double ptSysPy = lepton0.Py() + lepton1.Py() + jet.Py() + metPy;
	    double pxHT = lepton0.Px() + lepton1.Px() + jet.Px();
	    double pyHT = lepton0.Py() + lepton1.Py() + jet.Py();
	    double mht = sqrt(pxHT*pxHT+pyHT*pyHT);
	    double ptSystem = sqrt(ptSysPx*ptSysPx + ptSysPy*ptSysPy);
	    double ht = lepton0.Pt() + lepton1.Pt() + jet.Pt() + metPt; 
	    
	    histo_mht->Fill(mht,xlWeight);
	    histo_mht_met->Fill(mht-metPt,xlWeight);
	    histo_mht_div_met->Fill((mht-metPt)/metPt,xlWeight);
	    
	    histo_ptsys->Fill(ptSystem, xlWeight);
	    histo_ht->Fill(ht, xlWeight);
	    histo_btagHE->Fill(btTCHEJet->at(iJet), btSSVHEJet->at(iJet), xlWeight);
	    histo_btagHP->Fill(btTCHPJet->at(iJet), btSSVHPJet->at(iJet), xlWeight);
	    histo_met_bt->Fill(metPt, xlWeight);
	    
	    histo_nvertex->Fill(nvertex, xlWeight);
	    histo_npu->Fill(npu, xlWeight);
	    
	    if (nvertex > 5) {
	      histo_ptsys_high->Fill(ptSystem, xlWeight);
	      histo_ht_high->Fill(ht, xlWeight);
	    } else {
	      histo_ptsys_low->Fill(ptSystem, xlWeight);
	      histo_ht_low->Fill(ht, xlWeight);
	      
	    }
	    
	    if (ptSystem <= ptsysCut){
	      histo->Fill(8, xlWeight);
	      histo_ht_cut->Fill(ht, xlWeight);
	      if (ht > htMin || mode !=0){
		histo->Fill(9, xlWeight);

		histo_met_sqrtHT->Fill(	(metPt-mht_alljets)/sqrt(ht_alljets),xlWeight);	


	      }
	    }
	  }
	}
      }
    } // 2 leptons
  }// event loop.


  
  if (!silent){ 
   std::cout << "------------------------------------------" << std::endl;
   std::cout << "Results: " << plotName <<  std::endl;
   std::cout << "------------------------------------------" << std::endl;  
    for (int i = 1; i < 10; i++){
      std::cout << i << " "  <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << std::endl;

//       if (i == 2)std::cout << " leptons: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << std::endl;
      
//       if (i == 2+1)std::cout << " gamma veto: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << std::endl;
      
//       if (i == 3+1)std::cout << " inv. mass: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << std::endl;
//       if (i == 4+1)std::cout << " met: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << std::endl;
//       if (i == 5+1)std::cout << " jet: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << std::endl;
//       if (i == 6+1)std::cout << " jet_bt: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << std::endl;
//       if (i == 7+1)std::cout << " pt system: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << std::endl;
//       if (i == 8+1)std::cout << " ht: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << std::endl;
    }
   std::cout << "------------------------------------------" << std::endl; 
  }
  f_var.Write();
  f_var.Close();
}
