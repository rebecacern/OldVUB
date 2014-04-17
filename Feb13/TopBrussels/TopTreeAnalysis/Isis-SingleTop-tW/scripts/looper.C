#include "looper.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "inputs.h"
#include "Riostream.h"
#include <vector>
#include <string>
#include "TFile.h"
#include "TChain.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"

void looper::Loop(){
  // running default loop:
  myLoop(0,0,0);
}
void looper::myLoop(int nsel, int mode, bool silent)
{

  char plotName[300];
  sprintf(plotName,"test");
  
  if (nsel == 0)                	{sprintf(plotName,"tt");}
  else if (nsel == 1)   		{sprintf(plotName,"twdr");}
  else if (nsel == -1)   		{sprintf(plotName,"twds");}
  else if (nsel == 2)   		{sprintf(plotName,"zjets");}
  else if (nsel == 3)   		{sprintf(plotName,"di");}
  else if (nsel == 4)			{sprintf(plotName, "st");}
  else if (nsel == 5)   		{sprintf(plotName,"wjets");}
  else if (nsel == 6)   		{sprintf(plotName,"qcd_mu");}
  else if (nsel == 7)                	{sprintf(plotName,"others");}
  
  else if (nsel == 555)                	{sprintf(plotName,"mc");}
  
  else if (nsel == 666)                	{sprintf(plotName,"data");}
  
  bool nosf = false;
  
  
  char newRootFile[300];
  double lumi = luminosity; 
  if (mode == 0 )        lumi = 5085.246; 
  else if ( mode == 1)   lumi = 1000;
  else if ( mode == 2)   lumi = 5103.58;
  sprintf(newRootFile,"results/an_%dpb_%d.root", (int)lumi, mode);
 
  TFile f_var(newRootFile, "UPDATE");
  
  if(!silent){
    std::cout << "[Info:] results root file " << newRootFile << std::endl;
  }
  
  
  //////////
  char title[300];
  sprintf(title,"cuts_%s",plotName);
  TH1F* histo = new TH1F( title, " ", 10,  0, 10 );
  histo->Sumw2();
  
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
  
  sprintf(title,"mll_after_%s",plotName);
  TH1F* histo_mll_after = new TH1F( title, " ", 100,  0, 200 );
  histo_mll_after->Sumw2();
  
  sprintf(title,"njets_cut_%s",plotName);
  TH1F* histo_njets_cut = new TH1F( title, " ", 10,  -0.5, 9.5 );
  histo_njets_cut->Sumw2();
  
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
  
  sprintf(title,"ptsys_high_%s",plotName);
  TH1F* histo_ptsys_high = new TH1F( title, " ", 100,  0, 200 );
  histo_ptsys_high->Sumw2();
  
  sprintf(title,"ptsys_low_%s",plotName);
  TH1F* histo_ptsys_low = new TH1F( title, " ", 100,  0, 200 );
  histo_ptsys_low->Sumw2();
  
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
  
  sprintf(title,"npu_%s",plotName);
  TH1F* histo_npu = new TH1F( title, " ", 30,   -0.5, 29.5 );
  histo_npu->Sumw2();
  
  
  /// Classic plotmaker plots
  sprintf(title,"met_%s",plotName);
  TH1F* histo_met = new TH1F( title, " ", 100,  0, 200 );
  histo_met->Sumw2();
  
  sprintf(title,"mll_%s",plotName);
  TH1F* histo_mll = new TH1F( title, " ", 100,  0, 200 );
  histo_mll->Sumw2();
  
  sprintf(title,"njets_%s",plotName);
  TH1F* histo_njets = new TH1F( title, " ", 10,  0, 10 );
  histo_njets->Sumw2();
  
  sprintf(title,"njetsbt_%s",plotName);
  TH1F* histo_njetsbt = new TH1F( title, " ", 10,  -0.5, 9.5 );
  histo_njetsbt->Sumw2();
  
  sprintf(title,"ptsys_%s",plotName);
  TH1F* histo_ptsys = new TH1F( title, " ", 100,  0, 200 );
  histo_ptsys->Sumw2();
  
  sprintf(title,"ht_%s",plotName);
  TH1F* histo_ht = new TH1F( title, " ", 300,  0, 600 );
  histo_ht->Sumw2();
  
  sprintf(title,"pt_leading_%s",plotName);
  TH1F* histo_pt_leading = new TH1F( title, " ", 100,  0, 200 );
  histo_pt_leading->Sumw2();
  
  sprintf(title,"nvertex_%s",plotName);
  TH1F* histo_nvertex = new TH1F( title, " ", 70,   -0.5, 69.5 );
  histo_nvertex->Sumw2();
 
  
  // 1 jet level
  /// Classic plotmaker plots
  sprintf(title,"met_1j_%s",plotName);
  TH1F* histo_met_1j = new TH1F( title, " ", 100,  0, 200 );
  histo_met_1j->Sumw2();
  
  sprintf(title,"mll_1j_%s",plotName);
  TH1F* histo_mll_1j = new TH1F( title, " ", 100,  0, 200 );
  histo_mll_1j->Sumw2();
  
  sprintf(title,"ptsys_1j_%s",plotName);
  TH1F* histo_ptsys_1j = new TH1F( title, " ", 100,  0, 200 );
  histo_ptsys_1j->Sumw2();
  
  sprintf(title,"ht_1j_%s",plotName);
  TH1F* histo_ht_1j = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_1j->Sumw2();
  
  sprintf(title,"pt_leading_1j_%s",plotName);
  TH1F* histo_pt_leading_1j = new TH1F( title, " ", 100,  0, 200 );
  histo_pt_leading_1j->Sumw2();
  
  sprintf(title,"nvertex_1j_%s",plotName);
  TH1F* histo_nvertex_1j = new TH1F( title, " ", 30,   -0.5, 29.5 );
  histo_nvertex_1j->Sumw2();
  
  // 1 jet level
  /// Classic plotmaker plots
  sprintf(title,"met_1j1t_%s",plotName);
  TH1F* histo_met_1j1t = new TH1F( title, " ", 100,  0, 200 );
  histo_met_1j1t->Sumw2();
  
  sprintf(title,"mll_1j1t_%s",plotName);
  TH1F* histo_mll_1j1t = new TH1F( title, " ", 100,  0, 200 );
  histo_mll_1j1t->Sumw2();
  
  sprintf(title,"ptsys_1j1t_%s",plotName);
  TH1F* histo_ptsys_1j1t = new TH1F( title, " ", 100,  0, 200 );
  histo_ptsys_1j1t->Sumw2();
  
  sprintf(title,"ht_1j1t_%s",plotName);
  TH1F* histo_ht_1j1t = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_1j1t->Sumw2();
  
  sprintf(title,"pt_leading_1j1t_%s",plotName);
  TH1F* histo_pt_leading_1j1t = new TH1F( title, " ", 100,  0, 200 );
  histo_pt_leading_1j1t->Sumw2();
  
  sprintf(title,"nvertex_1j1t_%s",plotName);
  TH1F* histo_nvertex_1j1t = new TH1F( title, " ", 30,   -0.5, 29.5 );
  histo_nvertex_1j1t->Sumw2();
  
  
  // 2j1t
  /// Classic plotmaker plots
  sprintf(title,"met_2j1t_%s",plotName);
  TH1F* histo_met_2j1t = new TH1F( title, " ", 100,  0, 200 );
  histo_met_2j1t->Sumw2();
  
  sprintf(title,"mll_2j1t_%s",plotName);
  TH1F* histo_mll_2j1t = new TH1F( title, " ", 100,  0, 200 );
  histo_mll_2j1t->Sumw2();
  
  sprintf(title,"ptsys_2j1t_%s",plotName);
  TH1F* histo_ptsys_2j1t = new TH1F( title, " ", 100,  0, 200 );
  histo_ptsys_2j1t->Sumw2();
  
  sprintf(title,"ht_2j1t_%s",plotName);
  TH1F* histo_ht_2j1t = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_2j1t->Sumw2();
  
  sprintf(title,"pt_leading_2j1t_%s",plotName);
  TH1F* histo_pt_leading_2j1t = new TH1F( title, " ", 100,  0, 200 );
  histo_pt_leading_2j1t->Sumw2();
  
  sprintf(title,"nvertex_2j1t_%s",plotName);
  TH1F* histo_nvertex_2j1t = new TH1F( title, " ", 30,   -0.5, 29.5 );
  histo_nvertex_2j1t->Sumw2();
  
  
  // 2j2t
  /// Classic plotmaker plots
  sprintf(title,"met_2j2t_%s",plotName);
  TH1F* histo_met_2j2t = new TH1F( title, " ", 100,  0, 200 );
  histo_met_2j2t->Sumw2();
  
  sprintf(title,"mll_2j2t_%s",plotName);
  TH1F* histo_mll_2j2t = new TH1F( title, " ", 100,  0, 200 );
  histo_mll_2j2t->Sumw2();
  
  sprintf(title,"ptsys_2j2t_%s",plotName);
  TH1F* histo_ptsys_2j2t = new TH1F( title, " ", 100,  0, 200 );
  histo_ptsys_2j2t->Sumw2();
  
  sprintf(title,"ht_2j2t_%s",plotName);
  TH1F* histo_ht_2j2t = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_2j2t->Sumw2();
  
  sprintf(title,"pt_leading_2j2t_%s",plotName);
  TH1F* histo_pt_leading_2j2t = new TH1F( title, " ", 100,  0, 200 );
  histo_pt_leading_2j2t->Sumw2();
  
  sprintf(title,"nvertex_2j2t_%s",plotName);
  TH1F* histo_nvertex_2j2t = new TH1F( title, " ", 30,   -0.5, 29.5 );
  histo_nvertex_2j2t->Sumw2();

  
  // all regions
  sprintf(title,"R_%s",plotName);
  TH1F* histo_R = new TH1F( title, " ", 40,  0, 40 );
  histo_R->Sumw2();
  
  
  // checking pu reweighting
  sprintf(title,"nvertex_final_%s",plotName);
  TH1F* histo_nvertex_final = new TH1F( title, " ", 70,   -0.5, 69.5 );
  histo_nvertex_final->Sumw2();
  
  sprintf(title,"nvertex_final_3D_%s",plotName);
  TH1F* histo_nvertex_final_3D = new TH1F( title, " ", 70,   -0.5, 69.5 );
  histo_nvertex_final_3D->Sumw2();
  
  sprintf(title,"nvertex_final_purw_%s",plotName);
  TH1F* histo_nvertex_final_purw = new TH1F( title, " ", 70,   -0.5, 69.5 );
  histo_nvertex_final_purw->Sumw2();
  


  sprintf(title,"nvertex_2lep_%s",plotName);
  TH1F* histo_nvertex_2lep = new TH1F( title, " ", 70,   -0.5, 69.5 );
  histo_nvertex_2lep->Sumw2();

  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    
    if (lumi != lum && nsel != 666 && mode !=3){
      if (jentry == 0)std::cout << "[Warning:] This tree was made with a different luminosity (" << lum << ") than " << lumi << std::endl;
      //xlWeight*=(lumi/lum);
    }
    
    if(ptLepton->size() != 2){
      std::cout << "[Warning:] Something is wrong, your Tree is not correctly filled" << std::endl;
      break;
    } else {
      histo->Fill(0.,xlWeight);
      histo_nvertex_2lep->Fill(nvertex,xlWeight);
      TLorentzVector lepton0(pxLepton->at(0),pyLepton->at(0), pzLepton->at(0), eLepton->at(0));
      TLorentzVector lepton1(pxLepton->at(1),pyLepton->at(1), pzLepton->at(1), eLepton->at(1)); 
      TLorentzVector pair = lepton0+lepton1;
      
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
      
      if (pair.M() > 20){
	histo->Fill(1, xlWeight);

	double SFval, SFerror;
	if ( nsel == 666 || nosf){
	  SFval = 1;
	  SFerror = 0;
	} else if (nsel == 0){
	  SFval = 0.95;
	  SFerror = 0.03;
	} else {
	  SFval = 0.97;
	  SFerror = 0.03;
	}
	
	
	int nJetsBT = 0;
	int nTightJetsBT = 0;
	int nJets = 0;
	bool bTagged = false;
	int iJet = -5;
	int iSF;
	double tempSF = SFval;
	
	int SFvalue = int(tempSF*100);
	
	for (unsigned int i =0; i < ptJet->size(); i ++){ 
	  TLorentzVector tempJet(pxJet->at(i),pyJet->at(i), pzJet->at(i), eJet->at(i));
	  if (ptJet->at(i) > 30 && fabs(tempJet.Eta()) < 2.5 && TMath::Min(fabs(lepton0.DeltaR(tempJet)), fabs(lepton1.DeltaR(tempJet))) > 0.3) {
	    nJets++;
	    iJet = i;
	    if (btCSVBJet->at(i) > 0.679){
	      iSF = rand() % 100;
	      if (iSF < SFvalue ){
		bTagged = true;
		nJetsBT++;
		nTightJetsBT++;
	      } 
	    } 
	  } else if (btCSVBJet->at(i) > 0.679 && fabs(tempJet.Eta()) < 2.5){
	    iSF = rand() % 100;
	    if (iSF < SFvalue ) nJetsBT++;
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
	  histo->Fill(2, xlWeight);
	  histo_mll_after->Fill(pair.M(),  xlWeight);
	  histo_met_cut->Fill(metPt,  xlWeight);
	  if (metPt >= metCut || mode ==0){
	   
	    histo->Fill(3, xlWeight);
	    histo_njets_cut->Fill(nJets, xlWeight);
	    if (nJets == 1){
	      histo->Fill(4, xlWeight);
	      histo_njetsbt_cut->Fill(nJetsBT, xlWeight);
	      TLorentzVector jet(pxJet->at(iJet),pyJet->at(iJet), pzJet->at(iJet), eJet->at(iJet));
		
	      double ptSysPx = lepton0.Px() + lepton1.Px() + jet.Px() + metPx;
	      double ptSysPy = lepton0.Py() + lepton1.Py() + jet.Py() + metPy;
	      double ptSystem = sqrt(ptSysPx*ptSysPx + ptSysPy*ptSysPy);
	      double ht = lepton0.Pt() + lepton1.Pt() + jet.Pt() + metPt; 
	      
	      histo_mll_1j->Fill(pair.M(),  xlWeight);
	      histo_met_1j->Fill(metPt,  xlWeight);
	      histo_ptsys_1j->Fill(ptSystem, xlWeight);
	      histo_ht_1j->Fill(ht, xlWeight);
	      histo_pt_leading_1j->Fill(jet.Pt(), xlWeight);
	      
	      
	      if (nJets == 1 && nTightJetsBT == 1 && bTagged && nJetsBT == 1){
		histo->Fill(5, xlWeight);
		
	      
	        histo_mll_1j1t->Fill(pair.M(),  xlWeight);
	        histo_met_1j1t->Fill(metPt,  xlWeight);
	        histo_ptsys_1j1t->Fill(ptSystem, xlWeight);
	        histo_ht_1j1t->Fill(ht, xlWeight);
	        histo_pt_leading_1j1t->Fill(jet.Pt(), xlWeight);
	      
		histo_ptsys->Fill(ptSystem, xlWeight);
		histo_ht->Fill(ht, xlWeight);
		
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
		
		if (ht > htMin || mode !=0){
		  histo->Fill(6, xlWeight);
		  histo_ht_cut->Fill(ht, xlWeight);
		  
		  //Example to access the pu reweighting!
		  histo_nvertex_final->Fill(nvertex, rawWeight);
		  histo_nvertex_final_3D->Fill(nvertex, rawWeight*puweight3D);
		  histo_nvertex_final_purw->Fill(nvertex, rawWeight*puweight);
		  
		}
	      }
	    }
	  }
	  
	  
	  //Filling of all region from here
	  if (metPt >= metCut || mode ==0){
	
	    if (nJets != 0){
	      
	      TLorentzVector jet(pxJet->at(iJet),pyJet->at(iJet), pzJet->at(iJet), eJet->at(iJet));
	      
	      double ptSysPx = lepton0.Px() + lepton1.Px() + jet.Px() + metPx;
	      double ptSysPy = lepton0.Py() + lepton1.Py() + jet.Py() + metPy;
	      double ptSystem = sqrt(ptSysPx*ptSysPx + ptSysPy*ptSysPy);
	      double ht = lepton0.Pt() + lepton1.Pt() + jet.Pt() + metPt; 
	      
	      if (nJets == 2 && nTightJetsBT == 1 && nJetsBT == 1) {
	        histo_mll_2j1t->Fill(pair.M(),  xlWeight);
	        histo_met_2j1t->Fill(metPt,  xlWeight);
	        histo_ptsys_2j1t->Fill(ptSystem, xlWeight);
	        histo_ht_2j1t->Fill(ht, xlWeight);
	        histo_pt_leading_2j1t->Fill(jet.Pt(), xlWeight);
	      } else if (nJets == 2 && nTightJetsBT == 2 && nJetsBT == 2)  {
	        histo_mll_2j2t->Fill(pair.M(),  xlWeight);
	        histo_met_2j2t->Fill(metPt,  xlWeight);
	        histo_ptsys_2j2t->Fill(ptSystem, xlWeight);
	        histo_ht_2j2t->Fill(ht, xlWeight);
	        histo_pt_leading_2j2t->Fill(jet.Pt(), xlWeight);
	      }
	      
	      
	   
		  
	      //All possible regions
	      if (nJets == 1 && nTightJetsBT == 1 && bTagged && nJetsBT == 1 && ptSystem <= ptsysCut && (ht > htMin || mode !=0))histo_R->Fill(1, xlWeight); //signal
	      if (nJets == 1 && nTightJetsBT == 2)  histo_R->Fill(2, xlWeight);
	      if (nJets == 1 && nTightJetsBT > 0)  histo_R->Fill(3, xlWeight);
	      if (nJets == 1 && nTightJetsBT > 1)  histo_R->Fill(4, xlWeight);
	      if (nJets == 2 && nTightJetsBT == 0)  histo_R->Fill(5, xlWeight);
	      if (nJets == 2 && nTightJetsBT == 1)  histo_R->Fill(6, xlWeight); //CR1 no ht no ptsys
	      if (nJets == 2 && nTightJetsBT == 2)  histo_R->Fill(7, xlWeight); //CR2 no ht no ptsys
	      if (nJets == 2 && nTightJetsBT > 0)  histo_R->Fill(8, xlWeight);
	      if (nJets == 2 && nTightJetsBT > 1)  histo_R->Fill(9, xlWeight);
	      if (nJets > 1 && nTightJetsBT == 0)  histo_R->Fill(10, xlWeight);
	      if (nJets > 1 && nTightJetsBT == 1)  histo_R->Fill(11, xlWeight);
	      if (nJets > 1 && nTightJetsBT == 2)  histo_R->Fill(12, xlWeight);
	      if (nJets > 1 && nTightJetsBT !=0 )  histo_R->Fill(13, xlWeight);
	      if (nJets > 1 && nTightJetsBT > 1 )  histo_R->Fill(14, xlWeight);
	      if (nJets == 3 && nTightJetsBT ==3 )  histo_R->Fill(15, xlWeight);
	      if (nJets == 1 && nTightJetsBT ==1 && bTagged && nJetsBT == 1)  histo_R->Fill(16, xlWeight);
	      if (nJets == 2 && nTightJetsBT == 1 && ptSystem <= ptsysCut && (ht > htMin || mode !=0))  histo_R->Fill(17, xlWeight); //CR 1 regular
	      if (nJets == 2 && nTightJetsBT == 2 && ptSystem <= ptsysCut && (ht > htMin || mode !=0))  histo_R->Fill(18, xlWeight); //CR 2 regular
	      if (nJets == 2 && nTightJetsBT == 1 && nJetsBT == 1 && ptSystem <= ptsysCut && (ht > htMin || mode !=0))  histo_R->Fill(19, xlWeight);
	      if (nJets == 2 && nTightJetsBT == 2 && nJetsBT == 2 && ptSystem <= ptsysCut && (ht > htMin || mode !=0))  histo_R->Fill(20, xlWeight);
	      if (nJets == 2 && nJetsBT == 1 && ptSystem <= ptsysCut && (ht > htMin || mode !=0))  histo_R->Fill(21, xlWeight); //CR 1 another way
	      if (nJets == 2 && nJetsBT == 2 && ptSystem <= ptsysCut && (ht > htMin || mode !=0))  histo_R->Fill(22, xlWeight); //CR 2 another way
	      if (nJets == 2 && nTightJetsBT == 1 && nJetsBT == 1)  histo_R->Fill(23, xlWeight); //CR1 no ht no ptsys tighter
	      if (nJets == 2 && nTightJetsBT == 2 && nJetsBT == 2)  histo_R->Fill(24, xlWeight); //CR2 no ht no ptsys tighter
	      if (nJets == 2 && nJetsBT == 1)  histo_R->Fill(25, xlWeight); //CR1 no ht no ptsys another flavor
	      if (nJets == 2 && nJetsBT == 2)  histo_R->Fill(26, xlWeight); //CR2 no ht no ptsys another flavor
	      if (nJets == 1 && nTightJetsBT == 1 && bTagged && nJetsBT == 1 && (ht > htMin || mode !=0))histo_R->Fill(27, xlWeight); //signal no ptsys
	      if (nJets == 1 && nTightJetsBT == 1 && bTagged && nJetsBT == 1 && ptSystem <= ptsysCut)histo_R->Fill(28, xlWeight); //signal no ht
	      if (nJets == 2 && nTightJetsBT == 1 &&  (ht > htMin || mode !=0))  histo_R->Fill(29, xlWeight); //CR 1 
	      if (nJets == 2 && nTightJetsBT == 2 &&  (ht > htMin || mode !=0))  histo_R->Fill(30, xlWeight); //CR 2 
				
	      
	    } //jets in the event
	  } //all CR
	  
	} // mll
      } //mll pre
    } // 2 leptons
  }// event loop.
  
  
  
  if (!silent){ 
    cout << "------------------------------------------" << endl;
    cout << "[Results:] " << plotName <<  endl;
    cout << "------------------------------------------" << endl;  
    for (int i = 2; i < 9; i++){
      if (i == 2) cout << " leptons: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << endl;
      if (i == 3) cout << " inv. mass: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << endl;
      if (i == 4) cout << " met: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << endl;
      if (i == 5) cout << " jet: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << endl;
      if (i == 6) cout << " jet_bt: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << endl;
      if (i == 7) cout << " ht: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << endl;
    }
    cout << "------------------------------------------" << endl; 
  }
  f_var.Write();
  f_var.Close();
}
