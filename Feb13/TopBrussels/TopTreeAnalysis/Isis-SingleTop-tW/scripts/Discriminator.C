#include "Discriminator.h"
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

void Discriminator::Loop(){
  // running default loop:
  myLoop(0,0,0);
}
void Discriminator::myLoop(int nsel, int mode, bool silent)
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
  double lumi = 4399; 

  sprintf(newRootFile,"results/dis_an_%dpb_%d.root", (int)lumi, mode);
 
  TFile f_var(newRootFile, "UPDATE");
  
  if(!silent){
    std::cout << "[Info:] results root file " << newRootFile << std::endl;
  }
  
  
  //////////
  char title[300];
  
  
  sprintf(title,"promet_%s",plotName);
  TH1F* histo_promet = new TH1F( title, " ", 100,  0, 200 );
  histo_promet->Sumw2();
  
 
  
  
  sprintf(title,"pt_max_%s",plotName);
  TH1F* histo_pt_max = new TH1F( title, " ", 100,  0, 200 );
  histo_pt_max->Sumw2();
  
  sprintf(title,"pt_min_%s",plotName);
  TH1F* histo_pt_min = new TH1F( title, " ", 100,  0, 200 );
  histo_pt_min->Sumw2();
  

  
 
  sprintf(title,"met_%s",plotName);
  TH1F* histo_met = new TH1F( title, " ", 100,  0, 200 );
  histo_met->Sumw2();
  
  sprintf(title,"mll_%s",plotName);
  TH1F* histo_mll = new TH1F( title, " ", 100,  0, 200 );
  histo_mll->Sumw2();
  
  sprintf(title,"njets_%s",plotName);
  TH1F* histo_njets = new TH1F( title, " ", 10,  0, 10 );
  histo_njets->Sumw2();
  
  sprintf(title,"ptsys_%s",plotName);
  TH1F* histo_ptsys = new TH1F( title, " ", 100,  0, 200 );
  histo_ptsys->Sumw2();
  
  sprintf(title,"ht_%s",plotName);
  TH1F* histo_ht = new TH1F( title, " ", 300,  0, 600 );
  histo_ht->Sumw2();
  
  sprintf(title,"pt_leading_%s",plotName);
  TH1F* histo_pt_leading = new TH1F( title, " ", 100,  0, 200 );
  histo_pt_leading->Sumw2();
  
    sprintf(title,"pt_all_%s",plotName);
  TH1F* histo_pt_all = new TH1F( title, " ", 100,  0, 200 );
  histo_pt_all->Sumw2();
  
  
  sprintf(title,"eta_leading_%s",plotName);
  TH1F* histo_eta_leading = new TH1F( title, " ", 101,  -3, 3);
  histo_eta_leading->Sumw2();
 
  sprintf(title,"nvertex_2lep_%s",plotName);
  TH1F* histo_nvertex_2lep = new TH1F( title, " ", 70,   -0.5, 69.5 );
  histo_nvertex_2lep->Sumw2();


//__________________________________________________ END HISTO DEF _________________________________________________________


  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    


    
    if(ptLepton->size() != 2){
      std::cout << "[Warning:] Something is wrong, your Tree is not correctly filled" << std::endl;
      break;
    } 
    else {
    	
      
        // Fill the number of vertices for events with at least 2 leptons
     	 histo_nvertex_2lep->Fill(nvertex,xlWeight);
      
      	// Take the lepton with highest, and second highest pt
      	TLorentzVector lepton0(pxLepton->at(0),pyLepton->at(0), pzLepton->at(0), eLepton->at(0));
      	TLorentzVector lepton1(pxLepton->at(1),pyLepton->at(1), pzLepton->at(1), eLepton->at(1)); 
      
      	// Define the lepton pair
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
      
      	// Demand that mll is higher than 20 GeV
      	if (pair.M() > 20){
		
		int nJets = 0;
		int iJetn[5]={-1, -1,-1,-1,-1};
		int iJet = -5;
	        // define number of jets 
		for (unsigned int i =0; i < ptJet->size(); i ++){ 
			TLorentzVector tempJet(pxJet->at(i),pyJet->at(i), pzJet->at(i), eJet->at(i));
			histo_pt_all->Fill(tempJet.Pt(),xlWeight);
		
	  		
	  		if (ptJet->at(i) > 30 && TMath::Min(fabs(lepton0.DeltaR(tempJet)), fabs(lepton1.DeltaR(tempJet))) > 0.3) {
	        		iJetn[nJets] = i;
	    			iJet = i;
	    			nJets++;
	     		}  // end if loop
		} // end for loop
	

	
	
		histo_pt_max->Fill(TMath::Max(lepton0.Pt(), lepton1.Pt()), xlWeight);
		histo_pt_min->Fill(TMath::Min(lepton0.Pt(), lepton1.Pt()), xlWeight);
		histo_njets->Fill(nJets,  xlWeight);
		//histo_njetsbt->Fill(nJetsBT,  xlWeight);
		histo_mll->Fill(pair.M(),  xlWeight);
		histo_met->Fill(metPt,  xlWeight);
		histo_promet->Fill(promet, xlWeight);
	


	
		if (nJets) {
		          		
			histo_pt_leading->Fill(ptJet->at(0), xlWeight);
			
			TLorentzVector jet_aux(pxJet->at(0),pyJet->at(0), pzJet->at(0), eJet->at(0));
	  		histo_eta_leading->Fill(jet_aux.Eta(), xlWeight);
			
			
		}
	
	
	

	
	  
	  
	  	//Filling of all region from here
	  
	
	    	if (nJets != 0){
	      
	      		TLorentzVector jet(pxJet->at(iJet),pyJet->at(iJet), pzJet->at(iJet), eJet->at(iJet));
	      
	      		double ptSysPx = lepton0.Px() + lepton1.Px() + jet.Px() + metPx;
	      		double ptSysPy = lepton0.Py() + lepton1.Py() + jet.Py() + metPy;
        		double ptSystem = sqrt(ptSysPx*ptSysPx + ptSysPy*ptSysPy);
	      		double ht = lepton0.Pt() + lepton1.Pt() + jet.Pt() + metPt;
			
			histo_ptsys->Fill(ptSystem, xlWeight);
	  		histo_ht->Fill(ht, xlWeight);
	    	} //jets in the event
      	} //mll pre
    } // 2 leptons
  }// event loop.
  
  
  
  
  f_var.Write();
  f_var.Close();
}
