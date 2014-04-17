// rebeca@cern.ch
#include "TLorentzVector.h"
#include "TVector3.h"
#include "inputs.h"

void dycontrol(int nsel, int mode = 0, bool silent = false){  
  

  bool SFplus = false;
  bool SFminus = false;
  // samples used
  double x_sec = 0.;
  char plotName[300];
  sprintf(plotName,"test");
  
  if (nsel == 0)                	{sprintf(plotName,"tt");}
  else if (nsel == 1)   		{sprintf(plotName,"twdr");}
  else if (nsel == -1)   		{sprintf(plotName,"twds");}
  else if (nsel == 2)			{sprintf(plotName, "st");}
  else if (nsel == 3)   		{sprintf(plotName,"wjets");}
  else if (nsel == 4)   		{sprintf(plotName,"zjetsall");}
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
    if      (mode == 0) 	cout << " Electron-Muon Mixed channel " << endl;
    else if (mode == 1) 	cout << " Di-Muon channel " << endl;
    else if (mode == 2) 	cout << " Di-Electron channel " << endl;
    cout << "*************************************************" << endl;
  }
  
  bool nosf = true;
  char myRootFile[300];
  sprintf(myRootFile,"outputs/out_%d_%s.root", mode, plotName);
  TFile *input = TFile::Open(myRootFile);
  
  // tree variables
  /////
  double xlWeight; 
  double lum;
  
  double metPt;
  double metPx;
  double metPy;
  
  int npu;
  int nvertex;
      
  std::vector<double> *ptLepton;
  std::vector<double> *pxLepton;
  std::vector<double> *pyLepton;
  std::vector<double> *pzLepton;
  std::vector<double> *eLepton;
  std::vector<double> *qLepton;
  
  std::vector<double> *ptJet;
  std::vector<double> *pxJet;
  std::vector<double> *pyJet;
  std::vector<double> *pzJet;
  std::vector<double> *eJet;
  std::vector<double> *qJet;
  std::vector<double> *btSSVHEJet;
  std::vector<double> *btSSVHPJet;
  std::vector<double> *btTCHEJet;
  std::vector<double> *btTCHPJet;
  
  TTree*Tree = (TTree*) gROOT->FindObject("myTree");
  
  Tree->SetBranchAddress("xlWeight", &xlWeight);
  Tree->SetBranchAddress("lum", &lum);
  
  Tree->SetBranchAddress("npu", &npu);
  Tree->SetBranchAddress("nvertex", &nvertex);
  
  Tree->SetBranchAddress("metPt", &metPt);
  Tree->SetBranchAddress("metPx", &metPx);
  Tree->SetBranchAddress("metPy", &metPy);
  
  Tree->SetBranchAddress("ptLepton",&ptLepton);
  Tree->SetBranchAddress("pxLepton",&pxLepton);
  Tree->SetBranchAddress("pyLepton",&pyLepton);
  Tree->SetBranchAddress("pzLepton",&pzLepton);
  Tree->SetBranchAddress("eLepton",&eLepton);
  Tree->SetBranchAddress("qLepton",&qLepton);
  
  Tree->SetBranchAddress("ptJet",&ptJet);
  Tree->SetBranchAddress("pxJet",&pxJet);
  Tree->SetBranchAddress("pyJet",&pyJet);
  Tree->SetBranchAddress("pzJet",&pzJet);
  Tree->SetBranchAddress("eJet",&eJet);
  Tree->SetBranchAddress("qJet",&qJet);
  Tree->SetBranchAddress("btSSVHEJet",&btSSVHEJet);
  Tree->SetBranchAddress("btSSVHPJet",&btSSVHPJet);
  Tree->SetBranchAddress("btTCHEJet",&btTCHEJet);
  Tree->SetBranchAddress("btTCHPJet",&btTCHPJet);
  
  int nEvents = Tree->GetEntries();
  if(!silent){
    cout << endl;
    cout << "******************************************" << endl;
    cout << "------------------------------------------" << endl;
    cout << "Starting the analysis: " << plotName <<  endl;
    cout << "------------------------------------------" << endl;
    cout << "Number of Raw events: " <<  nEvents << endl;
    cout << "------------------------------------------" << endl;
    cout << "******************************************" << endl;
  }
  
  char newRootFile[300];
  double lumi = luminosity;
  if (mode == 0 )        lumi = 4626.297;
  else if ( mode == 1)   lumi = 4534.871;
  else if ( mode == 2)   lumi = 4593.348;
  sprintf(newRootFile,"results/crdy_%dpb_%d.root", lumi, mode);
  
  TFile f_var(newRootFile, "UPDATE");
  
  //////////
  char title[300];
  
  sprintf(title,"R_%s",plotName);
  TH1F* histo_R = new TH1F( title, " ", 40,  0, 40 );
  histo_R->Sumw2();
  
  //////////
  for(int event = 0; event<nEvents; event++){
    
    Tree->GetEntry(event);
    
    if (lumi != lum && nsel != 666 && mode !=3){
      if (event == 0) cout << "Warning: This tree was made with a different luminosity (" << lum << ") than " << lumi << endl;
      xlWeight*=(lumi/lum);
    }
    
    if(ptLepton->size() != 2) cout << "Something is wrong, your Tree is not correctly filled" << endl;
    else {
      
      TLorentzVector lepton0(pxLepton->at(0),pyLepton->at(0), pzLepton->at(0), eLepton->at(0));
      TLorentzVector lepton1(pxLepton->at(1),pyLepton->at(1), pzLepton->at(1), eLepton->at(1)); 
      TLorentzVector pair = lepton0+lepton1;
      
      if (pair.M() > 20){
        if (pair.M() > invMax || pair.M() < invMin) histo_R->Fill(1, xlWeight);
        else histo_R->Fill(2, xlWeight);

	// double SFval = 0.95;  //Summer11 version
	double SFval, SFerror;
	if ( nsel == 666 || !nosf){
	  SFval = 1;
	  SFerror = 0;
	} else if (nsel == 0){
	  SFval = 0.956;
	  SFerror = 0.030;
	} else {
	  SFval = 0.96;
	  SFerror = 0.04;
	}
	
	
	int nJetsBT = 0;
	int nTightJetsBT = 0;
	int nJets = 0;
	bool bTagged = false;
	int iJet = -5;
	int iSF;
	double tempSF = SFval;
	if (SFminus) 	tempSF = SFval - SFerror;
	if (SFplus) 	tempSF = SFval + SFerror;
	int SFvalue = int(tempSF*100);

	for (int i =0; i < ptJet->size(); i ++){ 
	  TLorentzVector tempJet(pxJet->at(i),pyJet->at(i), pzJet->at(i), eJet->at(i));
	  if (ptJet->at(i) > 30 && TMath::Min(fabs(lepton0.DeltaR(tempJet)), fabs(lepton1.DeltaR(tempJet))) > 0.3) {
	    nJets++;
	    iJet = i;
	    if (btSSVHEJet->at(i) > 1.74){
	      iSF = rand() % 100;
	      if (iSF < SFvalue ){
		bTagged = true;
		nJetsBT++;
		nTightJetsBT++;
	      } 
	    } 
	  } else if (btSSVHEJet->at(i) > 1.74){
	    iSF = rand() % 100;
	    if (iSF < SFvalue ) nJetsBT++;
	  }
	}
	//
      

	//
	
	
	bool invMass = true;
	
	if (invMass){
	  if (metPt >= metCut || mode ==0){
	    
	    if (nJets == 1 && nJetsBT == 1 && bTagged ){
	      
	      TLorentzVector jet(pxJet->at(iJet),pyJet->at(iJet), pzJet->at(iJet), eJet->at(iJet));
	      
	      double ptSysPx = lepton0.Px() + lepton1.Px() + jet.Px() + metPx;
	      double ptSysPy = lepton0.Py() + lepton1.Py() + jet.Py() + metPy;
	      double ptSystem = sqrt(ptSysPx*ptSysPx + ptSysPy*ptSysPy);
	      double ht = lepton0.Pt() + lepton1.Pt() + jet.Pt() + metPt; 
	      
	      if (ptSystem <= ptsysCut){
		
		if (ht > htMin || mode !=0){
		  if (pair.M() > invMax || pair.M() < invMin) histo_R->Fill(3, xlWeight);
                  else histo_R->Fill(4, xlWeight);  
		  
		}
	      }
	      
	    }
	  }
	}
      } //mll
    } // 2 leptons
  } // events
  
  
  if (!silent){ 
    cout << plotName ;
    cout << " ---> dy cr done" << endl;
   
    cout << "Lepton selection: " << endl;
    cout <<  "In: " << histo_R->GetBinContent(3) << "+-" <<  histo_R->GetBinError(3) << ", Out: " <<  histo_R->GetBinContent(2) << "+-" <<  histo_R->GetBinError(2) << endl;
    cout << "Final cut: " << endl;
    cout <<  "In: " << histo_R->GetBinContent(5) << "+-" <<  histo_R->GetBinError(5) << ", Out: " <<  histo_R->GetBinContent(4) << "+-" <<  histo_R->GetBinError(4) << endl;
  }
  f_var.Write();
  f_var.Close();
  
  
}
