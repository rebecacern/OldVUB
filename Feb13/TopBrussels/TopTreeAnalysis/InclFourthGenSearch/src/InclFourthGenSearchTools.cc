#include "../interface/InclFourthGenSearchTools.h"

InclFourthGenSearchTools::InclFourthGenSearchTools(bool semiMuon, bool semiElectron, vector<Dataset* > datasets, float Luminosity, bool doKinematicFit)
{ 
  semiMuon_ = semiMuon;
  semiElectron_ = semiElectron;
  for(unsigned int d=0;d<datasets.size();d++)
  {
      datasets_.push_back(datasets[d]);
  }
  Luminosity_ = Luminosity;
  doKinematicFit_ = doKinematicFit;
  resFitLightJets_ = 0; 
  
  string histoName_HT,MShistoName_HT,MShistoName,histoHT_dataset,histoName_nEvents,MShistoName_nEvents,histonEvents_dataset,histoName_Mtop,histoMtop_dataset;
  string histoPreFixes_SemiMu[4] = {"JetsPt","MET","MuonsPt","JetsHt"};
  string histoPreFixes_SemiEl[4] = {"JetsPt","MET","ElectronsPt","JetsHt"};
  for(int B = 1; B<3; B++)
  {
		for(int W = 1; W<=4; W++)
		{
			histoName_HT = "HT_"+IntToStr(B)+"B_"+IntToStr(W)+"W";
			MShistoName_HT = "MS_"+histoName_HT;
			histoName_nEvents = "nEvents_"+IntToStr(B)+"B_"+IntToStr(W)+"W";
			MShistoName_nEvents = "MS_"+histoName_nEvents;
			histoName_Mtop = "Mtop_MVA_"+IntToStr(B)+"B_"+IntToStr(W)+"W";
			if(B==1 && W==1) {
				MSPlot[MShistoName_HT.c_str()] = new MultiSamplePlot(datasets,MShistoName_HT.c_str(), 9, 150, 600, "HT (GeV)");
			}
			else if(B==2 && W==1)
			{
				MSPlot[MShistoName_HT.c_str()] = new MultiSamplePlot(datasets,MShistoName_HT.c_str(), 15, 150, 900, "HT (GeV)");
			}
			else if((B==1 || B==2) && W==2)
			{
				MSPlot[MShistoName_HT.c_str()] = new MultiSamplePlot(datasets,MShistoName_HT.c_str(), 17, 150, 1000, "HT (GeV)");
			}
			else if((B==1 || B==2) && W==3)
			{
				MSPlot[MShistoName_HT.c_str()] = new MultiSamplePlot(datasets,MShistoName_HT.c_str(), 8, 200, 1000, "HT (GeV)");
			}
			else if((B==1 || B==2) && W==4)
			{
				MSPlot[MShistoName_HT.c_str()] = new MultiSamplePlot(datasets,MShistoName_HT.c_str(), 2, 400, 1200, "HT (GeV)");
				MSPlot[MShistoName_nEvents.c_str()] = new MultiSamplePlot(datasets,MShistoName_nEvents.c_str(), 1, 0.5, 1.5, "");
			}
			
			for(int prefix=0; prefix<4; prefix++)
			{
				if(semiMuon)
				{
				  MShistoName = "MS_"+histoPreFixes_SemiMu[prefix]+IntToStr(B)+"B_"+IntToStr(W)+"W";
				  MSPlot[MShistoName.c_str()] = new MultiSamplePlot(datasets,MShistoName.c_str(), 50, 0, 300, histoPreFixes_SemiMu[prefix]);			  
				  if(prefix==3) MSPlot[MShistoName.c_str()] = new MultiSamplePlot(datasets,MShistoName.c_str(), 80, 0, 800, histoPreFixes_SemiMu[prefix]);			  
				}
				else if(semiElectron)
				{
				  MShistoName = "MS_"+histoPreFixes_SemiEl[prefix]+IntToStr(B)+"B_"+IntToStr(W)+"W";
				  MSPlot[MShistoName.c_str()] = new MultiSamplePlot(datasets,MShistoName.c_str(), 50, 0, 300, histoPreFixes_SemiEl[prefix]);			        
					if(prefix==3) MSPlot[MShistoName.c_str()] = new MultiSamplePlot(datasets,MShistoName.c_str(), 80, 0, 800, histoPreFixes_SemiEl[prefix]);			  
				}
			}
			
			for(unsigned int d = 0; d < datasets.size (); d++){
				histoHT_dataset = histoName_HT+(datasets[d]->Name()).c_str();
				histonEvents_dataset = histoName_nEvents+(datasets[d]->Name()).c_str();
				histoMtop_dataset = histoName_Mtop+(datasets[d]->Name()).c_str();
				if(B==1 && W==1) histo1D[histoHT_dataset.c_str()] = new TH1F(histoHT_dataset.c_str(),histoHT_dataset.c_str(), 9, 150, 600);
				else if(B==2 && W==1) histo1D[histoHT_dataset.c_str()] = new TH1F(histoHT_dataset.c_str(),histoHT_dataset.c_str(), 15, 150, 900);
				else if((B==1 || B==2) && W==2)
				{
					histo1D[histoHT_dataset.c_str()] = new TH1F(histoHT_dataset.c_str(),histoHT_dataset.c_str(), 17, 150, 1000);
 					histo1D[histoMtop_dataset.c_str()] = new TH1F(histoMtop_dataset.c_str(),histoMtop_dataset.c_str(), 80, 0, 700);
				}
				else if((B==1 || B==2) && W==3) histo1D[histoHT_dataset.c_str()] = new TH1F(histoHT_dataset.c_str(),histoHT_dataset.c_str(), 8, 200, 1000);
				else if((B==1 || B==2) && W==4)
				{ 
					histo1D[histoHT_dataset.c_str()] = new TH1F(histoHT_dataset.c_str(),histoHT_dataset.c_str(), 2, 400, 1200);
				  histo1D[histonEvents_dataset.c_str()] = new TH1F(histonEvents_dataset.c_str(),histonEvents_dataset.c_str(), 1, 0.5, 1.5);
				}
			}
		}
  }

  //Reconstructed top mass plots for the 1B_2W and 2B_2W boxes... Names should match these in the FillMassPlots function
  MSPlot["MS_MTop_MVA_1B_2W"] = new MultiSamplePlot(datasets,"MTop_1B_2W (MVA jet combination)", 80, 0, 700, "m_{top} (GeV)");
  MSPlot["MS_MTop_MVA_2B_2W"] = new MultiSamplePlot(datasets,"MTop_2B_2W (MVA jet combination)", 80, 0, 700, "m_{top} (GeV)");
//  MSPlot["MS_MTop_MVAkinfit_1B_2W"] = new MultiSamplePlot(datasets,"MTop_1B_2W (MVA jet combination + kinfit)", 80, 0, 700, "m_{top} (GeV)");
//  MSPlot["MS_MTop_MVAkinfit_2B_2W"] = new MultiSamplePlot(datasets,"MTop_2B_2W (MVA jet combination + kinfit)", 80, 0, 700, "m_{top} (GeV)");
  
//  MSPlot["MS_MTop_largestMass_WJetsFixed_2B_2W"] = new MultiSamplePlot(datasets,"MTop_2B_2W (jet combination with largest Top mass, W-jets fixed)", 80, 0, 700, "m_{top} (GeV)");
//  MSPlot["MS_MTop_smallestMass_WJetsFixed_2B_2W"] = new MultiSamplePlot(datasets,"MTop_2B_2W (jet combination with smallest Top mass, W-jets fixed)", 80, 0, 700, "m_{top} (GeV)");
//  MSPlot["MS_MTop_largestPt_WJetsFixed_2B_2W"] = new MultiSamplePlot(datasets,"MTop_2B_2W (jet combination with largest Top Pt, W-jets fixed)", 80, 0, 700, "m_{top} (GeV)");
//  MSPlot["MS_MTop_smallestPt_WJetsFixed_2B_2W"] = new MultiSamplePlot(datasets,"MTop_2B_2W (jet combination with smallest Top Pt, W-jets fixed)", 80, 0, 700, "m_{top} (GeV)");
		

  HT_ = -9999.;
  Mtop_ = -9999.;
  Mtop_kinfit_ = -9999.;
  
  counter_all4JetsMatched_MCdef = 0;
  counter_hadronictopJetsMatched_MCdef = 0;
  counter_hadronictopJetsMatched_largestMass = 0;
  counter_hadronictopJetsMatched_smallestMass = 0;
  counter_hadronictopJetsMatched_largestPt = 0;
  counter_hadronictopJetsMatched_smallestPt = 0;
  counter_hadronictopJetsMatched_largestMass_WJetsFixed = 0;
  counter_hadronictopJetsMatched_smallestMass_WJetsFixed = 0;
  counter_hadronictopJetsMatched_largestPt_WJetsFixed = 0;
  counter_hadronictopJetsMatched_smallestPt_WJetsFixed = 0;
	
	//Wjets plots in the different boxes, all with the same range
	for(int B = 1; B<3; B++)
  {
		for(int W = 1; W<=4; W++)
		{
			string CatName_LF = "SameRange_HT_WJets_LF_"+IntToStr(B)+"B_"+IntToStr(W)+"W";			
			string CatName_HF = "SameRange_HT_WJets_HF_"+IntToStr(B)+"B_"+IntToStr(W)+"W";			
			histo1D[CatName_LF.c_str()] = new TH1F(CatName_LF.c_str(),"HT;HT;#events",150,0,1000);
			histo1D[CatName_HF.c_str()] = new TH1F(CatName_HF.c_str(),"HT;HT;#events",150,0,1000);
		}
	}
}

InclFourthGenSearchTools::InclFourthGenSearchTools(const InclFourthGenSearchTools & i)
{
  //copy constructor to be done
}

InclFourthGenSearchTools::~InclFourthGenSearchTools()
{
}

void InclFourthGenSearchTools::SetResolutionFit(ResolutionFit* resFitLightJets)
{
	resFitLightJets_ = resFitLightJets;
}

void InclFourthGenSearchTools::FillPlots(int d, int nbOfBtags, int nbOfWs, float HT, vector<TRootMuon*> selectedMuons, vector<TRootElectron*> selectedElectrons, vector<TRootMET*> mets, vector<TRootJet*> selectedJets, float scaleFactor)
{
  vector<TLorentzVector> selectedMuonsTLV,selectedElectronsTLV,selectedJetsTLV;
	for(unsigned int i=0; i<selectedMuons.size(); i++)
  	selectedMuonsTLV.push_back((TLorentzVector) (*selectedMuons[i]));
	for(unsigned int i=0; i<selectedElectrons.size(); i++)
  	selectedElectronsTLV.push_back((TLorentzVector) (*selectedElectrons[i]));
	for(unsigned int i=0; i<selectedJets.size(); i++)
  	selectedJetsTLV.push_back((TLorentzVector) (*selectedJets[i]));
		
	FillPlots(d, nbOfBtags, nbOfWs, HT, selectedMuonsTLV, selectedElectronsTLV, mets[0]->Et(), selectedJetsTLV, scaleFactor);
}

void InclFourthGenSearchTools::FillPlots(int d, int nbOfBtags, int nbOfWs, float HT, vector<TLorentzVector> selectedMuons, vector<TLorentzVector> selectedElectrons, float met, vector<TLorentzVector> selectedJets, float scaleFactor) 
{
	string histoName_HT,MShistoName_HT,histoName_HT_dataset,MShistoName,histoName_nEvents,MShistoName_nEvents,histoName_nEvents_dataset;
	
	histoName_HT = "HT_"+IntToStr(nbOfBtags)+"B_"+IntToStr(nbOfWs)+"W";
	MShistoName_HT = "MS_"+histoName_HT; //MSPlot for HT
	histoName_HT_dataset = histoName_HT+(datasets_[d]->Name()).c_str(); //1D plot HT for each dataset
	
	histoName_nEvents = "nEvents_"+IntToStr(nbOfBtags)+"B_"+IntToStr(nbOfWs)+"W";
	MShistoName_nEvents = "MS_"+histoName_nEvents; //MSPlot for nEvents for 4W boxes
	histoName_nEvents_dataset = histoName_nEvents+(datasets_[d]->Name()).c_str(); //1D plot HT for each dataset
	
	
	string histoPreFixes_SemiMu[4] = {"JetsPt","MET","MuonsPt","JetsHt"};
  string histoPreFixes_SemiEl[4] = {"JetsPt","MET","ElectronsPt","JetsHt"};
  
	if(semiMuon_) HT_ = HT + selectedMuons[0].Pt();
	else if(semiElectron_) HT_ = HT + selectedElectrons[0].Pt();									
	MSPlot[MShistoName_HT.c_str()] -> Fill(HT_, datasets_[d], true, Luminosity_*scaleFactor);
	histo1D[histoName_HT_dataset.c_str()] -> Fill(HT_,datasets_[d]->NormFactor()*Luminosity_*scaleFactor);
  
	if((nbOfBtags==1 || nbOfBtags==2) && nbOfWs==4)
	{
		MSPlot[MShistoName_nEvents.c_str()] -> Fill(1, datasets_[d], true, Luminosity_*scaleFactor);
		histo1D[histoName_nEvents_dataset.c_str()] -> Fill(1,datasets_[d]->NormFactor()*Luminosity_*scaleFactor);  
	}
	
	for(int prefix=0; prefix<4; prefix++)
	{
		if(semiMuon_) MShistoName = "MS_"+histoPreFixes_SemiMu[prefix]+IntToStr(nbOfBtags)+"B_"+IntToStr(nbOfWs)+"W";
		else if(semiElectron_) MShistoName = "MS_"+histoPreFixes_SemiEl[prefix]+IntToStr(nbOfBtags)+"B_"+IntToStr(nbOfWs)+"W";
						
		if(prefix==0)
		{
			for(unsigned int j =0; j < selectedJets.size(); j++)
			{
				MSPlot[MShistoName.c_str()]->Fill(selectedJets[j].Pt(), datasets_[d], true, Luminosity_*scaleFactor);
			}
		}
		else if(prefix==1)
		{
			MSPlot[MShistoName.c_str()]->Fill(met, datasets_[d], true, Luminosity_*scaleFactor);
		}
		else if(prefix==2)
		{
			if(semiMuon_) MSPlot[MShistoName.c_str()]->Fill(selectedMuons[0].Pt(), datasets_[d], true, Luminosity_*scaleFactor);
			else if(semiElectron_) MSPlot[MShistoName.c_str()]->Fill(selectedElectrons[0].Pt(), datasets_[d], true, Luminosity_*scaleFactor);
		}
		else if(prefix==3)
		{
		  float JetsHT = 0;
		  for(unsigned int j =0; j < selectedJets.size(); j++)
			{
			  JetsHT = JetsHT + selectedJets[j].Pt();
			}
			MSPlot[MShistoName.c_str()]->Fill(JetsHT, datasets_[d], true, Luminosity_*scaleFactor);
		}
							
	}
	
	//Wjets plots in the different boxes, all with the same range
	if(datasets_[d]->Name() == "WJets_LF")
	{
		string CatName_LF = "SameRange_HT_WJets_LF_"+IntToStr(nbOfBtags)+"B_"+IntToStr(nbOfWs)+"W";	
		histo1D[CatName_LF.c_str()]-> Fill(HT_,datasets_[d]->NormFactor()*Luminosity_*scaleFactor);	
	}
	if(datasets_[d]->Name() == "WJets_HF")
	{
		string CatName_HF = "SameRange_HT_WJets_HF_"+IntToStr(nbOfBtags)+"B_"+IntToStr(nbOfWs)+"W";	
		histo1D[CatName_HF.c_str()]-> Fill(HT_,datasets_[d]->NormFactor()*Luminosity_*scaleFactor);
	}
}

void InclFourthGenSearchTools::CalculateTopMass(TRootJet* WJet1, TRootJet* WJet2, TRootJet* HadBJet)
{
  TLorentzVector Top;
	Top = (*WJet1 + *WJet2 + *HadBJet);
	Mtop_ = Top.M();
	
	if(doKinematicFit_)
	{	
		///////////////////
		//for kinematic fit
		///////////////////
		//coutObjectsFourVector(init_muons,init_electrons,selectedJets_MVAinput,mets,"*** Before kinematic fit ***");
		
		bool KinFitConverged = true;
		const TLorentzVector* afterfitLight1 = 0; //const? A variable pointer to a constant TLV... (otherwise doesn't compile)
		const TLorentzVector* afterfitLight2 = 0;
		
		TMatrixD Ml1(3,3), Ml2(3,3), Mb(3,3);
    		Ml1.Zero(); Ml2.Zero(); Mb.Zero();
		TLorentzVector lightJet1 = *WJet1;
		TLorentzVector lightJet2 = *WJet2;
    		//cout<<"  lightJet1->Pt() = "<<lightJet1.Pt()<<endl;
		//cout<<"  lightJet2->Pt() = "<<lightJet2.Pt()<<endl;
		Ml1(0,0) = pow(resFitLightJets_->EtResolution(&lightJet1), 2);
    		Ml1(1,1) = pow(resFitLightJets_->ThetaResolution(&lightJet1), 2);
    		Ml1(2,2) = pow(resFitLightJets_->PhiResolution(&lightJet1), 2);
    		Ml2(0,0) = pow(resFitLightJets_->EtResolution(&lightJet2), 2);
    		Ml2(1,1) = pow(resFitLightJets_->ThetaResolution(&lightJet2), 2);
    		Ml2(2,2) = pow(resFitLightJets_->PhiResolution(&lightJet2), 2); 
		TKinFitter *theFitter = new TKinFitter("fit", "fit");
	        theFitter->setVerbosity(0);
	        TFitParticleEtThetaPhi *fitLight1 = new TFitParticleEtThetaPhi("lightJet1", "lightJet1", &lightJet1, &Ml1);
	        TFitParticleEtThetaPhi *fitLight2 = new TFitParticleEtThetaPhi("lightJet2", "lightJet2", &lightJet2, &Ml2);

	        theFitter->addMeasParticles(fitLight1,fitLight2);
	        
 	        TFitConstraintM *consW = new TFitConstraintM("WMass", "MassConstraint", 0, 0, 80.4 );
	        consW->addParticles1(fitLight1,fitLight2);

	        theFitter->addConstraint(consW);
	        theFitter->setMaxNbIter(30);
	        theFitter->setMaxDeltaS(5e-5);
	        theFitter->setMaxF(1e-4);
						
		//do the fit!
	        theFitter->fit();
	        if (theFitter->getStatus() == 0) // if the fitter converged
		{
			//TMath::Prob(theFitter->getS(), theFitter->getNDF()) ); //this is the method to get the probchi2 (getS returns chi2 and getNDF returns number of degrees of freedom)
	        	afterfitLight1 = fitLight1->getCurr4Vec(); //get the 4 vector after the fit for the first ...
	        	afterfitLight2 = fitLight2->getCurr4Vec(); // ... and the second jet
			//as a cross-check: check if the LorentzVectors are indeed different after the fit!
			//cout<<"  afterfitLight1->Pt() = "<<afterfitLight1->Pt()<<endl;
			//cout<<"  afterfitLight2->Pt() = "<<afterfitLight2->Pt()<<endl;
		}
		else
		{
	        	cout << "FIT NOT CONVERGED... Should discard the kinematic fit and just calculate the top mass with the MVA (?)" << endl;
			KinFitConverged = false;
		}
	        
	        delete theFitter;
	        delete fitLight1;
	        delete fitLight2;
	        delete consW;  

	        if(KinFitConverged)
	        {	   
		    Top = (*afterfitLight1 + *afterfitLight2 + *HadBJet);
	        }
		else
		{
		    Top = (*WJet1 + *WJet2 + *HadBJet);
		}
		Mtop_kinfit_ = Top.M();		
	}	
}

void InclFourthGenSearchTools::CalculateTopMass(TLorentzVector WJet1, TLorentzVector WJet2, TLorentzVector HadBJet)
{
  TLorentzVector Top;
	Top = (WJet1 + WJet2 + HadBJet);
	Mtop_ = Top.M();	
}

void InclFourthGenSearchTools::FillMassPlots(int d, int nbOfBtags, int nbOfWs, float scaleFactor)
{
    string histoName_MTop,MShistoName_MTop;    
    histoName_MTop = "MTop_MVA_"+IntToStr(nbOfBtags)+"B_"+IntToStr(nbOfWs)+"W";
    MShistoName_MTop = "MS_"+histoName_MTop; //MSPlot for MTop
    MSPlot[MShistoName_MTop]->Fill(Mtop_, datasets_[d], true, Luminosity_*scaleFactor);
		
		string histoName_Mtop,histoMtop_dataset;
		histoName_Mtop = "Mtop_MVA_"+IntToStr(nbOfBtags)+"B_"+IntToStr(nbOfWs)+"W";
		histoMtop_dataset = histoName_Mtop+(datasets_[d]->Name()).c_str();
		histo1D[histoMtop_dataset.c_str()]-> Fill(Mtop_,datasets_[d]->NormFactor()*Luminosity_*scaleFactor);
    
    if(doKinematicFit_)
    {
      histoName_MTop = "MTop_MVAkinfit_"+IntToStr(nbOfBtags)+"B_"+IntToStr(nbOfWs)+"W";
      MShistoName_MTop = "MS_"+histoName_MTop; //MSPlot for MTop
      MSPlot[MShistoName_MTop]->Fill(Mtop_kinfit_, datasets_[d], true, Luminosity_*scaleFactor);
    }
					   	
}

void InclFourthGenSearchTools::WritePlots(TFile* fout, TDirectory* th1dir, bool savePNG, string pathPNG)
{
    cout<<" ... Writings MS plots of InclFourthGenSearchTools"<<endl;
    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
        MultiSamplePlot *temp = it->second;
        string name = it->first;
        temp->Draw(false, name, true, true, true, true, true,5,false,true,true);//(bool addRandomPseudoData, string label, bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST,int scaleNPsignal,bool addRatio, bool mergeVV, bool mergeTTV)
        temp->Write(fout, name, savePNG, pathPNG+"MSPlot/");//bool savePNG
    }
    
    //1D
    fout->cd();
    th1dir->cd();
    cout<<" ... Writing 1D plots of InclFourthGenSearchTools"<<endl;
    for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
    {
			TH1F *temp = it->second;
			temp->Write();
			TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
			if(savePNG) tempCanvas->SaveAs( (pathPNG+it->first+".pdf").c_str() );
    } 

}

void InclFourthGenSearchTools::TestPurityGoodCombinations(int d, int nbOfBtags, int nbOfWs, TRootGenEvent* genEvt, vector<TRootMCParticle*> mcParticles, vector<TRootJet*> selectedJets_MVAinput, bool TprimeEvaluation, float scaleFactor)
{     
	////////////////based on ProcessEvent in JESMeasurement/src/JetCombiner.cc
	////////////////
	///block to test 'purity' (aka 'efficiency') of good jet combinations with a 'simple' method (not MVA):
	//take the b-jet which gives the highest top mass as the hadronic top
	//NEED PROTECTION TO EXECUTE ONLY FOR TTBAR (or tprime eventually...)??
	//BTW: this entire block screws up the kinfitter result (when the kinfitter block comes before), apparently?? --> TO BE IVESTIGATED WHY... e.g. sometimes the afterfitLight1 has before this block normal values, and after this block totally crazy... 
    
    if(nbOfBtags==2 && nbOfWs==2) //not tested in 1B_2W case, modifications needed (name of the plots etc)
    {			
	////initialize stuff for each event
  	//bool all4PartonsMatched = false; // True if the 4 ttbar semi-lep partons are matched to 4 jets (not necessarily the 4 highest pt jets)
  	//bool all4JetsMatched_MCdef_ = false; // True if the 4 highest pt jets are matched to the 4 ttbar semi-lep partons
  	bool hadronictopJetsMatched_MCdef_ = false;
	pair<unsigned int, unsigned int> leptonicBJet_, hadronicBJet_, hadronicWJet1_, hadronicWJet2_; //First index is the JET number, second one is the parton  
  	leptonicBJet_ = pair<unsigned int,unsigned int>(9999,9999);
  	hadronicBJet_ = pair<unsigned int,unsigned int>(9999,9999);
  	hadronicWJet1_ = pair<unsigned int,unsigned int>(9999,9999);
  	hadronicWJet2_ = pair<unsigned int,unsigned int>(9999,9999);
	vector<TRootMCParticle*> mcParticlesMatching_; // MCParticles used for the matching
	int pdgID_top = 6;
	if(TprimeEvaluation)
	   pdgID_top = 8;
	if( (genEvt && (genEvt->isSemiLeptonic( TRootGenEvent::kMuon ) || genEvt->isSemiLeptonic( TRootGenEvent::kElec )) ) || TprimeEvaluation )
  	{
					        vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV;
						TLorentzVector topQuark, antiTopQuark;    
   					        bool muPlusFromTop = false, muMinusFromTop = false, elPlusFromTop = false, elMinusFromTop = false;
   						int nTTbarQuarks = 0;
						for(unsigned int i=0; i<mcParticles.size(); i++)
    						{
//      						cout << i << ":  status: " << mcParticles[i]->status() << "  pdgId: " << mcParticles[i]->type()<< "  motherPdgId: " << mcParticles[i]->motherType() << "  grannyPdgId: " << mcParticles[i]->grannyType() << endl;
      							if( mcParticles[i]->status() != 3) continue;

      							if( mcParticles[i]->type() == pdgID_top )
      							  topQuark = *mcParticles[i];
      							else if( mcParticles[i]->type() == -pdgID_top )
      							  antiTopQuark = *mcParticles[i];
          
      							if( mcParticles[i]->type() == 13 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -pdgID_top )
      							  muMinusFromTop = true;
      							if( mcParticles[i]->type() == -13 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == pdgID_top )
      							  muPlusFromTop = true;
      							if( mcParticles[i]->type() == 11 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -pdgID_top )
      							  elMinusFromTop = true;
      							if( mcParticles[i]->type() == -11 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == pdgID_top )
      							  elPlusFromTop = true;
      
      							if( abs(mcParticles[i]->type()) < 6 || abs(mcParticles[i]->type()) == 21 ) //light quarks, 6 should stay hardcoded
      							{
      							  mcParticlesTLV.push_back(*mcParticles[i]);
      							  mcParticlesMatching_.push_back(mcParticles[i]);
            
      							  if( fabs(mcParticles[i]->motherType()) == pdgID_top || fabs(mcParticles[i]->grannyType()) == pdgID_top )
       							  {
      							    nTTbarQuarks++;
      							  }
      							}
    						}
						// take all the selectedJets_ to study the radiation stuff, selectedJets_ are already ordened in decreasing Pt() 
						//-> WAIT, not true!!! PROBLEM???
						for(unsigned int i=0; i<selectedJets_MVAinput.size(); i++)
      							selectedJetsTLV.push_back(*selectedJets_MVAinput[i]);
      
						JetPartonMatching matching = JetPartonMatching(mcParticlesTLV, selectedJetsTLV, 2, true, true, 0.3);
						if(matching.getNumberOfAvailableCombinations() != 1)
 						     cerr << "matching.getNumberOfAvailableCombinations() = "<<matching.getNumberOfAvailableCombinations()<<"  This should be equal to 1 !!!"<<endl;
						vector< pair<unsigned int, unsigned int> > JetPartonPair; // First one is jet number, second one is mcParticle number
        
    						for(unsigned int i=0; i<mcParticlesTLV.size(); i++)
    						{
    						  int matchedJetNumber = matching.getMatchForParton(i, 0);
						  //cout<<"matchedJetNumber = "<<matchedJetNumber<<", i = "<<i<<endl;
				 		  if(matchedJetNumber != -1)
        			 		    JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
				 		}
						
						for(unsigned int i=0; i<JetPartonPair.size(); i++)
    						{
    						  unsigned int j = JetPartonPair[i].second;
    						  if( fabs(mcParticlesMatching_[j]->type()) < 6 ) //light quarks, 6 should stay hardcoded
    						  {
    						    if( ( ( muPlusFromTop || elPlusFromTop ) && mcParticlesMatching_[j]->motherType() == -24 && mcParticlesMatching_[j]->grannyType() == -pdgID_top )
   						       || ( ( muMinusFromTop || elMinusFromTop ) && mcParticlesMatching_[j]->motherType() == 24 && mcParticlesMatching_[j]->grannyType() == pdgID_top ) )
   						    {
   						       if(hadronicWJet1_.first == 9999) 
							   hadronicWJet1_ = JetPartonPair[i];
    						       else if(hadronicWJet2_.first == 9999) 
							   hadronicWJet2_ = JetPartonPair[i];
    						       else cerr<<"Found a third jet coming from a W boson which comes from a top quark..."<<endl;
  						     }
  						  }
    						  if( fabs(mcParticlesMatching_[j]->type()) == 5 )
    						  {
  						       if( ( ( muPlusFromTop || elPlusFromTop ) && mcParticlesMatching_[j]->motherType() == -pdgID_top )
    						      || ( ( muMinusFromTop || elMinusFromTop ) && mcParticlesMatching_[j]->motherType() == pdgID_top ) )
    						        hadronicBJet_ = JetPartonPair[i];
    						       else if( ( ( muPlusFromTop || elPlusFromTop ) && mcParticlesMatching_[j]->motherType() == pdgID_top )
    						      || ( ( muMinusFromTop || elMinusFromTop ) && mcParticlesMatching_[j]->motherType() == -pdgID_top ) )
    						        leptonicBJet_ = JetPartonPair[i];
    						  }					  						  
						}
						
						//if(hadronicWJet1_.first != 9999 && hadronicWJet2_.first != 9999 && hadronicBJet_.first != 9999 && leptonicBJet_.first != 9999)
    						//{ 
   						//   all4PartonsMatched = true;
    						//   if(hadronicWJet1_.first < 4 && hadronicWJet2_.first < 4 && hadronicBJet_.first < 4 && leptonicBJet_.first < 4)
    						//     all4JetsMatched_MCdef_ = true;
    						//}
					////	cout<<"hadronicWJet1_.first = "<<hadronicWJet1_.first<<", hadronicWJet2_.first = "<<hadronicWJet2_.first<<", hadronicBJet_.first = "<<hadronicBJet_.first<<endl;
    						if(hadronicWJet1_.first < 4 && hadronicWJet2_.first < 4 && hadronicBJet_.first < 4)
						{
    						  hadronictopJetsMatched_MCdef_ = true;
						  counter_hadronictopJetsMatched_MCdef++;//only the jets from the hadronic top should be matched
					////	  cout<<"   counter_hadronictopJetsMatched_MCdef = "<<counter_hadronictopJetsMatched_MCdef<<endl;
						}
    						
						//if(all4PartonsMatched)
    						//{      
   						//     if(all4JetsMatched_MCdef_)
   						//     {
   						//	     // all 4 jets found and matched, now do something with them!
						//	     //so this means: for this event, "a good jet combination exists": this is the denominator
						//	     counter_all4JetsMatched_MCdef++;
      						//     }
   						//}
						
						//in principle, with 4 input jets: 4 * 3 * 2 / 2 = 12 combinations (2 because the W-jets can be interchanged)
						TLorentzVector HadTop_choice;
						float Mtop_largestMass = -99999,Mtop_smallestMass = 99999,Mtop_largestPt = -99999,Mtop_smallestPt = 99999;
						unsigned int jetindexWJet1_largestMass = 99999,jetindexWJet2_largestMass = 99999,jetindexHadB_largestMass = 99999;
						unsigned int jetindexWJet1_smallestMass = 99999,jetindexWJet2_smallestMass = 99999,jetindexHadB_smallestMass = 99999;
						unsigned int jetindexWJet1_largestPt = 99999,jetindexWJet2_largestPt = 99999,jetindexHadB_largestPt = 99999;
						unsigned int jetindexWJet1_smallestPt = 99999,jetindexWJet2_smallestPt = 99999,jetindexHadB_smallestPt = 99999;
						float Mtop_largestMass_WJetsFixed = -99999,Mtop_smallestMass_WJetsFixed = 99999,Mtop_largestPt_WJetsFixed = -99999,Mtop_smallestPt_WJetsFixed = 99999;
						float Mtop_WJetsFixed_choice1 = -99999, Mtop_WJetsFixed_choice2 = -99999; //when WJets are fixed, only 2 choices
						float Pttop_WJetsFixed_choice1 = -99999, Pttop_WJetsFixed_choice2 = -99999; //when WJets are fixed, only 2 choices
						unsigned int jetindexWJet1_WJetsFixed = 99999,jetindexWJet2_WJetsFixed = 99999;
						unsigned int jetindexHadB_largestMass_WJetsFixed = 99999, jetindexHadB_smallestMass_WJetsFixed = 99999,jetindexHadB_largestPt_WJetsFixed = 99999,jetindexHadB_smallestPt_WJetsFixed = 99999;
						for(unsigned int i=0; i<4; i++)
						{
						   for(unsigned int j=i+1; j<4; j++)
						   {
						      for(unsigned int k=j+1; k<4; k++)
						      {
							
							   HadTop_choice = (*selectedJets_MVAinput[i] + *selectedJets_MVAinput[j] + *selectedJets_MVAinput[k]);
							   if(HadTop_choice.M() > Mtop_largestMass)
							   {
							     Mtop_largestMass = HadTop_choice.M();
							     jetindexWJet1_largestMass = i;
							     jetindexWJet2_largestMass = j;
							     jetindexHadB_largestMass = k;
							   }
							   if(HadTop_choice.M() < Mtop_smallestMass)
							   {
							     Mtop_smallestMass = HadTop_choice.M();
							     jetindexWJet1_smallestMass = i;
							     jetindexWJet2_smallestMass = j;
							     jetindexHadB_smallestMass = k;
							   }
							   if(HadTop_choice.Pt() > Mtop_largestPt)
							   {
							     Mtop_largestPt = HadTop_choice.M();
							     jetindexWJet1_largestPt = i;
							     jetindexWJet2_largestPt = j;
							     jetindexHadB_largestPt = k;
							   }
							   if(HadTop_choice.Pt() < Mtop_smallestPt)
							   {
							     Mtop_smallestPt = HadTop_choice.M();
							     jetindexWJet1_smallestPt = i;
							     jetindexWJet2_smallestPt = j;
							     jetindexHadB_smallestPt = k;
							   }
							   
							   //now for the case when you take the jets from the W counting as WJets
							   if(j==2 && k==3) // 2 and 3 are the jets from the W counting of before, only 0 and 1 (the 'b-jets') are left free to choose...
							   {
							     jetindexWJet1_WJetsFixed = 2; //j
							     jetindexWJet2_WJetsFixed = 3; //k
							     if(i==0)
							     {
							        Mtop_WJetsFixed_choice1 = HadTop_choice.M();
								Pttop_WJetsFixed_choice1 = HadTop_choice.Pt();
								
							     }
							     else if(i==1) // 2 and 3 are the jets from the W counting of before...
							     {
							        Mtop_WJetsFixed_choice2 = HadTop_choice.M();
								Pttop_WJetsFixed_choice2 = HadTop_choice.Pt();
							     }
							     
							     if(Mtop_WJetsFixed_choice1 > Mtop_WJetsFixed_choice2)
							     {
							       Mtop_largestMass_WJetsFixed = Mtop_WJetsFixed_choice1;
							       Mtop_smallestMass_WJetsFixed = Mtop_WJetsFixed_choice2;
							       jetindexHadB_largestMass_WJetsFixed = 0;
							       jetindexHadB_smallestMass_WJetsFixed = 1;
							     }
							     else
							     {
							       Mtop_largestMass_WJetsFixed = Mtop_WJetsFixed_choice2;
							       Mtop_smallestMass_WJetsFixed = Mtop_WJetsFixed_choice1;
							       jetindexHadB_largestMass_WJetsFixed = 1;
							       jetindexHadB_smallestMass_WJetsFixed = 0;
							     }
							     if(Pttop_WJetsFixed_choice1 > Pttop_WJetsFixed_choice2)
							     {
							       Mtop_largestPt_WJetsFixed = Mtop_WJetsFixed_choice1;
							       Mtop_smallestPt_WJetsFixed = Mtop_WJetsFixed_choice2;
							       jetindexHadB_largestPt_WJetsFixed = 0;
							       jetindexHadB_smallestPt_WJetsFixed = 1;
							     }
							     else
							     {
							       Mtop_largestPt_WJetsFixed = Mtop_WJetsFixed_choice2;
							       Mtop_smallestPt_WJetsFixed = Mtop_WJetsFixed_choice1;
							       jetindexHadB_largestPt_WJetsFixed = 1;
							       jetindexHadB_smallestPt_WJetsFixed = 0;
							     }
							   }
						      }
						   }
						}
						
						//only the jets from the hadronic top should be matched (minimum...)
						if(hadronictopJetsMatched_MCdef_) //necessary??
						{
      						  if ( ((jetindexWJet1_largestMass == hadronicWJet1_.first && jetindexWJet2_largestMass ==  hadronicWJet2_.first) 
      						    || (jetindexWJet1_largestMass == hadronicWJet2_.first && jetindexWJet2_largestMass == hadronicWJet1_.first)) 
      						    && jetindexHadB_largestMass == hadronicBJet_.first)
						  {
							counter_hadronictopJetsMatched_largestMass++;
							////cout<<"    counter_hadronictopJetsMatched_largestMass = "<<counter_hadronictopJetsMatched_largestMass<<endl;
						  }						  
						  if ( ((jetindexWJet1_smallestMass == hadronicWJet1_.first && jetindexWJet2_smallestMass ==  hadronicWJet2_.first) 
      						    || (jetindexWJet1_smallestMass == hadronicWJet2_.first && jetindexWJet2_smallestMass == hadronicWJet1_.first)) 
      						    && jetindexHadB_smallestMass == hadronicBJet_.first)
						  {
						  	counter_hadronictopJetsMatched_smallestMass++;
							////cout<<"    counter_hadronictopJetsMatched_smallestMass = "<<counter_hadronictopJetsMatched_smallestMass<<endl;
						  }	  
						  if ( ((jetindexWJet1_largestPt == hadronicWJet1_.first && jetindexWJet2_largestPt ==  hadronicWJet2_.first) 
      						    || (jetindexWJet1_largestPt == hadronicWJet2_.first && jetindexWJet2_largestPt == hadronicWJet1_.first)) 
      						    && jetindexHadB_largestPt == hadronicBJet_.first)
						  {
						        counter_hadronictopJetsMatched_largestPt++;
							////cout<<"    counter_hadronictopJetsMatched_largestPt = "<<counter_hadronictopJetsMatched_largestPt<<endl;
						  
						  }
						  if ( ((jetindexWJet1_smallestPt == hadronicWJet1_.first && jetindexWJet2_smallestPt ==  hadronicWJet2_.first) 
      						    || (jetindexWJet1_smallestPt == hadronicWJet2_.first && jetindexWJet2_smallestPt == hadronicWJet1_.first)) 
      						    && jetindexHadB_smallestPt == hadronicBJet_.first)
						  {
						  	counter_hadronictopJetsMatched_smallestPt++;
							////cout<<"    counter_hadronictopJetsMatched_smallestPt = "<<counter_hadronictopJetsMatched_smallestPt<<endl;
						  }
						  
						  //now for the case when you take the jets from the W counting as WJets
						  if ( ((jetindexWJet1_WJetsFixed == hadronicWJet1_.first && jetindexWJet2_WJetsFixed ==  hadronicWJet2_.first) 
      						    || (jetindexWJet1_WJetsFixed == hadronicWJet2_.first && jetindexWJet2_WJetsFixed == hadronicWJet1_.first)) 
      						    && jetindexHadB_largestMass_WJetsFixed == hadronicBJet_.first)
						  {
							counter_hadronictopJetsMatched_largestMass_WJetsFixed++;
							////cout<<"    counter_hadronictopJetsMatched_largestMass_WJetsFixed = "<<counter_hadronictopJetsMatched_largestMass_WJetsFixed<<endl;
						  }
						  ////cout<<"jetindexWJet1_WJetsFixed = "<<jetindexWJet1_WJetsFixed<<", jetindexWJet2_WJetsFixed = "<<jetindexWJet2_WJetsFixed<<", jetindexHadB_smallestMass_WJetsFixed = "<<jetindexHadB_smallestMass_WJetsFixed<<endl;			  
						  if ( ((jetindexWJet1_WJetsFixed == hadronicWJet1_.first && jetindexWJet2_WJetsFixed ==  hadronicWJet2_.first) 
      						    || (jetindexWJet1_WJetsFixed == hadronicWJet2_.first && jetindexWJet2_WJetsFixed == hadronicWJet1_.first)) 
      						    && jetindexHadB_smallestMass_WJetsFixed == hadronicBJet_.first)
						  {
						  	counter_hadronictopJetsMatched_smallestMass_WJetsFixed++;
							////cout<<"    counter_hadronictopJetsMatched_smallestMass_WJetsFixed = "<<counter_hadronictopJetsMatched_smallestMass_WJetsFixed<<endl;
						  }
						  if ( ((jetindexWJet1_WJetsFixed == hadronicWJet1_.first && jetindexWJet2_WJetsFixed ==  hadronicWJet2_.first) 
      						    || (jetindexWJet1_WJetsFixed == hadronicWJet2_.first && jetindexWJet2_WJetsFixed == hadronicWJet1_.first)) 
      						    && jetindexHadB_largestPt_WJetsFixed == hadronicBJet_.first)
						  {
						        counter_hadronictopJetsMatched_largestPt_WJetsFixed++;
							////cout<<"    counter_hadronictopJetsMatched_largestPt_WJetsFixed = "<<counter_hadronictopJetsMatched_largestPt_WJetsFixed<<endl;
						  
						  }
						  if ( ((jetindexWJet1_WJetsFixed == hadronicWJet1_.first && jetindexWJet2_WJetsFixed ==  hadronicWJet2_.first) 
      						    || (jetindexWJet1_WJetsFixed == hadronicWJet2_.first && jetindexWJet2_WJetsFixed == hadronicWJet1_.first)) 
      						    && jetindexHadB_smallestPt_WJetsFixed == hadronicBJet_.first)
						  {
						  	counter_hadronictopJetsMatched_smallestPt_WJetsFixed++;
							////cout<<"    counter_hadronictopJetsMatched_smallestPt_WJetsFixed = "<<counter_hadronictopJetsMatched_smallestPt_WJetsFixed<<endl;
						  }
						}
						//cout<<"Mtop_largestMass_WJetsFixed = "<<Mtop_largestMass_WJetsFixed<<endl;
						MSPlot["MS_MTop_largestMass_WJetsFixed_2B_2W"]->Fill(Mtop_largestMass_WJetsFixed, datasets_[d], true, Luminosity_*scaleFactor);
						MSPlot["MS_MTop_smallestMass_WJetsFixed_2B_2W"]->Fill(Mtop_smallestMass_WJetsFixed, datasets_[d], true, Luminosity_*scaleFactor);				
						MSPlot["MS_MTop_largestPt_WJetsFixed_2B_2W"]->Fill(Mtop_largestPt_WJetsFixed, datasets_[d], true, Luminosity_*scaleFactor);
						MSPlot["MS_MTop_smallestPt_WJetsFixed_2B_2W"]->Fill(Mtop_smallestPt_WJetsFixed, datasets_[d], true, Luminosity_*scaleFactor);				
							
	}//if Semi mu or semi el ttbar or Tprime...
   }
   else
     cout<<"Modify the InclFourthGenSearchTools::TestPurityGoodCombinations() function first if you want to test for other boxes than 2B_2W !!"<<endl;
	
	///end block to test 'purity' of good jet combinations with a 'simple' method (not MVA)
	////////////////
}

void InclFourthGenSearchTools::PrintPurityGoodCombinations()
{
    cout<<endl;
    
    float eff_goodjetcomb = -9999;
    cout<<"**************************  All 4 input jets free  **********************"<<endl;
    cout<<" ----- For simple jet choice largestMass -----"<<endl;
    cout<<"    number of times where a good jet combination exists = "<<counter_hadronictopJetsMatched_MCdef<<endl;
    cout<<"    number of times the good jet combination is chosen = "<<counter_hadronictopJetsMatched_largestMass<<endl;
    eff_goodjetcomb = float(counter_hadronictopJetsMatched_largestMass)/float(counter_hadronictopJetsMatched_MCdef);
    cout<<"      => 'efficiency' of taking the good jet combination = "<<eff_goodjetcomb<<endl;
    cout<<" ----- For simple jet choice smallestMass -----"<<endl;
    cout<<"    number of times where a good jet combination exists = "<<counter_hadronictopJetsMatched_MCdef<<endl;
    cout<<"    number of times the good jet combination is chosen = "<<counter_hadronictopJetsMatched_smallestMass<<endl;
    eff_goodjetcomb = float(counter_hadronictopJetsMatched_smallestMass)/float(counter_hadronictopJetsMatched_MCdef);
    cout<<"      => 'efficiency' of taking the good jet combination = "<<eff_goodjetcomb<<endl;
    cout<<" ----- For simple jet choice largestPt -----"<<endl;
    cout<<"    number of times where a good jet combination exists = "<<counter_hadronictopJetsMatched_MCdef<<endl;
    cout<<"    number of times the good jet combination is chosen = "<<counter_hadronictopJetsMatched_largestPt<<endl;
    eff_goodjetcomb = float(counter_hadronictopJetsMatched_largestPt)/float(counter_hadronictopJetsMatched_MCdef);
    cout<<"      => 'efficiency' of taking the good jet combination = "<<eff_goodjetcomb<<endl;
    cout<<" ----- For simple jet choice smallestPt -----"<<endl;
    cout<<"    number of times where a good jet combination exists = "<<counter_hadronictopJetsMatched_MCdef<<endl;
    cout<<"    number of times the good jet combination is chosen = "<<counter_hadronictopJetsMatched_smallestPt<<endl;
    eff_goodjetcomb = float(counter_hadronictopJetsMatched_smallestPt)/float(counter_hadronictopJetsMatched_MCdef);
    cout<<"      => 'efficiency' of taking the good jet combination = "<<eff_goodjetcomb<<endl;
    
    cout<<endl;
    
    cout<<"**************************  Only 2 input jets free (the 'b-jets') **********************"<<endl;
    cout<<" ----- For simple jet choice largestMass -----"<<endl;
    cout<<"    number of times where a good jet combination exists = "<<counter_hadronictopJetsMatched_MCdef<<endl;
    cout<<"    number of times the good jet combination is chosen = "<<counter_hadronictopJetsMatched_largestMass_WJetsFixed<<endl;
    eff_goodjetcomb = float(counter_hadronictopJetsMatched_largestMass_WJetsFixed)/float(counter_hadronictopJetsMatched_MCdef);
    cout<<"      => 'efficiency' of taking the good jet combination = "<<eff_goodjetcomb<<endl;
    cout<<" ----- For simple jet choice smallestMass -----"<<endl;
    cout<<"    number of times where a good jet combination exists = "<<counter_hadronictopJetsMatched_MCdef<<endl;
    cout<<"    number of times the good jet combination is chosen = "<<counter_hadronictopJetsMatched_smallestMass_WJetsFixed<<endl;
    eff_goodjetcomb = float(counter_hadronictopJetsMatched_smallestMass_WJetsFixed)/float(counter_hadronictopJetsMatched_MCdef);
    cout<<"      => 'efficiency' of taking the good jet combination = "<<eff_goodjetcomb<<endl;
    cout<<" ----- For simple jet choice largestPt -----"<<endl;
    cout<<"    number of times where a good jet combination exists = "<<counter_hadronictopJetsMatched_MCdef<<endl;
    cout<<"    number of times the good jet combination is chosen = "<<counter_hadronictopJetsMatched_largestPt_WJetsFixed<<endl;
    eff_goodjetcomb = float(counter_hadronictopJetsMatched_largestPt_WJetsFixed)/float(counter_hadronictopJetsMatched_MCdef);
    cout<<"      => 'efficiency' of taking the good jet combination = "<<eff_goodjetcomb<<endl;
    cout<<" ----- For simple jet choice smallestPt -----"<<endl;
    cout<<"    number of times where a good jet combination exists = "<<counter_hadronictopJetsMatched_MCdef<<endl;
    cout<<"    number of times the good jet combination is chosen = "<<counter_hadronictopJetsMatched_smallestPt_WJetsFixed<<endl;
    eff_goodjetcomb = float(counter_hadronictopJetsMatched_smallestPt_WJetsFixed)/float(counter_hadronictopJetsMatched_MCdef);
    cout<<"      => 'efficiency' of taking the good jet combination = "<<eff_goodjetcomb<<endl;
}
