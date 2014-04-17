#include "iomanip.h"
#include <cmath>
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TBranch.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TList.h"
#include <stdio.h>
//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeProducer/interface/TRootGenEvent.h"

#include "../config/Datasets.cc"
#include "../Selection/interface/SelectionTable.h"
//#include "../BkgEstimationMethods/interface/WJetEstimation.h"
//#include "../BkgEstimationMethods/interface/WJetEstPseudoExp.h"
#include "../BkgEstimationMethods/interface/TtJetEstimation.h"
#include "../BkgEstimationMethods/interface/ABCDEstimation.h"
#include "../BkgEstimationMethods/interface/QCDShapeEstimation.h"
//#include "../BkgEstimationMethods/interface/WJetShapeEstimation.h"
#include "../BkgEstimationMethods/interface/BkgEstimationSummary.h"
#include "../Tools/interface/PlottingTools.h"
#include "../Tools/interface/MultiSamplePlot.h"
#include "../Tools/interface/MultiCutPlot.h"
#include "../Tools/interface/CutImpactEvaluation.h"
#include "../Tools/interface/TemplateComparator.h"
#include "../Reconstruction/interface/MEzCalculator.h"
#include "../Reconstruction/interface/Observables.h"
#include "../Reconstruction/interface/TtSemiMuJetCombination.h"
#include "../Reconstruction/interface/PlotObservables.h"
#include "../Reconstruction/interface/TTreeObservables.h"
#include "../Reconstruction/interface/ObservablesRanker.h"

#include "../Reconstruction/interface/Combination.h"

#include "Style.C"

using namespace std;
using namespace TopTree;



void Scale(vector<Dataset> datasets, double*& numbers, float Luminosity ){
        for(unsigned int i=0;i<datasets.size();i++){
	     numbers[i] = numbers[i]*datasets[i].NormFactor()*Luminosity;
	}
}




void    JetID(std::vector<TRootJet>&, Int_t &, Float_t &, Float_t &, Bool_t &);
void    JetSelection(std::vector<TRootJet>&, Float_t &, Float_t &, Bool_t &);
std::vector<TRootJet*> BJetSelection(std::vector<TRootJet*> &, Int_t &, Float_t &, Bool_t &);


Float_t ChiSquareMatching( TRootMuon &, std::vector<TRootJet> &, UInt_t &, UInt_t *, Float_t &, Float_t &, Float_t &, Float_t &, Float_t &, Float_t &, Bool_t &);
			

//main function




int main(int argc, char*argv[]){

  /////////////////////
  // Configuration
  /////////////////////
  Bool_t DEBUG = (argc>2 ?(atoi(argv[2]) == 1) : 0);
  int verbose = 3;
  float Luminosity = 100; // in 1/fb
  cout<<"The results will be obtained for a luminosity of "<<Luminosity<<" x 1/pb"<<endl;
  
  ////////////////////////
  //Define list of variable
  /////////////////////////



  /////////////////////
 
  // Load Datasets
  /////////////////////
  if(verbose>0) cout<<" - Load datasets ..."<<endl;
  vector<Dataset> datasets;
  LoadDatasets(datasets);
  /////////////////////
 	Observables obss;


    TList *FileList;
    //TFile *f1 = new TFile("plots.root", "recreate");
 
  

 for (int j=0;j<datasets.size();j++){

   
   string ds=datasets[j].Name();
 
  
   cout <<ds<<endl;

  
     }

  cout<<"**********************************************************************"<<endl;
  cout<<"Begining of the program for bkg estimation in lepton+jets channels !" <<endl;
  cout<<"**********************************************************************"<<endl;

  //SetStyle if needed
  //setTDRStyle(); 
  //setMyStyle();


 
  
  if(verbose>0) cout<<" - Variable declaration ..."<<endl;

  //vector of objects
  vector<TRootMuon> init_muons;
  vector<TRootElectron> init_electrons;
  vector<TRootJet> init_jets;
  vector<TRootMET> mets;
	TRootJet qFromW1_;  
	TRootJet qFromW2_;  
	TRootJet bFromHadTop_;  
	TRootJet bFromLepTop_;
  const  TRootMET cmets;
  TBranch *jets_br = 0;
  TClonesArray *tcjets;
  TBranch *muons_br = 0;
  TClonesArray *tcmuons;
  TBranch *electrons_br = 0;
  TClonesArray *tcelectrons;
  TBranch *mets_br = 0;
  TClonesArray *tcmets;
  TBranch *genEvt_br = 0;
  TClonesArray *tcgenEvt;
  TBranch* run_br = 0;
  TRootRun* runInfos = 0;
  TBranch* event_br = 0;
  TRootEvent* event = 0;

 ////////////////////////////////////
 /// Conditions: b-tagging 
 ////////////////////////////////////

 

 ////////////////////////////////////
 /// Selection 
 ////////////////////////////////////
  
   float ElectronPtCut = 10;
   float ElectronEtaCut = 2.0;
   float ElectronRelIsoCut = 0.1;
   //SR
   float MuonPtCutSR = 30;
   float MuonEtaCutSR = 2.1;
   float MuonRelIsoCutSR = 0.1;
   float JetsPtCutSR = 30;
   float JetsEtaCutSR = 2.4;
   //CR: ttjetEstimation
   float MuonPtCutCR = 30;
   float MuonEtaCutCR = 2.1;
   float MuonRelIsoCutCR = 0.1;
   float JetsPtCutCR = 30;
   float JetsEtaCutCR = 2.4;
   float JetshcalFraction = 0.1;


  // Jet Chi2 combination
   UInt_t Permutation[4] = {0,0,0,0};
   Float_t WMass   = 90.0; // Derived from the SanityChecker (Truth module)
   Float_t TopMass = 174.5;
   Float_t blMass  = 91.1;
   //Float_t LepTopMassResol = 15.6;
   Float_t LepblMassResol  = 36.0;
   Float_t HadTopMassResol = 21.1;
   Float_t HadWMassResol   = 12.7;
   UInt_t  NjetsForComb    = 7;
   Float_t ChiSq           = 0;


	Float_t ElectronPtThreshold  = 10;
	Float_t ElectronEtaThreshold = 2.1;
	Float_t MuonPtThreshold      = 20;
	Float_t MuonEtaThreshold     = 2.1;
	Float_t RelCaloIso  = 9999; // Relative Calorimeter isolation
	Float_t RelTrackIso = 9999; // Relative Tracker isolation
	Float_t RelIso      = 0.1;  // Relative Calorimeter and Tracker isolation
	Bool_t  VetoIso     = true; // Make use of the veto iso cone (VetoHad<6 and VetoEm<4) otherwise use the cuts below :
	Float_t VetoEm      = 9999; // Electromagnetic energy in 0.07 Cone
	Float_t VetoHad     = 9999; // Hadronic energy in 0.1 Cone
	Float_t DRMuJets    = -1;   // (Eta,Phi) DeltaR between the muon and its closest jet
	Float_t d0Sig       = 3;    // IP significance
	
	Float_t d0value = 0;
	Float_t relIsovalue = 0;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Jet ID requirements
	Int_t   NbOfConsts  = 1;  // Min Number of jet constituents
	Float_t EcalNRJFrac = 0.01 ; // Lowest jet energy fraction of ecal energy (-1 : no cut, usually 0.05)
	Float_t HcalNRJFrac = -1; // Lowest jet energy fraction of hcal energy (-1 : no cut, usually 0.05)
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Jet selection
	Float_t JetPtThreshold = 30;    // Jet pt threshold [GeV/c]
	Float_t JetEtaThreshold = 2.4;  // Jet eta threshold [GeV/c]
	UInt_t Njets  = 4; // Min number of jets with pt>JetPtThreshold
	Bool_t ExactNbOfJets  = false;
	UInt_t Nbjets = 0; // Nin number of b-jets with pt>JetPtThreshold
	Bool_t ExactNbOfBJets = false;
	const int NbOfJetBins = 3; // 0 : 4 jets, 1 : 5 jets and 2 : >=6jets
	// B-tag discriminator :
	TString BtagAlgoNames[6] = {
	"TrackCountingHighEff",  // - 1 == TrackCountingHighEff    : loose/medium/tight = 2.030/4.380/14.200
	"TrackCountingHighPur",  // - 2 == TrackCountingHighPur    : loose/medium/tight = 1.470/2.360/ 5.360
	"JetProbability",        // - 3 == JetProbability          : loose/medium/tight = 0.241/0.490/ 0.795
	"JetBProbability",       // - 4 == JetBProbability         : loose/medium/tight = 1.100/2.370/ 3.390
	"simpleSecondaryVertex", // - 5 == simpleSecondaryVertex   : loose/medium/tight = 1.250/2.050/ 4.070
	"combinedSecondaryVertex"// - 6 == combinedSecondaryVertex : loose/medium/tight = 0.387/0.838/ 0.940
	};

  //////////////////////////////////////////////
  //loop on datasets
  if(verbose>0) cout<<" - Loop over datasets ... "<<datasets.size()<<" datasets !"<<endl;
  for(unsigned int d=0;d<datasets.size();d++){
  	if(verbose>1) cout<<"   Dataset "<<d<<": "<<datasets[d].Name()<<endl;
  	//open files and load
	TChain* eventTree = datasets[d].eventTree();
	TChain* runTree = datasets[d].runTree();


	string ds = datasets[d].Name();
  	///////////////////////////////////////
  	//Branches and TCLonesArray
  	///////////////////////////////////////
  	//jets
  	jets_br = (TBranch *) eventTree->GetBranch ("Jets");
  	tcjets = new TClonesArray ("TopTree::TRootJet", 0);
  	jets_br->SetAddress (&tcjets);
  	//muons
  	muons_br = (TBranch *) eventTree->GetBranch ("Muons");
  	tcmuons = new TClonesArray ("TopTree::TRootMuon", 0);
  	muons_br->SetAddress (&tcmuons);
  	//electrons
  	electrons_br = (TBranch *) eventTree->GetBranch ("Electrons");
  	tcelectrons = new TClonesArray ("TopTree::TRootElectron", 0);
  	electrons_br->SetAddress (&tcelectrons);
  	//mets
  	mets_br = (TBranch *) eventTree->GetBranch ("MET");
  	tcmets = new TClonesArray ("TopTree::TRootMET", 0);
  	mets_br->SetAddress (&tcmets);
	//GenEvent
	genEvt_br = (TBranch *) eventTree->GetBranch ("GenEvent");
  	tcgenEvt = new TClonesArray ("TopTree::TRootGenEvent", 0);
  	if(genEvt_br) genEvt_br->SetAddress (&tcgenEvt);
	//HLT info
	//run_br = (TBranch *) runTree->GetBranch("runInfos");
	//run_br->SetAddress(&runInfos);
	//event_br = (TBranch *) eventTree->GetBranch("Event");
	//event_br->SetAddress(&event);
	////////////////////////////////////////
	


	//loop on events
        Observables MyObss;
	Observables a;

	//	ObservablesRanker ObsRank;

	PlotObservables Plots(a.ListOfVariables(), a.RangeVariables()); 
	TTreeObservables TtreeObs(a.ListOfVariables());



        //ObservablesRanker obsr;
	// ObservablesRanker (fake);
        //ObservablesRanker ObsRank(fake.ListOfVa;
  	if(verbose>1) cout<<"	Loop over events "<<endl;
	std::vector<TRootJet>        SelectedJets;
	std::vector<TRootJet>        selectedJetss;	
        std::vector<TRootMuon>        selectedMuons;
	//	std::vector<TRootMET*>        SelectedMets;
	std::vector<TRootJet>        CombinedJets;
	std::vector<TRootJet>        ChiJets;
	int nevt = eventTree->GetEntries();
	nevt = 25000;
	for(unsigned int ievt=0;ievt<nevt;ievt++){//eventTree->GetEntries();ievt++){
	  // Observables a;
	 

	  if (ievt%10000 == 0)cout<< " NOW WORKING ON "<<ievt<<flush<<"\r";//endl;
	  if (ievt == nevt-1)cout<< " Finishing in this dataset "<<nevt<<" events"<<endl;;//endl;
		TLorentzVector JESbalance(0,0,0,0);
  		
		//if(verbose>2) cout<<"Event "<<ievt<<endl;
	        eventTree->GetEvent(ievt);
		
		// clear vectors	
		init_jets.clear();
		init_muons.clear();
		init_electrons.clear();
		mets.clear();
	  	
		TRootGenEvent* genEvt;
	       
	//std::vector<TRootElectron*>   SelectedElectrons;

	//
	SelectedJets.clear();	
        selectedJetss.clear();;

      

        CombinedJets.clear();
        ChiJets.clear();
        selectedMuons.clear();
        SelectedJets.clear();
	// qFromW1_  = TRootJet ();  qFromW2_  = TRootJet ();   bFromHadTop_  = TRootJet ();  bFromLepTop_  = TRootJet ();
		//fill vectors
		for(int i=0;i<tcjets->GetEntriesFast();i++) init_jets.push_back(* (TRootJet*) tcjets->At(i));
		 for(Int_t i = 0; i <tcjets->GetEntriesFast(); i++) SelectedJets.push_back(*(TRootJet*) tcjets->At(i));
		for(int i=0;i<tcmuons->GetEntriesFast();i++) init_muons.push_back(* (TRootMuon*) tcmuons->At(i));

	
		for(int i=0;i<tcelectrons->GetEntriesFast();i++) init_electrons.push_back(* (TRootElectron*) tcelectrons->At(i));
		for(int i=0;i<tcmets->GetEntriesFast();i++) mets.push_back(* (TRootMET*) tcmets->At(i));
		//	 for(Int_t i = 0; i <tcmets->GetEntriesFast(); i++) SelectedMets.push_back((TRootMET*) tcmets->At(i));
  	        if(tcgenEvt->GetEntriesFast()==1) genEvt = (TRootGenEvent*) tcgenEvt->At(0);

		  
	
		if(verbose>3) cout<<"Nof muons "<<init_muons.size()<<endl;
		if(verbose>3) cout<<"Nof electrons "<<init_electrons.size()<<endl;
		if(verbose>3) cout<<"Nof jets "<<init_jets.size()<<endl;
		
		/////////////////////////////
		//   Selection
		/////////////////////////////
		//Declare selection instance	
		 Selection selection(init_jets, init_muons, init_electrons, mets);
               
	
		 //SelectedJets = selection.GetSelectedJets(JetsPtCutSR,JetsEtaCutSR, JetshcalFraction);
		
		selectedJetss = selection.GetSelectedJets(JetsPtCutSR,JetsEtaCutSR, JetshcalFraction);


	   	 ChiJets = selection.GetSelectedJets(JetsPtCutSR,JetsEtaCutSR, JetshcalFraction);
	       	 selectedMuons = selection.GetSelectedMuons(MuonPtCutSR,MuonEtaCutSR,MuonRelIsoCutSR);
		 // TtSemiMuJetEstimation SemiSelection(init_jets , selectedMuons[0] , NjetsForComb , JetPtThreshold, JetEtaThreshold);

		//vector<TRootMuon*> SelectedMuons = selection.GetSelectedMuons(MuonPtCutSR,MuonEtaCutSR,MuonRelIsoCutSR);
		vector<TRootMET> selectedMET = mets;
                
	       
	


	// Jets :	

		JetID(SelectedJets, NbOfConsts, EcalNRJFrac, HcalNRJFrac, DEBUG);
		JetSelection(SelectedJets, JetPtThreshold, JetEtaThreshold, DEBUG);

		

		
	///////////////////////////////////////////////////////////////////
	// Jet pairing 
	///////////////////////////////////////////////////////////////////
		
	
 if ( selectedMuons.size()>0 && SelectedJets.size()>3){
	     
   //	cout <<"  Muon "<<selectedMuons.size()<<"  "<<SelectedJets.size()<<"  "<<selectedJetss.size()<<"  "<<ievt<<endl;
 
   	ChiSq = ChiSquareMatching(selectedMuons[0], SelectedJets, NjetsForComb, Permutation, WMass, TopMass, blMass, HadWMassResol, HadTopMassResol, LepblMassResol, DEBUG);
	CombinedJets.clear();
        
	for(Int_t i = 0; i < 4; i++) {CombinedJets.push_back(SelectedJets[Permutation[i]]);}
	//CombinedJets.push_back(ChiSq);	
    

 
   // cout << selectedMuons.size()<< "   "<<selectedJetss.size()<< "  "<<ChiJets.size()<<endl;
   //DEBUG = true;

		//MyObs  = new Observables(); 
		//Observables a(selectedMuons[0],CombinedJets,SelectedJets,mets[0],ds); 

      
		//Observables a(selectedMuons[0],selectedJetss,mets[0],ds); 

		//	Observables * a = new Observables();
              
		 Observables a(selectedMuons[0],CombinedJets,SelectedJets,mets[0],ds , ChiSq); 

	

		 //a.ComputeVariables(selectedMuons[0], CombinedJets[0], CombinedJets[1], CombinedJets[2], CombinedJets[3],SelectedJets, mets[0],ds , ChiSq);
	
                
 
		//cout << " Combined Jets" << CombinedJets.size()<<"  "<<SelectedJets.size()<<"  "<<selectedJetss.size()<<" blepquark  "<<CombinedJets[0].Pt()<<endl;

		 /*	 if ( a.Variable(string("TransMassLepTop")) <0 ||  a.Variable(string("TransMassTtbar"))<0 ){

		   cout << a.Variable(string("TransMassLepTop")) << " TransMassLepTop  "<<a.Variable(string("TransMassTtbar"))<<" TransMassTtbar "<<" inside Obs  "<<"  transMassTtbar "<<a.Variable(string("TransMassTtbar2"))<<a.ListOfVariables().size()<<endl;
		   }*/
		
			Plots.Fill(a);



			TtreeObs.Fill(a);

 }


     }////end of loop in events



		Plots.Write(ds,true);
		cout<<"  ended  plots loop in events for "<<ds<<endl;
		string ds1="TTJets"; 
                string ds2="SUSY";

	      string dsroot=ds+".root";


	      TtreeObs.Write(ds,true);
   cout<<"  ended write loop in events "<<endl;
	      //ObservablesRanker obsr;obsr.
	   
	



   //delete eventTree;
   //	delete runTree;
	

	}//loop on datasets
 

  
 string fS, fBkg,ff1,ff2,ff3;

            //fS = "SUSY";
             fBkg = datasets[0].Name()+"_tree.root"; 
	     fS = datasets[0].Name() +"_tree.root";
             ff3=datasets[0].Name() ;//should be the TTjet
	     
	     char name[100];   char name2[100];
          
	     sprintf (name, "%s_tree.root",datasets[0].Name().c_str() ); 

	     ff1="Merged_Plots_"+datasets[0].Name()   ;


	     TFile *hfileS, *hfileB; TTree *tSignal, *tBkg;
 for (unsigned int j=1;j<datasets.size();j++){

   cout<<" Will merge and compute overlap for "<< datasets[0].Name()<<"  and  "<<datasets[j].Name()<<endl;

        
          
              ff2=datasets[j].Name();
	      
	 
	      sprintf (name2, "%s_tree.root",datasets[j].Name().c_str() ); 
  
	      hfileS = new TFile(name);
	      hfileB = new TFile(name2);
  

	      cout<<name<<"  xixixi "<<name2<<endl;
   bool SameFiles = false;
    if (fS == fBkg) (SameFiles=1);


   tSignal = (TTree*)hfileS->Get("OBS");
   tBkg = (TTree*)hfileB->Get("OBS");

  

    ObservablesRanker ObsR(tSignal,tBkg,ff1,ff2,ff3,SameFiles);

    // delete hfileS;delete hfileB;delete tSignal;delete tBkg;
    
    }
 //delete hfileS;delete hfileB; delete tSignal; delete tBkg;
  //Once everything is filled ...
  if(verbose>0) cout<<" We ran over all the data ;-)"<<endl;
 
  cout<<"**********************************************************************"<<endl;
  cout<<"           End of the program !!" <<endl;
  cout<<" 		doesn't crashed yet ;-) "<<endl;
  cout<<"**********************************************************************"<<endl;
  return 0;

 
}



Float_t ChiSquareMatching(TRootMuon &muon, std::vector<TRootJet> &jets, UInt_t &n, UInt_t *Permutation, Float_t &WMass_, Float_t &TopMass_, Float_t &blMass_, Float_t &HadWMassResol_, Float_t &HadTopMassResol_, Float_t &LepblMassResol_, Bool_t &DEBUG)
{
	
	if(jets.size()<4) return -1;
	UInt_t NbOfJets = 0;
	(jets.size() < n ? NbOfJets = jets.size() : NbOfJets = n );
	UInt_t *numbers = new UInt_t[NbOfJets];
   	for(UInt_t i=0;i<NbOfJets;i++) numbers[i]=i;
   	UInt_t *comb = new UInt_t[4];
   	for(UInt_t i=0;i<4;i++) comb[i]=i;
	UInt_t NbOfPermutations = 0;
	Float_t ChiSquare    =  0, ChiSquare_tmp   =  0;
	Float_t HadWChiSquare=  0, HadTopChiSquare =  0, LepblChiSquare =  0;
	Float_t HadWCandMass_tmp = -1, HadTopCandMass_tmp  = -1, LepBlCandMass_tmp  = -1;
	Float_t HadWCandMass = -1, HadTopCandMass  = -1, LepBlCandMass  = -1;
	//DEBUG = true;
	ChiSquare = 9999999;
	do
	{
		if(DEBUG)
		{
			std::cout<<"-- Enter loop for jet combination "<<std::endl;
			std::cout<<"-- Considered Jet combination : "<<comb[0]<<"/"<<comb[1]<<"/"<<comb[2]<<"/"<<comb[3]<<std::endl;
		}
	   	NbOfPermutations = 0;
		do
		{
			NbOfPermutations++;
			if(DEBUG) std::cout<<"--- Current permutations : "<<comb[0]<<"/"<<comb[1]<<"/"<<comb[2]<<"/"<<comb[3]<<std::endl;
			HadWCandMass_tmp   = (jets[comb[0]]+jets[comb[1]]).M();
			HadTopCandMass_tmp = (jets[comb[0]]+jets[comb[1]]+jets[comb[2]]).M();
			LepBlCandMass_tmp  = (jets[comb[3]]+muon).M();

			HadWChiSquare   = pow((HadWCandMass_tmp-WMass_),2)/pow(HadWMassResol_,2);
			HadTopChiSquare = pow((HadTopCandMass_tmp-TopMass_),2)/pow(HadTopMassResol_,2);
			LepblChiSquare  = pow((LepBlCandMass_tmp-blMass_),2)/pow(LepblMassResol_,2);
	
			ChiSquare_tmp   = HadWChiSquare+HadTopChiSquare+LepblChiSquare;
			if(ChiSquare_tmp<ChiSquare)
			{
				//for(unsigned int j=0;j<4;j++) Permutation_tmp[i][j] = comb[j];
				//ChiSquare[i] = ChiSquare_tmp[i];
				for(UInt_t i=0;i<4;i++) Permutation[i] = comb[i];
				ChiSquare = ChiSquare_tmp;
				HadWCandMass   = HadWCandMass_tmp;
				HadTopCandMass = HadTopCandMass_tmp;
				LepBlCandMass  = LepBlCandMass_tmp;
			}
		}
		while(next_permutation(comb,comb+4));
	   	if(DEBUG) 
		{
			std::cout<<"--- ChiSquare matching : Nb of permutations found = "<<NbOfPermutations<<std::endl;
			std::cout<<"--- ChiSquare matching : Min ChiSquare found = "<<ChiSquare<<std::endl;
			std::cout<<"--- ChiSquare matching : associated to the permutation : "<<Permutation[0]<<"/"<<Permutation[1]<<"/"<<Permutation[2]<<"/"<<Permutation[3]<<std::endl;
		}
	}
	while(stdcomb::next_combination(numbers,numbers+NbOfJets,comb,comb+4));
	if(DEBUG)
	{
		std::cout<<"ChiSquare = "<<ChiSquare<<std::endl;
		std::cout<<"Associated had W mass   = "<<HadWCandMass<<std::endl;
		std::cout<<"Associated had top mass = "<<HadTopCandMass<<std::endl;
		std::cout<<"Associated lep b+l mass = "<<LepBlCandMass<<std::endl;
	}
	return ChiSquare;delete comb;delete numbers;
}

       
  




void JetID(std::vector<TRootJet>& jets, Int_t &NbOfConst, Float_t &EcalEnergyFrac, Float_t &HcalEnergyFrac, Bool_t &DEBUG)
{
	// This function removes the "bad ID" jets. 
	if(DEBUG) std::cout<<"JetID : Initial/Final nb of jets : "<<jets.size();
	for(Int_t i = (jets.size()-1); i>=0; i--)
	{
		if(jets[i].nConstituents()      <= NbOfConst)      {jets.erase(jets.begin()+i); continue;}
		if(jets[i].ecalEnergyFraction() <= EcalEnergyFrac) {jets.erase(jets.begin()+i); continue;}
		if(jets[i].hcalEnergyFraction() <= HcalEnergyFrac) {jets.erase(jets.begin()+i); continue;}
	}
	if(DEBUG) std::cout<<"/"<<jets.size()<<std::endl;
}

void JetSelection(std::vector<TRootJet> &jets, Float_t &ptThreshold, Float_t &etaThreshold, Bool_t &DEBUG)
{
	// This function removes jets not passing the selection cuts
	if(DEBUG) std::cout<<"JetSelection : Initial/Final nb of jets : "<<jets.size();
	for(Int_t i = (jets.size()-1); i>=0; i--)
	{
		if(jets[i].Pt()<ptThreshold || fabs(jets[i].Eta())>etaThreshold) {jets.erase(jets.begin()+i); continue;}
	}
	if(DEBUG) std::cout<<"/"<<jets.size()<<std::endl;
}
