
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
#include "TBranch.h"
#include "TList.h"
//user code
#include "TopTreeProducer/interface/TRootGenEvent.h"
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeProducer/interface/TRootJet.h"

#include "../config/Datasets.cc"
#include "../interface/SelectionTable.h"
#include "../interface/WJetEstimation.h"
#include "../interface/WJetEstPseudoExp.h"
#include "../interface/TtJetEstimation.h"
#include "../interface/ABCDEstimation.h"
#include "../interface/QCDShapeEstimation.h"
#include "../interface/IterativeShapeExtrapolator.h"
#include "../interface/PlottingTools.h"
#include "../interface/BkgEstimationSummary.h"
#include "../interface/MultiSamplePlot.h"
#include "../interface/MultiCutPlot.h"
#include "../interface/CutImpactEvaluation.h"
//#include "../interface/EstimationPerJets_BJetsMultiplicity.h"
#include "../interface/MEzCalculator.h"
#include "../interface/Observables.h"
#include "../interface/PlotObservables.h"
#include "../interface/TTreeObservables.h"
#include "../interface/ObservablesRanker.h"



//#include "Style.C"

using namespace std;
using namespace TopTree;



void Scale(vector<Dataset> datasets, double*& numbers, float Luminosity ){
        for(unsigned int i=0;i<datasets.size();i++){
	     numbers[i] = numbers[i]*datasets[i].NormFactor()*Luminosity;
	}
}
				

//main function

 

int main(int argc, char*argv[]){

  /////////////////////
  // Configuration
  /////////////////////

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
 
    
TH1F *h1= new TH1F("h1", "My First ", 100, 0, 1000); 


 TH1F *hPtJets[15][datasets.size()] ;

 char hPt[150]; 


 for (int j=0;j<datasets.size();j++){

   
   string ds=datasets[j].Name();
   //f1->mkdir(ds);
  

   for (int i=0; i<15; i++) {
   sprintf(hPt,"PtJets_%d for dataset %s",i,ds.c_str());

     hPtJets[i][j] = new TH1F(hPt, "Pt of Jets",100,0,600);}

     }

  cout<<"**********************************************************************"<<endl;
  cout<<"Begining of the program for bkg estimation in lepton+jets channels !" <<endl;
  cout<<"**********************************************************************"<<endl;

  //SetStyle if needed
  //setTDRStyle(); 
  //setMyStyle();


 
  
  if(verbose>0) cout<<" - Variable decleration ..."<<endl;

  //vector of objects
  vector<TRootMuon> init_muons;
  vector<TRootElectron> init_electrons;
  vector<TRootJet> init_jets;
  vector<TRootMET> mets;

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
        Observables fake;
	TTreeObservables ttre;
	//	ObservablesRanker ObsRank;
	PlotObservables plot(fake.ListOfVariables(), fake.RangeVariables()); 
	TTreeObservables TtreeObs(fake.ListOfVariables());
        //ObservablesRanker obsr;
	// ObservablesRanker (fake);
        //ObservablesRanker ObsRank(fake.ListOfVa;
  	if(verbose>1) cout<<"	Loop over events "<<endl;
	for(unsigned int ievt=0;ievt<500000;ievt++){//eventTree->GetEntries();ievt++){
		TLorentzVector JESbalance(0,0,0,0);
  		
		//if(verbose>2) cout<<"Event "<<ievt<<endl;
	        eventTree->GetEvent(ievt);
		
		// clear vectors	
		init_jets.clear();
		init_muons.clear();
		init_electrons.clear();
		mets.clear();
	  	
		TRootGenEvent* genEvt;
		
		//fill vectors
		for(int i=0;i<tcjets->GetEntriesFast();i++) init_jets.push_back(* (TRootJet*) tcjets->At(i));
		for(int i=0;i<tcmuons->GetEntriesFast();i++) init_muons.push_back(* (TRootMuon*) tcmuons->At(i));
		for(int i=0;i<tcelectrons->GetEntriesFast();i++) init_electrons.push_back(* (TRootElectron*) tcelectrons->At(i));
		for(int i=0;i<tcmets->GetEntriesFast();i++) mets.push_back(* (TRootMET*) tcmets->At(i));
  	        if(tcgenEvt->GetEntriesFast()==1) genEvt = (TRootGenEvent*) tcgenEvt->At(0);
	
		if(verbose>3) cout<<"Nof muons "<<init_muons.size()<<endl;
		if(verbose>3) cout<<"Nof electrons "<<init_electrons.size()<<endl;
		if(verbose>3) cout<<"Nof jets "<<init_jets.size()<<endl;
		
		/////////////////////////////
		//   Selection
		/////////////////////////////
		//Declare selection instance	
		Selection selection(init_jets, init_muons, init_electrons, mets);

             

		
		vector<TRootJet> selectedJetss = selection.GetSelectedJets(JetsPtCutSR,JetsEtaCutSR, JetshcalFraction);
	       	vector<TRootMuon> selectedMuons = selection.GetSelectedMuons(MuonPtCutSR,MuonEtaCutSR,MuonRelIsoCutSR);
	
		vector<TRootMET> selectedMET = mets;
                
		vector<TRootJet> jets = init_jets;

		//for (int i=0; i<selectedJetss.size(); i++) {
    
   //  if (datasets[d].Name()=="SUSY" || datasets[d].Name()=="QCD" )
		   // hPtJets[i][d]->Fill(selectedJetss[i].Pt());
		   //cout<<datasets[d].Name()<<"   "<<selectedJetss[i].Pt()<<endl;
		//}
 

 if (selectedMuons.size()>0){
   //  for (int i=0; i<selectedMuons.size(); i++) {
   //cout<<datasets[d].Name()<<"   "<<selectedMuons[i].Pt()<<endl;

   //}

 Observables a(selectedMuons[0],selectedJetss,mets[0],ds);

      
	     
 //cout << a.Variable(string("MET")) << " MET  "<<a.Variable(string("ET3"))<<" ET3 "<<" inside Obs  "<<a.ListOfVariables().size()<<endl;
 
    
        plot.Fill(a);
 
        TtreeObs.Fill(a);


	
 }

     }////end of loop in events

		plot.Write(ds,true);
		string ds1="TTJets"; string ds2="SUSY";
	      string dsroot=ds+".root";


	      TtreeObs.Write(ds,true);

	       //ObservablesRanker obsr;obsr.
		cout<<" NOW I AM IN THE DS  "<<ds.c_str()<<endl;
	
	    
	}//loop on datasets
 
 
 string fS, fBkg,ff1,ff2,ff3;

            //fS = "SUSY";
             fBkg = "TTJets_tree.root";
	     fS = "TTJets_tree.root";
             ff1="Merged_Plots";ff2="SUSY";ff3="TTJets";

   TFile * hfileS = new TFile("TTJets_tree.root");
   //TFile * hfileB = new TFile("SUSY_tree.root");
   TFile * hfileB = new TFile("TTJets_tree.root");
   bool SameFiles = false;
    if (fS == fBkg) (SameFiles=1);


   TTree *tSignal = (TTree*)hfileS->Get("OBS");
   TTree *tBkg = (TTree*)hfileB->Get("OBS");

  


  ObservablesRanker ObsR(tSignal,tBkg,ff1,ff2,ff3,SameFiles);
  
  //Once everything is filled ...
  if(verbose>0) cout<<" We ran over all the data ;-)"<<endl;
 
  cout<<"**********************************************************************"<<endl;
  cout<<"           End of the program !!" <<endl;
  cout<<" 		doesn't crashed yet ;-) "<<endl;
  cout<<"**********************************************************************"<<endl;

}
