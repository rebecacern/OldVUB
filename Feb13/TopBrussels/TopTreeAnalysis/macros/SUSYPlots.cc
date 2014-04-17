#include "TStyle.h"
#include <cmath>
#include <fstream>

#include <iostream>		/*added */
#include <sstream>		/*added */
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include "TF1.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TBranch.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h" /**/
#include "TROOT.h" /**/
//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeProducer/interface/TRootGenEvent.h"
#include "../Selection/interface/SelectionTable.h"
#include "../BkgEstimationMethods/interface/VJetEstimation.h"
#include "../BkgEstimationMethods/interface/VJetEstPseudoExp.h"
#include "../BkgEstimationMethods/interface/ABCDEstimation.h"
#include "../BkgEstimationMethods/interface/ABCDPseudoExp.h"
#include "../BkgEstimationMethods/interface/CrossSectionPseudoExp.h"
#include "../Tools/interface/PlottingTools.h"
#include "../Tools/interface/MultiSamplePlot.h"
#include "../Tools/interface/MultiCutPlot.h"
#include "../Tools/interface/CutImpactEvaluation.h"
#include "../Tools/interface/TemplateComparator.h"
#include "../Tools/interface/TTreeLoader.h"
#include "../Tools/interface/AnalysisEnvironmentLoader.h"
#include "../Content/interface/AnalysisEnvironment.h"
#include "../Content/interface/MCExpectation.h"
#include "../Content/interface/MCObsExpectation.h"
#include "../Content/interface/Dataset.h"
#include "../StatProcedure/interface/MCPseudoExp.h"
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



void JetSelection (std::vector < TRootJet > &, Float_t &, Float_t &, Bool_t &);
std::vector < TRootJet * >BJetSelection (std::vector < TRootJet * >&, Int_t &, Float_t &, Bool_t &);

void DoObservablesRanker();

Float_t ChiSquareMatching (TRootMuon &, std::vector < TRootJet > &, UInt_t &, UInt_t *, Float_t &, Float_t &, Float_t &, Float_t &, Float_t &, Float_t &, Bool_t &);



void
Scale (vector < Dataset > datasets, double *&numbers, float Luminosity)
{
  for (unsigned int i = 0; i < datasets.size (); i++) {
    numbers[i] = numbers[i] * datasets[i].NormFactor () * Luminosity;
  }
}

// the following function return a vector of TRoot*Jet with P4 changed according to JES and with the main members required fill
// Improvement: create a method in TRootJet which scale the P4 directly

vector < TRootJet > JESRescale (vector < TRootJet > &jets, float factor)
{
  vector < TRootJet > output;
  for (unsigned int i = 0; i < jets.size (); i++) {
    TRootJet jet (TLorentzVector (jets[i].Px () * factor, jets[i].Py () * factor, jets[i].Pz () * factor, jets[i].Energy () * factor));
    jet.setJetType (jets[i].jetType ());
    jet.setPartonFlavour (jets[i].partonFlavour ());
    jet.setBtag_trackCountingHighEffBJetTags (jets[i].btag_trackCountingHighEffBJetTags ());
    jet.setBtag_trackCountingHighPurBJetTags (jets[i].btag_trackCountingHighPurBJetTags ());
    jet.setBtag_jetProbabilityBJetTags (jets[i].btag_jetProbabilityBJetTags ());
    jet.setBtag_jetBProbabilityBJetTags (jets[i].btag_jetBProbabilityBJetTags ());
    jet.setBtag_simpleSecondaryVertexBJetTags (jets[i].btag_simpleSecondaryVertexBJetTags ());
    output.push_back (jet);
  }
  return output;
}

vector < TRootPFJet > JESRescale (vector < TRootPFJet > &jets, float factor)
{
  vector < TRootPFJet > output;
  for (unsigned int i = 0; i < jets.size (); i++) {
    TRootPFJet jet (TLorentzVector (jets[i].Px () * factor, jets[i].Py () * factor, jets[i].Pz () * factor, jets[i].Energy () * factor));
    jet.setJetType (jets[i].jetType ());
    jet.setPartonFlavour (jets[i].partonFlavour ());
    jet.setBtag_trackCountingHighEffBJetTags (jets[i].btag_trackCountingHighEffBJetTags ());
    jet.setBtag_trackCountingHighPurBJetTags (jets[i].btag_trackCountingHighPurBJetTags ());
    jet.setBtag_jetProbabilityBJetTags (jets[i].btag_jetProbabilityBJetTags ());
    jet.setBtag_jetBProbabilityBJetTags (jets[i].btag_jetBProbabilityBJetTags ());
    jet.setBtag_simpleSecondaryVertexBJetTags (jets[i].btag_simpleSecondaryVertexBJetTags ());
    output.push_back (jet);
  }
  return output;
}

vector < TRootCaloJet > JESRescale (vector < TRootCaloJet > &jets, float factor)
{
  vector < TRootCaloJet > output;
  for (unsigned int i = 0; i < jets.size (); i++) {
    TRootCaloJet jet (TLorentzVector (jets[i].Px () * factor, jets[i].Py () * factor, jets[i].Pz () * factor, jets[i].Energy () * factor));
    jet.setJetType (jets[i].jetType ());
    jet.setPartonFlavour (jets[i].partonFlavour ());
    jet.setBtag_trackCountingHighEffBJetTags (jets[i].btag_trackCountingHighEffBJetTags ());
    jet.setBtag_trackCountingHighPurBJetTags (jets[i].btag_trackCountingHighPurBJetTags ());
    jet.setBtag_jetProbabilityBJetTags (jets[i].btag_jetProbabilityBJetTags ());
    jet.setBtag_jetBProbabilityBJetTags (jets[i].btag_jetBProbabilityBJetTags ());
    jet.setBtag_simpleSecondaryVertexBJetTags (jets[i].btag_simpleSecondaryVertexBJetTags ());
    jet.setEcalEnergyFraction (jets[i].ecalEnergyFraction ());
    jet.setHcalEnergyFraction (jets[i].hcalEnergyFraction ());
    jet.setn90Hits (jets[i].n90Hits ());
    jet.setfHPD (jets[i].fHPD ());
    output.push_back (jet);
  }
  return output;
}

/*
Aim: for each SUSY samples, providing a list of information in a text file with the name of the sample
and reuse this information later to produce 2-D plots in the m0,m1/2 plan
List of variables:
 - Xsection
 - #evts for a given lumi
 - #evts expected after selection
 - efficiency for a given selection
 
To add:
 - differentiate the several cuts of the selection
 - plots to caracterise the events (mean of #jets etc ...)
*/

int
main (int argc, char *argv[])
{
  Bool_t DEBUG = (argc > 2 ? (atoi (argv[2]) == 1) : 0);
  cout << "**********************************************************************" << endl;
  cout << "Begining of the program for bkg estimation in lepton+jets channels !" << endl;
  cout << "**********************************************************************" << endl;

  //SetStyle if needed
  //setTDRStyle(); 
  //setMyStyle();

  /////////////////////
  // Configuration
  /////////////////////
  ////////////////////////////////////////

  //xml file
  string xmlFileName = "../config/mySUSYconfig.xml";
  if (argc >= 2)
    xmlFileName = string (argv[1]);
  const char *xmlfile = xmlFileName.c_str ();
  //Output files
  string rootFileName ("SUSYPlots.root");
 

  // Jet Chi2 combination
  UInt_t Permutation[4] = { 0, 0, 0, 0 };
  Float_t WMass = 90.0;		// Derived from the SanityChecker (Truth module)
  Float_t TopMass = 173.1;
  Float_t blMass = 91.1;
  //Float_t LepTopMassResol = 15.6;
  Float_t LepblMassResol = 36.0;
  Float_t HadTopMassResol = 21.1;
  Float_t HadWMassResol = 12.7;
  UInt_t NjetsForComb = 7;
  Float_t ChiSq = 0;

  
  //Configuration output format

 
  TTree *configTree = new TTree ("configTree", "configuration Tree");
  TClonesArray *tcdatasets = new TClonesArray ("Dataset", 1000);
  configTree->Branch ("Datasets", "TClonesArray", &tcdatasets);
  TClonesArray *tcAnaEnv = new TClonesArray ("AnalysisEnvironment", 1000);
  configTree->Branch ("AnaEnv", "TClonesArray", &tcAnaEnv);
  
  ////////////////////////////////////
  /// AnalysisEnvironment  
  ////////////////////////////////////

  AnalysisEnvironment anaEnv;
  cout << "Loading environment ..." << endl;
  AnalysisEnvironmentLoader anaLoad (anaEnv, xmlfile);
  new ((*tcAnaEnv)[0]) AnalysisEnvironment (anaEnv);
  int verbose = anaEnv.Verbose;
  float Luminosity = anaEnv.Luminosity;	// in 1/pb
  cout << "The results will be obtained for a luminosity of " << Luminosity << " x 1/pb" << endl;

  /////////////////////
  // Load Datasets
  /////////////////////

  TTreeLoader treeLoader;
  cout << " - Load datasets ..." << endl;
  vector < Dataset > datasets;
  treeLoader.LoadDatasets (datasets, xmlfile);
  for (unsigned int i = 0; i < datasets.size (); i++)
    new ((*tcdatasets)[i]) Dataset (datasets[i]);
  /////////////////////

  //vector of objects
  cout << " - Variable declaration ..." << endl;
  vector < TRootVertex > vertex;
  vector < TRootMuon > init_muons;
  vector < TRootElectron > init_electrons;
  vector < TRootJet > init_jets;
  vector < TRootCaloJet > init_Calojets;
  vector < TRootPFJet > init_PFjets;
  vector < TRootMET > mets;
  vector < TRootJet > CombinedJets;
  //after JES applied
  vector < TRootJet > jets;
  vector < TRootCaloJet > Calojets;
  vector < TRootPFJet > PFjets;
  ////// jets needed for Chi-Square sorting


  //vector < TRootMuon > selectedMuons;
  vector < TRootMuon > vetoMuons;
  vector < TRootElectron > vetoElectrons;


  std::vector < TRootJet > SelectedJets;
  std::vector < TRootJet > selectedJetss;
  std::vector < TRootMuon > selectedMuons;
  vector < TRootJet > selectedJets;


  vector < TRootJet > ChiJets;

  //Global variable
  TRootEvent *event = 0;

  init_jets.clear ();
  init_muons.clear ();
  init_electrons.clear ();
  mets.clear ();

  SelectedJets.clear ();
  selectedJetss.clear ();;

  CombinedJets.clear ();
  ChiJets.clear ();
  selectedMuons.clear ();
  SelectedJets.clear ();
  ////////////////////////////////////////

  //Global variable

  int ievtmin = 0;

  char name[100];

  //nof selected events
  int NEvtsData = 0;

  ///////////////
  //
  //////////////
  //  TTree *tree_TtJets = 0;

  ////////////////////////////////////
  /// List of output information
  ////////////////////////////////////
    TTree *tree = new TTree ("SUSYInfo", "Info for SUSY samples");
  Double_t m0 = 0;
  Double_t m12 = 0;
  Double_t Xsection = 0;
  float TotalNofEvts = 0;
  float NofSelEvts = 0;
  float Efficiency = 0;
  float MuonSelEfficiency = 0;
  float SecondLeptonVetoEfficiency = 0;
  float JetSelEfficiency = 0;
  float MeanMET = 0;
  
  tree->Branch ("m0", &m0, "m0/F");
  tree->Branch ("m12", &m12, "m12/F");
  tree->Branch ("Xsection", &Xsection, "Xsection/F");
  tree->Branch ("TotalNofEvts", &TotalNofEvts, "TotalNofEvts/F");
  tree->Branch ("NofSelEvts", &NofSelEvts, "NofSelEvts/F");
  tree->Branch ("Efficiency", &Efficiency, "Efficiency/F");
  tree->Branch ("MuonEfficiency", &Efficiency, "MuonEfficiency/F");
  tree->Branch ("SecondLeptonVetoEfficiency", &Efficiency, "SecondLeptonVetoEfficiency/F");
  tree->Branch ("JetSelEfficiency", &JetSelEfficiency, "JetSelEfficiency/F");
  tree->Branch ("MeanMET", &MeanMET, "MeanMET/F");
  
  ////////////////////////////////////
  /// Plots
  ////////////////////////////////////
  float m0min = 0;
  float m0max = 1000;
  float m0Nbins = 10;
  float m12min = 0;
  float m12max = 1000;
  
  float m12Nbins = 10;
  TH2F *hXsection = new TH2F ("hXsection", "#sigma (pb)", m0Nbins, m0min, m0max, m12Nbins, m12min, m12max);
  hXsection->SetXTitle("m_0");
   hXsection->SetXTitle("m_12");
  TH2F *hGenLuminosity = new TH2F ("hGenLuminosity", "Generated Luminosity", m0Nbins, m0min, m0max, m12Nbins, m12min, m12max);
  //Effiencies;
  TH2F *hMuonEfficiency = new TH2F ("hMuonEfficiency", "Muon Selection Efficiency", m0Nbins, m0min, m0max, m12Nbins, m12min, m12max);
  TH2F *hSecondLeptonVetoEfficiency = new TH2F ("hSecondLeptonVetoEfficiency", "Second Lepton Veto Efficiency", m0Nbins, m0min, m0max, m12Nbins, m12min, m12max);
  TH2F *hJetSelEfficiency = new TH2F ("hJetSelEfficiency", "Jet Selection Efficiency", m0Nbins, m0min, m0max, m12Nbins, m12min, m12max);
  TH2F *hEfficiency = new TH2F ("hEfficiency", "Total Efficiency", m0Nbins, m0min, m0max, m12Nbins, m12min, m12max);
  TH2F *hNofSelEvents = new TH2F ("hNofSelEvents", "# selected events)", m0Nbins, m0min, m0max, m12Nbins, m12min, m12max);
  //Variables
  TH2F *hMeanMET = new TH2F ("hMeanMET", "Mean Value of MET", m0Nbins, m0min, m0max, m12Nbins, m12min, m12max);
  
  
  ////////////////////////////////////
  /// Selection table
  ////////////////////////////////////
  vector < string > CutsSelecTable;
  CutsSelecTable.push_back (string ("all"));
  CutsSelecTable.push_back (string ("trigged"));
  CutsSelecTable.push_back (string ("Good PV"));
  CutsSelecTable.push_back (string ("$\\geq$ 1 selected muon"));
  CutsSelecTable.push_back (string ("Veto 2nd muon"));
  CutsSelecTable.push_back (string ("Veto electron"));
  char LabelNJets[100];
  sprintf (LabelNJets, "$\\geq$ %d jets", anaEnv.NofJets - 2);
  CutsSelecTable.push_back (string (LabelNJets));
  sprintf (LabelNJets, "$\\geq$ %d jets", anaEnv.NofJets - 1);
  CutsSelecTable.push_back (string (LabelNJets));
  sprintf (LabelNJets, "$\\geq$ %d jets", anaEnv.NofJets);
  CutsSelecTable.push_back (string (LabelNJets));
  sprintf (LabelNJets, "$\\geq$ %d jets", anaEnv.NofJets + 1);
  CutsSelecTable.push_back (string (LabelNJets));
  SelectionTable selecTable (CutsSelecTable, datasets);
  selecTable.SetLuminosity (Luminosity);

  ////////////////////////////////////
  //    Loop on datasets
  ////////////////////////////////////
  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
  for (unsigned int d = 0; d < datasets.size (); d++) {
    if (verbose > 1)
      cout << "   Dataset " << d << ": " << datasets[d].Name () << endl;

    //open files and load
    treeLoader.LoadDataset (datasets[d], anaEnv);
    int itrigger = treeLoader.iTrigger (datasets[d], string ("HLT_Mu9"));
    string ds = datasets[d].Title ();
    /////////////////////////////
    // initialising values
    /////////////////////////////
       NEvtsData = 0;
    m0 = 0;
    m12 = 0;
    Xsection = 0;
    TotalNofEvts = 0;
    NofSelEvts = 0;
    Efficiency = 0;
    MuonSelEfficiency = 0;
    SecondLeptonVetoEfficiency = 0;
    JetSelEfficiency = 0;

    //Create plots used in the loop on Datasets
      TH1F hMET ("hMET", "MET", 100, 0, 2000);
    
    
    ////////////////////////////////////////////////
    // read info from input tree: m0,m12,Xsection
    ////////////////////////////////////////////////
    string filename;
    if (datasets[d].Filenames ().size () > 0)
      filename = datasets[d].Filenames ()[0];

    if (d>0){
    TFile ifile (filename.c_str (), "READ");
    TTree *m0Tree = (TTree *) ifile.Get ("m0");
    TTree *m12Tree = (TTree *) ifile.Get ("m12");
    TTree *XsectionTree = (TTree *) ifile.Get ("CrossX");
    TBranch *m0Branch = (TBranch *) m0Tree->GetBranch ("m0_input");
    TBranch *m12Branch = (TBranch *) m12Tree->GetBranch ("m12_input");
    TBranch *XsectionBranch = (TBranch *) XsectionTree->GetBranch ("x_sec");
    m0Branch->SetAddress (&m0);
    m12Branch->SetAddress (&m12);
    XsectionBranch->SetAddress (&Xsection);
    m0Tree->GetEvent (0);
    m12Tree->GetEvent (0);
    XsectionTree->GetEvent (0);
    }

  //cout<<"General information: Xsection = "<<Xsection<<" m0:"<<m0<<" m12:"<<m12<<endl;
    if (ds == "TTJets") {
           cout << " Sample is a SM...." << ds << endl;
      cout << "General information: Xsection = " << Xsection << endl;
    }
    if (ds != "SUSY") {
      cout << " Sample is not a SM...." << ds << endl;
      cout << "General information: Xsection = " << Xsection << " m0:" << m0 << " m12:" << m12 << endl;
    }
    
    //string currentTreeFileName = ds+"_tree.root";
    //TFile* currentTreeFile = new TFile(currentTreeFileName.c_str(),"RECREATE");

	  ////////////////////////////////////////
	  // Observables
	  Observables a; 
	  PlotObservables Plots(a.ListOfVariables(), a.RangeVariables()); 
	  TTreeObservables TtreeObs (a.ListOfVariables ());


	    ////////////////////////////////////
	    //  Loop on events
	    ////////////////////////////////////
	    if (verbose > 1)
	      cout << "	Loop over events " << endl;
	    int nevt;
	    nevt = 2000;
            nevt = datasets[d].NofEvtsToRunOver ();

	    //  for (unsigned int ievt = 0; ievt < nevt; ievt++) {
	    //    for (int ievt = 0; ievt < datasets[d].NofEvtsToRunOver (); ievt++) {

	    for (unsigned int ievt = 0; ievt < nevt; ievt++) {
	      if (ievt % 5000 == 0)
		std::cout << "Processing the " << ievt << "th event" << flush << "\r";
                  if (ievt == nevt-1)cout<< " Finishing in this dataset "<<nevt<<" events"<<endl;;//endl;

	      TLorentzVector JESbalance (0, 0, 0, 0);

	      //load event
	      if (anaEnv.JetType == 0)
		event = treeLoader.LoadEvent (ievt, datasets[d], vertex, init_muons, init_electrons, init_jets, mets);
	      if (anaEnv.JetType == 1)
		event = treeLoader.LoadEvent (ievt, datasets[d], vertex, init_muons, init_electrons, init_Calojets, mets);
	      if (anaEnv.JetType == 2)
		event = treeLoader.LoadEvent (ievt, datasets[d], vertex, init_muons, init_electrons, init_PFjets, mets);
	      ////////////////////////////////////////

	      //apply JES
	      jets = JESRescale (init_jets, anaEnv.JES);
	      Calojets = JESRescale (init_Calojets, anaEnv.JES);
	      PFjets = JESRescale (init_PFjets, anaEnv.JES);

	      /////////////////////////////
	      //   Selection
	      /////////////////////////////

	      selectedMuons.clear ();
	      vetoMuons.clear ();
	      vetoElectrons.clear ();
	      selectedJets.clear ();
		CombinedJets.clear ();
	      //Declare selection instance    
	      Selection selection (anaEnv.JetType, jets, Calojets, PFjets, init_muons, init_electrons, mets);
	      selection.SetConfiguration (anaEnv);




	      bool trigged = treeLoader.EventTrigged (itrigger);

	      /////////////////////////////
	      //   Selection table
	      /////////////////////////////

	      bool isInSR = false;
	      selectedMuons = selection.GetSelectedMuons (anaEnv.MuonPtCutSR, anaEnv.MuonEtaCutSR, anaEnv.MuonRelIsoCutSR);
	      vetoMuons = selection.GetSelectedLooseMuons (anaEnv.MuonPtCutVetoSR, anaEnv.MuonEtaCutVetoSR, anaEnv.MuonRelIsoCutVetoSR);
	      vetoElectrons = selection.GetSelectedLooseElectrons (anaEnv.ElectronPtCut, anaEnv.ElectronEtaCut, anaEnv.ElectronRelIsoCut);
	      selectedJets = selection.GetSelectedJets (anaEnv.JetsPtCutSR, anaEnv.JetsEtaCutSR, anaEnv.applyJetID);




	      //ChiJets =selection.GetSelectedJets (anaEnv.JetsPtCutSR, anaEnv.JetsEtaCutSR, anaEnv.applyJetID);

	      bool isGoodPV = selection.isPVSelected (vertex, anaEnv.PVertexNdofCut, anaEnv.PVertexZCut, anaEnv.PVertexRhoCut);


	      selecTable.Fill (d, 0, 1.);
	      if (trigged) {
		selecTable.Fill (d, 1, 1.);
		if (isGoodPV) {
		  selecTable.Fill (d, 2, 1.);
		  if (selectedMuons.size () == 1) {
		    MuonSelEfficiency++;
		    selecTable.Fill (d, 3, 1.);
		    if (vetoMuons.size () == 1) {

		      selecTable.Fill (d, 4, 1.);
		      if (vetoElectrons.size () == 0) {
			SecondLeptonVetoEfficiency++;
			selecTable.Fill (d, 5, 1.);
			if (selectedJets.size () >= (unsigned int) anaEnv.NofJets - 2) {;
			  selecTable.Fill (d, 6, 1.);
			  if (selectedJets.size () >= (unsigned int) anaEnv.NofJets - 1) {
			    selecTable.Fill (d, 7, 1.);
			    if (selectedJets.size () >= (unsigned int) anaEnv.NofJets) {

		      JetSelEfficiency++;
		      selecTable.Fill (d, 8, 1.);
		      isInSR = true;
		      if (selectedJets.size () >= (unsigned int) anaEnv.NofJets + 1)
			selecTable.Fill (d, 9, 1.);
		      if (verbose > 2)
			cout << event->runId () << ":" << event->eventId () << ":" << event->lumiBlockId () << ":" << setprecision (8) << selectedMuons[0].Pt () << endl;

		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      //////////////////////////////

      if (!trigged)
	continue;


      //////////////////////////
      //the different regions
      //////////////////////////
      float RelIso = 100.;
      //////////////////////////

      ///////////////////////////////////////////
      //fill booleans for the different region
      ///////////////////////////////////////////

      if (isInSR) {
	NofSelEvts = NofSelEvts + datasets[d].NormFactor () * Luminosity;
	NEvtsData++;

	//Fill plots
	//fout->cd();
	hMET.Fill (mets[0].Et ());
	


	ChiSq = ChiSquareMatching (selectedMuons[0], selectedJets, NjetsForComb, Permutation, WMass, TopMass, blMass, HadWMassResol, HadTopMassResol, LepblMassResol, DEBUG);



	for (Int_t i = 0; i < 4; i++) {
	  CombinedJets.push_back (selectedJets[Permutation[i]]);
	  if (Permutation[i] >= selectedJets.size ())
	    cout << " perm " << Permutation[i] << endl;
	}


	// if ( selectedMuons.size()>0 && selectedJets.size()>3){
	Observables a (selectedMuons[0], CombinedJets, selectedJets, mets[0], ds, ChiSq);
       

		Plots.Fill(a);
		TtreeObs.Fill (a);
	
      }

      //delete selection;
    }				//loop on events



    //comput efficiencies
        if (datasets[d].NofEvtsToRunOver () > 0) {
         Efficiency = (float) NEvtsData / datasets[d].NofEvtsToRunOver ();
      MuonSelEfficiency = MuonSelEfficiency / datasets[d].NofEvtsToRunOver ();
      SecondLeptonVetoEfficiency = SecondLeptonVetoEfficiency / datasets[d].NofEvtsToRunOver ();
      JetSelEfficiency = JetSelEfficiency / datasets[d].NofEvtsToRunOver ();
    }
      
    //Fill other quantities used in the TTree
    //fout->cd();
    MeanMET = hMET.GetMean ();
    //Fill the TTree
    tree->Fill ();
    cout<<""<<endl;
    cout << "NofSelEvents: " << NEvtsData << endl;
      
       //fill the plots

    cout << "  ended  plots loop in events for " << ds << endl;

  
    string dsroot = ds + "_tree.root";
    Plots.Write(ds,true);
    TtreeObs.Write(ds,true);


      //fout->cd();
    hXsection->Fill (m0, m12, Xsection);
    hGenLuminosity->Fill (m0, m12, (float) datasets[d].eventTree ()->GetEntries () / datasets[d].Xsection ());
    //Efficiencies
    hMuonEfficiency->Fill (m0, m12, MuonSelEfficiency);
    hSecondLeptonVetoEfficiency->Fill (m0, m12, SecondLeptonVetoEfficiency);
    hJetSelEfficiency->Fill (m0, m12, JetSelEfficiency);
    hEfficiency->Fill (m0, m12, Efficiency);
    hNofSelEvents->Fill (m0, m12, NofSelEvts);
    //variables
    hMeanMET->Fill (m0, m12, MeanMET);
    
     
    cout<<"delete TtreeObs"<<endl;
    cout<<"written"<<endl;
    //currentTreeFile->Write();
    cout<<"..."<<endl;
    //delete TtreeObs;
    cout<<"..."<<endl;
    //currentTreeFile->Close();
    cout<<"closed"<<endl;
    //currentTreeFile->Delete();
    cout<<"deleted"<<endl;
    // delete currentTreeFile;
    
    //important: free memory
     treeLoader.UnLoadDataset ();
     
  }				//loop on datasets


 
  //Once everything is filled ...
  if (verbose > 0)
    cout << " We ran over all the data ;-)" << endl;



  ///////////////////
  // Writting
  //////////////////
  if (verbose > 1)
    cout << " - Writting outputs on files ..." << endl;

  //add configuration info
  TFile *fout = new TFile (rootFileName.c_str (), "UPDATE");
   fout->cd ();

  configTree->Fill ();
  configTree->Write ();
  tree->Write ();
  
  ///////////////////
  //write the plots
  //////////////////
  
  hXsection->Write ();
  hGenLuminosity->Write ();
  //efficiencies
  hMuonEfficiency->Write ();
  hSecondLeptonVetoEfficiency->Write ();
  hJetSelEfficiency->Write ();
  hEfficiency->Write ();
  hNofSelEvents->Write ();
  delete hXsection; delete hGenLuminosity;delete hMuonEfficiency;delete hSecondLeptonVetoEfficiency;delete hJetSelEfficiency;delete hNofSelEvents;delete hEfficiency;
 //Variables
  if (verbose > 1)
    cout << " - Done with writing the module outputs in the ouput file ..." << endl;
  cout << " - Closing the output file now..." << endl;
  fout->Write ();
 
   fout->Close();



   //delete
  delete fout;
  delete event;
  delete tcdatasets;
  delete tcAnaEnv;
  delete configTree;
 
   DoObservablesRanker();
  cout << "**********************************************************************" << endl;
  cout << "           End of the program !!" << endl;
  cout << " 		hasn't crashed yet ;-) " << endl;
  cout << "**********************************************************************" << endl;

 
  //return 0;
}


void DoObservablesRanker(){

  cout<< "  entering ObservablesRanker.... "<<endl;
 string xmlFileName = "../config/mySUSYconfig.xml";
  const char *xmlfile = xmlFileName.c_str ();
 TTreeLoader treeLoader;
  cout << " - Load datasets ..." << endl;
  vector < Dataset > datasets;
  treeLoader.LoadDatasets (datasets, xmlfile);
string fS, fBkg,merged_file_,np_file_,sm_file_;
   TFile *hfileS, *hfileB; TTree *tSignal, *tBkg;


            //fS = "SUSY";


             fBkg = datasets[0].Title()+"_tree.root"; 
             fS = datasets[0].Title() +"_tree.root";
             sm_file_=datasets[0].Title() ;//should be the TTjet
             
             char name_[100];   char name2[100];
          
             sprintf (name_, "%s_tree.root",datasets[0].Title().c_str() ); 

 for (unsigned int j=1;j<datasets.size();j++){
             merged_file_="Merged_Plots_"+datasets[0].Title()+"_"+datasets[j].Title()+".root"   ;


	     cout<<" Will merge and compute overlap for "<< datasets[0].Title()<<"  and  "<<datasets[j].Title()<<"  file will be "<<merged_file_.c_str()<<endl;

        
              np_file_=datasets[j].Title();
              
         
              sprintf (name2, "%s_tree.root",datasets[j].Title().c_str() ); 
  
              hfileS = new TFile(name_);
              hfileB = new TFile(name2);
  


              cout<<name_<<"  xixixi "<<name2<<endl;



//ffff->Print();
  
 bool SameFiles = false;
    if (fS == fBkg) (SameFiles=1);


   tSignal = (TTree*)hfileS->Get("OBS");
   tBkg = (TTree*)hfileB->Get("OBS");

  
cout<<"  GOT OBS TREE "<<endl;

    ObservablesRanker ObsR(tSignal,tBkg,merged_file_,np_file_,sm_file_,SameFiles);

    cout<<"  EXITED OBSERVABLESRANKER "<<endl;
    // 
    if (j==datasets.size()-1) cout<<"Processed all datasets...."<<endl;
 }
     treeLoader.UnLoadDataset ();
 
       delete hfileS;
        cout<< " deleted 1 "<<endl;
	delete hfileB;  cout<< " deleted 2 "<<endl;
	//	  delete tSignal;  cout<< " deleted 3 "<<endl;
	//delete tBkg;   cout<< " deleted 4 "<<endl;
 


}





Float_t
ChiSquareMatching (TRootMuon & muon, std::vector < TRootJet > &jets, UInt_t & n, UInt_t * Permutation, Float_t & WMass_, Float_t & TopMass_, Float_t & blMass_, Float_t & HadWMassResol_, Float_t & HadTopMassResol_, Float_t & LepblMassResol_,
		   Bool_t & DEBUG)
{

  if (jets.size () < 4)
    return -1;
  UInt_t NbOfJets = 0;
  (jets.size () < n ? NbOfJets = jets.size () : NbOfJets = n);
  UInt_t *numbers = new UInt_t[NbOfJets];
  for (UInt_t i = 0; i < NbOfJets; i++)
    numbers[i] = i;
  UInt_t *comb = new UInt_t[4];
  for (UInt_t i = 0; i < 4; i++)
    comb[i] = i;
  UInt_t NbOfPermutations = 0;
  Float_t ChiSquare = 0, ChiSquare_tmp = 0;
  Float_t HadWChiSquare = 0, HadTopChiSquare = 0, LepblChiSquare = 0;
  Float_t HadWCandMass_tmp = -1, HadTopCandMass_tmp = -1, LepBlCandMass_tmp = -1;
  Float_t HadWCandMass = -1, HadTopCandMass = -1, LepBlCandMass = -1;
  //DEBUG = true;
  ChiSquare = 9999999;
  do {
    if (DEBUG) {
      std::cout << "-- Enter loop for jet combination " << std::endl;
      std::cout << "-- Considered Jet combination : " << comb[0] << "/" << comb[1] << "/" << comb[2] << "/" << comb[3] << std::endl;
    }
    NbOfPermutations = 0;
    do {
      NbOfPermutations++;
      if (DEBUG)
	std::cout << "--- Current permutations : " << comb[0] << "/" << comb[1] << "/" << comb[2] << "/" << comb[3] << std::endl;
      HadWCandMass_tmp = (jets[comb[0]] + jets[comb[1]]).M ();
      HadTopCandMass_tmp = (jets[comb[0]] + jets[comb[1]] + jets[comb[2]]).M ();
      LepBlCandMass_tmp = (jets[comb[3]] + muon).M ();

      HadWChiSquare = pow ((HadWCandMass_tmp - WMass_), 2) / pow (HadWMassResol_, 2);
      HadTopChiSquare = pow ((HadTopCandMass_tmp - TopMass_), 2) / pow (HadTopMassResol_, 2);
      LepblChiSquare = pow ((LepBlCandMass_tmp - blMass_), 2) / pow (LepblMassResol_, 2);

      ChiSquare_tmp = HadWChiSquare + HadTopChiSquare + LepblChiSquare;
      if (ChiSquare_tmp < ChiSquare) {
	//for(unsigned int j=0;j<4;j++) Permutation_tmp[i][j] = comb[j];
	//ChiSquare[i] = ChiSquare_tmp[i];
	for (UInt_t i = 0; i < 4; i++)
	  Permutation[i] = comb[i];
	ChiSquare = ChiSquare_tmp;
	HadWCandMass = HadWCandMass_tmp;
	HadTopCandMass = HadTopCandMass_tmp;
	LepBlCandMass = LepBlCandMass_tmp;
      }
    }
    while (next_permutation (comb, comb + 4));
    if (DEBUG) {
      std::cout << "--- ChiSquare matching : Nb of permutations found = " << NbOfPermutations << std::endl;
      std::cout << "--- ChiSquare matching : Min ChiSquare found = " << ChiSquare << std::endl;
      std::cout << "--- ChiSquare matching : associated to the permutation : " << Permutation[0] << "/" << Permutation[1] << "/" << Permutation[2] << "/" << Permutation[3] << std::endl;
    }
  }
  while (stdcomb::next_combination (numbers, numbers + NbOfJets, comb, comb + 4));
  if (DEBUG) {
    std::cout << "ChiSquare = " << ChiSquare << std::endl;
    std::cout << "Associated had W mass   = " << HadWCandMass << std::endl;
    std::cout << "Associated had top mass = " << HadTopCandMass << std::endl;
    std::cout << "Associated lep b+l mass = " << LepBlCandMass << std::endl;
  }
  return ChiSquare;
  delete comb;
  delete numbers;
}
