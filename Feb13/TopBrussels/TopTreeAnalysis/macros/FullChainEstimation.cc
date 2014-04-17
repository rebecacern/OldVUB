#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
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
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TChainElement.h"
#include <cmath>
#include <fstream>
 
//user code
#include "TopTreeProducer/interface/TRootGenEvent.h"
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "../Selection/interface/SelectionTable.h"
#include "../BkgEstimationMethods/interface/TtJetEstimation.h"
#include "../BkgEstimationMethods/interface/TtJetEstPseudoExp.h"
#include "../BkgEstimationMethods/interface/BkgEstimationSummary.h"
#include "../Tools/interface/PlottingTools.h"
#include "../Tools/interface/MultiSamplePlot.h"
#include "../Tools/interface/MultiCutPlot.h"
#include "../Tools/interface/CutImpactEvaluation.h"
#include "../Tools/interface/TemplateComparator.h"
#include "../Tools/interface/TTreeLoader.h"
#include "../Tools/interface/AnalysisEnvironmentLoader.h"
#include "../Reconstruction/interface/MEzCalculator.h"
#include "../Reconstruction/interface/Observables.h"
#include "../Reconstruction/interface/PlotObservables.h"
#include "../Reconstruction/interface/TTreeObservables.h"
#include "../Reconstruction/interface/ObservablesRanker.h"
#include "../Reconstruction/interface/MakeBinning.h"
#include "../Content/interface/AnalysisEnvironment.h"
#include "../Content/interface/Dataset.h"
#include "../Content/interface/MCExpectation.h"
#include "../Content/interface/MCObsExpectation.h"
#include "../Content/interface/Container.h"
#include "../StatProcedure/interface/EventCombinedWeightCalculator.h"
#include "../StatProcedure/interface/SampleCombinedWeightCalculator.h"
#include "../StatProcedure/interface/MCPseudoExp.h"
#include "../StatProcedure/interface/WeightProbaCalculator.h"
#include "Style.C"

using namespace std;
using namespace TopTree;

void
Scale (vector < Dataset > datasets, double *&numbers, float Luminosity)
{
  for (unsigned int i = 1; i < datasets.size (); i++) {
    numbers[i] = numbers[i] * datasets[i].NormFactor () * Luminosity;
  }
}
Float_t ChiSquareMatching (TRootMuon* &, std::vector < TRootJet* > &, UInt_t &, UInt_t *, Float_t &, Float_t &, Float_t &, Float_t &, Float_t &, Float_t &, Bool_t &);


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

  //xml file
  const char *xmlfile = "../config/myNPconfig.xml";
  //Output files
  string rootFileName ("FullChainEstimation.root");
  string rootStatOutputFileName ("StatOutput.root");
  string rootStatVtreeFileName ("Vtree_StatOutput.root");
  // Jet Chi2 combination
  UInt_t Permutation[4] = { 0, 0, 0, 0 };
  Float_t WMass = 90.0;		// Derived from the SanityChecker (Truth module)
  Float_t TopMass = 174.5;
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
  TClonesArray *tcMCObsExp = new TClonesArray ("MCObsExpectation", 1000);
  configTree->Branch ("MCObsExp", "TClonesArray", &tcMCObsExp);
  TClonesArray* tcMCExpBjetMult = new TClonesArray("MCObsExpectation",1000);
  configTree->Branch("BjetMultiplicity","TClonesArray",&tcMCExpBjetMult);



  ////////////////////////////////////
  /// AnalysisEnvironment 
  ////////////////////////////////////
  AnalysisEnvironment anaEnv;
  cout << "Loading environment ..." << endl;
  AnalysisEnvironmentLoader anaLoad (anaEnv, xmlfile);
  new ((*tcAnaEnv)[0]) AnalysisEnvironment (anaEnv);
  int verbose = anaEnv.Verbose;
  float Luminosity = anaEnv.Luminosity;	// in 1/pb
 
  int nbins = anaEnv.nbins;
  int EventsPerBin=anaEnv.eventsperbin;
  float epsilon = anaEnv.EpsilonValue;
  float Correl_Cut = anaEnv.Correl_cut;
  int nPseudos = anaEnv.nPseudoSession;
  cout << "The results will be obtained for a luminosity of " << Luminosity << " x 1/pb" <<" and the EventsPerBin will be "<<EventsPerBin<< endl;
  cout << "Nb of jet bins considered : " << anaEnv.NofJetBins << " , starting with events with at least " << anaEnv.NofJets << " jets." << endl;



  //////////////////////////////
  //for statistical procedure
  //////////////////////////////
  string SignificanceCombination("SUM"); // to combine squared bin significances to get a combined weight for the event
  TString Combination("sum");  // to combine weights to get V value
  vector<float> FractionHWEvts = anaEnv.FractionHWEvts;
  unsigned int nofFractionsHWEvts = FractionHWEvts.size();
  //vector<unsigned int> sensitive_Obs; //for historical reasons; this is in principle not needed, but this vector should be kept empty and passed as an argument to one of the StatProcedure classes (should be removed/adapted in future)
  vector< pair< float , float > > Vvect; //first entry: fraction x, second entry: V value (vector without adding uncertainty to estimation)
  vector< pair< float , float > > FracNPEventsInHWSubsampleVect; //first entry: fraction x, second entry: fraction NP events in highest weight subsample
  ////////////////////////////////////////
  //for the statistical procedure: write a root-file per dataset containing all the info per pseudo-exp
  TTree* PseudoExpDumpTree = new TTree("PseudoExp","Info per pseudo-experiment");
  Container container;
  int itContainer = 0;
  TClonesArray* tcContainer = new TClonesArray("Container",1000);
  PseudoExpDumpTree->Branch("SelEvents","TClonesArray",&tcContainer);
  

  /////////////////////
  // Load Datasets
  /////////////////////
  TTreeLoader treeLoader;
  if (verbose > 0)
    cout << " - Load datasets ..." << endl;
  vector < Dataset* > datasets;
  treeLoader.LoadDatasets (datasets, xmlfile);
  for (unsigned int i = 0; i < datasets.size (); i++)
    new ((*tcdatasets)[i]) Dataset (*datasets[i]);

  TH1F* hRelIsoData = new TH1F("hRelIso_Data","hRelIso_Data",50,0.,1.);
  MultiSamplePlot myMultiSamplePlot(datasets,"hRelIso",50, 0.,1., "RelIso", "#Events");

  /////////////////////
  if (verbose > 0)
    cout << " - Variable declaration ..." << endl;

  //vector of objects
  vector < TRootVertex* > vertex;
  vector < TRootMuon* > init_muons;
  vector < TRootElectron* > init_electrons;
  vector < TRootJet* > init_jets;
  vector < TRootMET* > mets;
  vector < TRootJet* > CombinedJets;
  //after JES applied
  vector < TRootJet* > jets;

  //Global variable
  TRootEvent *event = 0;
  TTree *tree = new TTree ("OBS", "list of observables");


	
	//nof selected events
  float NEvtsData = 0;
  Double_t *nEvents = new Double_t[datasets.size()];

  vector<int > nofSelEvts_vector; //(z) added


  ////////////////////////////////////////
  // Observables
  Observables obs;
 // MakeBinning NewBins;
 //
 
  //MakeBinning NewBins (obs.ListOfVariables (), obs.RangeVariables ());
 //
  ////////////////////////////////////
  /// Selection table
  ////////////////////////////////////

  vector<string> CutsSelecTable;
  CutsSelecTable.push_back(string("initial"));
  CutsSelecTable.push_back(string("preselected"));
  CutsSelecTable.push_back(string("trigged"));
  CutsSelecTable.push_back(string("Good PV"));
  CutsSelecTable.push_back(string("1 selected muon"));
  CutsSelecTable.push_back(string("Veto 2nd muon"));
  CutsSelecTable.push_back(string("Veto electron"));
  char LabelNJets[100];
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-3);
  CutsSelecTable.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-4);
  CutsSelecTable.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-1);
  CutsSelecTable.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets);
  CutsSelecTable.push_back(string("Veto electron"));
  CutsSelecTable.push_back(string("maxprob"));
  CutsSelecTable.push_back(string("probnocorr"));
  CutsSelecTable.push_back(string("X=0.0"));
  CutsSelecTable.push_back(string("X=0.1"));
  CutsSelecTable.push_back(string("X=0.2"));
  CutsSelecTable.push_back(string("X=0.3"));
  CutsSelecTable.push_back(string("X=0.4"));
  CutsSelecTable.push_back(string("X=0.5"));
  CutsSelecTable.push_back(string("X=0.6"));
  CutsSelecTable.push_back(string("X=0.7"));
  CutsSelecTable.push_back(string("X=0.8"));
  CutsSelecTable.push_back(string("X=0.9"));
  CutsSelecTable.push_back(string("X=1.0"));
  if (verbose > 0)
    cout << " - CutsSelectionTable instantiated ..." << endl;
  SelectionTable selecTable(CutsSelecTable, datasets);
  selecTable.SetLuminosity(Luminosity);
  if (verbose > 0)
    cout << " - SelectionTable instantiated ..." << endl;


  ////////////////////////////////////////
  anaEnv.ChecklistOfObs ();	//necessary to be called before all the loops over the variables!
  if (!anaEnv.allObsDefined) {
    cout << "Aborting" << endl;
    return 1;
  }
  vector <string > lstVar;
  lstVar.clear();
  vector <string > VarsListString;
  vector <string > VarsListByInt;
  //vector < pair < string, int > > VarsListString_All_; 
  //VarsListString_All_.clear();////first dim the datasets, second dim the rankded variables;
  VarsListByInt.clear();
  VarsListString.clear();
  lstVar = obs.ListOfVariables ();
  vector <string > lstVarTest;
  lstVarTest.clear();
	    //Loading of the binning
  map < string, TAxis * >mapAxis;
  TFile fBinning (anaEnv.binningFile.c_str (), "READ");
  //Binning for the given observable
  for (vector < string >::iterator iter = lstVar.begin (); iter != lstVar.end (); iter++) {
    TAxis *axis = NULL;
    char tAxisName[150];

    sprintf (tAxisName, "Binning_%s_SM", iter->c_str ());
    fBinning.GetObject (tAxisName, axis);
  
    if (axis == NULL)
      {
      cout<<" While reading the  "<<tAxisName<<"  the axis is NULL ..... "<<endl;
      mapAxis.insert (make_pair (*iter, (TAxis *) NULL));}
    else{
      mapAxis.insert (make_pair (*iter, new TAxis (*axis))); 
	/*for( map<string, TAxis*>::iterator ii=mapAxis.begin(); ii!=mapAxis.end(); ++ii)
		{ 
		cout << (*ii).first << ": " <<(*ii).second->GetNbins()<< " Title  "<<(*ii).second->GetName()<<" min "<<(*ii).second->GetXmin()<<"  max   "<<(*ii).second->GetXmax()<<endl;
		}
	*/
    }
      	//delete axis;
  }  fBinning.Close ();

  /////////////this is the vector that it will be used to calculate the Correlance factors.....
  vector <float> Correl_Cut_;
  Correl_Cut_.clear();
  for (unsigned int j=0;j<5;j++){
     float v_cut=0.5 + 0.1*j;
     Correl_Cut_.push_back(v_cut);
  }
  
  vector<pair<string,float> > BinContentObsAll_;
  vector<pair<string,float> > BinContentObsPerSample_;
  BinContentObsAll_.clear();
  BinContentObsPerSample_.clear();
  TFile *fout = new TFile (rootFileName.c_str (), "RECREATE");
  TFile* foutStat = new TFile(rootStatOutputFileName.c_str(),"RECREATE");

  //now creating directories in the Stat root file
  TDirectory * mytdirFraction;
  TDirectory * mytdirVariable;
  for(unsigned int k=0;k<nofFractionsHWEvts;k++){
  	ostringstream Fractionstrstream;
        Fractionstrstream << FractionHWEvts[k];
  	mytdirFraction = foutStat->mkdir("x="+TString(Fractionstrstream.str()));
	for(unsigned int v=0;v<lstVar.size (); v++){
		if (!anaEnv.runOnObs ((int) v))
      			continue;
		mytdirFraction->cd();
		char varname[50];
		sprintf(varname,lstVar[v].c_str());
		mytdirVariable = mytdirFraction->mkdir(varname);
	}
  }
  TDirectory * mytdirall = foutStat->mkdir("AllEvents");
  for(unsigned int v=0;v<lstVar.size (); v++){
		if (!anaEnv.runOnObs ((int) v))
      			continue;
		mytdirall->cd();
		char varname[50];
		sprintf(varname,lstVar[v].c_str());
		mytdirVariable = mytdirall->mkdir(varname);
   }

  int* ievtmin = new int[datasets.size()];
  for(unsigned int i=0;i<datasets.size();i++) ievtmin[i] = 0;

  char name[100];
  


  ////////////////////////////////////
  //  MCExpectation (if MC)
  ////////////////////////////////////
  MCExpectation *MCExp = NULL;
  vector < MCObsExpectation * >MCObsExpBjetMult;
  //size must be [anaEnv.NofJetBins*anaEnv.NofBtagWorkingPoint_vjEst];
  //store in a array b-jet multiplicity for #jets per working point (ex: 3jWP0,4jWP0,3jWP1,4jWP1)
  vector<MCObsExpectation *> MCObsExp;
//temporary commenting because of problem with TStreamerinfo and classdef...
/*  
  if (anaEnv.isMC && anaEnv.nPseudoExp<1) {
    //create as first step
    MCExp = new MCExpectation ();
    MCExp->SetLuminosity (Luminosity);
    MCExp->SetLabel (string ("CrossSection"));
  }
  else {  
    //read one
    TFile fin (anaEnv.MCExpFilename.c_str (), "READ");
    TClonesArray *tcMCExpTemp = new TClonesArray ("MCExpectation", 0);
    TClonesArray *tcMCObsExpTemp = new TClonesArray ("MCObsExpectation", 0);
    TTree *tree = (TTree *) fin.Get ("configTree");
    TBranch *b1 = 0;
    TBranch *b3 = 0;
    if (tree) {
      b1 = tree->GetBranch ("MCExp");
      if (b1)
	b1->SetAddress (&tcMCExpTemp);
      b3 = tree->GetBranch ("MCObsExp");
      if (b3)
	b3->SetAddress (&tcMCObsExpTemp);
      tree->GetEntry (0);
    }
    MCExpectation *temp = 0;
    if (tcMCExpTemp->GetEntries ()) {
      temp = (MCExpectation *) tcMCExpTemp->At (0);
      MCExp = (MCExpectation *) temp->Clone ();
    }
    MCObsExpectation *temp3 = 0;
    for (unsigned x = 0; x < tcMCObsExpTemp->GetEntries (); x++) {
      temp3 = (MCObsExpectation *) tcMCObsExpTemp->At (x);
      MCObsExp.push_back ((MCObsExpectation *) temp3->Clone ());
    }
    delete temp;
    delete temp3;
    tcMCExpTemp->Delete ();
    tcMCObsExpTemp->Delete ();
    tree->Delete ();
  }
*/

  ////////////////////////////////////
  //    Loop on  Pseudo-experiments
  ////////////////////////////////////
  int nRuns = 1;
  if (anaEnv.nPseudoExp != 0)
    nRuns = anaEnv.nPseudoExp;
  for (int iPseudoExp = 0; iPseudoExp < nRuns; iPseudoExp++) {

    cout << " ** PseudoExp no: " << iPseudoExp << "/" << nRuns-1 << " **" << endl;
    //////////////////////////////////////// 
    // For the StatProcedure
    //////////////////////////////////////// 
    vector < vector < float > > VarValues; //first dim: events ; second dim: variables
    vector < unsigned int > NPbooleans;  //to indicate if an event is NP (1) or not (0)
    vector < vector < pair < string, float > > > VarValues_WithString; //first dim: events ; second dim: variables
    
    unsigned int isNP = 0;
    vector < TH1F * > ObsEstimated;
    vector < TH1F * > ObsData;
    vector < TH1F * > ObsSUSY;
    vector < TH1F * > ObsSUSYByVar;
    vector < TH1F * > ObsTtjetByVar;

    MCExpectation MCExp_PS;
    MCExp_PS.SetLuminosity (anaEnv.Luminosity);
    MCExp_PS.SetLabel (string ("CrossSection"));

    Double_t *nEvents = new Double_t[datasets.size ()];
    ////////////////////////////////////
    /// TtJetEstimation
    ////////////////////////////////////
    //Variables used to define the Control Region
    vector < Dataset > datasetsBkg;
    vector < Dataset > datasetsTTJets;
    vector < Dataset > datasetsNP;
    for (unsigned int d = 0; d < datasets.size (); d++) {
      if (datasets[d]->Name () == "QCD" || datasets[d]->Name () == "WJets" || datasets[d]->Name () == "ZJets" || datasets[d]->Name () == "SingleTop")
	datasetsBkg.push_back (*datasets[d]);
      if (datasets[d]->Name () == "TTJets")
	datasetsTTJets.push_back (*datasets[d]);
      if (datasets[d]->Name ().find ("LM") < datasets[d]->Name ().size ())
	datasetsNP.push_back (*datasets[d]);
    }


///////////////////// Loading the datasets for filling which Variables to run...

    if (anaEnv.Vars_ByFile){   
      cout << " - Loop over datasets for creating the VarsList... " << datasets.size () << " datasets !" << endl;
      for (unsigned int d = 0; d < datasets.size (); d++) {
        if(datasets[d]->Name().find("LM")<=datasets[d]->Name().size() || datasets[d]->Name().find("Zp")<=datasets[d]->Name().size() || datasets[d]->Name().find("NP")<=datasets[d]->Name().size() ||  datasets[d]->Name().find("SUSY")<=datasets[d]->Name().size() ){
	     isNP=1;
        }
        else {
	     isNP=0;
        }
        
        if (!isNP ) {
         string  filename = datasets[d]->Filenames ()[0];
  	
         if (d>0 && datasets[d]->Title() == datasets[d-1]->Title()) continue;

         //string  merged_f_ = "Merged_Plots_" + datasets[0]->Title () + "_" + datasets[d]->Title () + ".root"; 
         string  merged_f_ = "Merged_Plots_TTJets_m0_100_m12_100.root";// + datasets[0]->Title () + "_" + datasets[d]->Title () + ".root"; 
	 TFile *f = TFile::Open(merged_f_.c_str(),"read");
	   
	 if (f->IsZombie() || !f) {		    
		   lstVar =obs.ListOfVariables ();
                   cout<<"  Will load the Observables from Observables class, as file  "<<merged_f_<<" does not exist..."<<endl;
		   break;
	 }

         float cor = anaEnv.Correl_cut;
	 //cor=0.5;
	 char tree_name[100];
	 sprintf(tree_name,"ListVector with Cor of %f",cor);
	 TTree *t = 0;
	 f->GetObject(tree_name,t);
	 TBranch *br_list=0;
	 vector<pair<string, float> > overlaps_sorted_cleaned_;
	 vector<pair<string, float> >* VarsList=&overlaps_sorted_cleaned_;
	 t->SetBranchAddress("overlaps_sorted_cleaned_",&VarsList,&br_list);
	 //cout<<"  OK "<<endl;
	 Long64_t tentry = t->LoadTree(0);
	 br_list->GetEntry(tentry);
	 float best_var = overlaps_sorted_cleaned_[0].second;
	 for (UInt_t j=0;j<overlaps_sorted_cleaned_.size();++j){
   	   //cout<<" NOW READ OUT THE VECTOR   "<<overlaps_sorted_cleaned_[j].first<<"    "<<overlaps_sorted_cleaned_[j].second<<"   "<<datasets[d]->Title()<<"   "<<j<<"   "<<epsilon<<"   "<<Correl_Cut<<endl;
           if ( overlaps_sorted_cleaned_[j].second <= best_var*epsilon) 
	      VarsListString.push_back(overlaps_sorted_cleaned_[j].first);
	//   VarsListString_All_.push_back(pair < string , int > (overlaps_sorted_cleaned_[j].first , d));
	 }
	//anaEnv.IntToCut=4;
	 if (anaEnv.IntToCut>0){
	   for (UInt_t j=0;j<anaEnv.IntToCut;++j){
              VarsListByInt.push_back(overlaps_sorted_cleaned_[j].first);
	      cout<< " The Variables will be defined by the anaEnv.IntToCut paramaeter....---> Now got the variable "<<j+1<< " "<<VarsListByInt[j]<<endl;
	   }
	 }
	 delete t;
	 f->Close();
	 delete f;

         lstVar.clear(); 
         if (anaEnv.IntToCut>0) { 
	   lstVar = VarsListByInt;
           cout<<" The Variables will be loaded from file up to the "<<anaEnv.IntToCut<< " variable...."<<endl;	     
     	 }     
         else lstVar=VarsListString;	   

        }///end of isNP

	if (anaEnv.isMC){ 
	 cerr<< " You cant run MC and Pseudos at the same time....check the config file "<<endl;
	 return -1 ;
	}//!Vars_ByFile_;///in case the accidentanly you put MC and Pseudos then,
    
      }////end of loop in dataset for setting the vars according to NP
      
      //  cout<<" Conditions are  --> isMC  "<<anaEnv.isMC<<" Get Variables from file  "<<anaEnv.Vars_ByFile<<" mSUGRA point  "<<anaEnv.VarsFile<<"  size of variables lstVar "<<lstVar.size()<<" VarListString  "<<VarsListString.size()<<endl;

      lstVar.clear();
      lstVar.push_back("ET1oET4");
      //lstVar.push_back("HToHZ");
      lstVar.push_back("AllJetsPtMET");
      lstVar.push_back("MET");
      lstVar.push_back("PtMuon");
    }////end of if(anaEnv.VarsBy_File)

     /*lstVar.clear();
      lstVar.push_back("ET1oET4");
      //lstVar.push_back("HToHZ");
      lstVar.push_back("AllJetsPtMET");
      lstVar.push_back("MET");
      lstVar.push_back("PtMuon");*/
       //if you do this there will be problems if the TAxes of the other variables are NULL etc. Basically, using FullChainEstimation now, you never need to do VarsBy_File = 1 (since pseudoexp are done by GOFanalysis.cc)!!

    //For Background summary
    TH1F **hData = new TH1F *[lstVar.size ()];
    TH1F **hTtJetMC = new TH1F *[lstVar.size ()];
    TH1F **hSUSYMC = new TH1F *[lstVar.size ()];
    TH1F **hSUSYMCByVar = new TH1F *[lstVar.size ()];
    TH1F **hTtjetMCByVar = new TH1F *[lstVar.size ()];
    ///////////////////////////////
    //Declare Shape estimation methods objects as an array of variable
    ///////////////////////////////
/*    TtJetEstimation **ttjEstimation = new TtJetEstimation *[lstVar.size ()];
    TtJetEstPseudoExp **ttjEstPseudoExp = new TtJetEstPseudoExp *[lstVar.size ()];
*/    ///////////////////////////////


    //LOOP OVER VARIABLES ... only for shape estimation
    if (verbose > 0)
      cout << "Run on " << lstVar.size () << " variables " << endl;
    for (unsigned int i = 0; i < lstVar.size (); i++) {
      if (!anaEnv.runOnObs ((int) i)){
	if(anaEnv.isMC && anaEnv.nPseudoExp<1) MCObsExp.push_back(new MCObsExpectation (nbins, obs.RangeVariables ()[i].first, obs.RangeVariables ()[i].second, lstVar[i], Luminosity));
	continue;
      }
      if (verbose > 0)
	cout << "For the " << i << " st variable " << lstVar[i] << " with nbins read out from config file... " <<nbins<<" and mapaxis.size....." <<mapAxis.size()<< endl;

      TAxis *axis = mapAxis[lstVar[i]];
      // cout<<"  AND FOR EVERY VARIABLE THE AXIS RANGE  "<<lstVar[i]<<"   "<<axis->GetNbins()<<"  "<<axis->GetXbins ()->fArray<<"  "<<axis->GetXbins ()<<"  "<<axis->GetName()<<endl;
  
      string decision="Bins";
      ////////////////////////////////////
      /// TtJetEstimation
      ////////////////////////////////////
/*      if (verbose > 0)
	cout << " - Configuration of TTJetEstimation ..." << endl;
      ttjEstimation[i] = new TtJetEstimation (anaEnv.isMC);
      //////////////////////////////
      //follow that order !!
      ttjEstimation[i]->SetListOfDatasets (datasetsBkg, datasetsTTJets, datasetsNP);
      //if (axis == NULL)
	ttjEstimation[i]->ConfigHistos (nbins, obs.RangeVariables ()[i].first, obs.RangeVariables ()[i].second, lstVar[i], lstVar[i], true);
*/     
     
      /*else {
	if (axis->GetNbins () == 0)
	  ttjEstimation[i]->ConfigHistos (nbins, obs.RangeVariables ()[i].first, obs.RangeVariables ()[i].second, lstVar[i], lstVar[i], true);
	else {
		if (decision !="Bins"){
	  Double_t *fD = axis->GetXbins ()->fArray;
	  float *fF = new float[axis->GetNbins () +1] ();
	  for (unsigned int j = 0; j <= (unsigned int) axis->GetNbins (); j++)
	    fF[j] = (float) fD[j];
	  ttjEstimation[i]->ConfigHistos (axis->GetNbins (), fF, lstVar[i], lstVar[i], true);
	  delete[]fF;
	}
	
     	else  {	if (decision =="Bins")
	  ttjEstimation[i]->ConfigHistos (nbins, axis->GetXmin(), axis->GetXmax(), lstVar[i], lstVar[i], true);
	 // ttjEstimation[i]->ConfigHistos (nbins, obs.RangeVariables ()[i].first, obs.RangeVariables ()[i].second, lstVar[i], lstVar[i], true);
	}
		
	}


      }
      */
     

      ////////////////////////////////////
      /// Plots used by BkgEstimationSummary
      ////////////////////////////////////
      if (verbose > 0)
	cout << " - Create histos used by BkgEstimationSummary ..." << endl;
      if (axis == NULL) {
	if(anaEnv.isMC && anaEnv.nPseudoExp<1) MCObsExp.push_back(new MCObsExpectation (nbins, obs.RangeVariables ()[i].first, obs.RangeVariables ()[i].second, lstVar[i], Luminosity));
	hData[i] = new TH1F (TString ("hData_") + lstVar[i], "Data", nbins, obs.RangeVariables ()[i].first, obs.RangeVariables ()[i].second);
	hTtJetMC[i] = new TH1F (TString ("TtJetMC_") + lstVar[i], "TtJets - MC", nbins, obs.RangeVariables ()[i].first, obs.RangeVariables ()[i].second);
	hSUSYMC[i] = new TH1F (TString ("SUSYMC_") + lstVar[i], "SUSY - MC", nbins, obs.RangeVariables ()[i].first, obs.RangeVariables ()[i].second);
	hSUSYMCByVar[i] = new TH1F (lstVar[i].c_str(),lstVar[i].c_str() , nbins, obs.RangeVariables ()[i].first, obs.RangeVariables ()[i].second);
        hTtjetMCByVar[i] = new TH1F (lstVar[i].c_str(),lstVar[i].c_str() , nbins, obs.RangeVariables ()[i].first, obs.RangeVariables ()[i].second);

      }
      else {
	if (axis->GetNbins () == 0) {
	  if(anaEnv.isMC && anaEnv.nPseudoExp<1) MCObsExp.push_back(new MCObsExpectation (axis, lstVar[i], Luminosity));
	  hData[i] = new TH1F (TString ("hData_") + lstVar[i], "Data", nbins, obs.RangeVariables ()[i].first, obs.RangeVariables ()[i].second);
	  hTtJetMC[i] = new TH1F (TString ("TtJetMC_") + lstVar[i], "TtJets - MC", nbins, obs.RangeVariables ()[i].first, obs.RangeVariables ()[i].second);
	  hSUSYMC[i] = new TH1F (TString ("SUSYMC_") + lstVar[i], "SUSY - MC", nbins, obs.RangeVariables ()[i].first, obs.RangeVariables ()[i].second);
	  hSUSYMCByVar[i] = new TH1F (lstVar[i].c_str(),lstVar[i].c_str() , nbins, obs.RangeVariables ()[i].first, obs.RangeVariables ()[i].second);
          hTtjetMCByVar[i] = new TH1F (lstVar[i].c_str(),lstVar[i].c_str() , nbins, obs.RangeVariables ()[i].first, obs.RangeVariables ()[i].second);

	}
	else {

	 if(anaEnv.isMC && anaEnv.nPseudoExp<1)
          //MCObsExp.push_back(new MCObsExpectation (axis->GetNbins(),axis->GetXmin(),axis->GetXmax(), lstVar[i], Luminosity));
          MCObsExp.push_back(new MCObsExpectation (axis,lstVar[i], Luminosity));
	  
	 hData[i] = new TH1F (TString ("hData_") + lstVar[i], "Data", axis->GetNbins (), axis->GetXbins ()->fArray);
	 //hData[i] = new TH1F (TString ( "hData_") + lstVar[i], "Data", axis->GetNbins (), axis->GetXmin(), axis->GetXmax());
	 hTtJetMC[i] = new TH1F (TString ("TtJetMC_") + lstVar[i], "TtJets - MC", axis->GetNbins (), axis->GetXbins ()->fArray);
	 hSUSYMC[i] = new TH1F (TString ("SUSYMC_") + lstVar[i], "SUSY - MC", axis->GetNbins (), axis->GetXbins ()->fArray );
	 hSUSYMCByVar[i] = new TH1F (lstVar[i].c_str(),lstVar[i].c_str() ,  axis->GetNbins (), axis->GetXbins ()->fArray);
         hTtjetMCByVar[i] = new TH1F (lstVar[i].c_str(),lstVar[i].c_str() ,  axis->GetNbins (), axis->GetXbins ()->fArray);
	//}
        }
      }
    }				//end of loop over variables

   
    float fillweight;float fillweight_TTJets;
 // MakeBinning NewBins;//(obs.ListOfVariables (), obs.RangeVariables ());
    //////////////////////////////////////////////
    //loop on datasets
    if (verbose > 0)
      cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
    for (unsigned int d = 0; d < datasets.size (); d++) {
      if (verbose > 1)
	cout << "   Dataset " << d << ": " << datasets[d]->Name () <<" with events  "<< datasets[d]->NofEvtsToRunOver ()<<endl;

     
      //open files and load
      ///////////////////////////////////////
      BinContentObsPerSample_.clear();
      treeLoader.LoadDataset (datasets[d], anaEnv); //called second time? Is this needed?

      string dataSetName = datasets[d]->Name();
    /*  
      //regarding the issue with the trigger...
      TObjArray *fileElements= datasets[d]->runTree()->GetListOfFiles();
      TIter next(fileElements);
      TChainElement *chEl=0;
      while (( chEl=(TChainElement*)next() )) {
         //TFile f(chEl->GetTitle());
	 cout<<"               This is... "<<chEl->GetTitle()<<endl;
      }    
     */ 
    
      selecTable.Fill(d,0, datasets[d]->Xsection() * datasets[d]->EquivalentLumi() );

    // int itrigger = treeLoader.iTrigger (datasets[d], string ("HLT_Mu9"));
      ////////////////////////////////////////
      cout << "datasets NormFactor: " << datasets[d]->NormFactor () << endl;

      int nofSelEvts = 0;

      if(anaEnv.nPseudoExp == 0)
      	fillweight = datasets[d]->NormFactor () * Luminosity; //this should be Lumi_data / eqLumi_MCdataset (to be clear, 'Lumi_data' is 'Luminosity' here)
      else
      	fillweight = 1;
      cout<<"for dataset  "<<datasets[d]->Title()<<"  fillweight = "<<fillweight<<endl;
      if (datasets[d]->Name()=="TTJets") fillweight_TTJets = datasets[d]->NormFactor()*Luminosity;

      MCPseudoExp pseudoExp;
      int nPseudos=anaEnv.nPseudoExp;
      int PseudoSessions=anaEnv.nPseudoSession;
      bool manyPseudoExp;
      int itrigger = -1, previousRun = -1;
      string previousFilename("");
      int iFile = -1;
      if (nPseudos <1 && PseudoSessions==0) manyPseudoExp=0;
      if (nPseudos <1 && PseudoSessions!=0) manyPseudoExp=1;
      if (nPseudos >0 && PseudoSessions==0) manyPseudoExp=0;
      if (nPseudos >0 && PseudoSessions!=0) manyPseudoExp=1;
      if (nPseudos >0 && PseudoSessions!=0 && datasets[d]->Name()=="Data") manyPseudoExp=0;

     //to be used when runonTTrees is true
      string filename;
      TFile *f1;
      TTree* seleventtree = 0;
      if(anaEnv.runonTTrees){
        filename = datasets[d]->Title() + "_tree.root";
        f1 = new TFile (filename.c_str (), "READ");
      	seleventtree = (TTree*) f1->Get("OBS");
      }
      //float met = 0;
      //seleventtree->SetBranchAddress("MET",&met);
     //
     
      if(!anaEnv.runonTTrees){
        cout<<"Running on toptrees!"<<endl;
      	pseudoExp.CalculatePseudoExp(datasets[d], Luminosity, ievtmin[d], manyPseudoExp);
      }
      else{
        cout<<"Running on TTrees!"<<endl;
	float eqlumi = datasets[d]->EquivalentLumi();
	pseudoExp.CalculatePseudoExp(seleventtree, Luminosity, eqlumi);
      }

      string dsTitle = datasets[d]->Title (); //now placed here

      PlotObservables Plots (obs.ListOfVariables (), obs.RangeVariables ()); 
      TTreeObservables TtreeObs(obs.ListOfVariables (),dsTitle);//Old way: TTreeObservables TtreeObs(obs.ListOfVariables ());
    
    //MakeBinning NewBins(obs.ListOfVariables (), obs.RangeVariables ());
    
      vector<int> indices;
      if (anaEnv.nPseudoExp==0 || datasets[d]->Name()=="Data" || datasets[d]->Name()=="Data1" || datasets[d]->Name()=="Data2"){
        indices.clear();
     	nEvents[d] = 0;
   
    	if (verbose > 1) cout << "	Loop over events " << endl;

    	int  nevt = datasets[d]->NofEvtsToRunOver ();

  	//cout<<" --------------------------------------------->>>>>>>>>>>>>>>>Will loop for "<<datasets[d]->Name()<<" for   "<<nevt<<"  events...... "<<endl;
  	// if(d==0) nevt=1;  
   	//if (d>0) nevt=2000;
    
    	for (int ii=0; ii<nevt; ii++){ //ii<nevt
        	indices.push_back(ii);
        }
      }
      else
	indices = pseudoExp.GetRandomNumbers();
         
  
  /////////////////////////////////////////////loop on events
      nEvents[d] = 0;
      cout << "	Loop over events, indices.size() =  "<< indices.size() << endl;
      //if (indices.size() >10*Luminosity){
   //  if (datasets[d]->Title()=="TTJets") ievtmin[d]=100000;
   //cout<<" --------------------------------------------->>>>>>>>>>>>>>>>Will loop for "<<datasets[d]->Name()<<" for   "<<ievtmin[d]<<"  events...... "<<endl;
      int mycounter = 0;
      for (vector<int>::iterator iter=indices.begin(); iter!=indices.end(); iter++) {
        int ievt = *iter;
	nEvents[d]++;
	//Reset container
	container.Clear();

 //  TtreeObs (obs.ListOfVariables ());

	TLorentzVector JESbalance (0, 0, 0, 0);

        if (anaEnv.nPseudoExp!=0 && ievt % 5000 == 0)
          cout << "Processing the " << ievt << "th event" << " from "<<indices.size()<<" for Pseudo "<<iPseudoExp<<flush << "\r";
       //  if (ievt %1 ==0)
       //  cout<< "Starting from the "<<ievt<<"th event"<<" from "<<indices.size()<<flush<<"\r";
        if (anaEnv.nPseudoExp==0 && ievt % 5000 == 0)
	  cout << "Processing the " << ievt << "th event" << " from "<<indices.size()<<flush << "\r";
	//cout<<"event "<<ievt<<endl;

	//vectors needed to store variables
        vector < float > ValuesInSR;
	vector < pair  < string , float > > ValuesInSR_WithString;
	ValuesInSR_WithString.clear();
	ValuesInSR.clear();
	
	vector < float > ValuesInNPSamples;
	ValuesInNPSamples.clear();
	vector < pair  < string , float > > ValuesInNPSamples_WithString;
	ValuesInNPSamples_WithString.clear();
	//lstVar=VarsListString;
	
        int Index_;
	Index_ = -1;

	//load event... NOTE: different when running on the toptrees, and on TTrees of variables of selected events
        if(!anaEnv.runonTTrees){
	
      	  event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);
      
	  //cout<<"  datasets[d]->eventTree()->GetFile()->GetName() = "<<datasets[d]->eventTree()->GetFile()->GetName()<<endl;
	  string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
	  if(previousFilename != currentFilename){
	     previousFilename = currentFilename;
	     iFile++;
	     cout<<"   iFile = "<<iFile<<endl;
	  }
	  
	  int currentRun = event->runId();
	  //cout<<", currentRun "<<currentRun<<endl;	  
      	  if(previousRun != currentRun){	  
	     cout<<" currentRun "<<currentRun<<endl;
             previousRun = currentRun;
	     if(currentRun<147196){
	       //cout<<"HLT_Mu9 should be used"<<endl;
	       itrigger = treeLoader.iTrigger (string ("HLT_Mu9"), currentRun, iFile);
	       cout<<" itrigger "<<itrigger<<endl;
	       }
	     else if(currentRun>=147196){
               itrigger = treeLoader.iTrigger (string ("HLT_Mu15_v1"), currentRun, iFile); // HLT_Mu9 is prescaled from a certain run in data...
	       cout<<" itrigger "<<itrigger<<endl;
	     }
	  }

	  //int eventId = event->eventId();
	  ////////////////////////////////////////

	  //apply JES
	  //jets = JESRescale (init_jets, anaEnv.JES);
	  //Calojets = JESRescale (init_Calojets, anaEnv.JES);
	  //PFjets = JESRescale (init_PFjets, anaEnv.JES);

	  /////////////////////////////
	  //   Selection
	  /////////////////////////////
	  //Declare selection instance    
	  //Selection selection (anaEnv.JetType, jets, Calojets, PFjets, init_muons, init_electrons, mets);
	  //selection.SetConfiguration (anaEnv);
          ///////////Selection selection(anaEnv.JetType, init_jets, init_Calojets, init_PFjets, init_muons, init_electrons, mets); //old way
          
	  Selection selection(init_jets, init_muons, init_electrons, mets);
	  //selection.setMuonCuts(20,2.1,0.05,10,0.02,0.3,0,9999,0);//RefSelV3 (like this, the last 3 cuts changed in RefSelV4: selection.setMuonCuts(20,2.1,0.05,10,0.02,0.3,1,1,1) is RefSelV4)
	  //selection.setMuonCuts(20,2.1,0.05,10,0.02,0.3,1,1,1); //RefSelV4
	  //selection.setMuonCuts(20,2.1,9999,10,0.02,0.3,1,1,1); //RefSelV4 without muon reliso cut REMEMBER THAT YOU HAVE PUT THIS HERE!!!!

          // apply preselection to be sure that the same is applied on all the samples...
          unsigned int nPreSelMuons = 0, nPreSelJets = 0;
          /*for(unsigned int i=0; i<init_jets.size(); i++){
             if(init_jets[i]->Pt() > 0 && fabs(init_jets[i]->Eta()) < 99) //pt 20, eta 2.4?? (10 in stead of 99??)
                nPreSelJets++;
          }      */
          for(unsigned int i=0; i<init_muons.size(); i++){
	 //	cout<<"  init muons Pt  "<<init_muons[i].Pt()<<endl;
             if(init_muons[i]->Pt() > 20 && fabs(init_muons[i]->Eta()) < 2.1){ //pt 15 (Stijn, only this...)? pt 20? (me)
                nPreSelMuons++;//	cout<<"  init muons Pt  "<<init_muons[i].Pt()<<endl;
	     }            
          }      
          /*if(nPreSelJets < 4 || nPreSelMuons < 1)*/
	  if(nPreSelMuons < 1)
             continue; // event not preselected!
      
  //      if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
  //      {
  //        bool passedJSON = treeLoader.EventPassedJSON(datasets[d], event->runId(), event->lumiBlockId());
  //        if(!passedJSON) continue;
  //      }


	//from TTbarJES.cc
        if(datasets[d]->Name() == "Data" || datasets[d]->Name() == "data" || datasets[d]->Name() == "DATA"){
           // Apply the scraping veto
           bool isBeamBG = true;
           if(event->nTracks() > 10){
              if( ( (float) event->nHighPurityTracks() ) / ( (float) event->nTracks() ) > 0.25 )
                isBeamBG = false;
           }
           if(isBeamBG) continue;       
           // Apply the JSON
           // bool passedJSON = treeLoader.EventPassedJSON(datasets[d], event->runId(), event->lumiBlockId());
          // if(!passedJSON) continue;
        }





//uncomment this below to temporary not use trigger!!!, because problem with trigger when you comma-seperate data input filenames (now +- solved, maybe not in the best way)
          bool trigged = treeLoader.EventTrigged (itrigger);
	  ////bool trigged = true;
	  
          bool isGoodPV = false;
          /*if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") // temporary fix since there is a different cut on Data and MC!
             isGoodPV = selection.isPVSelected(vertex, anaEnv.PVertexNdofCut, 24, anaEnv.PVertexRhoCut);
          else
             isGoodPV = selection.isPVSelected(vertex, anaEnv.PVertexNdofCut, anaEnv.PVertexZCut, anaEnv.PVertexRhoCut);
          */
	  isGoodPV = selection.isPVSelected(vertex, 4,24,2.); //this?? from SyncRefSel.cc
	                      
	  
	  vector<TRootJet*> selectedJets;
	  vector<TRootMuon*> selectedMuons;
	   
	  bool doPF2PAT = false; // not supported/synced atm...

	  if (init_jets.size() > 0) {
	      if (init_jets[0]->jetType() == 1 || doPF2PAT) { // calojets
		//cout << "Selecting for caloJets" << endl;
		selectedJets = selection.GetSelectedJets(true);
		selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
	      }	else {
		//cout << "Selecting for PF/JPT jets" << endl;
		vector<TRootMuon*> overlapMuons = selection.GetSelectedMuons(vertex[0]);
		//selection.setJetCuts(30.,2.4,0.01,1.,0.98,0.3,0.1); // refSelV4 values
		selectedJets = selection.GetSelectedJets(overlapMuons,true);
		selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);	
	      }
	  }
	  vector<TRootMuon*> vetoMuons = selection.GetSelectedLooseMuons();
	  vector<TRootElectron*> vetoElectrons = selection.GetSelectedLooseElectrons(false);

          bool eventSelected = false;
	  //////////////////////////
	  //the different regions
	  //////////////////////////
	  bool isInCR = false;
	  bool isInSR = false;
	  float RelIso = 100.;
	  //////////////////////////


          if(trigged){       
	     ///////////////////////////////////////////
	     //fill booleans for the different region
	     ///////////////////////////////////////////
	
	     if (isGoodPV && selectedMuons.size () == 1 && vetoMuons.size () == 1 && vetoElectrons.size () == 0 && selectedJets.size () >= (unsigned int) anaEnv.NofJets) {
	//      if(mets[0]->Et()>=40){ //just extra cut... REMEMBER THAT YOU HAVE PUT THIS HERE!!!!
	  //////////     RelIso = 1 / selectedMuons[0]->relativeIso03 (); //why?
	       RelIso = selectedMuons[0]->relativeIso03 ();
	       isInSR = true;
	       // cout<<"selectedMuons[0].Pt() = "<<selectedMuons[0].Pt()<<endl;	      
	       mycounter++;	       
	//      }
	     }
	     /*if (isGoodPV && selection.GetSelectedMuonsInvIso (anaEnv.MuonPtCutSR, anaEnv.MuonEtaCutSR, anaEnv.MuonRelIsoCutSR).size () == 1 && vetoElectrons.size () == 0 && selectedJets.size () >= (unsigned int) anaEnv.NofJets) {
	       if (!isInSR)
	         RelIso = 1 / selection.GetSelectedMuonsInvIso (anaEnv.MuonPtCutSR, anaEnv.MuonEtaCutSR, anaEnv.MuonRelIsoCutSR)[0].relativeIso03 ();
	     } */     
          }
	   	
	  ///////////////////////
	  //LOOP OVER VARIABLES
	  for (unsigned int v = 0; v < lstVar.size (); v++) {
	     if (!anaEnv.runOnObs ((int) v))
	        continue;
	     //variable on which on want to make bgk estimation  & search for New Physics
	     float variable = -999;
	     string variable_string;
	     //fill the variable  (take into account if it's in SR of CR(s) ...)
	     //to be true, a jet combination should be performed first. Here jets are ranked by Pt
	     if (isInSR) {	        
	        ChiSq = ChiSquareMatching (selectedMuons[0], selectedJets, NjetsForComb, Permutation, WMass, TopMass, blMass, HadWMassResol, HadTopMassResol, LepblMassResol, DEBUG);
	        CombinedJets.clear ();
	   
 	        if (selectedJets.size () >= 4)
	          for (Int_t i = 0; i < 4; i++) {
		     CombinedJets.push_back (selectedJets[Permutation[i]]);
	          }
	        else
	          CombinedJets = selectedJets;
	        string a ("isInSR");
		
	        //Observables obsEvt (selectedMuons[0], CombinedJets, selectedJets, mets[0], a, ChiSq);
		
		//now other arguments needed... this may be a temporary workaround in stead of adapting the whole observables class...
		vector<TRootJet> CombinedJets_nopointers;
		vector<TRootJet> selectedJets_nopointers;	
		for(unsigned int mmm=0;mmm<CombinedJets.size();mmm++){
			CombinedJets_nopointers.push_back(*(CombinedJets[mmm]));
		}
		for(unsigned int mmm=0;mmm<CombinedJets.size();mmm++){
			selectedJets_nopointers.push_back(*(selectedJets[mmm]));
		}	
		Observables obsEvt (*(selectedMuons[0]), CombinedJets_nopointers, selectedJets_nopointers, *(mets[0]), a, ChiSq);

	        for (unsigned int kkk=0;kkk<obsEvt.Variables().size();kkk++){	  
      		  if ( lstVar[v]== obsEvt.Variables()[kkk].first) { 	     
		     //cout<<"FOUND THE CORRECT ONEEE!!!!!!!!!    "<<lstVar[v]<< "------------->>>>>>>>>>>>>>>>  "<< obsEvt.Variables ()[kkk].second<<"  " <<obsEvt.Variables ()[kkk].first<<"   "<<kkk<<endl;
		     Index_=kkk;
		  }
	  	}
	        //Index_=v;

		if (anaEnv.MCRound==0)
		  BinContentObsAll_.push_back( pair< string,float > (lstVar[v],obsEvt.Variables ()[Index_].second));
		  
	        BinContentObsPerSample_.push_back( pair< string,float > (lstVar[v],obsEvt.Variables ()[Index_].second));
	////        if(lstVar[v]=="MET"){}
		
		bool fill;fill=0;
	        if (v==lstVar.size()-1) 
		  fill=1;
                TtreeObs.Fill (obsEvt,fill);
                Plots.Fill(obsEvt,fillweight,fill);
	        variable = obsEvt.Variables ()[Index_].second;
                variable_string = lstVar[v];
	        if (anaEnv.nPseudoExp!=0) {
	    	  container.Add(variable);
	  	}	
             }//end if(isInSR)

	     ///////////////////////
	     //  TtJetEstimation
	     ///////////////////////
	////     ttjEstimation[v]->Fill (variable, trigged, selection, anaEnv, isInSR, datasets[d]->Name ());
	     ///////////////////////

	     if (isInSR) {       
	        ValuesInSR.push_back (variable);	    
	        ValuesInSR_WithString.push_back (pair < string, float > (variable_string,variable));
	     
               //NewBins.Fill(variable_string, variable);
	       //cout<<" WILL STORE NOW-------------> "<<variable_string<<"  "<<variable<<"    "<<num_<<endl;
	        if (anaEnv.nPseudoExp>0) {
                  for (int k=0;k<VarsListString.size();++k){
                     if ( VarsListString[k]==lstVar[v])
                       ValuesInNPSamples.push_back(variable);
		     ValuesInNPSamples_WithString.push_back(pair <string,float> (variable_string,variable));
                       break;
		  }
	        }
	    
	        if (anaEnv.isMC)
	          if(anaEnv.isMC && anaEnv.nPseudoExp<1){
		     MCObsExp[v]->Fill (datasets[d], variable, fillweight);
	          }
	       //Data = all datasets
	        hData[v]->Fill (variable, fillweight);
	        if (datasets[d]->Name () == "TTJets"){
	           hTtJetMC[v]->Fill (variable, fillweight);
	           hTtjetMCByVar[v]->Fill (variable, fillweight);
	        }
	        if (datasets[d]->Name () == "SUSY"){
	           hSUSYMC[v]->Fill (variable, fillweight);  hSUSYMCByVar[v]->Fill (variable, fillweight);
                   TH1F *htemp ;   htemp= hData[v];
                   tree->Fill();
	        }
  
  	        string dstitle=datasets[d]->Title();
		if(anaEnv.isMC == 0){
		   hRelIsoData->Fill(RelIso);
		}
		if(anaEnv.isMC == 1){
		   myMultiSamplePlot.Fill(RelIso,datasets[d],true,Luminosity);
		}

	     }/////end of isInSR	
    
          /////////
	  }			//loop over variables

          //TtreeObs.FillTtree();

          // TtreeObs.Fill (BinContentObsPerSample_);

	  if(isInSR && anaEnv.nPseudoExp==0) {
	     VarValues.push_back (ValuesInSR);
             VarValues_WithString.push_back(ValuesInSR_WithString);
	     //ValuesInSR_WithString.push_back (pair < string, float > (variable_string,variable));
	     //push back once per event (in SR); actually, one such 'ValuesInSR' vector IS an event (in SR) in some sense
      	     NPbooleans.push_back(isNP); 
	  }
          if (isInSR && anaEnv.nPseudoExp>0){ 
	     VarValues.push_back(ValuesInNPSamples);  
             VarValues_WithString.push_back(ValuesInNPSamples_WithString);
             NPbooleans.push_back(isNP);
	  }	
          //fill the TClonesArray
	  if(isInSR){
	     new ((*tcContainer)[itContainer]) Container (container);
	     nofSelEvts++;
	  }
	} //end if(!anaEnv.runonTTrees)
	
	//when running on TTrees rather than toptrees...


        else if(anaEnv.runonTTrees){	  
	  ///////////////////////
	  //LOOP OVER VARIABLES
	  int nvar = lstVar.size ();
	  float variable[nvar]; //TEMPORARY?
	  string variable_string;
	  for (unsigned int v = 0; v < lstVar.size (); v++) {	     
	               variable_string = lstVar[v].c_str();
		  seleventtree->SetBranchAddress(variable_string.c_str(),&variable[v]); //is this the way to do it? doesn't seem to work	     
	  }
	  
	  seleventtree->GetEntry(ievt);
	  
	  for (unsigned int v = 0; v < lstVar.size (); v++) {
	     variable_string = lstVar[v];
              BinContentObsPerSample_.push_back( pair< string,float > (lstVar[v],variable[v]));
	  
	    
	     ValuesInSR_WithString.push_back(pair < string, float > (variable_string,variable[v]));
	  

	     hData[v]->Fill (variable[v], fillweight);
    
	  
	  }/////////////end of looping for variables

	  
	  VarValues_WithString.push_back(ValuesInSR_WithString);
	  isNP = 0;
	  NPbooleans.push_back(isNP);
	  nofSelEvts++;   
	
	  
	}/////////////////////end of runonTTrees



	
      }	////////////loop on 'indices' ('events')
      cout<<"mycounter = "<<mycounter<<endl;
      
      cout<<"XS/PreSelEff/nEvents : " << datasets[d]->Xsection () << "/" << datasets[d]->PreSelEfficiency () << "/" << nEvents[d] << endl;
      //cout<<"pseudoExp.GetRandomNumbers().size() = "<<pseudoExp.GetRandomNumbers().size()<<endl;
      cout<<" VarValues.size() = "<<VarValues.size()<<endl;
      cout<<" NPbooleans.size() = "<<NPbooleans.size()<<" and total BinContentObs size  "<<BinContentObsAll_.size()<<"  and for this dataset  "<<BinContentObsPerSample_.size()<<endl;
      //string dsTitle = datasets[d]->Title ();
   
      //if( anaEnv.nPseudoExp==0)      
      if( anaEnv.nPseudoExp==0 || dsTitle=="Data" || dsTitle=="Data1" || dsTitle=="Data2") //(z)!!
        TtreeObs.Write (true); ////////////////write all tree OBS for each dataset
    	
      sort( BinContentObsAll_.begin(), BinContentObsAll_.end());	
      sort( BinContentObsPerSample_.begin(), BinContentObsPerSample_.end());
	
      string bin_decision = "Bins";
      if(d>0 && anaEnv.nPseudoExp==0){   
        cout<<"  Will try to create and fill .... "<<dsTitle<<"   "<<BinContentObsPerSample_.size()<<endl;
        TString dss= datasets[d]->Title();
        if (anaEnv.MCRound==1 && bin_decision == "Bins")   Plots.WriteNBins(dss,BinContentObsPerSample_,anaEnv.nbins,fillweight,true);
        if (anaEnv.MCRound==1 && bin_decision == "Events")     Plots.WriteEvtsPerBin(dss,BinContentObsPerSample_,anaEnv.eventsperbin,fillweight,true);
      }

      BinContentObsPerSample_.clear();

      indices.clear();
      cout<<"   DONE FOR DATASET........................ "<<dsTitle<<"   and the total selected events.... " <<BinContentObsAll_.size()<<endl;
      cout<<endl;
 //  treeLoader.UnLoadDataset();
 
      nofSelEvts_vector.push_back(nofSelEvts);
      cout<<"   real number of selected events (inSR) of this dataset: "<<nofSelEvts_vector[d]<<endl<<endl;
    } //////loop on datasets
      
      
   //regarding binning   
    if (anaEnv.MCRound==0){	
      cout<<"  Will try to create the Binning  .... "<<"   "<<BinContentObsAll_.size()<<"  "<<EventsPerBin<<"  "<<endl;
      string bin_decision = "Bins";	
      MakeBinning NewBins;
      cout<<"HERE 1"<<endl;
      cout<<"BinContentObsAll_.size() = "<<BinContentObsAll_.size()<<endl;
      NewBins.Binning(bin_decision,BinContentObsAll_,anaEnv.nbins,EventsPerBin,fillweight_TTJets);
      cout<<"HERE 2"<<endl;
    } 
    BinContentObsAll_.clear();
    BinContentObsPerSample_.clear();
   
   //Once everything is filled ...
    if (verbose > 0)
      cout << " We ran over all the data ;-)" << endl;
 
   //Selection tables
    selecTable.TableCalculator(false, true, true, true);
    string selectiontable = "SelectionTable_JES";
    if (argc >= 3){
      string sample=string(argv[2]);
      selectiontable = selectiontable +"_"+sample;
    }
    selectiontable = selectiontable +".tex"; 	
    selecTable.Write(selectiontable.c_str());
     
     
   string lstChosenVar[lstVar.size ()]; //needed for the variable names in a later loop over ObsData; only related to the subdirectory structure of the Stat output file
   int kk =0;
  //LOOP OVER VARIABLES
   cout<<"lstVar.size () = "<<lstVar.size ()<<endl;
   for (unsigned int v = 0; v < lstVar.size (); v++) {
      if (!anaEnv.runOnObs ((int) v))
	continue;

      /////////////////////////////////
      //  Summary plot
      /////////////////////////////////
      if (verbose > 0)
	cout << " - Summary plots ..." << endl;

      //       cout << "  CREATED THE OBSESTIMATED !!!!!!!!!!!!! --------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> "<<endl;
      cout<<"MCObsExp.size() = "<<MCObsExp.size()<<endl;
      for (int l=0;l<MCObsExp.size();l++){
	if (MCObsExp[l]->GetHistoSMProcesses()->GetName()=="SMProcess_"+lstVar[v]){
	   cout<<"  found MCObsExp plots....... "<<MCObsExp[l]->GetHistoSMProcesses()->GetName()<<"   "<<l<<endl;
	   ObsEstimated.push_back(MCObsExp[l]->GetHistoSMProcesses());
	   break;
	}
      }
      //      ObsEstimated.push_back(MCObsExp[v]->GetHistoSMProcesses());
      
      ObsData.push_back (hData[v]);
      ObsSUSY.push_back (hSUSYMC[v]);  
      ObsSUSYByVar.push_back (hSUSYMCByVar[v]);
      ObsTtjetByVar.push_back(hTtjetMCByVar[v]);
      lstChosenVar[kk] = lstVar[v];
    
      kk++;
      /////////////////////////////////////////
   }	//end of loop on variables
   
   //writing observable plots in Stat Procedure output file
   if(anaEnv.nPseudoExp > 0){
      foutStat->cd();
      for(unsigned int i=0;i<ObsData.size();i++){
	 mytdirall->cd(lstChosenVar[i].c_str());
	 TH1F hdatatemp;
	 TH1F hestimtemp;
         // TH1F htemp;
         // htemp = *ObsSUSY[i];
         //htemp.Write();
	 hdatatemp = *ObsData[i];
	 hdatatemp.Write();  
	 hestimtemp = *ObsEstimated[i];
	 hestimtemp.Write(); //}// always the same when MC expectation as estimation, but always different when data-driven estimation...
       
         /*TH1F hestimtemp_StackSM; // SM processes stacked
         hestimtemp_StackSM = MCObsExp[i]->GetTHStackSMProcesses();
	 hestimtemp_StackSM.Write();
	 TH1F hestimtemp_hW_LF; // W_LF
	 hestimtemp_hW_LF = *(MCObsExp[i]->Get_hW_LF());
	 hestimtemp_hW_LF.Write();*/
       }
   }
 
 
  ///Here, the ObservablesRanker will be intialized and run...
   bool over_table;bool correl_table;bool over_cleaned_sorted_table,TTjets_first_loop;
   over_table=true;
   correl_table=true;
   over_cleaned_sorted_table= true;
   TTjets_first_loop=false;
   string merged_file_ ;
   string sm_file_ ;
   string np_file_;
   char sgnl_file_char_[100];
   char bkg_file_char_[100];
  
   TTree *tSignal, *tBkg; TFile *hfileS, *hfileB;
   if(anaEnv.nPseudoExp == 0){
     //ObservablesRanker ObsR;    
      if (anaEnv.MCRound==2){
	 for (unsigned int j = 0; j < datasets.size (); j++) {
           if ( j>0 ) TTjets_first_loop=true;
 	   if (datasets[j]->Name().find("LM1")<=datasets[j]->Name().size() || datasets[j]->Name().find("Zp")<=datasets[j]->Name().size() || datasets[j]->Name().find("NP")<=datasets[j]->Name().size() ||  datasets[j]->Name().find("SUSY")<=datasets[j]->Name().size() ||  datasets[j]->Name().find("Data")<=datasets[j]->Name().size() ||  datasets[j]->Name().find("Data1")<=datasets[j]->Name().size() ||  datasets[j]->Name().find("Data2")<=datasets[j]->Name().size()  ){

	     merged_file_ = "Merged_Plots_" + datasets[1]->Title () + "_" + datasets[j]->Title () + ".root";
    	     cout << " Will merge and compute overlap for " << datasets[1]->Title () << "  and  " << datasets[j]->Title () << "  file will be " << merged_file_.c_str () << endl;

             np_file_ = datasets[j]->Title ();  
    	     sm_file_ = datasets[1]->Title();
             sprintf (bkg_file_char_, "%s_tree.root", datasets[1]->Title ().c_str ());  
             sprintf (sgnl_file_char_, "%s_tree.root", datasets[j]->Title ().c_str ());  
             hfileS = new TFile (sgnl_file_char_);
             hfileB = new TFile (bkg_file_char_);   
             cout << "Filename for signal: " << sgnl_file_char_ << "  - Filename for background: " << bkg_file_char_ << endl;  
             tSignal = (TTree *) hfileS->Get ("OBS");
             tBkg = (TTree *) hfileB->Get ("OBS");
    
             ObservablesRanker ObsR(tSignal, tBkg, merged_file_, np_file_, sm_file_, over_table,correl_table,over_cleaned_sorted_table, TTjets_first_loop,Correl_Cut_);
 	   }//end of specifically looking for NP datasets
         }///end of looping inside datasets for rankings
      }///end of MCRound requirement
   }///end of Pseudo==0
  /////////end of ObservablesRanker block.....


  ///////////////////
  // Writing
  //////////////////
   if (verbose > 1)
      cout << " - Start writing the module outputs in the ouput file ..." << endl;
   char nom[200] = "";
   sprintf (nom, "_exp%d_", iPseudoExp);

   //loop on variables
   int it = 0;
   for (unsigned int v = 0; v < lstVar.size (); v++) {
      if (anaEnv.runOnObs ((int) v)) {
	 sprintf (nom, "%s_exp%d_", lstVar[v].c_str (), iPseudoExp);
	 if (anaEnv.isMC) {
	   MCObsExp[v]->SetColors (datasets);
	   MCObsExp[v]->Compute ();
	   new ((*tcMCObsExp)[it]) MCObsExpectation (*MCObsExp[v]);
	 }
	 it++;
      }
   }


    //Run the statistical procedure
    // -> Loop over the variables
    // -> use only the selected one
    // -> Inputs: "2 histo" 
    //            - "data" = sum of all process  "MCObsExp"
    //            - "estimation" = sum of bkg estimated  "bkgSum"

  ///////////////////////////////
  ///  STATISTICAL PROCEDURE
  //////////////////////////////
   if(anaEnv.nPseudoExp > 0){
     if (verbose > 0)
        cout << " - Statistical procedure ..." << endl;
     vector<float> Svect;
     SampleCombinedWeightCalculator mySampleCombinedWeightCalculator;////// For the V-value
     if(verbose>0) 
        cout << endl << "Starting calculation of weights..." << endl;
     // LOOP OVER SELECTED EVENTS

     int nVar = VarValues_WithString.size();
     //cout<<"lstVar size "<<lstVar.size()<<" VarValues_WithString.size() = "<<VarValues_WithString.size()<<" VarValues "<<VarValues.size()<<"   obsData  "<<ObsData.size()<<" obsEstimated  "<<ObsEstimated.size()<<endl;
     for (int k=0;k<ObsData.size();++k){
  	cout<<ObsData[k]->GetName()<<"  "<<ObsEstimated[k]->GetName()<<"   "<<k<<endl;
     }

     VarValues.clear();

     vector <float> dummy;
     dummy.clear();

     sort(VarValues_WithString.begin(),VarValues_WithString.end());

     for(unsigned int i=0;i<VarValues_WithString.size();i++){
	for ( int ii=0;ii<VarValues_WithString[i].size();ii++){
//cout<<"  Varvalues......"<<i<<"   "<<VarValues_WithString[i][ii].second<<endl;
	   
	    dummy.push_back(VarValues_WithString[i][ii].second);
        }
        VarValues.push_back(dummy);
     }
	 
     for(unsigned int i=0;i<VarValues_WithString.size();i++){
	/*for(unsigned int j=0; j<nVar;j++){ 
		cout << " VarValues["<<i<<"]["<<j<<"]: "<< VarValues[i][j];
	}*/
  //cout<<"lstVar size "<<lstVar.size()<<" VarValues_WithString.size() = "<<VarValues_WithString.size()<<" VarValues "<<VarValues.size()<<"   obsData  "<<ObsData.size()<<" obsEstimated  "<<ObsEstimated.size()<<endl;
        EventCombinedWeightCalculator myEventCombinedWeightCalculator(SignificanceCombination);////// For the squared - significance and S-product
   	myEventCombinedWeightCalculator.CalculateWeight(VarValues_WithString[i],ObsData,ObsEstimated);
    //	myEventCombinedWeightCalculator.CalculateWeight(VarValues[i],ObsData,ObsEstimated,sensitive_Obs); //recall, the sensitive_Obs is redundant (but was implemented in the beginning) and should be an empty vector

	if(!myEventCombinedWeightCalculator.isUnderflowWarning()){
	    Svect.push_back(myEventCombinedWeightCalculator.GetCombinedWeight());
	    if(verbose>3) cout << "Weight event " << i << ": S = " << Svect[i] << endl;
	}
	else{
	    cout<<"WARNING: myEventCombinedWeightCalculator.isUnderflowWarning() = true. This problem has to be solved!"<<endl;
	    Svect.push_back(-9999);  // or -999? By setting it to 0 or a negative number, these events will not end up in the highest-weight subsample 
	}
     }


     if(Svect.size()!=NPbooleans.size()) cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                  ---------------->>>>>>>>>>>> WARNING: vectors Svect and NPbooleans not equal size!"<<endl;
   
     vector< eventinfo > Svect_eventsinfos; //the eventinfo struct is defined in SampleCombinedWeightCalculator.h
     eventinfo w;

     for(unsigned int i=0;i<Svect.size();i++){
   	w.weight = Svect[i];
	w.isnp = NPbooleans[i];
	w.variableValues = &VarValues[i];
	Svect_eventsinfos.push_back(w);
     }
     foutStat->cd();
     for(unsigned int k=0;k<nofFractionsHWEvts;k++){
        mySampleCombinedWeightCalculator.CalculateV(Svect_eventsinfos, Combination, FractionHWEvts[k]);   // after this step, Svect_eventsinfos is sorted (according to the weigths), and the V value is calculated
        ostringstream Fractionstrstream;
        Fractionstrstream << FractionHWEvts[k];
        mytdirFraction = foutStat->GetDirectory("x="+TString(Fractionstrstream.str()));
        mytdirFraction->cd();
	
        //writing highest-weight events histograms			
        for(unsigned int i=0;i<ObsData.size();i++){  // why does this in an if(ips==0) loop give no such histograms written in the file
	
           TH1F *hdatatemp=new TH1F(ObsData[i]->GetName()+TString("_highestweights_fr")+TString(Fractionstrstream.str()), ObsData[i]->GetTitle()+TString(" (highest-weight subsample [x = ")+TString(Fractionstrstream.str())+TString("])"), ObsData[i]->GetNbinsX(), ObsData[i]->GetXaxis()->GetXbins ()->fArray);
           for(unsigned int j=0; j<mySampleCombinedWeightCalculator.GetnEventsInHighestWeightSubsample();j++){
	  	hdatatemp->Fill((*mySampleCombinedWeightCalculator.GetEventsInfos()[j].variableValues)[i]);
	   }
	   //cout<< "NAME OF CURRENT DIRECTORY: " << gDirectory->GetName() <<endl;
           mytdirFraction->cd(lstChosenVar[i].c_str());
	   hdatatemp->Write();
	   delete hdatatemp;
        }
        if(verbose>0) cout<<"FractionHWEvts["<<k<<"]: "<<FractionHWEvts[k]<<endl;
        if(verbose>0) cout << "  V = " << mySampleCombinedWeightCalculator.GetV() << endl << endl;
        Vvect.push_back(pair< float , float >(FractionHWEvts[k],mySampleCombinedWeightCalculator.GetV()));
        float fracNP = -9999;
        fracNP = mySampleCombinedWeightCalculator.GetFracNPEventsInHighestWeightSubsample();
        FracNPEventsInHWSubsampleVect.push_back(pair< float , float >(FractionHWEvts[k],fracNP));  
     }  //end loop on fraction
 
     if(verbose>0) cout<<"   Number of events in final sample = "<<mySampleCombinedWeightCalculator.GetnEventsInFinalSample()<<endl; //final sample = sample in SR
     if(verbose>0) cout<<"   Number of NP events in final sample = "<<mySampleCombinedWeightCalculator.GetnNPEventsInFinalSample()<<endl;
     
    // loop to free memory allocated to hData //worked (and was necessary) in 22X version, should be checked in 35X version
     if(verbose>0) cout<<"Releasing memory of hData..."<<endl; 
     for (unsigned int v = 0; v < lstVar.size();v++){//obs.ListOfVariables ().size (); v++) {
        if(!anaEnv.runOnObs((int)v)) continue; 
        if(hData[v]) delete hData[v];  // if(hData[v]!=NULL)?
        if(hTtJetMC[v]) delete hTtJetMC[v];  // if(hData[v]!=NULL)?
        if(hSUSYMC[v]) delete hSUSYMC[v];  // if(hData[v]!=NULL)?
        if(hSUSYMCByVar[v]) delete hSUSYMCByVar[v];  // if(hData[v]!=NULL)?
        if(hTtjetMCByVar[v]) delete hTtjetMCByVar[v];  // if(hData[v]!=NULL)?
     }
   } //end if(anaEnv.nPseudoExp > 0) statement
    
    
   if(anaEnv.doDumpPseudoExpInfoInTTree) PseudoExpDumpTree->Fill(); 
     
   delete MCExp;
   MCExp = NULL;
    

   ObsEstimated.clear();
   ObsData.clear();
   ObsSUSY.clear();
   ObsSUSYByVar.clear();
   ObsTtjetByVar.clear();
   cout << "** End pseudoexperiment  **" << endl << endl;
 }				// loop over pseudo-experiments
 
 ///////////////////////////////
 ///  STATISTICAL PROCEDURE (to write V values)
 //////////////////////////////
 if(anaEnv.nPseudoExp > 0){
  //in principle calculations related to V_alfa can or even should be (re)done in a seperate macro, because all the V values will be stored
   WeightProbaCalculator myWeightProbaCalculator;
   float alfa= 0.1;
  //writing this vector to a TFile... in TTree format for the moment; can be changed in the future
   cout<<"Writing V vector in TTree format to TFile..."<<endl;  // experimenting with trees to write and read a collection of numbers (the V values)!
   TFile treefile(rootStatVtreeFileName.c_str(),"RECREATE");
   float Vtemp, Valfatemp;
   for(int k=0;k<nofFractionsHWEvts;k++){
     ostringstream Fractionstrstream, alfastrstream;
     Fractionstrstream << FractionHWEvts[k];
     TTree *Vtree = new TTree("Vtree_fr"+TString(Fractionstrstream.str()),"Tree of V values [fr="+TString(Fractionstrstream.str())+"]");
     Vtree->Branch("V",&Vtemp,"Vleaf");
     vector<float> VvectCurrentFraction;
     for(int i=0;i<Vvect.size();i++){
	cout<<"Vvect  "<<Vvect[i].first<<"  "<<Vvect[i].second<<"   "<<Vvect.size()<<endl;
	if(Vvect[i].first == FractionHWEvts[k]){
	   VvectCurrentFraction.push_back(Vvect[i].second);
	}
     }

     sort(VvectCurrentFraction.begin(),VvectCurrentFraction.end());
	
     for(int i=0;i<VvectCurrentFraction.size();i++){
	//cout<<"   V values output  "<<VvectCurrentFraction[i]<<"  and fraction  "<<(k+1)*0.1<<endl;
	Vtemp = VvectCurrentFraction[i];
	Vtree->Fill();
     }	
     Vtree->Write();
     delete Vtree;
     //new "tree"... (maybe just one number per tree...)
     TTree *Valfatree = new TTree("Valfatree_fr"+TString(Fractionstrstream.str()),"Tree with Valfa values [fr="+TString(Fractionstrstream.str())+"]");
     alfastrstream << alfa;
     Valfatree->Branch("Valfa"+TString(alfastrstream.str()),&Valfatemp,"Valfa"+TString(alfastrstream.str())+"leaf");
     Valfatemp = myWeightProbaCalculator.GetProbaPercentile(VvectCurrentFraction, alfa*100);
     cout<<"Valfatemp [x="<<FractionHWEvts[k]<<"] = "<<Valfatemp<<endl;
     Valfatree->Fill();
     Valfatree->Write();
   }
  // duplicate each step for the AU = 15% on estimation case (new trees/branches, but same file)
  //was old way, but will be abandonded
   
   treefile.Close();
 
  //related to filling of histograms of fractions of NP events in highest weight subsample etc...
   cout<<"Updating output file"<<endl;
   foutStat->cd();
   for(int k=0;k<nofFractionsHWEvts;k++){
     ostringstream Fractionstrstream;
     Fractionstrstream << FractionHWEvts[k];
     TH1F* FracNPHisto = new TH1F("Fraction of events in highest weight subsample that are NP "+TString(Fractionstrstream.str()),"Fraction of events in highest weight subsample that are NP "+TString(Fractionstrstream.str()),100,0,1);
     FracNPHisto->GetXaxis()->SetTitle("fraction NP events");
     for(int i=0;i<FracNPEventsInHWSubsampleVect.size();i++){
	if(FracNPEventsInHWSubsampleVect[i].first == FractionHWEvts[k]){
	   FracNPHisto->Fill(FracNPEventsInHWSubsampleVect[i].second);
	}
     }
     mytdirFraction = foutStat->GetDirectory("x="+TString(Fractionstrstream.str()));
     mytdirFraction->cd();
     FracNPHisto->Write();
   }
 // foutStat->Close();
 } //end if(anaEnv.nPseudoExp > 0) statement
 
 ///////////////////
 //  Write the TTree    //not really used
 //////////////////
 if(anaEnv.doDumpPseudoExpInfoInTTree){
   TFile treeFile(anaEnv.DumpTreeName.c_str(),"RECREATE");
   treeFile.cd();
   PseudoExpDumpTree->Write();
   treeFile.Write();
   treeFile.Close();
 }
  
  /////////////////////////////////////


  ///////////////////
  // Writing
  //////////////////
  if(anaEnv.nPseudoExp > 0){
    if (verbose > 1)
      cout << " - Writing configTree in Stat output file ..." << endl;
    //add configuration
    foutStat->cd ();
    configTree->Fill ();
    configTree->Write ();
    foutStat->Close();
  }
  
  if (verbose > 1)
    cout << " - Writing outputs on files ..." << endl;
  //add configuration
  fout->cd ();
  configTree->Fill ();
  configTree->Write ();
  //
  if (verbose > 1)
    cout << " - Writing  the file ..." << endl;
  fout->Write ();
/*
  if(anaEnv.isMC == 0){
  	hRelIsoData->SaveAs("RelIso_Data.root");
  }
  if(anaEnv.isMC == 1){
  	TFile* Outputfile_ExtraVariables = new TFile("Outputfile_ExtraVariables.root","RECREATE");
	//  TFile *fData = TFile::Open("dummy.root");
	//  myMultiSamplePlot.AddDataHisto(DataHisto);
  	string label = "RelIso";
  	myMultiSamplePlot.Draw(false,label,false,false,false,false,false);
  	myMultiSamplePlot.Write(Outputfile_ExtraVariables,label,false,"");
  }
*/  
  //delete
  delete event;
  delete configTree;
  delete tcdatasets;
  delete tcAnaEnv;
  delete tcMCObsExp;
  delete tcMCExpBjetMult;


  cout << "**********************************************************************" << endl;
  cout << "           End of the program !!" << endl;
  cout << "**********************************************************************" << endl;
  return 0;
}




Float_t
ChiSquareMatching (TRootMuon* & muon, std::vector < TRootJet* > &jets, UInt_t & n, UInt_t * Permutation, Float_t & WMass_, Float_t & TopMass_, Float_t & blMass_, Float_t & HadWMassResol_, Float_t & HadTopMassResol_, Float_t & LepblMassResol_,
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
      HadWCandMass_tmp = (*jets[comb[0]] + *jets[comb[1]]).M ();
      HadTopCandMass_tmp = (*jets[comb[0]] + *jets[comb[1]] + *jets[comb[2]]).M ();
      LepBlCandMass_tmp = (*jets[comb[3]] + *muon).M ();

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




