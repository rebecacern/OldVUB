#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
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

#include "Style.C"

using namespace std;
using namespace TopTree;

void Scale (vector < Dataset > datasets, double *&numbers, float Luminosity)
{
  for (unsigned int i = 0; i < datasets.size (); i++) {
    numbers[i] = numbers[i] * datasets[i].NormFactor () * Luminosity;
  }
}

// the following function return a vector of TRoot*Jet with P4 changed according to JES and with the main members required fill
// Improvement: create a method in TRootJet which scale the P4 directly

vector<TRootJet> JESRescale(vector<TRootJet> &jets, float factor){
	vector<TRootJet> output;
	for(unsigned int i=0;i<jets.size();i++){
		TRootJet jet(TLorentzVector(jets[i].Px()*factor,jets[i].Py()*factor,jets[i].Pz()*factor,jets[i].Energy()*factor));
		jet.setJetType(jets[i].jetType());
		jet.setPartonFlavour(jets[i].partonFlavour());
		jet.setBtag_trackCountingHighEffBJetTags(jets[i].btag_trackCountingHighEffBJetTags());
		jet.setBtag_trackCountingHighPurBJetTags(jets[i].btag_trackCountingHighPurBJetTags());
		jet.setBtag_jetProbabilityBJetTags(jets[i].btag_jetProbabilityBJetTags());
		jet.setBtag_jetBProbabilityBJetTags(jets[i].btag_jetBProbabilityBJetTags());
		jet.setBtag_simpleSecondaryVertexBJetTags(jets[i].btag_simpleSecondaryVertexBJetTags());
		output.push_back(jet);
	}
	return output;
}

vector<TRootPFJet> JESRescale(vector<TRootPFJet> &jets, float factor){
	vector<TRootPFJet> output;
	for(unsigned int i=0;i<jets.size();i++){
		TRootPFJet jet(TLorentzVector(jets[i].Px()*factor,jets[i].Py()*factor,jets[i].Pz()*factor,jets[i].Energy()*factor));
		jet.setJetType(jets[i].jetType());
		jet.setPartonFlavour(jets[i].partonFlavour());
		jet.setBtag_trackCountingHighEffBJetTags(jets[i].btag_trackCountingHighEffBJetTags());
		jet.setBtag_trackCountingHighPurBJetTags(jets[i].btag_trackCountingHighPurBJetTags());
		jet.setBtag_jetProbabilityBJetTags(jets[i].btag_jetProbabilityBJetTags());
		jet.setBtag_jetBProbabilityBJetTags(jets[i].btag_jetBProbabilityBJetTags());
		jet.setBtag_simpleSecondaryVertexBJetTags(jets[i].btag_simpleSecondaryVertexBJetTags());
		output.push_back(jet);
	}
	return output;
}

vector<TRootCaloJet> JESRescale(vector<TRootCaloJet> &jets, float factor){
	vector<TRootCaloJet> output;
	for(unsigned int i=0;i<jets.size();i++){
		TRootCaloJet jet(TLorentzVector(jets[i].Px()*factor,jets[i].Py()*factor,jets[i].Pz()*factor,jets[i].Energy()*factor));
		jet.setJetType(jets[i].jetType());
		jet.setPartonFlavour(jets[i].partonFlavour());
		jet.setBtag_trackCountingHighEffBJetTags(jets[i].btag_trackCountingHighEffBJetTags());
		jet.setBtag_trackCountingHighPurBJetTags(jets[i].btag_trackCountingHighPurBJetTags());
		jet.setBtag_jetProbabilityBJetTags(jets[i].btag_jetProbabilityBJetTags());
		jet.setBtag_jetBProbabilityBJetTags(jets[i].btag_jetBProbabilityBJetTags());
		jet.setBtag_simpleSecondaryVertexBJetTags(jets[i].btag_simpleSecondaryVertexBJetTags());
		jet.setEcalEnergyFraction(jets[i].ecalEnergyFraction());
		jet.setHcalEnergyFraction(jets[i].hcalEnergyFraction());
		jet.setn90Hits(jets[i].n90Hits());
		jet.setfHPD(jets[i].fHPD());
		output.push_back(jet);
	}
	return output;
}

/*
What is missing:
- perform pseudo ex if isMC for ABCD, QCD, SingleTop (Thierry & Greg)
- pseudo MC ... (Eric) 
- do SingleTop subtraction ->OK
- compute cross-section error ->ok  -> break down in efficiencies missing -> done not yet tested + diviser error in stat & syst (not yet correct)
- table with results (latex) -> ok -> few info missing -> ok
- selection table (latex) -> implemented, to be tested
- class a not yet implemented class with plots of variables used for selection etc ... (Thierry)
- do a MC/data comparison using a not yet implement class which store numbers if isMC (loading or reading) -> OK but too be tested
	-> partialy done, not yet tested , btag to be added
*/


int main (int argc, char *argv[])
{
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
  string xmlFileName ="../config/myconfig.xml";
  if (argc >= 2)
    xmlFileName=string(argv[1]);
  const char *xmlfile = xmlFileName.c_str();
  //Output files
  string rootFileName ("CrossSectionMeasurement.root");

  //Configuration output format
  TTree* configTree = new TTree("configTree","configuration Tree");
  TClonesArray* tcdatasets = new TClonesArray("Dataset",1000);
  configTree->Branch("Datasets","TClonesArray",&tcdatasets);
  TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000);
  configTree->Branch("AnaEnv","TClonesArray",&tcAnaEnv);
  TClonesArray* tcMCExp = new TClonesArray("MCExpectation",1000);
  configTree->Branch("MCExp","TClonesArray",&tcMCExp);
  TClonesArray* tcMCExpBjetMult = new TClonesArray("MCObsExpectation",1000);
  configTree->Branch("BjetMultiplicity","TClonesArray",&tcMCExpBjetMult);

  ////////////////////////////////////
  /// AnalysisEnvironment  
  ////////////////////////////////////

  AnalysisEnvironment anaEnv;
  cout<<"Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);
  int verbose = anaEnv.Verbose;
  float Luminosity = anaEnv.Luminosity;	// in 1/pb
  cout << "The results will be obtained for a luminosity of " << Luminosity << " x 1/pb" << endl;
  cout << "Nb of jet bins considered : " << anaEnv.NofJetBins << " , starting with events with at least " << anaEnv.NofJets << " jets." <<endl;

  /////////////////////
  // Load Datasets
  /////////////////////

  TTreeLoader treeLoader;
  cout << " - Load datasets ..." << endl;
  vector < Dataset > datasets;
  treeLoader.LoadDatasets (datasets, xmlfile);
  for(unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(datasets[i]);
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
  //after JES applied
  vector < TRootJet > jets;
  vector < TRootCaloJet > Calojets;
  vector < TRootPFJet > PFjets;
  ////////////////////////////////////////
  
  TFile *fout = new TFile (rootFileName.c_str (), "RECREATE");
  //Global variable
  TRootEvent* event = 0;
  //Pseudo-Exp
  CrossSectionPseudoExp* XSMCPseudoExp = new CrossSectionPseudoExp(10,2*anaEnv.Luminosity,80*anaEnv.Luminosity,25*anaEnv.Luminosity);
  int* ievtmin = new int[datasets.size()];
  for(unsigned int i=0;i<datasets.size();i++) ievtmin[i] = 0;
  
  char name[100];
  ////////////////////////////////////
  //  MCExpectation (if MC)
  ////////////////////////////////////
  MCExpectation* MCExp = NULL;
  vector<MCObsExpectation*> MCObsExpBjetMult;
  //size must be [anaEnv.NofJetBins*anaEnv.NofBtagWorkingPoint_vjEst];
  //store in a array b-jet multiplicity for #jets per working point (ex: 3jWP0,4jWP0,3jWP1,4jWP1)
  if(anaEnv.isMC){
  	//create one
	MCExp = new MCExpectation();
  	MCExp->SetLuminosity(Luminosity);
	MCExp->SetLabel(string("CrossSection"));
  	for(int i=0;i<anaEnv.NofBtagWorkingPoint_vjEst;i++){
  		for(int j=0;j<anaEnv.NofJetBins;j++){
			sprintf(name,"bjet_multiplicity_%d_wp_%d_jets",i,j);
			MCObsExpBjetMult.push_back(new MCObsExpectation(4,0,4,name,1));
		}
  	}
  }
  else{
  	//read one
  	TFile fin(anaEnv.MCExpFilename.c_str(),"READ");
  	TClonesArray* tcMCExpTemp = new TClonesArray("MCExpectation",0);
  	TClonesArray* tcMCExpTBjetsTemp = new TClonesArray("MCObsExpectation",0);
	TTree* tree = (TTree*) fin.Get("configTree");
	TBranch* b1 = 0; 
	TBranch* b2 = 0; 
	if(tree){ 
		b1 = tree->GetBranch("MCExp");
		if(b1) b1->SetAddress(&tcMCExpTemp);
		b2 = tree->GetBranch("BjetMultiplicity");
		if(b2) b2->SetAddress(&tcMCExpTBjetsTemp);
		tree->GetEntry(0);
	}
	MCExpectation* temp = 0;
	if(tcMCExpTemp->GetEntries()){
		temp = (MCExpectation*) tcMCExpTemp->At(0);
		//MCExp = new MCExpectation((*temp)); 
		MCExp = (MCExpectation*) temp->Clone(); 
	}
	MCObsExpectation* temp2 = 0;
	for(unsigned x=0;x<tcMCExpTBjetsTemp->GetEntries();x++){
			temp2 = (MCObsExpectation*) tcMCExpTBjetsTemp->At(x);
			MCObsExpBjetMult.push_back((MCObsExpectation*) temp2->Clone());
	}
  	delete temp;
	delete temp2;
	tcMCExpTemp->Delete();
	tcMCExpTBjetsTemp->Delete();
	tree->Delete();
	//delete b1;//->Delete();
	//delete b2;//->Delete();
  }
  
  ////////////////////////////////////
  //	Loop on  Pseudo-experiments
  ////////////////////////////////////
  int nRuns = 1;
  if (anaEnv.nPseudoExp!=0)
    nRuns=anaEnv.nPseudoExp;
  for(int iPseudoExp=0; iPseudoExp<nRuns ; iPseudoExp++)
 {

  //Object which store info per pseudo-exp
  MCExpectation MCExp_PS;
  MCExp_PS.SetLuminosity (anaEnv.Luminosity);
  MCExp_PS.SetLabel (string ("CrossSection"));


  //nof selected events
  float NEvtsData = 0;
  float NEvtsTrigged = 0;
  Double_t *nEvents = new Double_t[datasets.size()];
  pair<float, float> MuonSelEff = pair<float,float>(0,0);
  pair<float, float> SecondLeptonVetoEff = pair<float,float>(0,0);
  pair<float, float> JetSelEff = pair<float,float>(0,0);

  ////////////////////////////////////
  /// Conditions: b-tagging 
  ////////////////////////////////////
  // for V+jets estimation
  int BtagAlgo_vjEst            = anaEnv.BtagAlgo_vjEst;
  int NofBtagWorkingPoint_vjEst = anaEnv.NofBtagWorkingPoint_vjEst;
  float *BtagWorkingPoint_vjEst= anaEnv.BtagWorkingPoint_vjEst;  
  //for(int i=0;i<NofBtagWorkingPoint_vjEst;i++) cout<<"Btagging wp "<<i<<" : "<<BtagWorkingPoint_vjEst[i]<<endl;

  // for tt+jets estimation
  float BtagDiscriCut_ttjEst  = anaEnv.BtagDiscriCut_ttjEst;
  int BtagAlgo_ttjEst         = anaEnv.BtagAlgo_ttjEst;


  ////////////////////////////////////
  /// V+jets Estimation
  ////////////////////////////////////
  if (verbose > 0)
    cout << " - Configuration of VJetEstimation ..." << endl;
  //Declare VJetEstimation
  //
  vector < int >iDTTLike;
  vector < int >iDWLike;
  vector < int >iDWbLike;
  for (unsigned int d = 0; d < datasets.size (); d++) {
    if (datasets[d].Name () == "TTJets")   iDTTLike.push_back(d);
    if (datasets[d].Name () == "WJets")    iDWLike.push_back(d);
    if (datasets[d].Name () == "ZJets")    iDWLike.push_back(d);
    if (datasets[d].Name () == "QCD")      iDWLike.push_back(d);
    if (datasets[d].Name () == "Wcc")      iDWLike.push_back(d);
    if (datasets[d].Name () == "Zcc")      iDWLike.push_back(d);
    if (datasets[d].Name () == "Vcc")      iDWLike.push_back(d);
  }
  vector < string > BckgdNames;
  BckgdNames.push_back ("SingleTop_s_ch");
  BckgdNames.push_back ("SingleTop_t_ch");
  BckgdNames.push_back ("SingleTop_tW_ch");
  
  VJetEstimation   *vjEstimation   = 0;
  VJetEstPseudoExp *vjEstPseudoExp = 0;

  if (anaEnv.doVJEstim){
	//n = 3 jets : e0bq/e1bq/e2bq = 0.0466/0.4450/0.5063/
	//n > 4 jets : e0bq/e1bq/e2bq = 0.0201/0.2687/0.6943/
	double** EffEbq = new double*[anaEnv.NofJetBins];
	for(int i=0; i<anaEnv.NofJetBins; i++) EffEbq[i] = new double[3];
	if(anaEnv.NofJetBins > anaEnv.EffEbsel.size()) {cout<<"Mismatch between NofJetBins ("<<anaEnv.NofJetBins<<") and EffEbsel.size("<<anaEnv.EffEbsel.size()<<")... Abort!"<<endl;return 0;}
	for(unsigned int i=0;i<anaEnv.NofJetBins;i++){
  		if(anaEnv.EffEbsel[i].size() != 3 ) {cout<<"EffEbsel["<<i<<"].size() != 3... Abort!"<<endl;return 0;}
  		for(unsigned int j=0;j<3;j++) EffEbq[i][j] = anaEnv.EffEbsel[i][j];
	}

  	//vjEstimation = new VJetEstimation (NofBtagWorkingPoint_vjEst,BtagWorkingPoint_vjEst,anaEnv.NofJets,anaEnv.NofJetBins,(int) datasets.size (),iDTTLike,iDWLike,iDWbLike);
	vjEstimation   = new VJetEstimation (NofBtagWorkingPoint_vjEst,BtagWorkingPoint_vjEst,anaEnv.NofJets,anaEnv.NofJetBins,EffEbq,(int) datasets.size (),iDTTLike,iDWLike,iDWbLike);
	if (anaEnv.doVJEstPE) vjEstPseudoExp = new VJetEstPseudoExp(anaEnv.NVJetPE,vjEstimation);
  }

  ////////////////////////////////////
  /// ABCD Estimation for QCD
  ////////////////////////////////////
  if (verbose > 0)
    cout << " - Configuration of ABCDEstimation ..." << endl;
  float empiettementX = 0.;// empiettement sur la region voisine en X, en % de la partie voisine (la bande supprimée comptant comme empiettement, les 100% échus avant la borne extrême atteinte)
  float empiettementY = 0.;// empiettement sur la region voisine en X, en % de la partie voisine (la bande supprimée comptant comme empiettement, les 100% échus avant la borne extrême atteinte)
  ABCDEstimation abcdEstim (TString ("ABCD"), TString ("ABCD"), TString ("RelIso"), TString ("IP Significance"), anaEnv.NXbinsABCD, anaEnv.XbinMinABCD, anaEnv.XbinMaxABCD, anaEnv.NYbinsABCD, anaEnv.YbinMinABCD, anaEnv.YbinMaxABCD);
  
  ABCDEstimation abcdEstimA (TString ("ABCD_regionA"), TString ("ABCD_regionA"), TString ("RelIso"), TString ("IP Significance"), anaEnv.NXbinsABCD, anaEnv.XbinMinABCD, anaEnv.cutX0+(anaEnv.XbinMaxABCD-anaEnv.cutX1)*empiettementX/100., anaEnv.NYbinsABCD, anaEnv.YbinMinABCD, anaEnv.cutY0+(anaEnv.YbinMaxABCD-anaEnv.cutY1)*empiettementY/100.);
  ABCDEstimation abcdEstimB (TString ("ABCD_regionB"), TString ("ABCD_regionB"), TString ("RelIso"), TString ("IP Significance"), anaEnv.NXbinsABCD, anaEnv.cutX1-(anaEnv.cutX0-anaEnv.XbinMinABCD)*empiettementX/100., anaEnv.XbinMaxABCD, anaEnv.NYbinsABCD, anaEnv.YbinMinABCD, anaEnv.cutY0+(anaEnv.YbinMaxABCD-anaEnv.cutY1)*empiettementY/100.);
  ABCDEstimation abcdEstimC (TString ("ABCD_regionC"), TString ("ABCD_regionC"), TString ("RelIso"), TString ("IP Significance"), anaEnv.NXbinsABCD, anaEnv.XbinMinABCD, anaEnv.cutX0+(anaEnv.XbinMaxABCD-anaEnv.cutX1)*empiettementX/100., anaEnv.NYbinsABCD, anaEnv.cutY1-(anaEnv.cutY0-anaEnv.YbinMinABCD)*empiettementY/100., anaEnv.YbinMaxABCD);
  ABCDEstimation abcdEstimD (TString ("ABCD_regionD"), TString ("ABCD_regionD"), TString ("RelIso"), TString ("IP Significance"), anaEnv.NXbinsABCD, anaEnv.cutX1-(anaEnv.cutX0-anaEnv.XbinMinABCD)*empiettementX/100., anaEnv.XbinMaxABCD, anaEnv.NYbinsABCD, anaEnv.cutY1-(anaEnv.cutY0-anaEnv.YbinMinABCD)*empiettementY/100., anaEnv.YbinMaxABCD);
  if (verbose > 0)
    cout << " - ABCDEstimation configured ..." << endl;
 
  ////////////////////////////////////
  /// Selection plots
  ////////////////////////////////////
  MultiSamplePlot ***multiplot = new MultiSamplePlot**[NofBtagWorkingPoint_vjEst];
  for(int i=0;i<NofBtagWorkingPoint_vjEst;i++){
  	multiplot[i] = new MultiSamplePlot*[anaEnv.NofJetBins];
  	for(int j=0;j<anaEnv.NofJetBins;j++){
  		sprintf(name,"bjet_multiplicity_%d_wp_%d_jets",i,j+anaEnv.NofJets);
  		multiplot[i][j] = new MultiSamplePlot(datasets,name,4,0,4,"Number of b-tagged jets","");
  	}
  }
  if (verbose > 0)
    cout << " - MultiSamplePlot instantiated ..." << endl;
   
  ////////////////////////////////////
  /// Selection table
  ////////////////////////////////////
  vector<string> CutsSelecTable;
  CutsSelecTable.push_back(string("all"));
  CutsSelecTable.push_back(string("trigged"));
  CutsSelecTable.push_back(string("Good PV"));
  CutsSelecTable.push_back(string("$\\geq$ 1 selected muon"));
  CutsSelecTable.push_back(string("Veto 2nd muon"));
  CutsSelecTable.push_back(string("Veto electron"));
  char LabelNJets[100];
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-2);
  CutsSelecTable.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-1);
  CutsSelecTable.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets);
  CutsSelecTable.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets+1);
  CutsSelecTable.push_back(string(LabelNJets));
  if (verbose > 0)
    cout << " - CutsSelectionTable instantiated ..." << endl;
  SelectionTable selecTable(CutsSelecTable, datasets);
  selecTable.SetLuminosity(Luminosity);
  if (verbose > 0)
    cout << " - SelectionTable instantiated ..." << endl;

  ////////////////////////////////////
  //	Loop on datasets
  ////////////////////////////////////
  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
  for (unsigned int d = 0; d < datasets.size (); d++) {
    if (verbose > 1)
      cout << "   Dataset " << d << ": " << datasets[d].Name () << "/ title : " << datasets[d].Title () << endl;

    //open files and load
    cout<<"LoadEvent"<<endl;
    treeLoader.LoadDataset (datasets[d], anaEnv);
    cout<<"LoadEvent"<<endl;
    ////////////////////////////////////////
    // Pseudo-exp
    MCPseudoExp pseudoExp;
    pseudoExp.CalculatePseudoExp(datasets[d], Luminosity, ievtmin[d]);
    vector<int> indices;
    if (anaEnv.nPseudoExp==0)
    {
     indices.clear();
     for (int ii=0; ii<datasets[d].NofEvtsToRunOver(); ii++)
       indices.push_back(ii);
    }
    else
      indices = pseudoExp.GetRandomNumbers();

    ////////////////////////////////////
    //	Loop on events
    ////////////////////////////////////
    nEvents[d] = 0;
    if (verbose > 1)
      cout << "	Loop over events " << endl;
    for (vector<int>::iterator iter=indices.begin(); iter!=indices.end(); iter++) {
      int ievt = *iter;
      //    for (unsigned int ievt = 0; ievt < datasets[d].NofEvtsToRunOver(); ievt++) { 
      nEvents[d]++;
      if(ievt%5000 == 0) std::cout<<"Processing the "<<ievt<<"th event"<<flush<<"\r";

      TLorentzVector JESbalance (0, 0, 0, 0);
     
      //load event
      if(anaEnv.JetType == 0 ) event = treeLoader.LoadEvent (ievt, datasets[d], vertex, init_muons, init_electrons, init_jets, mets);
      if(anaEnv.JetType == 1 ) event = treeLoader.LoadEvent (ievt, datasets[d], vertex, init_muons, init_electrons, init_Calojets, mets);
      if(anaEnv.JetType == 2 ) event = treeLoader.LoadEvent (ievt, datasets[d], vertex, init_muons, init_electrons, init_PFjets, mets);
      ////////////////////////////////////////
    
      int itrigger = treeLoader.iTrigger (datasets[d], string ("HLT_Mu9"), event->runId());
       
      //apply JES
      jets = JESRescale(init_jets, anaEnv.JES);
      Calojets = JESRescale(init_Calojets, anaEnv.JES);
      PFjets = JESRescale(init_PFjets, anaEnv.JES);

      /////////////////////////////
      //   Selection
      /////////////////////////////
      //Declare selection instance    
      Selection selection(anaEnv.JetType, jets, Calojets, PFjets, init_muons, init_electrons, mets);
      selection.SetConfiguration(anaEnv);
      bool trigged = treeLoader.EventTrigged (itrigger);
		
      /////////////////////////////
      //   Selection table
      /////////////////////////////
      
      bool isInSR = false;
      vector<TRootJet> selectedJets = selection.GetSelectedJets(anaEnv.JetsPtCutSR,anaEnv.JetsEtaCutSR, anaEnv.applyJetID);
      vector<TRootMuon> selectedMuons = selection.GetSelectedMuons(anaEnv.MuonPtCutSR,anaEnv.MuonEtaCutSR,anaEnv.MuonRelIsoCutSR, selectedJets);
      vector<TRootMuon> vetoMuons = selection.GetSelectedLooseMuons(anaEnv.MuonPtCutVetoSR,anaEnv.MuonEtaCutVetoSR,anaEnv.MuonRelIsoCutVetoSR);
      vector<TRootElectron> vetoElectrons = selection.GetSelectedLooseElectrons(anaEnv.ElectronPtCut,anaEnv.ElectronEtaCut,anaEnv.ElectronRelIsoCut);
      bool isGoodPV = selection.isPVSelected(vertex, anaEnv.PVertexNdofCut, anaEnv.PVertexZCut, anaEnv.PVertexRhoCut);

      selecTable.Fill(d,0,1.);
      if(trigged){
        selecTable.Fill(d,1,1.);
        if(isGoodPV){
          selecTable.Fill(d,2,1.);
          if(selectedMuons.size()==1){
		      selecTable.Fill(d,3,1.);
      		if(vetoMuons.size()==1){
              selecTable.Fill(d,4,1.);
              if(vetoElectrons.size()==0){
                selecTable.Fill(d,5,1.);
                if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-2){
                  selecTable.Fill(d,6,1.);
                  if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-1){
                    selecTable.Fill(d,7,1.);
                    if(selectedJets.size()>=(unsigned int)anaEnv.NofJets){
                      selecTable.Fill(d,8,1.);
                      isInSR = true;
                      if(selectedJets.size()>=(unsigned int)anaEnv.NofJets+1){
                        selecTable.Fill(d,9,1.);
                        if(verbose>2) cout << event->runId() << ":" << event->eventId() << ":" << event->lumiBlockId()  << ":" << setprecision(8) << selectedMuons[0].Pt() << endl;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      //////////////////////////////
      
      if (!trigged) continue;

      if(datasets[d].Name()==string("TTJets")) NEvtsTrigged++;

      //////////////////////////
      //the different regions
      //////////////////////////
      bool isQCD = false;
      if(datasets[d].Name()==string("QCD")) isQCD=true;
      bool isInABCD = false;
      float RelIso = 100.;
      //////////////////////////

      //////////////////////////
      //btag jets collection
      //////////////////////////
      vector < TRootJet > btagJets[NofBtagWorkingPoint_vjEst][anaEnv.NofJetBins];
      // = selection.GetSelectedBJets(selectedJets,BtagAlgo_vjEst,BtagWorkingPoint_vjEst[0]);
      if(isInSR){
      	for(int i=0;i<NofBtagWorkingPoint_vjEst;i++){
      	      for(int j=0;j<anaEnv.NofJetBins-1;j++){
	      	if(selectedJets.size()==(unsigned int)anaEnv.NofJets+i){
			btagJets[i][j] = selection.GetSelectedBJets(selectedJets,BtagAlgo_vjEst,BtagWorkingPoint_vjEst[i]);
			multiplot[i][j]->Fill(btagJets[i][j].size(),datasets[d],true,Luminosity);
		}
	      }
	      if(selectedJets.size()>=(unsigned int)(anaEnv.NofJets+anaEnv.NofJetBins-1)){
		btagJets[i][anaEnv.NofJetBins-1] = selection.GetSelectedBJets(selectedJets,BtagAlgo_vjEst,BtagWorkingPoint_vjEst[i]);
	      	multiplot[i][anaEnv.NofJetBins-1]->Fill(btagJets[i][anaEnv.NofJetBins-1].size(),datasets[d],true,Luminosity);
	      }
	}
	if(anaEnv.isMC){
	      	for(int i=0;i<NofBtagWorkingPoint_vjEst;i++){
	      	      for(int j=0;j<anaEnv.NofJetBins;j++) MCObsExpBjetMult[i*anaEnv.NofJetBins+j]->Fill(datasets[d],btagJets[i][j].size(),datasets[d].NormFactor());
		}
	}
      }
      ///////////////////////////////////////////
      //fill booleans for the different region
      ///////////////////////////////////////////

      if(anaEnv.isMC){
        bool isMuonSelected = false;
	bool isPassVeto = false;
	bool isJetSelected = false;
	if(selectedMuons.size()>0) isMuonSelected = true;
	if(vetoMuons.size() == 1  && vetoElectrons.size () == 0) isPassVeto = true;
	if(selectedJets.size()>= (unsigned int) anaEnv.NofJets) isJetSelected = true;
	MCExp->FillEfficiencies(datasets[d], isMuonSelected, isPassVeto, isJetSelected );
	MCExp_PS.FillEfficiencies(datasets[d], isMuonSelected, isPassVeto, isJetSelected );
	if(datasets[d].Name()==string("TTJets")){
      		if(isMuonSelected) MuonSelEff.first+=datasets[d].NormFactor () * Luminosity; 
		if(isPassVeto) SecondLeptonVetoEff.first+=datasets[d].NormFactor () * Luminosity;
		if(isJetSelected) JetSelEff.first+=datasets[d].NormFactor () * Luminosity;
	}
      }

      if (isInSR){
	RelIso = selectedMuons[0].relativeIso03 ();
	if(anaEnv.isMC){
		NEvtsData = NEvtsData + datasets[d].NormFactor () * Luminosity;
        	MCExp->Fill(datasets[d], datasets[d].NormFactor () * Luminosity);
        	MCExp_PS.Fill(datasets[d], datasets[d].NormFactor () * Luminosity);
	}
	else NEvtsData++;
      }
      if (trigged && isGoodPV && selection.GetSelectedMuonsQCDEstim (anaEnv.MuonPtCutSR, anaEnv.MuonEtaCutSR, selectedJets).size () >= 1 && selectedMuons.size()<=1 && vetoMuons.size()<=1 && vetoElectrons.size()==0 && selectedJets.size () >= (unsigned int) anaEnv.NofJets) {
	isInABCD = true;
      }

      ////////////////////////////////////
      //  VJetEstimation: norm 
      ///////////////////////////////////

      //Signal Region
      if (isInSR && anaEnv.doVJEstim) {
      	vjEstimation->FillInputs(selectedJets,d,BtagAlgo_vjEst);
	if(!anaEnv.isMC) vjEstimation->BckgdSubstraction(MCObsExpBjetMult, BckgdNames, Luminosity);
      }
      ///////////////////////
      //  ABCDEstimation
      ///////////////////////
      if (isInABCD) {
	float d0      = selection.GetSelectedMuonsQCDEstim (anaEnv.MuonPtCutSR, anaEnv.MuonEtaCutSR, selectedJets)[0].d0();
	float d0error = selection.GetSelectedMuonsQCDEstim (anaEnv.MuonPtCutSR, anaEnv.MuonEtaCutSR, selectedJets)[0].d0error();
	float d0Signvalue = 100.;
	if (d0error > 0)
	  d0Signvalue = d0 / d0error;
	RelIso = selection.GetSelectedMuonsQCDEstim (anaEnv.MuonPtCutSR, anaEnv.MuonEtaCutSR, selectedJets)[0].relativeIso03 ();
	abcdEstim.Fill  (RelIso, d0Signvalue, datasets[d].NormFactor () * Luminosity, isQCD);
	abcdEstimA.Fill (RelIso, d0Signvalue, datasets[d].NormFactor () * Luminosity, isQCD);
	abcdEstimB.Fill (RelIso, d0Signvalue, datasets[d].NormFactor () * Luminosity, isQCD);
	abcdEstimC.Fill (RelIso, d0Signvalue, datasets[d].NormFactor () * Luminosity, isQCD);
	abcdEstimD.Fill (RelIso, d0Signvalue, datasets[d].NormFactor () * Luminosity, isQCD);
      }
      ///////////////////////

      //delete selection;
    }				//loop on events
    cout<<endl;
    if(anaEnv.isMC && anaEnv.doVJEstim && anaEnv.nPseudoExp==0) vjEstimation->ReScaleInputs(d,nEvents[d]/datasets[d].PreSelEfficiency(),(datasets[d].Xsection()*datasets[d].PreSelEfficiency()*Luminosity)/nEvents[d]);
    cout<<"XS/PreSelEff/nEvents : "<<datasets[d].Xsection()<<"/"<<datasets[d].PreSelEfficiency()<<"/"<<nEvents[d]<<endl;
    
    //important: free memory
    treeLoader.UnLoadDataset();
  }				//loop on datasets

  //Once everything is filled ...
  if (verbose > 0)
    cout << " We ran over all the data ;-)" << endl;

  ofstream ofile("CrossSectionTable.tex");
  //report on table ...
  ofile<<"\\documentclass[8pt]{article}"<<endl;
  ofile<<"\\begin{document}"<<endl;
  
  ////////////////////////////////////
  /// ABCDEstimation
  ////////////////////////////////////
  if (verbose > 0)
    cout << " - Perform ABCD method ..." << endl;

  //Cuts in the middle of the regions
  float cutX0mid = (anaEnv.cutXmin+anaEnv.cutX0)/2;
  float cutX1mid = (anaEnv.cutX1+anaEnv.cutXmax)/2;
  float cutY0mid = (anaEnv.cutYmin+anaEnv.cutY0)/2;
  float cutY1mid = (anaEnv.cutY1+anaEnv.cutYmax)/2;
  abcdEstim.ComputeEstimate (anaEnv.cutXmin, anaEnv.cutX0, anaEnv.cutX1, anaEnv.cutXmax, anaEnv.cutYmin, anaEnv.cutY0, anaEnv.cutY1, anaEnv.cutYmax, anaEnv.region);
  abcdEstimA.ComputeEstimate (anaEnv.cutXmin, cutX0mid, cutX0mid, anaEnv.cutX0, anaEnv.cutYmin, cutY0mid, cutY0mid, anaEnv.cutY0, anaEnv.region);
  abcdEstimB.ComputeEstimate (anaEnv.cutX1, cutX1mid, cutX1mid, anaEnv.cutXmax, anaEnv.cutYmin, cutY0mid, cutY0mid, anaEnv.cutY0, anaEnv.region);
  abcdEstimC.ComputeEstimate (anaEnv.cutXmin, cutX0mid, cutX0mid, anaEnv.cutX0, anaEnv.cutY1, cutY1mid, cutY1mid, anaEnv.cutYmax, anaEnv.region);
  abcdEstimD.ComputeEstimate (anaEnv.cutX1, cutX1mid, cutX1mid, anaEnv.cutXmax, anaEnv.cutY1, cutY1mid, cutY1mid, anaEnv.cutYmax, anaEnv.region);
  
   //    pseudoExpABCD.Fill(abcdEstim.GetData(1));
   //    pseudoExpABCD.Fill(abcdEstim.GetDataError(1));
   //    pseudoExpABCD.Fill(abcdEstim.GetQCDEstimate(1));
   //    pseudoExpABCD.Fill(abcdEstim.GetQCDEstimateError(1));
  
  if (verbose > 1)
    abcdEstim.Print2DEstimate ();
/*  
  printf("{{1.,2.,3.,4.}, {%16.15lE, %16.15lE, %16.15lE, %16.15lE}, {0., 0., 0., 0.}, {%16.15lE, %16.15lE, %16.15lE, %16.15lE}};\n",
         abcdEstim.GetData(1), abcdEstim.GetData(2), abcdEstim.GetData(3), abcdEstim.GetData(4),
         abcdEstim.GetDataError(1), abcdEstim.GetDataError(2), abcdEstim.GetDataError(3), abcdEstim.GetDataError(4));
  printf("[[1.,2.,3.,4.], [%16.15lE, %16.15lE, %16.15lE, %16.15lE], [0., 0., 0., 0.], [%16.15lE, %16.15lE, %16.15lE, %16.15lE]];\n",
         abcdEstim.GetEstimate(1), abcdEstim.GetEstimate(2), abcdEstim.GetEstimate(3), abcdEstim.GetEstimate(4),
         abcdEstim.GetEstimateError(1), abcdEstim.GetEstimateError(2), abcdEstim.GetEstimateError(3), abcdEstim.GetEstimateError(4));
*/
  ////////////////////////////////////
  /// VJets Estimation = NofEvents
  ////////////////////////////////////
if (anaEnv.doVJEstim){
  if (verbose > 0)
    cout << " - Perform VJets estimation ..." << endl;

  vjEstimation->SumOverAllInputs();
  // Fill internal arrays of numbers ; needed for the PrintInput/PrintResults method when running on MC
  if(anaEnv.isMC){
	// Compute from MC the b-tagging/mis-tagging efficiencies. Needed if comparison between predicted and estimated eff is required.
  	vjEstimation->ComputeEffFromMC();
  	std::string algo;
	std::stringstream salgo;
	salgo << BtagAlgo_vjEst;
	algo = salgo.str();
	for(int iwp=0;iwp<NofBtagWorkingPoint_vjEst;iwp++){
		MCExp->AddBtagConditions(algo,BtagWorkingPoint_vjEst[iwp],pair<float,float>(vjEstimation->GetPredEb(iwp),0.),pair<float,float>(vjEstimation->GetPredEudsc(iwp),0.),pair<float,float>(vjEstimation->GetPredEuds(iwp),0.));
		MCExp_PS.AddBtagConditions(algo,BtagWorkingPoint_vjEst[iwp],pair<float,float>(vjEstimation->GetPredEb(iwp),0.),pair<float,float>(vjEstimation->GetPredEudsc(iwp),0.),pair<float,float>(vjEstimation->GetPredEuds(iwp),0.));
  	}
  }
  if(anaEnv.isMC) vjEstimation->ComputeCondProb();

  // Compute from MC the efficiency to select 0/1/2 b-quarks for tt+jets events
  vjEstimation->ComputeEffbqFromMC(iDTTLike[0]);

  if (verbose > 1) for(int i=0;i<vjEstimation->GetNbOfJetsBins();i++) vjEstimation->PrintInputs(i);

  double** VJEstimParamInit = new double*[anaEnv.TagEffInit.size()];
  for(unsigned int ts=0;ts<anaEnv.TagEffInit.size();ts++){
  	VJEstimParamInit[ts] = new double[2+anaEnv.TagEffInit[ts].size()*anaEnv.TagEffInit[ts][0].size()];
	int it = 0;
	VJEstimParamInit[ts][it] = anaEnv.NTTlikeInit[ts]*Luminosity; it++;
	VJEstimParamInit[ts][it] = anaEnv.NVlikeInit[ts]*Luminosity; it++;
	for(unsigned int ts2=0;ts2<anaEnv.TagEffInit[ts].size();ts2++){
		for(unsigned int ts3=0;ts3<anaEnv.TagEffInit[ts][ts2].size();ts3++){
			VJEstimParamInit[ts][it] = (double) anaEnv.TagEffInit[ts][ts2][ts3]; 
			it++;
		}
	}
  }
  vjEstimation->SetInitialValues(VJEstimParamInit);

  bool UseUnBinMLE = anaEnv.useUnBinMLE;
  bool UseMJLE = anaEnv.useMJLE;
  bool doMinos = false;

  bool VJetVerbose = false;
  if(verbose>0) VJetVerbose = true;
  if(!UseMJLE){
  	if(!UseUnBinMLE) vjEstimation->  BinnedMaximumLikelihoodEst(anaEnv.MinMethod, anaEnv.MinOption, anaEnv.VJEstFixParam, true, VJetVerbose);
	else             vjEstimation->UnBinnedMaximumLikelihoodEst(anaEnv.VJEstFixParam, doMinos, false, VJetVerbose);
  }
  else{
	if(!UseUnBinMLE) vjEstimation->  BinnedMaximumJointWPLikelihoodEst(anaEnv.MinMethod, anaEnv.MinOption, anaEnv.VJEstFixParam, true, VJetVerbose); 
	else		 vjEstimation->UnBinnedMaximumJointWPLikelihoodEst(anaEnv.VJEstFixParam, doMinos, VJetVerbose);
	//else		 vjEstimation->UnBinnedMaximumNjetsLikelihoodEst(0,anaEnv.VJEstFixParam, VJetVerbose);
  }

  if (verbose > 1){
  	vjEstimation->PrintResults();
  	vjEstimation->PrintResults_LatexFormat(ofile);
  }
  if(anaEnv.doVJEstPE){
  //compute the errors for VJetEstimation with Pseudo-exp
    vjEstPseudoExp->SetInitNbjets(vjEstimation->GetN());
    vjEstPseudoExp->RollTheDice(anaEnv.MinMethod, anaEnv.MinOption, anaEnv.VJEstFixParam, doMinos, UseUnBinMLE, UseMJLE, false);
    vjEstPseudoExp->PrintResults();
    vjEstPseudoExp->ResetInitialValues();
  }
  vjEstimation->FillSummaryHistos();
}
  //////////////////////////////
  //  Cross-section measurement
  //////////////////////////////
  float NEvtsQCD = 0;
  float NEvtsQCDError = 0;
  float NEvtsVJets = 0;
  float NEvtsVJetsError = 0;
  float NEvtsTTJets = 0;
  float NEvtsTTJetsError = 0;
  float NEvtsSingleTop = 0;
  float NEvtsSingleTopError = -0;
  pair < float, float >QCDEstim;
  pair < float, float >VJetEstim;
  pair < float, float >TTJetEstim;
  pair < float, float >XsectionEstim;
 
  //NEvtsQCD        = abcdEstim.GetEstimate(anaEnv.region);
  //NEvtsQCDError   = abcdEstim.GetEstimateError(anaEnv.region);
  //temporary solution: protect QCDEstim to have non-nan results
  NEvtsQCD = 0.01;
  NEvtsQCDError = 0.1;
  //"0" should not be hard-coded
  NEvtsVJets       = (vjEstimation != 0 ? vjEstimation->GetEstNv(0) : 0);
  NEvtsVJetsError  = (vjEstimation != 0 ? vjEstimation->GetEstNvErr(0) : 0);
  NEvtsTTJets      = (vjEstimation != 0 ? vjEstimation->GetEstNtt(0) : 0);
  NEvtsTTJetsError = (vjEstimation != 0 ? vjEstimation->GetEstNttErr(0) : 0);

  QCDEstim = pair < float, float >(NEvtsQCD, NEvtsQCDError);
  VJetEstim = pair < float, float >(NEvtsVJets, NEvtsVJetsError);
  TTJetEstim = pair < float, float >(NEvtsTTJets, NEvtsTTJetsError);

  float Efficiencies = -9999;
  float EfficienciesError = 0;
  //Efficiencies is the product of many term
  Efficiencies = anaEnv.TriggerEff*anaEnv.SkimEff*anaEnv.MuonSelEff*anaEnv.SecondLeptonVetoEff*anaEnv.JetSelEff;
  //pow is ~low but lines are more readable !!
  EfficienciesError+=pow((Efficiencies/anaEnv.TriggerEff)*anaEnv.TriggerEffError,2);
  EfficienciesError+=pow((Efficiencies/anaEnv.SkimEff)*anaEnv.SkimEffError,2);
  EfficienciesError+=pow((Efficiencies/anaEnv.MuonSelEff)*anaEnv.MuonSelEffError,2);
  EfficienciesError+=pow((Efficiencies/anaEnv.SecondLeptonVetoEff)*anaEnv.SecondLeptonVetoEffError,2);
  EfficienciesError+=pow((Efficiencies/anaEnv.JetSelEff)*anaEnv.JetSelEffError,2);
  EfficienciesError = sqrt(EfficienciesError);

  float CrossSection = -9999;
  //compute cross-section
  CrossSection = (NEvtsData-NEvtsQCD-NEvtsVJets-NEvtsSingleTop)/(Efficiencies*Luminosity);
 
  //Compute the errors
  float CrossSectionError = 0;
  float StatError = pow(sqrt(NEvtsData)/(Efficiencies*Luminosity),2);
  float NQCDEffect = pow(NEvtsQCDError/(Efficiencies*Luminosity),2);
  float NVJetEffect = pow(NEvtsVJetsError/(Efficiencies*Luminosity),2);
  float NSingleTopEffect = pow(NEvtsSingleTopError/(Efficiencies*Luminosity),2);
  //break down in efficiencies
  float TriggerEffEffect = pow(anaEnv.TriggerEffError/(pow(anaEnv.TriggerEff,2))*(CrossSection*anaEnv.TriggerEff),2);
  float SkimEffEffect = pow(anaEnv.SkimEffError/(pow(anaEnv.SkimEff,2))*(CrossSection*anaEnv.SkimEff),2);
  float MuonSelEffEffect = pow(anaEnv.MuonSelEffError/(pow(anaEnv.MuonSelEff,2))*(CrossSection*anaEnv.MuonSelEff),2);
  float SecondLeptonVetoEffEffect = pow(anaEnv.SecondLeptonVetoEffError/(pow(anaEnv.SecondLeptonVetoEff,2))*(CrossSection*anaEnv.SecondLeptonVetoEff),2);
  float JetSelEffEffect = pow(anaEnv.JetSelEffError/(pow(anaEnv.JetSelEff,2))*(CrossSection*anaEnv.JetSelEff),2);
  //product efficiencies
  float EfficienciesEffect = pow(EfficienciesError/(pow(Efficiencies,2))*((NEvtsData-NEvtsQCD-NEvtsVJets-NEvtsSingleTop)/Luminosity),2);
  //anaEnv.LuminosityError has to be relative [0-1], ex: 0.1 means 10% uncertainty on luminosity
  float LuminosityEffect = pow(anaEnv.LuminosityError*CrossSection,2);
  CrossSectionError =  sqrt(NQCDEffect + NVJetEffect + NSingleTopEffect + EfficienciesEffect + LuminosityEffect);
 
  XsectionEstim = pair < float, float >(CrossSection, CrossSectionError);
  
  //Selection plots
  for(int i=0;i<NofBtagWorkingPoint_vjEst;i++){
  	for(int j=0;j<anaEnv.NofJetBins;j++){
  		sprintf(name,"bjet_multiplicity_%d_wp_%d_jets",i,j+anaEnv.NofJets);
  		multiplot[i][j]->Draw(name);
  		multiplot[i][j]->Write(fout,name);
  	}
  }
  //Selection tables
  selecTable.TableCalculator();
  selecTable.Write(ofile);
  
  //Table of errors
  ofile<<"\\begin{table}"<<endl;
  ofile<<"\\centering"<<endl;
  ofile<<"\\begin{tabular}{|l|cc|}"<<endl;
  ofile<<"\\hline"<<endl;
  ofile<<"&  RelError (\\%) & cross-section error(\\%)"<<"\\\\"<<endl;
  ofile<<"\\hline"<<endl;
  ofile<<"$\\sqrt(N_{data})$"<<" & "<< (1/sqrt(NEvtsData))*100<<" \\% & "<< sqrt(StatError)/CrossSection*100 <<"\\%\\\\"<<endl;
  ofile<<"$\\sigma(N_{QCD})$"<<" & "<< (NEvtsQCDError/NEvtsQCD)*100<<" \\% & "<< sqrt(NQCDEffect)/CrossSection*100 <<"\\%\\\\"<<endl;
  ofile<<"$\\sigma(N_{VJets})$"<<" & "<< (NEvtsVJetsError/NEvtsVJets)*100<<"\\% & "<< sqrt(NVJetEffect)/CrossSection*100 <<"\\%\\\\"<<endl;
  if(NEvtsSingleTop>0)
  	ofile<<"$\\sigma(N_{SingleTop})$"<<" & "<< (NEvtsSingleTopError/NEvtsSingleTop)*100 <<"\\% & "<< sqrt(NSingleTopEffect)/CrossSection*100 <<"\\%\\\\"<<endl;
  ofile<<"$\\sigma(\\epsilon(trigger))$"<<" & "<< (anaEnv.TriggerEffError/anaEnv.TriggerEff)*100 <<"\\% & "<< sqrt(TriggerEffEffect)/CrossSection*100 <<"\\%\\\\"<<endl;
  ofile<<"$\\sigma(\\epsilon(skimming))$"<<" & "<< (anaEnv.SkimEffError/anaEnv.SkimEff)*100 <<"\\% & "<< sqrt(SkimEffEffect)/CrossSection*100 <<"\\%\\\\"<<endl;
  ofile<<"$\\sigma(\\epsilon(muon))$"<<" & "<< (anaEnv.MuonSelEffError/anaEnv.MuonSelEff)*100 <<"\\% & "<< sqrt(MuonSelEffEffect)/CrossSection*100 <<"\\%\\\\"<<endl;
  ofile<<"$\\sigma(\\epsilon(veto 2nd lepton)$"<<" & "<< (anaEnv.SecondLeptonVetoEffError/anaEnv.SecondLeptonVetoEff)*100 <<"\\% & "<< sqrt(SecondLeptonVetoEffEffect)/CrossSection*100 <<"\\%\\\\"<<endl;
  ofile<<"$\\sigma(\\epsilon(jets))$"<<" & "<< (anaEnv.JetSelEffError/anaEnv.JetSelEff)*100 <<"\\% & "<< sqrt(JetSelEffEffect)/CrossSection*100 <<"\\%\\\\"<<endl;
  ofile<<"Lumi"<<" & "<< anaEnv.LuminosityError*100<<"\\% & "<< sqrt(LuminosityEffect)/CrossSection*100 <<"\\%\\\\"<<endl;
  ofile<<"\\hline"<<endl;
  ofile<<"Total"<<" & "<<" & "<< (CrossSectionError/CrossSection)*100 <<"\\%\\\\"<<endl;
  ofile<<"\\hline"<<endl;
  ofile<<"\\end{tabular}"<<endl;
  ofile<<"\\caption{}"<<endl;
  ofile<<"\\label{tab:}"<<endl;
  ofile<<"\\end{table}"<<endl<<endl;
  ofile<<"Final result:"<<endl; 
  ofile<<"\\begin{equation}"<<endl;
  ofile<<"\\sigma(ttbar) = "<<CrossSection<<" \\pm "<< sqrt(StatError) <<"(stat) \\pm "<< CrossSectionError <<"(syst)"<<endl;
  ofile<<"\\end{equation}"<<endl;


  ////////////////////////////////////
  /// MCExpectation
  ////////////////////////////////////
  if(anaEnv.isMC){
  	MCExp->Compute();
  	MCExp_PS.Compute();
  	float a = 0;
	a = MuonSelEff.first; if(NEvtsTrigged>0) MuonSelEff.first = a/NEvtsTrigged; MuonSelEff.second = sqrt(MuonSelEff.first*(1-MuonSelEff.first)/NEvtsTrigged); 
	a = SecondLeptonVetoEff.first; if(NEvtsTrigged>0) SecondLeptonVetoEff.first = a/NEvtsTrigged; SecondLeptonVetoEff.second = sqrt(SecondLeptonVetoEff.first*(1-SecondLeptonVetoEff.first)/NEvtsTrigged); 
	a = JetSelEff.first; if(NEvtsTrigged>0) JetSelEff.first = a/NEvtsTrigged; JetSelEff.second = sqrt(JetSelEff.first*(1-JetSelEff.first)/NEvtsTrigged); 
	MCExp->SetEfficiencies(MuonSelEff, SecondLeptonVetoEff, JetSelEff);
  	//MCExp->AddBtagConditions(string algo, float bdisCut, pair<float,float> Eb, pair<float,float> Eusdc, pair<float,float> Eudsc_vj );//from Vj
  	new ((*tcMCExp)[0]) MCExpectation(*MCExp);
  	for(unsigned int i=0;i<MCObsExpBjetMult.size();i++ ) new ((*tcMCExpBjetMult)[i]) MCObsExpectation(*MCObsExpBjetMult[i]);
  }
  else{
  	//do MC/data comparison
  	//Eff
  	//
  	ofile<<"\\begin{table}"<<endl;
  	ofile<<"\\centering"<<endl;
  	ofile<<"\\begin{tabular}{|l|cc|}"<<endl;
  	ofile<<"\\hline"<<endl;
  	ofile<<"&  MC Expectation  & Data-driven estimation"<<"\\\\"<<endl;
	ofile<<"\\hline"<<endl;
	/*
	vector< pair<float,float> > Eb;
	vector< pair<float,float> > Eudsc;
	vector< pair<float,float> > Eudsc_vj;
	*/
	//add comparison for btag efficiencies ...
	ofile<<"\\hline"<<endl;
	ofile<<"QCD events &"<< MCExp->NQCD.first<<" $\\pm$ "<<MCExp->NQCD.second<<" & "<< NEvtsQCD<<" $\\pm$ "<<NEvtsQCDError <<"\\\\"<<endl;
	ofile<<"V+jets events &"<< MCExp->GetVJets().first<<" $\\pm$ "<<MCExp->GetVJets().second<<" & "<< NEvtsVJets<<" $\\pm$ "<<NEvtsVJetsError <<"\\\\"<<endl;
	ofile<<"tt+jets events &"<< MCExp->NTtJets.first<<" $\\pm$ "<<MCExp->NTtJets.second<<" & "<< (NEvtsData-NEvtsQCD-NEvtsVJets-NEvtsSingleTop)<<" $\\pm$ "<<0<<"\\\\"<<endl;
	ofile<<"\\hline"<<endl;
	ofile<<"$\\epsilon(trigger)$ &"<<MCExp->TriggerEff.first<<" $\\pm$ "<<MCExp->TriggerEff.second<<" & "<<anaEnv.TriggerEff<<" $\\pm$ "<<anaEnv.TriggerEffError<<"\\\\"<<endl;
	ofile<<"$\\epsilon(skim)$ &"<<MCExp->SkimEff.first<<" $\\pm$ "<<MCExp->SkimEff.second<<" & "<<anaEnv.SkimEff<<" $\\pm$ "<<anaEnv.SkimEffError<<"\\\\"<<endl;
	ofile<<"$\\epsilon(muon)$ &"<<MCExp->MuonSelEff.first<<" $\\pm$ "<<MCExp->MuonSelEff.second<<" & "<<anaEnv.MuonSelEff<<" $\\pm$ "<<anaEnv.MuonSelEffError<<"\\\\"<<endl;
	ofile<<"$\\epsilon(2nd lepton veto)$ &"<<MCExp->SecondLeptonVetoEff.first<<" $\\pm$ "<<MCExp->SecondLeptonVetoEff.second<<" & "<<anaEnv.SecondLeptonVetoEff<<" $\\pm$ "<<anaEnv.SecondLeptonVetoEffError<<"\\\\"<<endl;
	ofile<<"$\\epsilon(jets)$ &"<<MCExp->JetSelEff.first<<" $\\pm$ "<<MCExp->JetSelEff.second<<" & "<<anaEnv.JetSelEff<<" $\\pm$ "<<anaEnv.JetSelEffError<<"\\\\"<<endl;
	ofile<<"\\hline"<<endl;
  	ofile<<"\\end{tabular}"<<endl;
  	ofile<<"\\caption{}"<<endl;
  	ofile<<"\\label{tab:}"<<endl;
  	ofile<<"\\end{table}"<<endl<<endl;
  }
  ofile<<"\\end{document}"<<endl;
  //  system("pdflatex CrossSectionTable.txt");
 
  //Fill the class for Pseudo-Exp
  XSMCPseudoExp->Fill (MCExp_PS, QCDEstim, VJetEstim, TTJetEstim, XsectionEstim);

  ///////////////////
  // Writing
  //////////////////
  if (verbose > 1)
    cout << " - Start writing the module outputs in the ouput file ..." << endl;
  if (iPseudoExp==0 ) {
   char nom[200] = "";
   sprintf(nom, "_exp%d_", iPseudoExp);
   if (anaEnv.doVJEstim) vjEstimation  ->Write (fout, string(nom));
   if (anaEnv.doVJEstPE) vjEstPseudoExp->Write (fout, string(nom));
   abcdEstim.Write (fout, string(nom));
   ABCDPseudoExp ape(abcdEstim.GetData(1), abcdEstim.GetData(2), abcdEstim.GetData(3), abcdEstim.GetData(4), abcdEstim.GetDataError(1), abcdEstim.GetEstimateError(1));
   for (int i=0; i<5000; i++)
    ape.GenExp();
   ape.Write();
   sprintf(nom, "_region_A_exp%d_", iPseudoExp);
   abcdEstimA.Write (fout,string(nom));
   sprintf(nom, "_region_B_exp%d_", iPseudoExp);
   abcdEstimB.Write (fout,string(nom));
   sprintf(nom, "_region_C_exp%d_", iPseudoExp);
   abcdEstimC.Write (fout,string(nom));
   sprintf(nom, "_region_D_exp%d_", iPseudoExp);
   abcdEstimD.Write (fout,string(nom));
  }
  delete MCExp;
  MCExp=NULL;
  delete vjEstPseudoExp;
  delete vjEstimation;
  vjEstimation=NULL;
  vjEstPseudoExp=NULL;
 } // loop over pseudo-experiments
 
  ///////////////////
  // Writting
  //////////////////
  if (verbose > 1)
  	cout << " - Writting outputs on files ..." << endl;
  if (anaEnv.nPseudoExp!=0)
	XSMCPseudoExp->Write(fout);
  delete XSMCPseudoExp;

  //compute
  if(anaEnv.isMC) for(unsigned int x=0;x<MCObsExpBjetMult.size();x++) MCObsExpBjetMult[x]->Compute();
  //add configuration info
  fout->cd();
  configTree->Fill();
  configTree->Write();
  //
  if (verbose > 1)
    cout << " - Done with writing the module outputs in the ouput file ..." << endl;
    cout << " - Closing the output file now..." << endl;
  //fout->Write ();
  fout->Close();

  //delete
  delete fout;
  delete event;
  delete tcdatasets;
  delete tcAnaEnv;
  delete configTree;

  cout << "**********************************************************************" << endl;
  cout << "           End of the program !!" << endl;
  cout << " 		hasn't crashed yet ;-) " << endl;
  cout << "**********************************************************************" << endl;

  return 0;
}
