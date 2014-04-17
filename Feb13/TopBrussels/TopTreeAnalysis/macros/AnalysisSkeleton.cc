
////////////////////////////
///// TODO & COMMENTS /////
/////////////////////////// 

#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "TopTreeAnalysisBase/MCInformation/interface/MCWeighter.h"
#include "TopTreeAnalysisBase/Selection/interface/ElectronPlotter.h"
#include "TopTreeAnalysisBase/Selection/interface/MuonPlotter.h"
#include "TopTreeAnalysisBase/Selection/interface/JetPlotter.h"
#include "TopTreeAnalysisBase/Selection/interface/VertexPlotter.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"
#include "TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "Style.C"
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "TopTreeAnalysisBase/MCInformation/interface/Lumi3DReWeighting.h"

using namespace std;
using namespace reweight;
using namespace TopTree;

int main (int argc, char *argv[])
{

  string rootFileName = "Output.root";

  clock_t start = clock();

  cout << "********************************************************" << endl;
  cout << " Beginning of the program for creating the BTag Trees ! " << endl;
  cout << "********************************************************" << endl;

  //SetStyle if needed
  //setTDRStyle(); 
  setMyStyle();

  /////////////////////
  // Configuration
  /////////////////////

  //xml file
  string xmlFileName ="../config/mySampleConfig.xml";
  //xmlFileName ="../config/myBTAGconfig_newcalib.xml";
  //string xmlFileName ="../config/myBTAGconfig_fall11.xml";

  if (argc > 3)
    xmlFileName = (string)argv[3];
  
  const char *xmlfile = xmlFileName.c_str();

  cout << "used config file: " << xmlfile << endl;

  //Configuration output format
  TTree *configTree = new TTree("configTree","configuration Tree");
  TClonesArray* tcdatasets = new TClonesArray("Dataset",1000);
  configTree->Branch("Datasets","TClonesArray",&tcdatasets);
  TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000);
  configTree->Branch("AnaEnv","TClonesArray",&tcAnaEnv);

  ////////////////////////////////////
  /// AnalysisEnvironment  
  ////////////////////////////////////

  AnalysisEnvironment anaEnv;
  cout<<"Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);
  int verbose = anaEnv.Verbose;
  float oldLuminosity = anaEnv.Luminosity;	// in 1/pb
 
  cout << "analysis environment luminosity for rescaling "<< oldLuminosity << endl;
  
  /////////////////////
  // Load Datasets
  /////////////////////

  TTreeLoader treeLoader;
  cout << " - Load datasets ..." << endl;
  vector < Dataset* > datasets;
  vector < Dataset* > datasetsMu;
  vector < Dataset* > datasetsEl;
  vector < Dataset* > datasetsPlot;

  treeLoader.LoadDatasets (datasets, xmlfile);
  for(unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
  
  float Luminosity = oldLuminosity;

  float LuminosityMu = oldLuminosity;
  float LuminosityEl = oldLuminosity;

  bool isSemiMu = false;
  bool isSemiE = false;

  bool foundMu = false;
  bool foundEl = false;

  for (unsigned int d = 0; d < datasets.size (); d++) {

    if(Luminosity > datasets[d]->EquivalentLumi() ) Luminosity = datasets[d]->EquivalentLumi();
    string dataSetName = datasets[d]->Name();
    
    if(dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) {
      LuminosityMu = datasets[d]->EquivalentLumi();
      foundMu=true;
    }  
    if(dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) {
      LuminosityEl = datasets[d]->EquivalentLumi();
      foundEl=true;  
    }  

    if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) { // when only using one channel
      Luminosity = datasets[d]->EquivalentLumi();
    }   

    if( dataSetName.find("QCD") == 0 ) datasets[d]->SetColor(kYellow);
    if( dataSetName.find("TT") == 0 ) datasets[d]->SetColor(kRed+1);
    if( dataSetName.find("TTbarJets_Other") == 0 ) datasets[d]->SetColor(kRed-7);
    if( dataSetName.find("WJets") == 0 )
    {
      datasets[d]->SetTitle("W#rightarrowl#nu");
      datasets[d]->SetColor(kGreen-3);
    }
    if( dataSetName.find("ZJets") == 0 )
    {
      datasets[d]->SetTitle("Z/#gamma*#rightarrowl^{+}l^{-}");
      datasets[d]->SetColor(kAzure-2);
    }
    if( dataSetName.find("ST") == 0 || dataSetName.find("SingleTop") ==0 )
      datasets[d]->SetColor(kMagenta);
    
  }

  if(!foundMu && !foundEl && Luminosity != oldLuminosity) cout << "changed analysis environment luminosity to "<< Luminosity << endl;
  else {
    if(LuminosityMu != oldLuminosity) cout << "Muon PD: changed analysis environment luminosity to "<< LuminosityMu << endl;
    if(LuminosityEl != oldLuminosity) cout << "Electron PD: changed analysis environment luminosity to "<< LuminosityEl << endl;
  }

  // make a datasets vector only for 
  if (foundMu) {
      for (unsigned int d = 0; d < datasets.size (); d++) {
	string dataSetName = datasets[d]->Name();
	if ( ! (dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) ) 
	  datasetsMu.push_back(datasets[d]);
      }
  }
      
  if (foundEl) {
    for (unsigned int d = 0; d < datasets.size (); d++) {
      string dataSetName = datasets[d]->Name();
      if ( ! (dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) ) 
	datasetsEl.push_back(datasets[d]);
    }
  }
  
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");

  //Global variable
  //TRootEvent* event = 0;
  
  //nof selected events
  double NEvtsData = 0;
  Double_t *nEvents = new Double_t[datasets.size()];
  
  ////////////////////////////////////
  /// MultiSamplePlot
  ////////////////////////////////////

  map<string,MultiSamplePlot*> MSPlot;

  ////////////////////////////////////
  /// Selection table
  ////////////////////////////////////
      
  vector<string> CutsSelecTableSemiMu;
  CutsSelecTableSemiMu.push_back(string("preselected"));
  CutsSelecTableSemiMu.push_back(string("trigged"));
  CutsSelecTableSemiMu.push_back(string("Good PV"));
  CutsSelecTableSemiMu.push_back(string("1 selected muon"));
  CutsSelecTableSemiMu.push_back(string("Veto 2nd muon"));
  CutsSelecTableSemiMu.push_back(string("Veto electron"));
  
  vector<string> CutsSelecTableSemiEl;
  CutsSelecTableSemiEl.push_back(string("preselected"));
  CutsSelecTableSemiEl.push_back(string("trigged"));
  CutsSelecTableSemiEl.push_back(string("Good PV"));
  CutsSelecTableSemiEl.push_back(string("1 selected electron"));
  CutsSelecTableSemiEl.push_back(string("Veto muon"));
  CutsSelecTableSemiEl.push_back(string("Veto 2nd electron from Z-decay"));
  CutsSelecTableSemiEl.push_back(string("Conversion veto"));
  
  char LabelNJets[100];
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-3);
  CutsSelecTableSemiMu.push_back(string(LabelNJets));
  CutsSelecTableSemiEl.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-2);
  CutsSelecTableSemiMu.push_back(string(LabelNJets));
  CutsSelecTableSemiEl.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-1);
  CutsSelecTableSemiMu.push_back(string(LabelNJets));
  CutsSelecTableSemiEl.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets);
  CutsSelecTableSemiMu.push_back(string(LabelNJets));
  CutsSelecTableSemiEl.push_back(string(LabelNJets));

  if (verbose > 0)
    cout << " - CutsSelectionTable instantiated ..." << endl;
  SelectionTable selecTableSemiMu(CutsSelecTableSemiMu, datasets);
  selecTableSemiMu.SetLuminosity(LuminosityMu);
  SelectionTable selecTableSemiEl(CutsSelecTableSemiEl, datasets);
  selecTableSemiEl.SetLuminosity(LuminosityEl);

  if (verbose > 0)
    cout << " - SelectionTable instantiated ..." << endl;

  ////////////////////////
  // PileUp Reweighting //
  ////////////////////////

  //cout << Luminosity << endl;

  LumiReWeighting LumiWeights, LumiWeightsUp, LumiWeightsDown;

  LumiWeights = LumiReWeighting("PileUpReweighting/pileup_MC_Summer12_S10.root", "PileUpReweighting/pileup_2012Data53X_UpToRun208357/nominal.root", "pileup", "pileup");
  LumiWeightsUp = LumiReWeighting("PileUpReweighting/pileup_MC_Summer12_S10.root", "PileUpReweighting/pileup_2012Data53X_UpToRun208357/sys_up.root", "pileup", "pileup");
  LumiWeightsDown = LumiReWeighting("PileUpReweighting/pileup_MC_Summer12_S10.root", "PileUpReweighting/pileup_2012Data53X_UpToRun208357/sys_down.root", "pileup", "pileup");

  cout << " Initialized LumiReWeighting stuff" << endl;
  
  ////////////////////////////////////
  //	Loop on datasets
  ////////////////////////////////////

  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
  for (unsigned int d = 0; d < datasets.size (); d++) {

    string previousFilename = "";
    int iFile = -1;
    string dataSetName = datasets[d]->Name();

    int nSelectedMu=0;
    int nSelectedEl=0;
    if (verbose > 1)
      cout << "   Dataset " << d << ": " << datasets[d]->Name () << "/ title : " << datasets[d]->Title () << endl;
    if (verbose > 1)
      std::cout<<"      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
    
    //open files and load
    cout<<"LoadEvent"<<endl;
    treeLoader.LoadDataset (datasets[d], anaEnv);
    cout<<"LoadEvent"<<endl;
    
    
    
    /////////////////////////////////////
    /// Initialize JEC factors
    /////////////////////////////////////
   	    
    vector<JetCorrectorParameters> vCorrParam;

    /*JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters("JECFiles/Summer12_V3_MC_L3Absolute_AK5PFchs.txt");
    JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters("JECFiles/Summer12_V3_MC_L2Relative_AK5PFchs.txt");
    JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters("JECFiles/Summer12_V3_MC_L1FastJet_AK5PFchs.txt");
    
    //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
    vCorrParam.push_back(*L1JetPar);
    vCorrParam.push_back(*L2JetPar);
    vCorrParam.push_back(*L3JetPar);

    if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) { // DATA!
      JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("JECFiles/Summer12_V3_DATA_L2L3Residual_AK5PFchs.txt");
      vCorrParam.push_back(*ResJetCorPar);
      }*/
    
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(*(new JetCorrectorParameters("JECFiles/Fall12_V6_DATA_UncertaintySources_AK5PFchs.txt", "Total")));
    
    // true means redo also the L1
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);
    

    ////////////////////////////////////
    //	Loop on events
    ////////////////////////////////////
    
    nEvents[d] = 0;
    int itriggerSemiMu = -1,itriggerSemiEl = -1, previousRun = -1;
    if (verbose > 1)
      cout << "	Loop over events " << endl;

    for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
      //for (unsigned int ievt = 0; ievt < 20000; ievt++)
    {
      
      vector < TRootVertex* > vertex;
      vector < TRootMuon* > init_muons;
      vector < TRootElectron* > init_electrons;
      vector < TRootJet* > init_jets_corrected;
      vector < TRootJet* > init_jets;
      vector < TRootMET* > mets;
      vector < TRootGenJet* > genjets;
      
      nEvents[d]++;
            
      if(ievt%1000 == 0)
	std::cout<<"Processing the "<<ievt<<"th event (" << ((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100  << "%)" << " +> # selected: " << nSelectedMu << " (mu+jets) " << nSelectedEl << " (e+jets)" << flush<<"\r";

      ////////////////
      // LOAD EVENT //
      ////////////////

      TRootEvent* event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets_corrected, mets);  

      if(! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) ) {
        genjets = treeLoader.LoadGenJet(ievt,false);
        sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
      }

      // check with genEvent which ttbar channel it is
      if(dataSetName.find("TTbarJets") == 0)  {
	//cout << "LOADING GenEvent" << endl;
	TRootGenEvent* genEvt = treeLoader.LoadGenEvent(ievt,false);
	if( genEvt->isSemiLeptonic(TRootGenEvent::kMuon) ) {
	  isSemiMu=true;
	  isSemiE=false;
	}
	else if( genEvt->isSemiLeptonic(TRootGenEvent::kElec) ) {
	  isSemiMu=false;
	  isSemiE=true;
	}
	else {
	  isSemiMu=false;
	  isSemiE=false;
	}
      }


      /////////////////////////////////
      // DETERMINE EVENT SCALEFACTOR //
      /////////////////////////////////
      
      // scale factor for the event
      float scaleFactor = 1.;
      
      // Load the GenEvent and calculate the branching ratio correction
      /*if(dataSetName.find("TTbarJets") == 0)
	{
	  //cout << "LOADING GenEvent" << endl;
	  TRootGenEvent* genEvt = treeLoader.LoadGenEvent(ievt,false);
	  if( genEvt->isSemiLeptonic() )
	    scaleFactor *= (0.108*9.)*(0.676*1.5);
	  else if( genEvt->isFullHadronic() )
	    scaleFactor *= (0.676*1.5)*(0.676*1.5);
	  else if( genEvt->isFullLeptonic() )
	    scaleFactor *= (0.108*9.)*(0.108*9.);

      }*/
      
      //////////////////////////////////////
      // Apply Jet Corrections on-the-fly //   
      //////////////////////////////////////

      // not needed for now, GT contains good stuff
      /*if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) {
      	//jetTools->correctJets(init_jets_corrected,event->kt6PFJetsPF2PAT_rho(),true); //last boolean: isData (needed for L2L3Residual...)
      } else {
	jetTools->correctJets(init_jets_corrected,event->kt6PFJets_rho(),false); //last boolean: isData (needed for L2L3Residual...)
	}*/

      // PU reweighting

      // old method

      double lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );

      if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
	lumiWeight=1;
      
      // up syst -> lumiWeight = LumiWeightsUp.ITweight( (int)event->nTruePU() );
      // down syst -> lumiWeight = LumiWeightsDown.ITweight( (int)event->nTruePU() );

      scaleFactor = scaleFactor*lumiWeight;

      ///////////////////
      // TRIGGER SETUP //
      ///////////////////

      string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
      if(previousFilename != currentFilename){
	previousFilename = currentFilename;
	iFile++;
	cout<<"File changed!!! => iFile = "<<iFile<<endl;
      }

      int currentRun = event->runId();

      if(previousRun != currentRun) {
        previousRun = currentRun;
	
	//semi-mu
	if(dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) {
	  
	  // 2.7/fb recalib 
	  if( event->runId() <= 190738 )
	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v11"), currentRun, iFile);

	  else if( event->runId() >= 191043 && event->runId() <= 193621 )
	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v12"), currentRun, iFile);

	  else if( event->runId() >= 193834 && event->runId() <= 196531 )
	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);

	  else if( event->runId() >= 198049 && event->runId() <= 199608)
	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v14"), currentRun, iFile);

	  else if( event->runId() >= 199698 && event->runId() <= 208357)
	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v15"), currentRun, iFile);

	  else {
	    cout << "Unknown run for HLTpath selection: " << event->runId() << endl;
	    exit(1);
	  }

	  if( itriggerSemiMu == 9999 )
          {
            cout << "itriggerSemiMu == 9999 for SemiMu HLTpath selection: " << event->runId() << endl;
            exit(-1);
          }

	} else {
	  itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun); // Summer12 DR53X
	}

	// semi-el
	// semi-electron
        if(dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0 ) {
         
	  
	  // 2.7/fb recalib 
	  if( event->runId() <= 190738 )
	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v8"), currentRun, iFile);
	  
	  else if( event->runId() >= 191043 && event->runId() <= 191411 )
	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v9"), currentRun, iFile);

	  else if( event->runId() >= 191695 && event->runId() <= 196531)
	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v10"), currentRun, iFile);
	  
	  else if( event->runId() >= 198049 && event->runId() <= 208357)
	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v11"), currentRun, iFile);
  
	  else { 
            cout << "Unknown run for SemiEl HLTpath selection: " << event->runId() << endl;
	    exit(1);
	  }
	  if( itriggerSemiEl == 9999 )
	    {
	      cout << "itriggerSemiEl == 9999 for SemiEl HLTpath selection: " << event->runId() << endl;
	      exit(-1);
	    }
        }
        else
        {
	  itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v10"), currentRun); // Summer12 DR53X 
	}
	
      }

      if (itriggerSemiMu == 9999 && itriggerSemiEl == 9999) {
	  
	cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT IN RUN " << event->runId() << endl;

	exit(1);
	  
      }

      /////////////////////////////////////////////////////////////////////////////
      // JES SYSTEMATICS && SMEAR JET RESOLUTION TO MIMIC THE RESOLUTION IN DATA //
      /////////////////////////////////////////////////////////////////////////////

      if( ! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) ) {

	jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "nominal",false);

	//jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "minus",false); //false means don't use old numbers but newer ones...
	//jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "plus",false);

	// Example how to apply JES systematics
	
	 
	//jetTools->correctJetJESUnc(init_jets_corrected, "minus",1);
	//jetTools->correctJetJESUnc(init_jets_corrected, "plus",1);
	
      }

      /////////////////////
      // EVENT SELECTION //
      /////////////////////
      
      //Declare selection instance    
      Selection selection(init_jets_corrected, init_muons, init_electrons, mets, event->kt6PFJets_rho());
      selection.setJetCuts(40,2.5,0.01,1.,0.98,0.3,0.1);
      selection.setMuonCuts(25,2.1,0.12,0.2,0.3,1,0.5,5,0); 
      selection.setElectronCuts(32,2.5,0.1,0.02,0.5,0.3,0); 
      selection.setLooseMuonCuts(10,2.5,0.2);
      selection.setLooseElectronCuts(20,2.5,0.15,0.);

      bool triggedSemiMu = false;
      bool triggedSemiEl = false;

      if( ! (dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) )
        triggedSemiMu = treeLoader.EventTrigged (itriggerSemiMu);
      if( ! (dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) )
        triggedSemiEl = treeLoader.EventTrigged (itriggerSemiEl);

      bool isGoodPV = selection.isPVSelected(vertex, 4, 24, 2.);

      vector<TRootJet*> selectedJets, selectedJetsNoMu, selectedJetsNoEl;
      vector<TRootMuon*> selectedMuons;
      vector<TRootElectron*> selectedElectrons;
      vector<TRootMuon*> vetoMuons = selection.GetSelectedLooseMuons(10,2.5,0.2);
      vector<TRootElectron*> vetoElectronsSemiMu = selection.GetSelectedLooseElectrons(20,2.5,0.15);
      vector<TRootElectron*> vetoElectronsSemiEl = selection.GetSelectedLooseElectrons(20,2.5,0.15);

      if( dataSetName.find("InvIso") != string::npos )  { // event selection for special Data TopTrees for ARI QCD
        vector<TRootMuon*> overlapMuons = selection.GetSelectedMuonsInvIso(0.2, vertex[0]);
        vector<TRootElectron*> overlapElectrons = selection.GetSelectedElectronsInvIso(0.2);
        selectedJetsNoMu = selection.GetSelectedJets(overlapMuons,true);
        selectedJetsNoEl = selection.GetSelectedJets(overlapElectrons,true);

	/*if (selectedJetsNoMu.size() >= 4) {
	  //cout << "ol" << endl;
	  if (selectedJetsNoMu[0]->Pt() < 45) selectedJetsNoMu.clear();
	  if (selectedJetsNoMu[1]->Pt() < 45) selectedJetsNoMu.clear();
	  if (selectedJetsNoMu[2]->Pt() < 40) selectedJetsNoMu.clear();
	  if (selectedJetsNoMu[3]->Pt() < 40) selectedJetsNoMu.clear();
	}

	if (selectedJetsNoEl.size() >= 4) {
	  //cout << "ol" << endl;
	  if (selectedJetsNoEl[0]->Pt() < 45) selectedJetsNoEl.clear();
	  if (selectedJetsNoEl[1]->Pt() < 45) selectedJetsNoEl.clear();
	  if (selectedJetsNoEl[2]->Pt() < 40) selectedJetsNoEl.clear();
	  if (selectedJetsNoEl[3]->Pt() < 40) selectedJetsNoEl.clear();
	  }*/

	//selectedJetsNoMu = selection.GetSelectedJets(true);
        //selectedJetsNoEl = selection.GetSelectedJets(true);

	selectedMuons = selection.GetSelectedMuonsInvIso(0.2, vertex[0], selectedJetsNoMu);
        selectedElectrons = selection.GetSelectedElectronsInvIso(0.2,selectedJetsNoEl);
      }
      else { // Normal selection
	selectedJets = selection.GetSelectedJets(true);

	/*if (selectedJets.size() >= 4) {
	  //cout << "ol" << endl;
	  if (selectedJets[0]->Pt() < 45) selectedJets.clear();
	  if (selectedJets[1]->Pt() < 45) selectedJets.clear();
	  if (selectedJets[2]->Pt() < 40) selectedJets.clear();
	  if (selectedJets[3]->Pt() < 40) selectedJets.clear();
	  }*/

	//selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);	
	selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
	selectedElectrons = selection.GetSelectedElectrons(selectedJets);
      }

      vector<TRootMCParticle*> mcParticles;
      
      if(dataSetName.find("TTbarJets") == 0)
      {
        treeLoader.LoadMCEvent(ievt, 0, 0, mcParticles,false);  
        sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      }

      bool eventselectedSemiMu = false;
      bool eventselectedSemiEl = false;

            // semi-mu selection
      if (triggedSemiMu) {
	if (isGoodPV) {
	  if(dataSetName.find("InvIso") != string::npos) selectedJets = selectedJetsNoMu;
	  if (selectedMuons.size() == 1) {
	    if( vetoMuons.size() == 1 || ( dataSetName.find("InvIso") != string::npos && vetoMuons.size() == 0 ) ) { // if InvertedIso, selected muon not part of vetoMuons vector!
	      if (vetoElectronsSemiMu.size() == 0) {
		if (selectedJets.size() >= 4) {
		  eventselectedSemiMu = true;
		}
	      }
	    }
	  }
	}
      }

     selecTableSemiEl.Fill(d,0,scaleFactor*lumiWeight);

     if( triggedSemiEl) {
       if(dataSetName.find("InvIso") != string::npos) selectedJets = selectedJetsNoEl;
       selecTableSemiEl.Fill(d,1,scaleFactor*lumiWeight);
       if (isGoodPV ) {
	 selecTableSemiEl.Fill(d,2,scaleFactor*lumiWeight);
	 if( selectedElectrons.size() == 1 ) {
	   selecTableSemiEl.Fill(d,3,scaleFactor*lumiWeight);
	   if( vetoMuons.size() == 0 ) {
	     selecTableSemiEl.Fill(d,4,scaleFactor*lumiWeight);
	     if (vetoElectronsSemiEl.size() == 1) {
	       selecTableSemiEl.Fill(d,5,scaleFactor*lumiWeight);
	       selecTableSemiEl.Fill(d,6,scaleFactor*lumiWeight);
	       if( selectedJets.size()>=1 ) {
		 selecTableSemiEl.Fill(d,7,scaleFactor*lumiWeight);
		 if( selectedJets.size()>=2 ) {
		   selecTableSemiEl.Fill(d,8,scaleFactor*lumiWeight);
		   if( selectedJets.size()>=3 ) {
		     selecTableSemiEl.Fill(d,9,scaleFactor*lumiWeight);
		     if( selectedJets.size()>=4 ) {
		       selecTableSemiEl.Fill(d,10,scaleFactor*lumiWeight);
		       eventselectedSemiEl=true;
		     }
		   }
		 }
	       }
	     }
	   }
	 }
       }
     }
     
     if (!eventselectedSemiMu && !eventselectedSemiEl) continue;
     
     if (eventselectedSemiMu)
       nSelectedMu++;
     if (eventselectedSemiEl)
       nSelectedEl++;
     
     TLorentzVector* selectedLepton;
     
     if (eventselectedSemiMu)
       selectedLepton = (TLorentzVector*)selectedMuons[0];
     else if (eventselectedSemiEl)
       selectedLepton = (TLorentzVector*)selectedElectrons[0];
     


     //-----------------//
     // do some data-mc //
     //-----------------//
     
     // when running both electron and muon data, pick the right dataset vector and lumi for the MSPlots
     
     if (!foundMu && !foundEl)
       datasetsPlot = datasets;
     else if (eventselectedSemiMu) {
       datasetsPlot = datasetsMu;
       Luminosity = LuminosityMu;
     }
     else if (eventselectedSemiEl) {
       datasetsPlot = datasetsEl;
       Luminosity = LuminosityEl;
     }
     
     string leptonFlav="_other";
     
     if (eventselectedSemiMu)
       leptonFlav="_mu";
     else if (eventselectedSemiEl)
       leptonFlav="_el";
     
     if (MSPlot.find("Selected_Events_pT_jet1"+leptonFlav) == MSPlot.end()){
       MSPlot["Selected_Events_pT_jet1"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_jet1"+leptonFlav, 30, 0, 600, "p_{T} (GeV)");
       MSPlot["Selected_Events_pT_jet2"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_jet2"+leptonFlav, 30, 0, 600, "p_{T} (GeV)");
       MSPlot["Selected_Events_pT_jet3"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_jet3"+leptonFlav, 30, 0, 600, "p_{T} (GeV)");
       MSPlot["Selected_Events_pT_jet4"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_jet4"+leptonFlav, 30, 0, 600, "p_{T} (GeV)");
       //MSPlot["Selected_Events_pT_jet4"] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_jet4", 30, 0, 600, "p_{T} (GeV)");
       MSPlot["Selected_Events_pT_4leadingjets"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_4leadingjets"+leptonFlav, 30, 0, 600, "p_{T} (GeV)");
       MSPlot["Selected_Events_pT_alljets"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_alljets"+leptonFlav, 30, 0, 600, "p_{T} (GeV)");
     }
     
     MSPlot["Selected_Events_pT_jet1"+leptonFlav]->Fill(selectedJets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
     MSPlot["Selected_Events_pT_jet2"+leptonFlav]->Fill(selectedJets[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
     MSPlot["Selected_Events_pT_jet3"+leptonFlav]->Fill(selectedJets[2]->Pt(), datasets[d], true, Luminosity*scaleFactor);
     MSPlot["Selected_Events_pT_jet4"+leptonFlav]->Fill(selectedJets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
     
     for (unsigned int q=0; q<selectedJets.size(); q++) {
       
       MSPlot["Selected_Events_pT_alljets"+leptonFlav]->Fill(selectedJets[q]->Pt(), datasets[d], true, Luminosity*scaleFactor);
       
       if (q<4)
	 MSPlot["Selected_Events_pT_4leadingjets"+leptonFlav]->Fill(selectedJets[q]->Pt(), datasets[d], true, Luminosity*scaleFactor);
       
     }


     ///////////////////////////////////////
     // END OF EVENT REMOVING SOME STUFF //
     //////////////////////////////////////
     
    }			//loop on events

    cout<<endl;

    cout << "+> " << nSelectedMu << " mu+jets events where selected"<< endl;
    cout << "+> " << nSelectedEl << " e+jets events where selected"<< endl;

    
    //////////////
    // CLEANING //
    //////////////

    if (jecUnc) delete jecUnc;
    if (jetTools) delete jetTools;
  
    //important: free memory
    treeLoader.UnLoadDataset();
    
    }				//loop on datasets

    //Once everything is filled ...
    if (verbose > 0)
      cout << " We ran over all the data ;-)" << endl;
 
    /////////////////////////
    // Write out the plots //
    /////////////////////////
    
    mkdir("DemoPlots",0777);
   
    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
      {
	MultiSamplePlot *temp = it->second;
	string name = it->first;
	temp->Draw(false, name, true, true, true, true, true, 1, false);
	temp->Write(fout, name, true, "DemoPlots/");
      }
    

  
    //Selection tables
    selecTableSemiEl.TableCalculator(false, true, true, true, true);
    string selectiontableEl = "SelectionTable_BTAG_SEMIEL.tex";
    selecTableSemiEl.Write(selectiontableEl.c_str());
        
    // Do some special things with certain plots (normalize, BayesDivide, ... )
    if (verbose > 0)
      cout << "Treating the special plots." << endl;
    
    
    delete fout;
    delete tcdatasets;
    delete tcAnaEnv;
    delete configTree;
        
    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "           hasn't crashed yet ;-)           " << endl;
    cout << "********************************************" << endl;
    
    return 0;
}
