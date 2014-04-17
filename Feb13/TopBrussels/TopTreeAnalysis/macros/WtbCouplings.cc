////////////////////////////
///// TODO & COMMENTS /////
/////////////////////////// 

#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <Math/VectorUtil.h>

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

template <typename T> string tostr(const T& t) { ostringstream os; os<<t; return os.str(); }

int main (int argc, char *argv[])
{

  string rootFileName = "WtbCouplings_output.root";

  clock_t start = clock();

  cout << "********************************************************" << endl;
  cout << " Beginning of the programme for anomalous couplings ! " << endl;
  cout << "********************************************************" << endl;

  //SetStyle if needed
  //setTDRStyle(); 
  setMyStyle();

  /////////////////////
  // Configuration
  /////////////////////

	bool useMassesAndResolutions = true;
	bool calculateTransferFunctions = true; // true will have only an effect is useMassesAndResolutions = false

  //xml file
  string xmlFileName ="../config/myWtbCouplings.xml"; // can contain both electron and muon datasets at the same time!

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
  bool isSemiEl = false;
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

    if( dataSetName.find("QCD") == 0 ) datasets[d]->SetColor(kYellow);
    if( dataSetName.find("TT") == 0 ) datasets[d]->SetColor(kRed+1);
    if( dataSetName.find("TTbarJets_Other") == 0 ) datasets[d]->SetColor(kRed-7);
    if( dataSetName.find("W_Jets") == 0 )
    {
      datasets[d]->SetTitle("W#rightarrowl#nu");
      datasets[d]->SetColor(kGreen-3);
    }
    if( dataSetName.find("ZJets") == 0 )
    {
      datasets[d]->SetTitle("Z/#gamma*#rightarrowl^{+}l^{-}");
      datasets[d]->SetColor(kBlue-2);
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
			if ( ! (dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) ) {
				datasetsMu.push_back(datasets[d]);
				cout << "datasets[d]->Name() " << datasets[d]->Name() << endl;
			}
		}
	}
      
  if (foundEl) {
    for (unsigned int d = 0; d < datasets.size (); d++) {
      string dataSetName = datasets[d]->Name();
      if ( ! (dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) ) 
				datasetsEl.push_back(datasets[d]);
    }
  }
  

  cout << " - Recreate output file ..." << endl;
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");

  //Global variable
  //TRootEvent* event = 0;
  
  //nof selected events
  //double NEvtsData = 0;
  Double_t *nEvents = new Double_t[datasets.size()];
  
  ////////////////////////////////////
  /// Normal Plots (TH1F* and TH2F*)
  ////////////////////////////////////
  
  map<string,TH1F*> histo1D;
  map<string,TH2F*> histo2D;
  
  histo1D["hadronicGenTopMass"] = new TH1F("hadronicGenTopMass","Hadronic Top Mass, using the quarks", 100, 100, 300);
  histo1D["hadronicGenWMass"] = new TH1F("hadronicGenWMass","Hadronic W Mass, using the quarks",100,0,200);
  histo1D["hadronicRecoTopMass"] = new TH1F("hadronicRecoTopMass","Hadronic Top Mass, using the RecoJets", 100, 100, 300);
  histo1D["hadronicRecoWMass"] = new TH1F("hadronicRecoWMass","Hadronic W Mass, using the RecoJets",100,0,200);

	histo2D["Eparton_vs_Enonbjet"] = new TH2F("Eparton_vs_Enonbjet","Eparton_vs_Enonbjet", 50, 0, 400, 50, 0, 400);
	histo2D["Eparton_vs_Eparton-Enonbjet"] = new TH2F("Eparton_vs_Eparton-Enonbjet","Eparton_vs_Eparton-Enonbjet", 4, 0, 250, 50, -100, 100);
	histo2D["Eparton_vs_Ebjet"] = new TH2F("Eparton_vs_Ebjet","Eparton_vs_Ebjet", 50, 0, 400, 50, 0, 400);
	histo2D["Eparton_vs_Eparton-Ebjet"] = new TH2F("Eparton_vs_Eparton-Ebjet","Eparton_vs_Eparton-Ebjet", 4, 0, 300, 50, -100, 100);
	histo2D["EgenEl_vs_ErecEl"] = new TH2F("EgenEl_vs_ErecEl","EgenEl_vs_ErecEl", 50, 0, 400, 50, 0, 400);
	histo2D["EgenEl_vs_EgenEl-ErecEl"] = new TH2F("EgenEl_vs_EgenEl-ErecEl","EgenEl_vs_EgenEl-ErecEl", 50, 0, 400, 50, -50, 50);
	histo2D["InvPtgenMu_vs_InvPtrecMu"] = new TH2F("InvPtgenMu_vs_InvPtrecMu","InvPtgenMu_vs_InvPtrecMu", 50, 0, 0.05, 50, 0, 0.05);
	histo2D["InvPtgenMu_vs_InvPtgenMu-InvPtrecMu"] = new TH2F("InvPtgenMu_vs_InvPtgenMu-InvPtrecMu","InvPtgenMu_vs_InvPtgenMu-InvPtrecMu", 5, 0.002, 0.04, 50, -0.005, 0.005);
	
	histo2D["Thparton_vs_Thnonbjet"] = new TH2F("Thparton_vs_Thnonbjet","Thparton_vs_Thnonbjet", 60, 0, 3.15, 60, 0, 3.15);
	histo2D["Thparton_vs_Thparton-Thnonbjet"] = new TH2F("Thparton_vs_Thparton-Thnonbjet","Thparton_vs_Thparton-Thnonbjet", 5, 0.1, 3.0, 60, -0.15, 0.15);
	histo2D["Thparton_vs_Thbjet"] = new TH2F("Thparton_vs_Thbjet","Thparton_vs_Thbjet", 60, 0, 3.15, 60, 0, 3.15);
	histo2D["Thparton_vs_Thparton-Thbjet"] = new TH2F("Thparton_vs_Thparton-Thbjet","Thparton_vs_Thparton-Thbjet", 5, 0.1, 3.0, 60, -0.15, 0.15);
	histo2D["ThgenEl_vs_ThrecEl"] = new TH2F("ThgenEl_vs_ThrecEl","ThgenEl_vs_ThrecEl", 60, 0, 3.15, 60, 0, 3.15);
	histo2D["ThgenEl_vs_ThgenEl-ThrecEl"] = new TH2F("ThgenEl_vs_ThgenEl-ThrecEl","ThgenEl_vs_ThgenEl-ThrecEl", 5, 0.2, 3.0, 60, -0.05, 0.05);
	histo2D["ThgenMu_vs_ThrecMu"] = new TH2F("ThgenMu_vs_ThrecMu","ThgenMu_vs_ThrecMu", 60, 0, 3.15, 60, 0, 3.15);
	histo2D["ThgenMu_vs_ThgenMu-ThrecMu"] = new TH2F("ThgenMu_vs_ThgenMu-ThrecMu","ThgenMu_vs_ThgenMu-ThrecMu", 5, 0.2, 3.0, 60, -0.05, 0.05);
	histo2D["Eparton_vs_Thnonbjet"] = new TH2F("Eparton_vs_Thnonbjet","Eparton_vs_Thnonbjet", 50, 0, 400, 60, 0, 3.15);
	histo2D["Eparton_vs_Thparton-Thnonbjet"] = new TH2F("Eparton_vs_Thparton-Thnonbjet","Eparton_vs_Thparton-Thnonbjet", 4, 0, 250, 60, -0.15, 0.15);
	histo2D["Eparton_vs_Thbjet"] = new TH2F("Eparton_vs_Thbjet","Eparton_vs_Thbjet", 50, 0, 400, 60, 0, 3.15);
	histo2D["Eparton_vs_Thparton-Thbjet"] = new TH2F("Eparton_vs_Thparton-Thbjet","Eparton_vs_Thparton-Thbjet", 4, 0, 300, 60, -0.15, 0.15);
	histo2D["EgenEl_vs_ThrecEl"] = new TH2F("EgenEl_vs_ThrecEl","EgenEl_vs_ThrecEl", 50, 0, 400, 60, 0, 3.15);
	histo2D["EgenEl_vs_ThgenEl-ThrecEl"] = new TH2F("EgenEl_vs_ThgenEl-ThrecEl","EgenEl_vs_ThgenEl-ThrecEl", 5, 0, 400, 60, -0.05, 0.05);
	histo2D["InvPtgenMu_vs_ThrecMu"] = new TH2F("InvPtgenMu_vs_ThrecMu","InvPtgenMu_vs_ThrecMu", 50, 0, 0.05, 60, 0, 3.15);
	histo2D["InvPtgenMu_vs_ThgenMu-ThrecMu"] = new TH2F("InvPtgenMu_vs_ThgenMu-ThrecMu","InvPtgenMu_vs_ThgenMu-ThrecMu", 5, 0.002, 0.04, 60, -0.05, 0.05);

	histo2D["Phiparton_vs_Phinonbjet"] = new TH2F("Phiparton_vs_Phinonbjet","Phiparton_vs_Phinonbjet", 120, -3.2, 3.2, 120, -3.2, 3.2);
	histo2D["Phiparton_vs_Phiparton-Phinonbjet"] = new TH2F("Phiparton_vs_Phiparton-Phinonbjet","Phiparton_vs_Phiparton-Phinonbjet", 5, -3.2, 3.2, 120, -0.3, 0.3);
	histo2D["Phiparton_vs_Phibjet"] = new TH2F("Phiparton_vs_Phibjet","Phiparton_vs_Phibjet", 120, -3.2, 3.2, 120, -3.2, 3.2);
	histo2D["Phiparton_vs_Phiparton-Phibjet"] = new TH2F("Phiparton_vs_Phiparton-Phibjet","Phiparton_vs_Phiparton-Phibjet", 5, -3.2, 3.2, 120, -0.3, 0.3);
	histo2D["PhigenEl_vs_PhirecEl"] = new TH2F("PhigenEl_vs_PhirecEl","PhigenEl_vs_PhirecEl", 120, -3.2, 3.2, 120, -3.2, 3.2);
	histo2D["PhigenEl_vs_PhigenEl-PhirecEl"] = new TH2F("PhigenEl_vs_PhigenEl-PhirecEl","PhigenEl_vs_PhigenEl-PhirecEl", 5, -3.2, 3.2, 120, -0.05, 0.05);
	histo2D["PhigenMu_vs_PhirecMu"] = new TH2F("PhigenMu_vs_PhirecMu","PhigenMu_vs_PhirecMu", 120, -3.2, 3.2, 120, -3.2, 3.2);
	histo2D["PhigenMu_vs_PhigenMu-PhirecMu"] = new TH2F("PhigenMu_vs_PhigenMu-PhirecMu","PhigenMu_vs_PhigenMu-PhirecMu", 5, -3.2, 3.2, 120, -0.05, 0.05);
	histo2D["Eparton_vs_Phinonbjet"] = new TH2F("Eparton_vs_Phinonbjet","Eparton_vs_Phinonbjet", 50, 0, 400, 120, -3.2, 3.2);
	histo2D["Eparton_vs_Phiparton-Phinonbjet"] = new TH2F("Eparton_vs_Phiparton-Phinonbjet","Eparton_vs_Phiparton-Phinonbjet", 5, 0, 250, 120, -0.3, 0.3);
	histo2D["Eparton_vs_Phibjet"] = new TH2F("Eparton_vs_Phibjet","Eparton_vs_Phibjet", 50, 0, 400, 120, -3.2, 3.2);
	histo2D["Eparton_vs_Phiparton-Phibjet"] = new TH2F("Eparton_vs_Phiparton-Phibjet","Eparton_vs_Phiparton-Phibjet", 5, 0, 300, 120, -0.3, 0.3);
	histo2D["EgenEl_vs_PhirecEl"] = new TH2F("EgenEl_vs_PhirecEl","EgenEl_vs_PhirecEl", 50, 0, 400, 120, -3.2, 3.2);
	histo2D["EgenEl_vs_PhigenEl-PhirecEl"] = new TH2F("EgenEl_vs_PhigenEl-PhirecEl","EgenEl_vs_PhigenEl-PhirecEl", 5, 0, 400, 120, -0.05, 0.05);
	histo2D["InvPtgenMu_vs_PhirecMu"] = new TH2F("InvPtgenMu_vs_PhirecMu","InvPtgenMu_vs_PhirecMu", 50, 0, 0.05, 120, -3.2, 3.2);
	histo2D["InvPtgenMu_vs_PhigenMu-PhirecMu"] = new TH2F("InvPtgenMu_vs_PhigenMu-PhirecMu","InvPtgenMu_vs_PhigenMu-PhirecMu", 5, 0.002, 0.04, 120, -0.05, 0.05);

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
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets);
  CutsSelecTableSemiMu.push_back(string(LabelNJets));
  CutsSelecTableSemiEl.push_back(string(LabelNJets));

  CutsSelecTableSemiMu.push_back(">= 1 b-jet (CSVM)");
  CutsSelecTableSemiMu.push_back(">= 2 b-jets (CSVM)");
  CutsSelecTableSemiEl.push_back(">= 1 b-jet (CSVM)");
  CutsSelecTableSemiEl.push_back(">= 2 b-jets (CSVM)");

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
    cout<< "equivalent luminosity of dataset "<< datasets[d]->EquivalentLumi() << endl;
    string previousFilename = "";
    int iFile = -1;
    string dataSetName = datasets[d]->Name();

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

    JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters("JECFiles/Summer12_V3_MC_L3Absolute_AK5PFchs.txt");
    JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters("JECFiles/Summer12_V3_MC_L2Relative_AK5PFchs.txt");
    JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters("JECFiles/Summer12_V3_MC_L1FastJet_AK5PFchs.txt");
    
    //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
    vCorrParam.push_back(*L1JetPar);
    vCorrParam.push_back(*L2JetPar);
    vCorrParam.push_back(*L3JetPar);

    if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) { // DATA!
      JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("JECFiles/Summer12_V3_DATA_L2L3Residual_AK5PFchs.txt");
      vCorrParam.push_back(*ResJetCorPar);
    }
    
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(*(new JetCorrectorParameters("JECFiles/Fall12_V6_DATA_UncertaintySources_AK5PFchs.txt", "Total")));
    
    // true means redo also the L1
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);
    


		////////////////////////////////////////////////////////
		// GET MASSES AND RESOLUTIONS FOR CHI2 JETCOMBINATION //					
		////////////////////////////////////////////////////////

		// load the mass values and resolutions for chi2 jetcombination
		float Chi2Wmass = -9999;
		float SigmaChi2Wmass = -9999;
		float Chi2Topmass = -9999;
		float SigmaChi2Topmass = -9999;

		if (useMassesAndResolutions) {

			//cout << "Test" << endl;
			string filename = "MassPlotsForChi2JetCombiner.root";

			cout << "INFO: Using masses and widths for the Chi2 jetcombiner from: " << filename << endl;

			TFile* res = new TFile(filename.c_str(),"READ");
			res->cd();
	
			TF1* WmassFit = (TF1*)res->Get("hadronicRecoWMass_Fitted");
			if (WmassFit) {
				cout << "Fitted Wmass: " << WmassFit->GetParameter(1) << "+-" << WmassFit->GetParameter(2) << endl;
				Chi2Wmass = WmassFit->GetParameter(1); 
				SigmaChi2Wmass = WmassFit->GetParameter(2); 	
			}

			TF1* TopmassFit = (TF1*)res->Get("hadronicRecoTopMass_Fitted");
			if (TopmassFit) {
				cout << "Fitted Topmass: " << TopmassFit->GetParameter(1) << "+-" << TopmassFit->GetParameter(2) << endl;
				Chi2Topmass = TopmassFit->GetParameter(1); 
				SigmaChi2Topmass = TopmassFit->GetParameter(2);
			}
			res->Close();
			delete res;
		
		}


    ////////////////////////////////////
    //	Loop on events
    ////////////////////////////////////
    
    nEvents[d] = 0;
    int itriggerSemiMu = -1,itriggerSemiEl = -1, previousRun = -1;
    if (verbose > 1)
      cout << "	Loop over events " << endl;

    for (int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
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
	std::cout<<"Processing the "<<ievt<<"th event (" << ((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100  << "%)" << flush<<"\r";
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
	  			isSemiEl=false;
				}
				else if( genEvt->isSemiLeptonic(TRootGenEvent::kElec) ) {
	  			isSemiMu=false;
	  			isSemiEl=true;
				}
				else {
	  			isSemiMu=false;
	  			isSemiEl=false;
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
      //cout<< "scalefactor " << scaleFactor << endl; 
      double lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );

      if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
				lumiWeight=1;
      
      // up syst -> lumiWeight = LumiWeightsUp.ITweight( (int)event->nTruePU() );
      // down syst -> lumiWeight = LumiWeightsDown.ITweight( (int)event->nTruePU() );

      scaleFactor = scaleFactor*lumiWeight;
      //cout << "scalefactor after lumiweight " << scaleFactor << endl;
      //cout << "lumiweight " << lumiWeight << endl;
       
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
			 	 	if( itriggerSemiMu == 9999 ) {
						cout << "itriggerSemiMu == 9999 for SemiMu HLTpath selection: " << event->runId() << endl;
						exit(-1);
					}
				} else {
					itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun); // Summer12 DR53X
					if( itriggerSemiMu == 9999 )
						itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v11"), currentRun); // Summer12 DR53X
				}

				// semi-electron
  		  if(dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0 ) {
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
	  			if( itriggerSemiEl == 9999 ) {
	    		  cout << "itriggerSemiEl == 9999 for SemiEl HLTpath selection: " << event->runId() << endl;
	   	 	  exit(-1);
	   		 }
				} else {
					itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v10"), currentRun); // Summer12 DR53X
					if( itriggerSemiEl == 9999 )
						itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v8"), currentRun); // Summer12 DR53X  
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
				//cout<<"JER smeared!!! "<<endl;
	
      }

      /////////////////////
      // EVENT SELECTION //
      /////////////////////
      
      //Declare selection instance    
      Selection selection(init_jets_corrected, init_muons, init_electrons, mets, event->kt6PFJets_rho());

      selection.setJetCuts(30.,2.5,0.01,1.,0.98,0.3,0.1);
      selection.setMuonCuts(35,2.1,0.12,0.2,0.3,1,0.5,5,0); // DR mu-jets cleaning still needed?
      selection.setElectronCuts(32,2.5,0.1,0.02,0.5,0.3,0); // Ele25 uses 30 so Ele27 32? DR el-jets cleaning still needed?
      selection.setLooseMuonCuts(10,2.5,0.2);
      selection.setLooseElectronCuts(20,2.5,0.15,0.); // semiMu looseElectron cuts

      bool triggedSemiMu = false;
      bool triggedSemiEl = false;

      if( ! (dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) )
        triggedSemiMu = treeLoader.EventTrigged (itriggerSemiMu);
      if( ! (dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) )
        triggedSemiEl = treeLoader.EventTrigged (itriggerSemiEl);

      bool isGoodPV = false;
        
      isGoodPV = selection.isPVSelected(vertex, anaEnv.PVertexNdofCut, anaEnv.PVertexZCut, anaEnv.PVertexRhoCut);

      vector<TRootJet*> selectedJets, selectedJetsNoMu, selectedJetsNoEl;
      vector<TRootMuon*> selectedMuons;
      vector<TRootElectron*> selectedElectrons;
      vector<TRootMuon*> vetoMuons = selection.GetSelectedLooseMuons();
      vector<TRootElectron*> vetoElectronsSemiMu = selection.GetSelectedLooseElectrons(30,2.5,0.2);
      vector<TRootElectron*> vetoElectronsSemiEl = selection.GetSelectedLooseElectrons(30,2.5,0.2);
      
      selectedJets = selection.GetSelectedJets(true);

      if (selectedJets.size() >= 4)
				if (selectedJets[3]->Pt() < 30) selectedJets.clear();

      selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
      selectedElectrons = selection.GetSelectedElectrons(selectedJets);

			//cout<<"blabla "<<endl;
      vector<TRootMCParticle*> mcParticles;
      
      if(dataSetName.find("TTbarJets") == 0)
      {
        treeLoader.LoadMCEvent(ievt, 0, 0, mcParticles,false);  
        sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
			//cout<<"blabli "<<endl;

      bool eventselectedSemiMu = false; //2 CSVM btags
      bool eventselectedSemiEl = false; //2 CSVM btags
      bool eventselectedSemiMu_onebtag = false;
      bool eventselectedSemiEl_onebtag = false;
			int nb_btags = 0;
      
    	selecTableSemiMu.Fill(d,0,scaleFactor);
			// semi-mu selection
			if (triggedSemiMu){
    		selecTableSemiMu.Fill(d,1,scaleFactor);
				if(isGoodPV) {
    			selecTableSemiMu.Fill(d,2,scaleFactor);
					if (selectedMuons.size() == 1) {
    				selecTableSemiMu.Fill(d,3,scaleFactor);
	  				if (vetoMuons.size() == 1) {
    					selecTableSemiMu.Fill(d,4,scaleFactor);
	    				if (vetoElectronsSemiMu.size() == 0) {
    						selecTableSemiMu.Fill(d,5,scaleFactor);
		 						if( selectedJets.size()>=4 ) {
		   						selecTableSemiMu.Fill(d,6,scaleFactor);
									//cout<<"blabli2 "<<endl;
			 
			 						for (unsigned int i=0; i < selectedJets.size(); i++) {
			   						if (selectedJets[i]->btag_combinedSecondaryVertexBJetTags() > 0.679) 
			     						nb_btags++;
									}
		 							
									if( nb_btags>=1 ) {
		   							selecTableSemiMu.Fill(d,7,scaleFactor);
										
										eventselectedSemiMu_onebtag = true;
										
										if( nb_btags>=2 ) {
		   								selecTableSemiMu.Fill(d,8,scaleFactor);
											eventselectedSemiMu = true;
										} // 2 btags
									} // at least 1 btag
								} // at least 4 jets
							} //no electrons
						} //no loose muons
					} // one good muon
				} //good PV
			} //trigger semimu
				//cout<<"nb_btags: "<<nb_btags<<endl;
				//cout<<"eventselectedSemiMu_onebtag: "<<eventselectedSemiMu_onebtag<<endl;
				//cout<<"eventselectedSemiMu: "<<eventselectedSemiMu<<endl;
									

/*     selecTableSemiEl.Fill(d,0,scaleFactor);

     if( triggedSemiEl) {
       selecTableSemiEl.Fill(d,1,scaleFactor);
       if (isGoodPV ) {
	 selecTableSemiEl.Fill(d,2,scaleFactor);
	 if( selectedElectrons.size() == 1 ) {
	   selecTableSemiEl.Fill(d,3,scaleFactor);
	   if( vetoMuons.size() == 0 ) {
	     selecTableSemiEl.Fill(d,4,scaleFactor);
	     if (vetoElectronsSemiEl.size() == 1) {
	       //if( !selection.foundZCandidate(selectedElectrons[0], vetoElectronsSemiEl) ) {
	       selecTableSemiEl.Fill(d,5,scaleFactor);
	       if( selection.passConversionRejection(selectedElectrons[0]) ) {
		 selecTableSemiEl.Fill(d,6,scaleFactor);
		 if( selectedJets.size()>=1 ) {
		   selecTableSemiEl.Fill(d,7,scaleFactor);
		   if( selectedJets.size()>=2 ) {
		     selecTableSemiEl.Fill(d,8,scaleFactor);
		     if( selectedJets.size()>=3 ) {
		       selecTableSemiEl.Fill(d,9,scaleFactor);
		       if( selectedJets.size()>=4 ) {
			 selecTableSemiEl.Fill(d,10,scaleFactor);

			 //if (selectedJets[0]->Pt() > 45 && selectedJets[1]->Pt() > 45 && selectedJets[2]->Pt() > 45 && selectedJets[3]->Pt() > 45) {
			   eventselectedSemiEl=true;
			 //}
		       }
		     }
		   }
		 }
	       }
	     }
	   }
	 }
       }
     }
*/

     
			//if (!eventselectedSemiMu_onebtag && !eventselectedSemiEl_onebtag) continue;
			if (!eventselectedSemiMu && !eventselectedSemiEl) continue;
			//cout << "check - event is selected!" << endl;
     
			////////////////////////////////////
			// JET PARTON MATCHING FOR MASSES //					
			////////////////////////////////////

			// this is only necessary for the chi2 jetcomb input masses. if usemassesandresolutions == true this will not be run

			int MCPermutation[4]; 

			bool all4PartonsMatched = false; // True if the 4 ttbar semi-lep partons are matched to 4 jets (not necessarily the 4 highest pt jets)
			bool all4JetsMatched_MCdef_ = false; // True if the 4 highest pt jets are matched to the 4 ttbar semi-lep partons
			bool hadronictopJetsMatched_MCdef_ = false;

			pair<unsigned int, unsigned int> leptonicBJet_ = pair<unsigned int,unsigned int>(9999,9999);
			pair<unsigned int, unsigned int> hadronicBJet_ = pair<unsigned int,unsigned int>(9999,9999);
			pair<unsigned int, unsigned int> hadronicWJet1_ = pair<unsigned int,unsigned int>(9999,9999);
			pair<unsigned int, unsigned int> hadronicWJet2_ = pair<unsigned int,unsigned int>(9999,9999);

		//	double relDiffEJetParton_b_ = -9999;
		//	double relDiffEJetParton_l1_ = -9999;
		//	double relDiffEJetParton_l2_ = -9999;
	
			int pdgID_top = 6; //top quark

			vector<TRootMCParticle*> mcParticlesMatching_;
			int genmuon = -9999; int genelectron = -9999;
			bool muonmatched = false; bool electronmatched = false;
			if (!useMassesAndResolutions && dataSetName.find("TTbarJets") == 0 && (isSemiMu || isSemiEl) ) {
	
				sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class)     
	
				vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV;
				TLorentzVector topQuark, antiTopQuark;
	
				bool muPlusFromTop = false, muMinusFromTop = false, elPlusFromTop = false, elMinusFromTop = false;
			//	int nTTbarQuarks = 0;
				mcParticlesMatching_.clear();
		
				for(unsigned int i=0; i<mcParticles.size(); i++) {
					//cout << i << ":  status: " << mcParticles[i]->status() << "  pdgId: " << mcParticles[i]->type()
					//  << "  motherPdgId: " << mcParticles[i]->motherType() << "  grannyPdgId: " << mcParticles[i]->grannyType() << endl;
					if( mcParticles[i]->status() != 3) continue;		// 0: empty line; 1: undecayed particle, stable in the generator; 2: particle decayed in the generator; 3: documentation line.
		
					if( mcParticles[i]->type() == pdgID_top )
						topQuark = *mcParticles[i];
					else if( mcParticles[i]->type() == -pdgID_top )
						antiTopQuark = *mcParticles[i];
		
					if( mcParticles[i]->type() == 13 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -pdgID_top ){		// mu-, W-, tbar
						muMinusFromTop = true;
						genmuon = i;
					}
					if( mcParticles[i]->type() == -13 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == pdgID_top ){		// mu+, W+, t
						muPlusFromTop = true;
						genmuon = i;
					}
					if( mcParticles[i]->type() == 11 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -pdgID_top ){ 	// e-, W-, tbar
						elMinusFromTop = true;
						genelectron = i;
					}
					if( mcParticles[i]->type() == -11 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == pdgID_top ){		// e+, W+, t
						elPlusFromTop = true;
						genelectron = i;
					}
		
					if( abs(mcParticles[i]->type()) < 6 || abs(mcParticles[i]->type()) == 21 ) {  //light/b quarks, 6 should stay hardcoded, OR gluon
						mcParticlesTLV.push_back(*mcParticles[i]);
						mcParticlesMatching_.push_back(mcParticles[i]);
					}
				}
	
				// take all the selectedJets_ to study the radiation stuff, selectedJets_ are already ordened in decreasing Pt()
				for(unsigned int i=0; i<selectedJets.size(); i++)
					selectedJetsTLV.push_back(*selectedJets[i]);
		
				//cout << "will do the jet parton matching now" << endl;
				JetPartonMatching matching = JetPartonMatching(mcParticlesTLV, selectedJetsTLV, 2, true, true, 0.3);		// partons, jets, choose algorithm, use maxDist, use dR, set maxDist=0.3
	
				if(matching.getNumberOfAvailableCombinations() != 1)
					cerr << "matching.getNumberOfAvailableCombinations() = "<<matching.getNumberOfAvailableCombinations()<<" .  This should be equal to 1 !!!"<<endl;
	
				vector< pair<unsigned int, unsigned int> > JetPartonPair, ISRJetPartonPair; // First one is jet number, second one is mcParticle number
		
				for(unsigned int i=0; i<mcParticlesTLV.size(); i++) {
					int matchedJetNumber = matching.getMatchForParton(i, 0);
					if(matchedJetNumber > -1)
						JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
				}
		
				for(unsigned int i=0; i<JetPartonPair.size(); i++) {
					unsigned int j = JetPartonPair[i].second;
		
					if( fabs(mcParticlesMatching_[j]->type()) < 6 ) {//light/b quarks, 6 should stay hardcoded
						if( ( ( muPlusFromTop || elPlusFromTop ) && mcParticlesMatching_[j]->motherType() == -24 && mcParticlesMatching_[j]->grannyType() == -pdgID_top )
				 || ( ( muMinusFromTop || elMinusFromTop ) && mcParticlesMatching_[j]->motherType() == 24 && mcParticlesMatching_[j]->grannyType() == pdgID_top ) ) {		// if mu+, check if mother of particle is W- and granny tbar --> then it is a quark from W- decay
							if(hadronicWJet1_.first == 9999) {
								hadronicWJet1_ = JetPartonPair[i];
								MCPermutation[0] = JetPartonPair[i].first;
							} else if(hadronicWJet2_.first == 9999) {
								hadronicWJet2_ = JetPartonPair[i];
								MCPermutation[1] = JetPartonPair[i].first;
							} else {
								cerr<<"Found a third jet coming from a W boson which comes from a top quark..."<<endl;
								cerr<<" -- isSemiMu: " << isSemiMu << " isSemiEl: " << isSemiEl << endl;
								cerr<<" -- muMinusFromMtop: " << muMinusFromTop << " muPlusFromMtop: " << muPlusFromTop << endl;
								cerr<<" -- pdgId: " << mcParticlesMatching_[j]->type() << " mother: " << mcParticlesMatching_[j]->motherType() << " granny: " << mcParticlesMatching_[j]->grannyType() << " Pt: " << mcParticlesMatching_[j]->Pt()<< endl;
								exit(1);
							}
						}
					}
					if( fabs(mcParticlesMatching_[j]->type()) == 5 ) {
						if( ( ( muPlusFromTop || elPlusFromTop ) && mcParticlesMatching_[j]->motherType() == -pdgID_top )
						 || ( ( muMinusFromTop || elMinusFromTop ) && mcParticlesMatching_[j]->motherType() == pdgID_top ) ) {		// if mu+ (top decay leptonic) and mother is antitop ---> hadronic b
							hadronicBJet_ = JetPartonPair[i];
							MCPermutation[2] = JetPartonPair[i].first;
						}
						else if( ( ( muPlusFromTop || elPlusFromTop ) && mcParticlesMatching_[j]->motherType() == pdgID_top )
								|| ( ( muMinusFromTop || elMinusFromTop ) && mcParticlesMatching_[j]->motherType() == -pdgID_top ) ) {
							leptonicBJet_ = JetPartonPair[i];
							MCPermutation[3] = JetPartonPair[i].first;
						}
					}
		
/*					// look for ISR stuff
					if( fabs(mcParticlesMatching_[j]->type()) != pdgID_top && fabs(mcParticlesMatching_[j]->motherType()) != 24 && fabs(mcParticlesMatching_[j]->motherType()) != pdgID_top &&
						 fabs(mcParticlesMatching_[j]->grannyType()) != 24 && fabs(mcParticlesMatching_[j]->grannyType()) != pdgID_top )
					{		// not top & mother not W (no standard light jet) & mother not top (no standard b-jet) & granny not W (no radiation from light quarks) & granny not top (no radiation from b)
						ISRJetPartonPair.push_back(JetPartonPair[i]);
					}
*/	
				}
				if(hadronicWJet1_.first != 9999 && hadronicWJet2_.first != 9999 && hadronicBJet_.first != 9999 && leptonicBJet_.first != 9999) {
		
					all4PartonsMatched = true;
					if(hadronicWJet1_.first < 4 && hadronicWJet2_.first < 4 && hadronicBJet_.first < 4 && leptonicBJet_.first < 4)
						all4JetsMatched_MCdef_ = true;
				}
			////cout<<"   ------> according to JetCombiner: hadronicWJet1_.first = "<<hadronicWJet1_.first<<", hadronicWJet2_.first = "<<hadronicWJet2_.first<<", hadronicBJet_.first = "<<hadronicBJet_.first<<endl;
				if(hadronicWJet1_.first < 4 && hadronicWJet2_.first < 4 && hadronicBJet_.first < 4)
					hadronictopJetsMatched_MCdef_ = true;
				if(isSemiMu && ROOT::Math::VectorUtil::DeltaR( (TLorentzVector) *mcParticles[genmuon], (TLorentzVector)*selectedMuons[0]) < 0.3)
					muonmatched =  true;
				if(isSemiEl && ROOT::Math::VectorUtil::DeltaR( (TLorentzVector)*mcParticles[genelectron], (TLorentzVector)*selectedElectrons[0]) < 0.3)
					electronmatched = true;
			}

			/*cout << "isSemiMu: " << isSemiMu;
			 cout << " isSemiEl: " << isSemiEl;
			 cout << " SelSemiMu: "<<eventselectedSemiMu;
			 cout << " SelSemiEl: " << eventselectedSemiEl;*/
			 //cout << " All 4 partons matched: " << all4JetsMatched_MCdef_ << endl;

			//if (eventselectedSemiEl) 
			//	exit(0);
	
			if (all4JetsMatched_MCdef_ && !useMassesAndResolutions && dataSetName.find("TTbarJets") == 0) {
		
				float genWMass = (*mcParticlesMatching_[hadronicWJet1_.second]+*mcParticlesMatching_[hadronicWJet2_.second]).M();
				float genTopMass = (*mcParticlesMatching_[hadronicWJet1_.second]+*mcParticlesMatching_[hadronicWJet2_.second]+*mcParticlesMatching_[hadronicBJet_.second]).M();

				float WMass = (*selectedJets[hadronicWJet1_.first]+*selectedJets[hadronicWJet2_.first]).M();
				float TopMass = (*selectedJets[hadronicWJet1_.first]+*selectedJets[hadronicWJet2_.first]+*selectedJets[hadronicBJet_.first]).M();
	
				histo1D["hadronicGenWMass"]->Fill(genWMass);
				histo1D["hadronicGenTopMass"]->Fill(genTopMass);

				histo1D["hadronicRecoWMass"]->Fill(WMass);
				histo1D["hadronicRecoTopMass"]->Fill(TopMass);
		
				//cout << "WMass " << WMass << "; TopMass " << TopMass << endl;
			
				//if(calculateTransferFunctions){
					
					//energies
					histo2D["Eparton_vs_Enonbjet"]->Fill(mcParticlesMatching_[hadronicWJet1_.second]->E(),selectedJets[hadronicWJet1_.first]->E());
					histo2D["Eparton_vs_Eparton-Enonbjet"]->Fill(mcParticlesMatching_[hadronicWJet1_.second]->E(),mcParticlesMatching_[hadronicWJet1_.second]->E()-selectedJets[hadronicWJet1_.first]->E());
					histo2D["Eparton_vs_Enonbjet"]->Fill(mcParticlesMatching_[hadronicWJet2_.second]->E(),selectedJets[hadronicWJet2_.first]->E());
					histo2D["Eparton_vs_Eparton-Enonbjet"]->Fill(mcParticlesMatching_[hadronicWJet2_.second]->E(),mcParticlesMatching_[hadronicWJet2_.second]->E()-selectedJets[hadronicWJet2_.first]->E());
					histo2D["Eparton_vs_Ebjet"]->Fill(mcParticlesMatching_[hadronicBJet_.second]->E(),selectedJets[hadronicBJet_.first]->E());
					histo2D["Eparton_vs_Eparton-Ebjet"]->Fill(mcParticlesMatching_[hadronicBJet_.second]->E(),mcParticlesMatching_[hadronicBJet_.second]->E()-selectedJets[hadronicBJet_.first]->E());
					histo2D["Eparton_vs_Ebjet"]->Fill(mcParticlesMatching_[leptonicBJet_.second]->E(),selectedJets[leptonicBJet_.first]->E());
					histo2D["Eparton_vs_Eparton-Ebjet"]->Fill(mcParticlesMatching_[leptonicBJet_.second]->E(),mcParticlesMatching_[leptonicBJet_.second]->E()-selectedJets[leptonicBJet_.first]->E());
					if(isSemiEl && electronmatched) {
						histo2D["EgenEl_vs_ErecEl"]->Fill(mcParticles[genelectron]->E(),selectedElectrons[0]->E());
						histo2D["EgenEl_vs_EgenEl-ErecEl"]->Fill(mcParticles[genelectron]->E(),mcParticles[genelectron]->E()-selectedElectrons[0]->E());
					}
					if(isSemiMu && muonmatched) {
						float InvPtgenMu = 1./mcParticles[genmuon]->Pt();
						float InvPtrecMu = 1./selectedMuons[0]->Pt();
						histo2D["InvPtgenMu_vs_InvPtrecMu"]->Fill(InvPtgenMu,InvPtrecMu);
						histo2D["InvPtgenMu_vs_InvPtgenMu-InvPtrecMu"]->Fill(InvPtgenMu,InvPtgenMu-InvPtrecMu);
					}
					//cout << "what's up?" << endl;
					
					//angles
					histo2D["Thparton_vs_Thnonbjet"]->Fill(mcParticlesMatching_[hadronicWJet1_.second]->Theta(),selectedJets[hadronicWJet1_.first]->Theta());
					histo2D["Thparton_vs_Thparton-Thnonbjet"]->Fill(mcParticlesMatching_[hadronicWJet1_.second]->Theta(),mcParticlesMatching_[hadronicWJet1_.second]->Theta()-selectedJets[hadronicWJet1_.first]->Theta());
					histo2D["Thparton_vs_Thnonbjet"]->Fill(mcParticlesMatching_[hadronicWJet2_.second]->Theta(),selectedJets[hadronicWJet2_.first]->Theta());
					histo2D["Thparton_vs_Thparton-Thnonbjet"]->Fill(mcParticlesMatching_[hadronicWJet2_.second]->Theta(),mcParticlesMatching_[hadronicWJet2_.second]->Theta()-selectedJets[hadronicWJet2_.first]->Theta());
					histo2D["Thparton_vs_Thbjet"]->Fill(mcParticlesMatching_[hadronicBJet_.second]->Theta(),selectedJets[hadronicBJet_.first]->Theta());
					histo2D["Thparton_vs_Thparton-Thbjet"]->Fill(mcParticlesMatching_[hadronicBJet_.second]->Theta(),mcParticlesMatching_[hadronicBJet_.second]->Theta()-selectedJets[hadronicBJet_.first]->Theta());
					histo2D["Thparton_vs_Thbjet"]->Fill(mcParticlesMatching_[leptonicBJet_.second]->Theta(),selectedJets[leptonicBJet_.first]->Theta());
					histo2D["Thparton_vs_Thparton-Thbjet"]->Fill(mcParticlesMatching_[leptonicBJet_.second]->Theta(),mcParticlesMatching_[leptonicBJet_.second]->Theta()-selectedJets[leptonicBJet_.first]->Theta());
					histo2D["Eparton_vs_Thnonbjet"]->Fill(mcParticlesMatching_[hadronicWJet1_.second]->E(),selectedJets[hadronicWJet1_.first]->Theta());
					histo2D["Eparton_vs_Thparton-Thnonbjet"]->Fill(mcParticlesMatching_[hadronicWJet1_.second]->E(),mcParticlesMatching_[hadronicWJet1_.second]->Theta()-selectedJets[hadronicWJet1_.first]->Theta());
					histo2D["Eparton_vs_Thnonbjet"]->Fill(mcParticlesMatching_[hadronicWJet2_.second]->E(),selectedJets[hadronicWJet2_.first]->Theta());
					histo2D["Eparton_vs_Thparton-Thnonbjet"]->Fill(mcParticlesMatching_[hadronicWJet2_.second]->E(),mcParticlesMatching_[hadronicWJet2_.second]->Theta()-selectedJets[hadronicWJet2_.first]->Theta());
					histo2D["Eparton_vs_Thbjet"]->Fill(mcParticlesMatching_[hadronicBJet_.second]->E(),selectedJets[hadronicBJet_.first]->Theta());
					histo2D["Eparton_vs_Thparton-Thbjet"]->Fill(mcParticlesMatching_[hadronicBJet_.second]->E(),mcParticlesMatching_[hadronicBJet_.second]->Theta()-selectedJets[hadronicBJet_.first]->Theta());
					histo2D["Eparton_vs_Thbjet"]->Fill(mcParticlesMatching_[leptonicBJet_.second]->E(),selectedJets[leptonicBJet_.first]->Theta());
					histo2D["Eparton_vs_Thparton-Thbjet"]->Fill(mcParticlesMatching_[leptonicBJet_.second]->E(),mcParticlesMatching_[leptonicBJet_.second]->Theta()-selectedJets[leptonicBJet_.first]->Theta());
					if(isSemiEl && electronmatched) {
						histo2D["ThgenEl_vs_ThrecEl"]->Fill(mcParticles[genelectron]->Theta(),selectedElectrons[0]->Theta());
						histo2D["ThgenEl_vs_ThgenEl-ThrecEl"]->Fill(mcParticles[genelectron]->Theta(),mcParticles[genelectron]->Theta()-selectedElectrons[0]->Theta());
						histo2D["EgenEl_vs_ThrecEl"]->Fill(mcParticles[genelectron]->E(),selectedElectrons[0]->Theta());
						histo2D["EgenEl_vs_ThgenEl-ThrecEl"]->Fill(mcParticles[genelectron]->E(),mcParticles[genelectron]->Theta()-selectedElectrons[0]->Theta());
					}
					if(isSemiMu && muonmatched) {
						histo2D["ThgenMu_vs_ThrecMu"]->Fill(mcParticles[genmuon]->Theta(),selectedMuons[0]->Theta());
						histo2D["ThgenMu_vs_ThgenMu-ThrecMu"]->Fill(mcParticles[genmuon]->Theta(),mcParticles[genmuon]->Theta()-selectedMuons[0]->Theta());
						histo2D["InvPtgenMu_vs_ThrecMu"]->Fill(1./mcParticles[genmuon]->Pt(),selectedMuons[0]->Theta());
						histo2D["InvPtgenMu_vs_ThgenMu-ThrecMu"]->Fill(1./mcParticles[genmuon]->Pt(),mcParticles[genmuon]->Theta()-selectedMuons[0]->Theta());
					}
					
					histo2D["Phiparton_vs_Phinonbjet"]->Fill(mcParticlesMatching_[hadronicWJet1_.second]->Phi(),selectedJets[hadronicWJet1_.first]->Phi());
					histo2D["Phiparton_vs_Phinonbjet"]->Fill(mcParticlesMatching_[hadronicWJet2_.second]->Phi(),selectedJets[hadronicWJet2_.first]->Phi());
					histo2D["Eparton_vs_Phinonbjet"]->Fill(mcParticlesMatching_[hadronicWJet1_.second]->E(),selectedJets[hadronicWJet1_.first]->Phi());
					histo2D["Eparton_vs_Phinonbjet"]->Fill(mcParticlesMatching_[hadronicWJet2_.second]->E(),selectedJets[hadronicWJet2_.first]->Phi());
					float DeltaPhi_nonbjet1 = ROOT::Math::VectorUtil::DeltaPhi((TLorentzVector)*mcParticlesMatching_[hadronicWJet1_.second],(TLorentzVector)*selectedJets[hadronicWJet1_.first]);
					histo2D["Phiparton_vs_Phiparton-Phinonbjet"]->Fill(mcParticlesMatching_[hadronicWJet1_.second]->Phi(),DeltaPhi_nonbjet1);
					histo2D["Eparton_vs_Phiparton-Phinonbjet"]->Fill(mcParticlesMatching_[hadronicWJet1_.second]->E(),DeltaPhi_nonbjet1);
					float DeltaPhi_nonbjet2 = ROOT::Math::VectorUtil::DeltaPhi((TLorentzVector)*mcParticlesMatching_[hadronicWJet2_.second],(TLorentzVector)*selectedJets[hadronicWJet2_.first]);
					histo2D["Phiparton_vs_Phiparton-Phinonbjet"]->Fill(mcParticlesMatching_[hadronicWJet2_.second]->Phi(),DeltaPhi_nonbjet2);
					histo2D["Eparton_vs_Phiparton-Phinonbjet"]->Fill(mcParticlesMatching_[hadronicWJet2_.second]->E(),DeltaPhi_nonbjet2);
					
					histo2D["Phiparton_vs_Phibjet"]->Fill(mcParticlesMatching_[hadronicBJet_.second]->Phi(),selectedJets[hadronicBJet_.first]->Phi());
					histo2D["Phiparton_vs_Phibjet"]->Fill(mcParticlesMatching_[leptonicBJet_.second]->Phi(),selectedJets[leptonicBJet_.first]->Phi());
					histo2D["Eparton_vs_Phibjet"]->Fill(mcParticlesMatching_[hadronicBJet_.second]->E(),selectedJets[hadronicBJet_.first]->Phi());
					histo2D["Eparton_vs_Phibjet"]->Fill(mcParticlesMatching_[leptonicBJet_.second]->E(),selectedJets[leptonicBJet_.first]->Phi());
					float DeltaPhi_bjet1 = ROOT::Math::VectorUtil::DeltaPhi((TLorentzVector)*mcParticlesMatching_[hadronicBJet_.second],(TLorentzVector)*selectedJets[leptonicBJet_.first]);				
					histo2D["Phiparton_vs_Phiparton-Phibjet"]->Fill(mcParticlesMatching_[hadronicBJet_.second]->Phi(),DeltaPhi_bjet1);
					histo2D["Eparton_vs_Phiparton-Phibjet"]->Fill(mcParticlesMatching_[hadronicBJet_.second]->E(),DeltaPhi_bjet1);
					float DeltaPhi_bjet2 = ROOT::Math::VectorUtil::DeltaPhi((TLorentzVector)*mcParticlesMatching_[leptonicBJet_.second],(TLorentzVector)*selectedJets[leptonicBJet_.first]);				
					histo2D["Phiparton_vs_Phiparton-Phibjet"]->Fill(mcParticlesMatching_[leptonicBJet_.second]->Phi(),DeltaPhi_bjet2);
					histo2D["Eparton_vs_Phiparton-Phibjet"]->Fill(mcParticlesMatching_[leptonicBJet_.second]->E(),DeltaPhi_bjet2);

					if(isSemiEl && electronmatched) { 
						histo2D["PhigenEl_vs_PhirecEl"]->Fill(mcParticles[genelectron]->Phi(),selectedElectrons[0]->Phi());
						histo2D["EgenEl_vs_PhirecEl"]->Fill(mcParticles[genelectron]->E(),selectedElectrons[0]->Phi());
						float DeltaPhi = ROOT::Math::VectorUtil::DeltaPhi((TLorentzVector)*mcParticles[genelectron],(TLorentzVector)*selectedElectrons[0]);
						histo2D["PhigenEl_vs_PhigenEl-PhirecEl"]->Fill(mcParticles[genelectron]->Phi(),DeltaPhi);
						histo2D["EgenEl_vs_PhigenEl-PhirecEl"]->Fill(mcParticles[genelectron]->E(),DeltaPhi);
					}
					if(isSemiMu && muonmatched) {
						histo2D["PhigenMu_vs_PhirecMu"]->Fill(mcParticles[genmuon]->Phi(),selectedMuons[0]->Phi());
						histo2D["InvPtgenMu_vs_PhirecMu"]->Fill(1./mcParticles[genmuon]->Pt(),selectedMuons[0]->Phi());
						float DeltaPhi = ROOT::Math::VectorUtil::DeltaPhi((TLorentzVector)*mcParticles[genmuon],(TLorentzVector)*selectedMuons[0]);
						histo2D["PhigenMu_vs_PhigenMu-PhirecMu"]->Fill(mcParticles[genmuon]->Phi(),DeltaPhi);
						histo2D["InvPtgenMu_vs_PhigenMu-PhirecMu"]->Fill(1./mcParticles[genmuon]->Pt(),DeltaPhi);
					}
				//}
					//cout << "what's up2?" << endl;
				
			}
			
			
		

			TLorentzVector* selectedLepton;
			if (eventselectedSemiMu_onebtag)
				selectedLepton = (TLorentzVector*)selectedMuons[0];
			else if (eventselectedSemiEl_onebtag)
				selectedLepton = (TLorentzVector*)selectedElectrons[0];
     

			//-----------------//
			// do some data-mc //
			//-----------------//
     
		 	//cout << "make data-mc plots now" << endl;
			// when running both electron and muon data, pick the right dataset vector and lumi for the MSPlots
			if (!foundMu && !foundEl) {
				//cout << "datasetsPlot is independent of muon or electron "  << endl;
				datasetsPlot = datasets;
			} else if (eventselectedSemiMu_onebtag) { 
				//cout << "datasetsPlot should be muonlist " << endl;
				datasetsPlot = datasetsMu;
				Luminosity = LuminosityMu;
			} else if (eventselectedSemiEl_onebtag) {
				datasetsPlot = datasetsEl;
				Luminosity = LuminosityEl;
			}
     
			string Flav="_other";
			if (eventselectedSemiMu_onebtag || eventselectedSemiMu)
				Flav="_mu";
			else if (eventselectedSemiEl_onebtag || eventselectedSemiEl)
				Flav="_el";
     
		 	//declaring plots
			if (MSPlot.find("Selected_Events_pT_jet1"+Flav) == MSPlot.end()){
		 		//cout << "declaring the plots..." << endl;
				MSPlot["Selected_Events_pT_jet1"+Flav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_jet1"+Flav, 30, 0, 600, "p_{T} (GeV)");
				MSPlot["Selected_Events_pT_jet2"+Flav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_jet2"+Flav, 30, 0, 600, "p_{T} (GeV)");
				MSPlot["Selected_Events_pT_jet3"+Flav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_jet3"+Flav, 30, 0, 600, "p_{T} (GeV)");
				MSPlot["Selected_Events_pT_jet4"+Flav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_jet4"+Flav, 30, 0, 600, "p_{T} (GeV)");
				MSPlot["Selected_Events_pT_4leadingjets"+Flav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_4leadingjets"+Flav, 30, 0, 600, "p_{T} (GeV)");
				MSPlot["Selected_Events_pT_alljets"+Flav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_alljets"+Flav, 30, 0, 600, "p_{T} (GeV)");
				MSPlot["NofSelectedJets"+Flav] = new MultiSamplePlot(datasetsPlot, "NofSelectedJets"+Flav, 12, 2, 14, "Number of Jets");
				MSPlot["Selected_Events_Btag_Values"+Flav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_Btag_Values"+Flav, 30, -1, 1.1, "BTag value");
				MSPlot["Selected_Events_nb_Btags"+Flav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_nb_Btags"+Flav, 30, -1, 2, "BTag value");
				MSPlot["NofPV"] = new MultiSamplePlot(datasetsPlot, "NofPV", 31, -0.5, 30.5, "Nb. of primary vertices");
				MSPlot["NofPV_after_lumiWeight"] = new MultiSamplePlot(datasetsPlot, "NofPV_after_lumiWeight", 31, -0.5, 30.5, "Nb. of primary vertices");
				MSPlot["Pileup_Reweighting"+Flav] = new MultiSamplePlot(datasetsPlot,"Pileup_Reweighting"+Flav, 40, 0, 20, "lumiWeight");
				MSPlot["MET_Pt"+Flav] = new MultiSamplePlot(datasetsPlot,"MET_Pt"+Flav, 200, 0, 600, "p_{T} (GeV)");
				MSPlot["MET_Et"+Flav] = new MultiSamplePlot(datasetsPlot,"MET_Et"+Flav, 200, 0, 600, "E_{T} (GeV)");
				MSPlot["Ht_4leadingjets"+Flav] = new MultiSamplePlot(datasetsPlot,"Ht_4leadingjets"+Flav, 100, 0, 1000, "H_{T} (GeV)");
				MSPlot["Selected_Events_pT_lepton"+Flav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_lepton"+Flav, 22, 0, 440, "p_{T} (GeV)");
				MSPlot["Selected_Events_Eta_lepton"+Flav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_Eta_lepton"+Flav, 22, 0, 2.2, "Eta");

				//MSPlot["Chi2_1btag"+Flav] = new MultiSamplePlot(datasetsPlot,"Chi2_1btag"+Flav, 100, 0, 50, "{#chi}^2");
				//MSPlot["hadTop_Ht_1btag"+Flav] = new MultiSamplePlot(datasetsPlot,"hadTopHt_1btag"+Flav, 100, 0, 1000, "H_{T} (GeV)");
				//MSPlot["hadTop_Mass_1btag"+Flav] = new MultiSamplePlot(datasetsPlot,"hadTop_Mass_1btag"+Flav, 300, 0, 600, "Mass (GeV)");
				//MSPlot["hadTopPt_1btag"+Flav] = new MultiSamplePlot(datasetsPlot,"hadTop_Pt_1btag"+Flav, 300, 0, 600, "Mass (GeV)");
				//MSPlot["TTbar_Mass_1btag"+Flav] = new MultiSamplePlot(datasetsPlot,"TTbar_Mass_1btag"+Flav, 300, 0, 600, "Mass (GeV)");
				//MSPlot["TTbar_Angle_1btag"+Flav] = new MultiSamplePlot(datasetsPlot,"TTbar_Angle_1btag"+Flav, 300, 0, 600, "Phi");

				MSPlot["Chi2_2btags"+Flav] = new MultiSamplePlot(datasetsPlot,"Chi2_2btags"+Flav, 100, 0, 50, "{#chi}^2");
				MSPlot["hadTop_Ht_2btags"+Flav] = new MultiSamplePlot(datasetsPlot,"hadTopHt_2btags"+Flav, 100, 0, 1000, "H_{T} (GeV)");
				MSPlot["hadTop_Mass_2btags"+Flav] = new MultiSamplePlot(datasetsPlot,"hadTop_Mass_2btags"+Flav, 100, 0, 600, "Mass (GeV)");
				MSPlot["hadTop_Pt_2btags"+Flav] = new MultiSamplePlot(datasetsPlot,"hadTop_Pt_2btags"+Flav, 100, 0, 600, "Mass (GeV)");
				//MSPlot["TTbar_Mass_2btags"+Flav] = new MultiSamplePlot(datasetsPlot,"TTbar_Mass_2btags"+Flav, 100, 0, 600, "Mass (GeV)");
				//MSPlot["TTbar_Angle_2btags"+Flav] = new MultiSamplePlot(datasetsPlot,"TTbar_Angle_2btags"+Flav, 100, 0, 600, "Phi");
			
				if (eventselectedSemiMu_onebtag) {
					MSPlot["Selected_Events_d0_Muon"] = new MultiSamplePlot(datasetsPlot, "Selected_Events_d0_Muon", 20, 0, 0.1, "d0");
					MSPlot["Selected_Events_relIso_Muon"] = new MultiSamplePlot(datasetsPlot, "Selected_Events_relIso_Muon", 20, 0, 0.2, "relIso");
       	}
				if (eventselectedSemiEl_onebtag) {
					MSPlot["Selected_Events_relIso_Electron"] = new MultiSamplePlot(datasetsPlot, "Selected_Events_relIso_Electron", 20, 0, 0.2, "relIso");
       	}
		 		//cout << "all plots declared" << endl;
			}
     
		 	//filling plots
		 	//cout << "filling the plots..." << endl;
		 
			MSPlot["Selected_Events_pT_jet1"+Flav]->Fill(selectedJets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
			MSPlot["Selected_Events_pT_jet2"+Flav]->Fill(selectedJets[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
			MSPlot["Selected_Events_pT_jet3"+Flav]->Fill(selectedJets[2]->Pt(), datasets[d], true, Luminosity*scaleFactor);
			MSPlot["Selected_Events_pT_jet4"+Flav]->Fill(selectedJets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
			MSPlot["NofSelectedJets"+Flav]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["Selected_Events_nb_Btags"+Flav]->Fill(nb_btags, datasets[d], true, Luminosity*scaleFactor);
			
			for (unsigned int q=0; q<selectedJets.size(); q++) {
				MSPlot["Selected_Events_pT_alljets"+Flav]->Fill(selectedJets[q]->Pt(), datasets[d], true, Luminosity*scaleFactor);
				if (q<4) {
					MSPlot["Selected_Events_pT_4leadingjets"+Flav]->Fill(selectedJets[q]->Pt(), datasets[d], true, Luminosity*scaleFactor);
					MSPlot["Selected_Events_Btag_Values"+Flav]->Fill(selectedJets[q]->btag_combinedSecondaryVertexBJetTags(), datasets[d], true, Luminosity*scaleFactor);
				}
			}
     
			MSPlot["NofPV"]->Fill(vertex.size(),datasets[d], true, Luminosity*1);
			MSPlot["NofPV_after_lumiWeight"]->Fill(vertex.size(),datasets[d], true, Luminosity*scaleFactor);
			MSPlot["Pileup_Reweighting"+Flav]->Fill(lumiWeight, datasets[d], true, Luminosity*scaleFactor);
     
			MSPlot["MET_Pt"+Flav]->Fill(mets[0]->Pt(),datasets[d], true, Luminosity*scaleFactor);
			MSPlot["MET_Et"+Flav]->Fill(mets[0]->Et(),datasets[d], true, Luminosity*scaleFactor);
     
			float HT = selectedJets[0]->Pt()+selectedJets[1]->Pt()+selectedJets[2]->Pt()+selectedJets[3]->Pt();
			MSPlot["Ht_4leadingjets"+Flav]->Fill(HT, datasets[d], true, Luminosity*scaleFactor);

			MSPlot["Selected_Events_pT_lepton"+Flav]->Fill(selectedLepton->Pt(), datasets[d], true, Luminosity*scaleFactor);
			MSPlot["Selected_Events_Eta_lepton"+Flav]->Fill(selectedLepton->Eta(), datasets[d], true, Luminosity*scaleFactor);
			
			if (eventselectedSemiMu_onebtag) {
				MSPlot["Selected_Events_d0_Muon"]->Fill(selectedMuons[0]->d0(), datasets[d], true, Luminosity*scaleFactor);
				float muon_relIso = (selectedMuons[0]->chargedHadronIso() + max( 0.0, selectedMuons[0]->neutralHadronIso() + selectedMuons[0]->photonIso() - 0.5*selectedMuons[0]->puChargedHadronIso() ) ) / selectedMuons[0]->Pt(); // dBeta corrected
				MSPlot["Selected_Events_relIso_Muon"]->Fill(muon_relIso, datasets[d], true, Luminosity*scaleFactor);
			}
     
			if (eventselectedSemiEl_onebtag) {
				float electron_relIso = (selectedElectrons[0]->chargedHadronIso() + max( 0.0, selectedElectrons[0]->neutralHadronIso() + selectedElectrons[0]->photonIso() - 0.5*selectedElectrons[0]->puChargedHadronIso() ) ) / selectedElectrons[0]->Pt(); // dBeta corrected
				MSPlot["Selected_Events_relIso_Electron"]->Fill(electron_relIso, datasets[d], true, Luminosity*scaleFactor);
			}
     
			//--------------------------------------------//
			// find the b-tagged jets with the highest pt //					
			//--------------------------------------------//
		 	//cout << "find the b-tagged jets with the highest pt" << endl;
			int labelBtag1 = -9999;
			int labelBtag2 = -9999;
			float PtBtag1 = -9999.;
			float PtBtag2 = -9999.;
			for (unsigned int i=0; i < 4; i++) {
				//cout << "working with jet[" << i << "] with pt=" << selectedJets[i]->Pt() << " and btag value = " << selectedJets[i]->btag_combinedSecondaryVertexBJetTags() << endl;
				if (selectedJets[i]->btag_combinedSecondaryVertexBJetTags() > 0.679) {		// CSVM
					//cout << "jet[" << i << "] is btagged (CSVM)" << endl;
					//cout << "labelBtag1:" << labelBtag1 << " and PtBtag1:"<<PtBtag1 << endl;
					//cout << "labelBtag2:" << labelBtag2 << " and PtBtag2:"<<PtBtag2 << endl;
					if (selectedJets[i]->Pt() > PtBtag1) {
						// Save previous as second best
						if(labelBtag1 >= 0){
							labelBtag2 = labelBtag1;
							PtBtag2 = selectedJets[labelBtag2]->Pt();
						}
						// Keep new one
						labelBtag1 = i;
						PtBtag1 = selectedJets[i]->Pt();
					}
					else if (selectedJets[i]->Pt() > PtBtag2) {
						labelBtag2 = i;
						PtBtag2 = selectedJets[i]->Pt();
					}
				}
			}
			//if(labelBtag1>=0) cout << "first b-tagged jet is jet[" << labelBtag1 << "] with pt " << selectedJets[labelBtag1]->Pt() << endl;
			//if(labelBtag2>=0) cout << "second b-tagged jet is jet[" << labelBtag2 << "] with pt " << selectedJets[labelBtag2]->Pt() <<endl;

/* we require at least two b-tags in the event selection, hence this part does not make much sense
			/////////////////////
			// CHI2 FOR 1 BTAG //
			/////////////////////

			if (useMassesAndResolutions && labelBtag1 != -9999) {
				float RecoWMass, RecoTopMassB, RecoTopMassOther, WTerm, TopTermB, TopTermOther, chi2B, chi2Other;
				float smallestChi2 = 9999.;
				int labelsReco[4];		// 0 = leptonic b-jet, 1 = hadronic b-jet, 2,3 = light jets.
				for (unsigned int ijet=0; ijet<4; ijet++) {
					for (unsigned int jjet=0; jjet<4; jjet++) {
						for (unsigned int kjet=0; kjet<4; kjet++) {
							if (ijet < jjet && ijet != jjet && ijet != kjet && ijet != labelBtag1 && jjet != kjet && jjet != labelBtag1 && kjet != labelBtag1) {
								RecoWMass = (*selectedJets[ijet] + *selectedJets[jjet]).M();
								RecoTopMassB = (*selectedJets[ijet] + *selectedJets[jjet] + *selectedJets[labelBtag1]).M();
								RecoTopMassOther = (*selectedJets[ijet] + *selectedJets[jjet] + *selectedJets[kjet]).M();
						
								WTerm = pow( (RecoWMass - Chi2Wmass)/SigmaChi2Wmass, 2);
								TopTermB = pow( (RecoTopMassB - Chi2Topmass)/SigmaChi2Topmass, 2);
								TopTermOther = pow( (RecoTopMassOther - Chi2Topmass)/SigmaChi2Topmass, 2);

								chi2B = WTerm + TopTermB;
								chi2Other = WTerm + TopTermOther;
						
								if (chi2B < smallestChi2) {
									smallestChi2 = chi2B;
									labelsReco[0] = kjet;
									labelsReco[1] = labelBtag1;
									labelsReco[2] = ijet;
									labelsReco[3] = jjet;
								}
								if (chi2Other < smallestChi2) {
									smallestChi2 = chi2Other;
									labelsReco[0] = labelBtag1;
									labelsReco[1] = kjet;
									labelsReco[2] = ijet;
									labelsReco[3] = jjet;
								}
							}
						}
					}
				}
		
				// Fill histos
				MSPlot["Chi2_1btag"+Flav]->Fill(smallestChi2, datasets[d], true, Luminosity*scaleFactor);
				
				float HtTop_1btag = selectedJets[labelsReco[1]]->Pt() + selectedJets[labelsReco[2]]->Pt() + selectedJets[labelsReco[3]]->Pt();
				MSPlot["hadTop_Ht_1btag"+Flav]->Fill(HtTop_1btag, datasets[d], true, Luminosity*scaleFactor);
				float hadtopmass = ( *selectedJets[labelsReco[1]] + *selectedJets[labelsReco[2]] +  *selectedJets[labelsReco[3]]).M();
				float hadtoppt = ( *selectedJets[labelsReco[1]] + *selectedJets[labelsReco[2]] +  *selectedJets[labelsReco[3]]).Pt();
				MSPlot["hadTop_Mass_1btag"+Flav]->Fill(hadtopmass, datasets[d], true, Luminosity*scaleFactor);
				MSPlot["hadTop_Pt_1btag"+Flav]->Fill(hadtoppt, datasets[d], true, Luminosity*scaleFactor);
		
				// the following is not correct -> z-component of MET needs to be calculated

				//float leptonicTopMass_1btag = (selectedMuons[0] + mets[0] + selectedJets[labelsReco[0]]).M();
				//float TTbarMass_1btag = leptonicTopMass_1btag + RecoTopMassForHisto_1btag;
				//MSPlot["TTbar_Mass_1btag"+Flav]->Fill(TTbarMass_1btag, datasets[d], true, Luminosity*scaleFactor);
		
				//float leptonicTopAngle_1btag = (selectedMuons[0] + mets[0] + selectedJets[labelsReco[0]]).Phi();
				//float hadronicTopAngle_1btag = (selectedJets[labelsReco[1]] + selectedJets[labelsReco[2]] + selectedJets[labelsReco[3]]).Phi();
				//float TTbarAngle_1btag = leptonicTopAngle_1btag - hadronicTopAngle_1btag;
				//if (TTbarAngle_1btag < 0) {
				//	TTbarAngle_1btag = - TTbarAngle_1btag;
				//}
				//MSPlot["TTbar_Mass_1btag"+Flav]->Fill(TTbarAngle_1btag, datasets[d], true, Luminosity*scaleFactor);
	
			
			}
*/	
			/////////////////////
			// CHI2 FOR 2 BTAGS //
			//////////////////////
		 	//cout << "find the best jet combination with the chi2 method" << endl;
			if (useMassesAndResolutions && labelBtag1 != -9999 && labelBtag2 != -9999) {
				float RecoWMass, RecoTopMassB1, RecoTopMassB2, WTerm, TopTermB1, TopTermB2, chi2B1, chi2B2;
				float smallestChi2 = 9999.;
				int labelsReco[4] = {0,0,0,0};		// 0 = leptonic b-jet, 1 = hadronic b-jet, 2,3 = light jets.
				for (int ijet=0; ijet<4; ijet++) {
					for (int jjet=0; jjet<4; jjet++) {
						if (ijet < jjet && ijet != jjet && ijet != labelBtag1 && ijet != labelBtag2 && jjet != labelBtag1 && jjet != labelBtag2) {
							RecoWMass = ( *selectedJets[ijet] + *selectedJets[jjet]).M();
							RecoTopMassB1 = (*selectedJets[ijet] + *selectedJets[jjet] + *selectedJets[labelBtag1]).M();
							RecoTopMassB2 = (*selectedJets[ijet] + *selectedJets[jjet] + *selectedJets[labelBtag2]).M();
					
							WTerm = pow( (RecoWMass - Chi2Wmass)/SigmaChi2Wmass, 2);
							TopTermB1 = pow( (RecoTopMassB1 - Chi2Topmass)/SigmaChi2Topmass, 2);
							TopTermB2 = pow( (RecoTopMassB2 - Chi2Topmass)/SigmaChi2Topmass, 2);
					
							chi2B1 = WTerm + TopTermB1;
							chi2B2 = WTerm + TopTermB2;
					
							if (chi2B1 < smallestChi2) {
								smallestChi2 = chi2B1;
								labelsReco[0] = labelBtag2;
								labelsReco[1] = labelBtag1;
								labelsReco[2] = ijet;
								labelsReco[3] = jjet;
							}
							if (chi2B2 < smallestChi2) {
								smallestChi2 = chi2B2;
								labelsReco[0] = labelBtag1;
								labelsReco[1] = labelBtag2;
								labelsReco[2] = ijet;
								labelsReco[3] = jjet;
							}
						}
					}
				}
		
				// Fill histos
				MSPlot["Chi2_2btags"+Flav]->Fill(smallestChi2, datasets[d], true, Luminosity*scaleFactor);
		
				float HtTop_2btags = selectedJets[labelsReco[1]]->Pt() + selectedJets[labelsReco[2]]->Pt() + selectedJets[labelsReco[3]]->Pt();
				MSPlot["hadTop_Ht_2btags"+Flav]->Fill(HtTop_2btags, datasets[d], true, Luminosity*scaleFactor);
				float hadtopmass = ( *selectedJets[labelsReco[1]] + *selectedJets[labelsReco[2]] +  *selectedJets[labelsReco[3]]).M();
				float hadtoppt = ( *selectedJets[labelsReco[1]] + *selectedJets[labelsReco[2]] +  *selectedJets[labelsReco[3]]).Pt();
				MSPlot["hadTop_Mass_2btags"+Flav]->Fill(hadtopmass, datasets[d], true, Luminosity*scaleFactor);
				MSPlot["hadTop_Pt_2btags"+Flav]->Fill(hadtoppt, datasets[d], true, Luminosity*scaleFactor);

				// the following is not correct -> z-component of MET needs to be calculated
		
				//float leptonicTopMass_2btags = (selectedMuons[0] + mets[0] + selectedJets[labelsReco[0]]).M();
				//float TTbarMass_2btags = leptonicTopMass_2btags + RecoTopMassForHisto_2btags;
				//MSPlot["TTbar_Mass_2btags"+Flav]->Fill(TTbarMass_2btags, datasets[d], true, Luminosity*scaleFactor);
		
				//float leptonicTopAngle_2btags = (selectedMuons[0] + mets[0] + selectedJets[labelsReco[0]]).Phi();
				//float hadronicTopAngle_2btags = (selectedJets[labelsReco[1]] + selectedJets[labelsReco[2]] + selectedJets[labelsReco[3]]).Phi();
				//float TTbarAngle_2btags = leptonicTopAngle_2btags - hadronicTopAngle_2btags;
				//if (TTbarAngle_2btags < 0) {
				//	TTbarAngle_2btags = - TTbarAngle_2btags;
				//}
				//MSPlot["TTbar_Mass_1btag"+Flav]->Fill(TTbarAngle_1btag, datasets[d], true, Luminosity*scaleFactor);
			}




     //////////////////
     // END OF EVENT //
     //////////////////
     
    }			//loop on events

		if (!useMassesAndResolutions && dataSetName.find("TTbarJets") == 0) {
			
      string filename = "MassPlotsForChi2JetCombiner.root";

      cout << "INFO: Creating output file to store the masses for the Chi2 jetcombiner: " << filename << endl;    

      TFile* outfile = new TFile(filename.c_str(),"RECREATE");
      
      outfile->cd();
      
      histo1D["hadronicGenWMass"]->Write();
      histo1D["hadronicGenTopMass"]->Write();

      histo1D["hadronicRecoWMass"]->Write();
      histo1D["hadronicRecoTopMass"]->Write();
      
      // fit the distributions
			for (unsigned int f=0; f<2;f++) {

				TH1F* histo;

				if (f==0) histo=histo1D["hadronicRecoWMass"]; 
				if (f==1) histo=histo1D["hadronicRecoTopMass"];

				TF1 *fitfunc;
				string func_title = string(histo->GetName())+"_Fitted";
	
				double rms = histo->GetRMS();
				double maxbin =  histo->GetBinCenter(histo->GetMaximumBin());

				fitfunc = new TF1(func_title.c_str(),"gaus");
				fitfunc->SetRange(maxbin-rms,maxbin+rms);
				histo->Fit(fitfunc,"RQ");

				fitfunc->Write();

				delete fitfunc;
      
      }
            
      outfile->Close(); 			
			
			if(calculateTransferFunctions) {
				
      	string filename = "PlotsForTransferFunctions.root";

      	cout << "INFO: Creating output file to store the plots for the transfer functions " << filename << endl;    

      	TFile* outfile = new TFile(filename.c_str(),"RECREATE");
      
      	outfile->cd();
      
      	histo2D["Eparton_vs_Enonbjet"]->Write();
      	histo2D["Eparton_vs_Ebjet"]->Write();
      	histo2D["EgenEl_vs_ErecEl"]->Write();
      	histo2D["InvPtgenMu_vs_InvPtrecMu"]->Write();
      	histo2D["Thparton_vs_Thnonbjet"]->Write();
      	histo2D["Thparton_vs_Thbjet"]->Write();
      	histo2D["ThgenEl_vs_ThrecEl"]->Write();
      	histo2D["ThgenMu_vs_ThrecMu"]->Write();
      	histo2D["Phiparton_vs_Phinonbjet"]->Write();
      	histo2D["Phiparton_vs_Phibjet"]->Write();
      	histo2D["PhigenEl_vs_PhirecEl"]->Write();
      	histo2D["PhigenMu_vs_PhirecMu"]->Write();
      	histo2D["Eparton_vs_Thnonbjet"]->Write();
      	histo2D["Eparton_vs_Thbjet"]->Write();
      	histo2D["EgenEl_vs_ThrecEl"]->Write();
      	histo2D["InvPtgenMu_vs_ThrecMu"]->Write();
      	histo2D["Eparton_vs_Phinonbjet"]->Write();
      	histo2D["Eparton_vs_Phibjet"]->Write();
      	histo2D["EgenEl_vs_PhirecEl"]->Write();
      	histo2D["InvPtgenMu_vs_PhirecMu"]->Write();
      	
				histo2D["Eparton_vs_Eparton-Enonbjet"]->Write();
      	histo2D["Eparton_vs_Eparton-Ebjet"]->Write();
      	histo2D["EgenEl_vs_EgenEl-ErecEl"]->Write();
      	histo2D["InvPtgenMu_vs_InvPtgenMu-InvPtrecMu"]->Write();
      	histo2D["Eparton_vs_Thparton-Thnonbjet"]->Write();
      	histo2D["Eparton_vs_Thparton-Thbjet"]->Write();
      	histo2D["EgenEl_vs_ThgenEl-ThrecEl"]->Write();
      	histo2D["InvPtgenMu_vs_ThgenMu-ThrecMu"]->Write();
      	histo2D["Eparton_vs_Phiparton-Phinonbjet"]->Write();
      	histo2D["Eparton_vs_Phiparton-Phibjet"]->Write();
      	histo2D["EgenEl_vs_PhigenEl-PhirecEl"]->Write();
      	histo2D["InvPtgenMu_vs_PhigenMu-PhirecMu"]->Write();
				
      	// fit the distributions
				for (unsigned int f=0; f<12;f++) {
					
					if(f==2 || f==6 || f==10) continue; //electron plots not filled at the moment...
					
					TH2F* histo;

					if (f==0) histo=histo2D["Eparton_vs_Eparton-Enonbjet"]; 
					if (f==1) histo=histo2D["Eparton_vs_Eparton-Ebjet"]; 
					if (f==2) histo=histo2D["EgenEl_vs_EgenEl-ErecEl"]; 
					if (f==3) histo=histo2D["InvPtgenMu_vs_InvPtgenMu-InvPtrecMu"]; 
					if (f==4) histo=histo2D["Eparton_vs_Thparton-Thnonbjet"]; 
					if (f==5) histo=histo2D["Eparton_vs_Thparton-Thbjet"]; 
					if (f==6) histo=histo2D["EgenEl_vs_ThgenEl-ThrecEl"]; 
					if (f==7) histo=histo2D["InvPtgenMu_vs_ThgenMu-ThrecMu"]; 
					if (f==8) histo=histo2D["Eparton_vs_Phiparton-Phinonbjet"]; 
					if (f==9) histo=histo2D["Eparton_vs_Phiparton-Phibjet"]; 
					if (f==10) histo=histo2D["EgenEl_vs_PhigenEl-PhirecEl"]; 
					if (f==11) histo=histo2D["InvPtgenMu_vs_PhigenMu-PhirecMu"]; 

	 				if(f==0 || f==1 || f==3 || f==4 || f==5|| f==7 || f==8 || f==9 || f==11) {
						int nbins = histo->GetXaxis()->GetNbins();
						cout << "nbins: " << nbins << endl;
						int npar = 6;

						//Create one histogram for each function parameter -> 5 histograms for each 2D plot
						TH1D **hlist = new TH1D*[npar];
						string parnames[6]={"a1","a2","a3","a4","a5","a6"};
						string name="";string title="";
						const TArrayD *bins = histo->GetXaxis()->GetXbins();
						for (int ipar=0;ipar<npar;ipar++) {
							name = string(histo->GetName())+ "_" + parnames[ipar]; 
							title = string(histo->GetName())+ ": Fitted value of " + parnames[ipar] ;
       				hlist[ipar] = new TH1D(name.c_str(),title.c_str(), nbins,histo->GetXaxis()->GetXmin(),histo->GetXaxis()->GetXmax());
      				hlist[ipar]->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
   					}
						//name = histo->GetName()+ "_" + myfit->GetParName(ipar)+"_chi2"; 
   					//TH1D *hchi2 = 0;
     				//hchi2 = new TH1D(name,"chisquare", nbins, bins->fArray);
						//hchi2->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
			
						//Loop on all bins in X, generate a projection along Y
						int cut = 0; // require a minimum number of bins in the slice to be is filled
   					for (int bin=1;bin < nbins+1;bin ++) {
      				string projection_title = string(histo->GetName())+"_sliceXbin"+tostr(bin);
							TH1D *hp = histo->ProjectionY(projection_title.c_str(),bin,bin,"e");
							if(bin==nbins) hp = histo->ProjectionY(projection_title.c_str(),bin,bin+1,"e"); //include overflow in last bin
							if(bin==1) hp = histo->ProjectionY(projection_title.c_str(),bin-1,bin,"e"); //include underflow in first bin
							if (hp == 0) continue;
							float nentries = float(hp->GetEntries());
      				if (nentries == 0 || nentries < cut) {delete hp; continue;}
							//declare the fit function! its range depends on the jet/lepton energy range (hence, the Y-axis)
							TF1 *myfit = new TF1("myfit","[2]*(TMath::Exp(-TMath::Power((x-[0]),2)/(2*TMath::Power([1],2)))+[5]*TMath::Exp(-TMath::Power((x-[3]),2)/(2*TMath::Power([4],2))))");	
							//give names to the parameters
							myfit->SetParName(0,"a1");
							myfit->SetParName(1,"a2");
							myfit->SetParName(2,"a3");
							myfit->SetParName(3,"a4");
							myfit->SetParName(4,"a5");
							myfit->SetParName(5,"a6");
							//set initial values
							if(f==0 || f==1) { // for E transfer function of JETS
								myfit->SetParameter(0, -8); //central value of first, broad gaussian
								myfit->SetParameter(1, 18); //sigma value of first, broad gaussian
								myfit->SetParameter(2, 63); //constant value of first, broad gaussian
								myfit->SetParameter(3, 0);; //central value of second, narrow gaussian
								myfit->SetParameter(4, 8.6); //sigma value of second, narrow gaussian
								myfit->SetParameter(5, 4.1); //constant value of second, narrow gaussian
							}else if (f==4 || f==5 || f==8 || f==9) { //for theta and phi transfer functions of JETS
								myfit->SetParameter(0, 0); 		 //central value of first, broad gaussian
								myfit->SetParameter(1, 0.038); //sigma value of first, broad gaussian
								myfit->SetParameter(2, 77); 	 //constant value of first, broad gaussian
								myfit->SetParameter(3, 0.004); //central value of second, narrow gaussian
								myfit->SetParameter(4, 0.011); //sigma value of second, narrow gaussian
								myfit->SetParameter(5, 6.5);   //constant value of second, narrow gaussian
							}else if (f==3) { //for 1/pt transfer function of muons
								myfit->SetParameter(0, -0.0008); 		 //central value of first, broad gaussian
								myfit->SetParameter(1, 0.001); //sigma value of first, broad gaussian
								myfit->SetParameter(2, 24); 	 //constant value of first, broad gaussian
								myfit->SetParameter(3, -0.0001); //central value of second, narrow gaussian
								myfit->SetParameter(4, 0.0001); //sigma value of second, narrow gaussian
								myfit->SetParameter(5, 4);   //constant value of second, narrow gaussian
							}else if (f==7 || f==11) { //for theta, phi transfer function of muons
								myfit->SetParameter(0, 0); 		 //central value of first, broad gaussian
								myfit->SetParameter(1, 0.01); //sigma value of first, broad gaussian
								myfit->SetParameter(2, 24); 	 //constant value of first, broad gaussian
								myfit->SetParameter(3, 0); //central value of second, narrow gaussian
								myfit->SetParameter(4, 0.001); //sigma value of second, narrow gaussian
								myfit->SetParameter(5, 4);   //constant value of second, narrow gaussian
							}
							string func_title = string(histo->GetName())+"_sliceXbin"+tostr(bin)+"_Fitted";
							myfit->SetName(func_title.c_str()); 
      				hp->Fit(myfit);
      				int npfits = myfit->GetNumberFitPoints();
      				if (npfits > npar && npfits >= cut) {
         				int binOn = bin + 1/2;
         				for (int ipar=0;ipar<npar;ipar++) {
            			//cout << "histo->GetXaxis()->GetBinCenter(binOn): " << histo->GetXaxis()->GetBinCenter(binOn) << endl;
            			//cout << "myfit->GetParameter("<<ipar<<") " << myfit->GetParameter(ipar) << endl;
									hlist[ipar]->Fill(histo->GetXaxis()->GetBinCenter(binOn),myfit->GetParameter(ipar)); // fill histogram for parameter i
            			hlist[ipar]->SetBinError(histo->GetXaxis()->GetBinCenter(binOn),myfit->GetParError(ipar));
								}
         				//hchi2->Fill(histo->GetXaxis()->GetBinCenter(binOn),myfit->GetChisquare()/(npfits-npar));
      				}
      				hp->Write();
							myfit->Write();
							delete hp;
 							delete myfit;
  					}
						
         		
   					// define the fitfunction for all parameters (6):
						// ai = ai0 + ai1*Ep + ai2*sqrt(Ep) its range depends on the parton energy range (hence, the X-axis)
						TF1 *myfit2 = new TF1("myfit2", "[0]+[1]*sqrt(x)+[2]*x",histo->GetXaxis()->GetXmin(),histo->GetXaxis()->GetXmax() );	
						//give names to the parameters
						myfit2->SetParName(0,"ai0");
						myfit2->SetParName(1,"ai1");
						myfit2->SetParName(2,"ai2");
						
						for (int ipar=0;ipar<npar;ipar++){
							int paramname = ipar+1;
							string func_title2 = string(histo->GetName())+"_a"+tostr(paramname)+"_Fitted";
							myfit2->SetName(func_title2.c_str()); 
							hlist[ipar]->Fit(myfit2);
							hlist[ipar]->Write();
							myfit2->Write();
						}
						
   					delete [] hlist;
						delete myfit2;
      		}	
				}
				
				outfile->Close();
				
				cout << endl;
				cout << endl;
				cout << "--- INFO: THE TRANSFER FUNCTIONS ARE PUT IN A TABLE  ---" << endl;
		
	
				string histonames[9] = {"Eparton_vs_Eparton-Ebjet","Eparton_vs_Thparton-Thbjet","Eparton_vs_Phiparton-Phibjet","Eparton_vs_Eparton-Enonbjet","Eparton_vs_Thparton-Thnonbjet","Eparton_vs_Phiparton-Phinonbjet","InvPtgenMu_vs_InvPtgenMu-InvPtrecMu","InvPtgenMu_vs_ThgenMu-ThrecMu","InvPtgenMu_vs_PhigenMu-PhirecMu"};
			
				string histodescription[9] = {"b-jet energy", "b-jet theta", "b-jet phi", "non-b jet energy", "non-b jet theta", "non-b jet phi", "muon inv. pt", "muon theta", "muon phi",};

				ofstream myTransferFunctions;
				string myTransferFunctions_TABLE = "TransferFunctions_TABLE.txt";
				myTransferFunctions.open(myTransferFunctions_TABLE.c_str());
				ofstream myTransferFunctionsForMadWeight;
				string myTransferFunctionsForMadWeight_dat = "transfer_card_user.dat";
				myTransferFunctionsForMadWeight.open(myTransferFunctionsForMadWeight_dat.c_str());
				for(int i = 0; i<9; i++){
					TF1 *TF_par1,*TF_par2,*TF_par3,*TF_par4,*TF_par5,*TF_par6;
					string name1 = histonames[i]+"_a1_Fitted";
					string name2 = histonames[i]+"_a2_Fitted";
					string name3 = histonames[i]+"_a3_Fitted";
					string name4 = histonames[i]+"_a4_Fitted";
					string name5 = histonames[i]+"_a5_Fitted";
					string name6 = histonames[i]+"_a6_Fitted";
			
					string tf_filename = "PlotsForTransferFunctions.root";
					TFile* tf = new TFile(tf_filename.c_str(),"READ");
					tf->cd();
					TF_par1 = (TF1*)tf->Get(name1.c_str());
					TF_par2 = (TF1*)tf->Get(name2.c_str());
					TF_par3 = (TF1*)tf->Get(name3.c_str());
					TF_par4 = (TF1*)tf->Get(name4.c_str());
					TF_par5 = (TF1*)tf->Get(name5.c_str());
					TF_par6 = (TF1*)tf->Get(name6.c_str());
					if (TF_par1 && TF_par2 && TF_par3 && TF_par4 && TF_par5 && TF_par6) {
						myTransferFunctions<< endl;
						myTransferFunctions<<"\\begin{table}" << endl;
  					myTransferFunctions<<"\\caption{Parameters of the transfer function for the " << histodescription[i]  << "}" << endl;
  					myTransferFunctions<<"\\label{tab:}" << endl;
  					myTransferFunctions<<"\\centering" << endl;
  					myTransferFunctions<<"\\begin{tabular}{c|ccc}" << endl;
  					myTransferFunctions<<"\\hline" << endl;
						myTransferFunctions << "Type	& $a_{i0}$ & $a_{i1}$ ($\\sqrt{E}$) & $a_{i2}$ ($E$)" << "\\\\" << endl;
  					myTransferFunctions<<"\\hline" << endl;
						myTransferFunctions << "Mean broad gaussian & $a_{10}$ = " << TF_par1->GetParameter(0) << "$\\pm$" << TF_par1->GetParError(0) << " & $a_{11}$ = " << TF_par1->GetParameter(1) << "$\\pm$" << TF_par1->GetParError(1) << " & $a_{12}$ = " << TF_par1->GetParameter(2) << "$\\pm$" << TF_par1->GetParError(2) << "\\\\" << endl;
						myTransferFunctions << "Width broad gaussian & $a_{20}$ = " << TF_par2->GetParameter(0) << "$\\pm$" << TF_par2->GetParError(0) << " & $a_{21}$ = " << TF_par2->GetParameter(1) << "$\\pm$" << TF_par2->GetParError(1) << " & $a_{22}$ = " << TF_par2->GetParameter(2) << "$\\pm$" << TF_par2->GetParError(2) << "\\\\" << endl;
						myTransferFunctions << "Constant broad gaussian & $a_{30}$ = " << TF_par2->GetParameter(0) << "$\\pm$" << TF_par2->GetParError(0) << " & $a_{31}$ = " << TF_par2->GetParameter(1) << "$\\pm$" << TF_par2->GetParError(1) << " & $a_{32}$ = " << TF_par2->GetParameter(2) << "$\\pm$" << TF_par2->GetParError(2) << "\\\\" << endl;
						myTransferFunctions << "Mean narrow gaussian & $a_{40}$ = " << TF_par4->GetParameter(0) << "$\\pm$" << TF_par4->GetParError(0) << " & $a_{41}$ = " << TF_par4->GetParameter(1) << "$\\pm$" << TF_par4->GetParError(1) << " & $a_{42}$ = " << TF_par4->GetParameter(2) << "$\\pm$" << TF_par4->GetParError(2) << "\\\\" << endl;
						myTransferFunctions << "Width narrow gaussian & $a_{50}$ = " << TF_par5->GetParameter(0) << "$\\pm$" << TF_par5->GetParError(0) << " & $a_{51}$ = " << TF_par5->GetParameter(1) << "$\\pm$" << TF_par5->GetParError(1) << " & $a_{52}$ = " << TF_par5->GetParameter(2) << "$\\pm$" << TF_par5->GetParError(2) << "\\\\" << endl;
						myTransferFunctions << "Constant narrow gaussian & $a_{60}$ = " << TF_par6->GetParameter(0) << "$\\pm$" << TF_par6->GetParError(0) << " & $a_{61}$ = " << TF_par6->GetParameter(1) << "$\\pm$" << TF_par6->GetParError(1) << " & $a_{62}$ = " << TF_par6->GetParameter(2) << "$\\pm$" << TF_par6->GetParError(2) << "\\\\" << endl;
  					myTransferFunctions<<"\\hline" << endl;
  					myTransferFunctions<<"\\end{tabular}"<<endl;
  					myTransferFunctions<<"\\end{table}"<<endl;
						myTransferFunctions<< endl;

						if(i==0) myTransferFunctionsForMadWeight<<"#+-----------------------------------------------------------------------+" << endl;
						if(i==0) myTransferFunctionsForMadWeight<<"#|    Parameter for particles: b                                         |" << endl;
						if(i==0) myTransferFunctionsForMadWeight<<"#+-----------------------------------------------------------------------+" << endl;
						if(i==3) myTransferFunctionsForMadWeight<<"#+-----------------------------------------------------------------------+" << endl;
						if(i==3) myTransferFunctionsForMadWeight<<"#|    Parameter for particles: nonb                                         |" << endl;
						if(i==3) myTransferFunctionsForMadWeight<<"#+-----------------------------------------------------------------------+" << endl;
						if(i==6) myTransferFunctionsForMadWeight<<"#+-----------------------------------------------------------------------+" << endl;
						if(i==6) myTransferFunctionsForMadWeight<<"#|    Parameter for particles: muon                                         |" << endl;
						if(i==6) myTransferFunctionsForMadWeight<<"#+-----------------------------------------------------------------------+" << endl;
						if(i==0) myTransferFunctionsForMadWeight<<"BLOCK TF_bjet_E" << endl;
						if(i==1) myTransferFunctionsForMadWeight<<"BLOCK TF_bjet_THETA" << endl;
						if(i==2) myTransferFunctionsForMadWeight<<"BLOCK TF_bjet_PHI" << endl;
						if(i==3) myTransferFunctionsForMadWeight<<"BLOCK TF_nonbjet_E" << endl;
						if(i==4) myTransferFunctionsForMadWeight<<"BLOCK TF_nonbjet_THETA" << endl;
						if(i==5) myTransferFunctionsForMadWeight<<"BLOCK TF_nonbjet_PHI" << endl;
						if(i==6) myTransferFunctionsForMadWeight<<"BLOCK TF_muon_InvPt" << endl;
						if(i==7) myTransferFunctionsForMadWeight<<"BLOCK TF_muon_THETA" << endl;
						if(i==8) myTransferFunctionsForMadWeight<<"BLOCK TF_muon_PHI" << endl;
						if(i<6) { //for jets, should also work for electrons if we have them
							for(int j = 0; j<3; j++)
								myTransferFunctionsForMadWeight<<j+1 << " " << TF_par1->GetParameter(j) << "  # bias broad gaussian b1=#1+#2*sqrt(E)*#3*E" << endl;
							for(int j = 3; j<6; j++)
								myTransferFunctionsForMadWeight<<j+1 << " " << TF_par2->GetParameter(j-3) << "  # sigma broad gaussian s1=#4+#5*sqrt(E)*#6*E" << endl;
							for(int j = 6; j<9; j++)
								myTransferFunctionsForMadWeight<<j+1 << " " << TF_par3->GetParameter(j-6) << "  # constant broad gaussian c1=#7+#8*sqrt(E)*#9*E" << endl;
							for(int j = 9; j<12; j++)
								myTransferFunctionsForMadWeight<<j+1 << " " << TF_par4->GetParameter(j-9) << "  # bias narrow gaussian b2=#10+#11*sqrt(E)*#12*E" << endl;
							for(int j = 12; j<15; j++)
								myTransferFunctionsForMadWeight<<j+1 << " " << TF_par5->GetParameter(j-12) << "  # sigma narrow gaussian s2=#13+#14*sqrt(E)*#15*E" << endl;
							for(int j = 15; j<18; j++)
								myTransferFunctionsForMadWeight<<j+1 << " " << TF_par6->GetParameter(j-15) << "  # constant narrow gaussian c2=#16+#17*sqrt(E)*#18*E" << endl;
						}else{ //for muon (only difference: paramtrization according to InvPt instead of E -> something textual...)
							for(int j = 0; j<3; j++)
								myTransferFunctionsForMadWeight<<j+1 << " " << TF_par1->GetParameter(j) << "  # bias broad gaussian b1=#1+#2*sqrt(InvPt)*#3*InvPt" << endl;
							for(int j = 3; j<6; j++)
								myTransferFunctionsForMadWeight<<j+1 << " " << TF_par2->GetParameter(j-3) << "  # sigma broad gaussian s1=#4+#5*sqrt(InvPt)*#6*InvPt" << endl;
							for(int j = 6; j<9; j++)
								myTransferFunctionsForMadWeight<<j+1 << " " << TF_par3->GetParameter(j-6) << "  # constant broad gaussian c1=#7+#8*sqrt(InvPt)*#9*InvPt" << endl;
							for(int j = 9; j<12; j++)
								myTransferFunctionsForMadWeight<<j+1 << " " << TF_par4->GetParameter(j-9) << "  # bias narrow gaussian b2=#10+#11*sqrt(InvPt)*#12*InvPt" << endl;
							for(int j = 12; j<15; j++)
								myTransferFunctionsForMadWeight<<j+1 << " " << TF_par5->GetParameter(j-12) << "  # sigma narrow gaussian s2=#13+#14*sqrt(InvPt)*#15*InvPt" << endl;
							for(int j = 15; j<18; j++)
								myTransferFunctionsForMadWeight<<j+1 << " " << TF_par6->GetParameter(j-15) << "  # constant narrow gaussian c2=#16+#17*sqrt(InvPt)*#18*InvPt" << endl;									}
					}
					tf->Close();
				}
				myTransferFunctions.close();
				myTransferFunctionsForMadWeight.close();
			
			}	
			
	
			useMassesAndResolutions = true;
			d--; // loop again over the first dataset (i.e. ttbar semilep) using the masses and resolutions needed for the chi2 combination
		}
		cout<<endl;


    
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
    
    string pathPNG = "PlotsWtbCouplings/";    
		mkdir(pathPNG.c_str(),0777);
		mkdir((pathPNG+"MSPlot/").c_str(),0777);
    
		//Write histograms
		fout->cd();
		for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
		{
				//cout << "MSPlot: " << it->first << endl;
				MultiSamplePlot *temp = it->second;
				string name = it->first;
				temp->Draw(false, name, true, true, true, true, true, 1, false);
				temp->Write(fout, name, true, pathPNG+"MSPlot/");
		}
    
		// 1D 
		TDirectory* th1dir = fout->mkdir("1D_histograms");
		th1dir->cd();
		for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
		{
			TH1F *temp = it->second;
			int N = temp->GetNbinsX();
			temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
			temp->SetBinContent(N+1,0);
			temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
			temp->Write();
			TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
			tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
		}
    
		// 2D
		TDirectory* th2dir = fout->mkdir("2D_histograms_graphs");
		th2dir->cd();
		for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
		{
			TH2F *temp = it->second;
			temp->Write();
			TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
			tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
		}

  
		//Selection tables
		selecTableSemiMu.TableCalculator(false, true, true, true, true);
    string selectiontableMu = "SelectionTable_WtbCouplings_SemiMu.tex";
    selecTableSemiMu.Write(selectiontableMu.c_str());
        
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

