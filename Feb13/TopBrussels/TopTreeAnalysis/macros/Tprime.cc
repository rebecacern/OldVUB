///////////////////////////
///// TODO & COMMENTS /////
/////////////////////////// 

#include "TStyle.h"
#include "TF2.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

// RooFit librairies
#include "RooArgSet.h"
#include "RooAddition.h"
#include "RooCategory.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooHistPdf.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooGenericPdf.h"
#include "RooMCStudy.h"
#include "RooMinuit.h"
#include "RooPlot.h"
#include "RooFitResult.h"

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "../Selection/interface/SelectionTable.h"
#include "../Tools/interface/PlottingTools.h"
#include "../Tools/interface/MultiSamplePlot.h"
#include "../Tools/interface/TTreeLoader.h"
#include "../Tools/interface/AnalysisEnvironmentLoader.h"
#include "../Content/interface/AnalysisEnvironment.h"
#include "../Content/interface/Dataset.h"
#include "../MCInformation/interface/MCWeighter.h"
#include "../Selection/interface/ElectronPlotter.h"
#include "../Selection/interface/MuonPlotter.h"
#include "../Selection/interface/JetPlotter.h"
#include "../Selection/interface/VertexPlotter.h"
#include "../Tools/interface/MVATrainer.h"
#include "../Tools/interface/MVAComputer.h"
#include "../Tools/interface/JetTools.h"
#include "../JESMeasurement/interface/JetCombiner.h"
#include "../Reconstruction/interface/JetCorrectorParameters.h"
#include "../Reconstruction/interface/JetCorrectionUncertainty.h"
#include "../Reconstruction/interface/Observables.h"
//#include "../Reconstruction/interface/PlotObservables.h"
#include "../Reconstruction/interface/TTreeObservables.h"
#include "../MCInformation/interface/LumiReWeighting.h"
#include "Style.C"

using namespace std;
using namespace TopTree;
using namespace RooFit;
using namespace reweight;


/// TGraphAsymmErrors
map<string,TGraphAsymmErrors*> graphAsymmErr;
map<string,TGraphErrors*> graphErr;

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

void ModifyHist(TH1F* &h, Color_t lcolor);

int main (int argc, char *argv[])
{

  bool doPF2PAT = false; //ignore this...

  clock_t start = clock();

  cout << "************************************************" << endl;
  cout << " Beginning of the program for Tprime searches ! " << endl;
  cout << "************************************************" << endl;

  //SetStyle if needed
  setTDRStyle();
  //setMyStyle();

  /////////////////////
  // Configuration
  /////////////////////

  //xml file
  string xmlFileName ="../config/myTprimeconfig.xml";
  if (argc >= 2)
    xmlFileName=string(argv[1]);
  const char *xmlfile = xmlFileName.c_str();
  int doJESShift = 0; // 0: off 1: minus 2: plus
   cout << "doJESShift: " << doJESShift << endl;
  int doJERShift = 0; // 0: off (except nominal scalefactor for jer) 1: minus 2: plus
   cout << "doJERShift: " << doJERShift << endl;
  int dobTagEffShift = 0; //0: off (except nominal scalefactor for btag eff) 1: minus 2:plus
   cout << "dobTagEffShift: " << dobTagEffShift << endl;
  int domisTagEffShift = 0; //0: off (except nominal scalefactor for mistag eff) 1: minus 2: plus
   cout << "domisTagEffShift: " << domisTagEffShift << endl;
  
  //TRandom3 * scalefactor_btageff = new TRandom3(0); // data/MC scalefactor = 0.94 +- 0.09
  float scalefactorbtageff;
  if(dobTagEffShift == 0)
     scalefactorbtageff = 0.94;
  if(dobTagEffShift == 1)
     scalefactorbtageff = 0.85;
  if(dobTagEffShift == 2)
     scalefactorbtageff = 1.03;
     
  float mistagfactor;
  if(domisTagEffShift == 0) // data/MC scalefactor = 1.21 +- 0.17
      mistagfactor = 1.21;
  if(domisTagEffShift == 1)
      mistagfactor = 1.04;
  if(domisTagEffShift == 2)
      mistagfactor = 1.38;	
		
  bool WJetsSmoothing = false, doBtagScaling = false; //switches off b-tagging and scales the WJets in the plots down with b-tag efficiency for WJets (has to be determined first)
  
  cout << "used config file: " << xmlfile << endl;

  bool writeTTrees = false;
  //Output ROOT file
  string rootFileName ("Tprime.root");
  if (argc >= 3){
		string sample=string(argv[2]);
 		rootFileName = "Tprime_"+ sample + ".root";
 	}
 	const char *rootfile = rootFileName.c_str();

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
  treeLoader.LoadDatasets (datasets, xmlfile);
  for(unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
  
  float Luminosity = oldLuminosity;

  Observables obs;
  vector <string > lstVar;
  lstVar = obs.ListOfVariables ();
	
  //vector of objects
  cout << " - Variable declaration ..." << endl;
  vector < TRootVertex* > vertex;
  vector < TRootMuon* > init_muons;
  vector < TRootElectron* > init_electrons;
  vector < TRootJet* > init_jets;
  vector < TRootJet* > init_jets_corrected;
  vector < TRootMET* > mets;
  vector < TRootMET* > mets_corrected;

  TFile *fout = new TFile (rootfile, "RECREATE");
  //Global variable
  TRootEvent* event = 0;
  
  string pathPNG = "PlotsTprime_png";
  if (argc >= 3){
		string sample=string(argv[2]);
 		pathPNG = pathPNG+"_"+sample;
  }
  pathPNG = pathPNG +"/"; 	
  
  mkdir(pathPNG.c_str(),0777);
	
	//nof selected events
  Double_t *nEvents = new Double_t[datasets.size()];
  
  // Define the plots
  MSPlot["nEventsAfterCuts"] = new MultiSamplePlot(datasets, "nEventsAfterCuts", 23, -0.5, 22.5, "Nr. of events after each cut");

  MSPlot["nPrimaryVertices"] = new MultiSamplePlot(datasets, "nPrimaryVertices", 12, -0.5, 11.5, "Nr. of primary vertices");
  MSPlot["nGoodPrimaryVertices"] = new MultiSamplePlot(datasets, "nGoodPrimaryVertices", 12, -0.5, 11.5, "Nr. of good primary vertices");
  // Variables used during event selection
  MSPlot["globalMuPt"] = new MultiSamplePlot(datasets, "globalMuPt", 200, 0, 200, "Global Muon P_{T}");
  MSPlot["globalMuEta"] = new MultiSamplePlot(datasets, "globalMuEta", 50, -2.5, 2.5, "Global Muon #eta");
  MSPlot["globalMuRelIso"] = new MultiSamplePlot(datasets, "globalMuRelIso", 100, -0.0001, 5, "Global Muon RelIso");
  MSPlot["globalMuChi2"] = new MultiSamplePlot(datasets, "globalMuChi2", 100, 0, 25, "Global Muon #Chi^{2}/ndf");
  MSPlot["globalMuNValidHits"] = new MultiSamplePlot(datasets, "globalMuNValidHits", 35, -0.5, 34.5, "Global Muon Number of Valid Hits");
  MSPlot["globalMuMinDRMuJet"] = new MultiSamplePlot(datasets, "globalMuMinDRMuJet", 110, 0, 5.5, "Global Muon Minimal #DeltaR(#mu, jet)");
  MSPlot["globalMud0"] = new MultiSamplePlot(datasets, "globalMud0", 100, -0.1, 0.1, "Global Muon d0");

  MSPlot["nSelectedMuons"] = new MultiSamplePlot(datasets, "nSelectedMuons", 5, -0.5, 4.5, "Nr. of selected Muons");

  MSPlot["otherGlobalMuPt"] = new MultiSamplePlot(datasets, "otherGlobalMuPt", 10, 0, 100, "Other Global Muon P_{T}");
  MSPlot["otherGlobalMuEta"] = new MultiSamplePlot(datasets, "otherGlobalMuEta", 7, -2.8, 2.8, "Other Global Muon #eta");
  MSPlot["otherGlobalMuRelIso"] = new MultiSamplePlot(datasets, "otherGlobalMuRelIso", 25, 0, 5, "Other Global Muon RelIso");
  
  MSPlot["nLooseOtherMuons"] = new MultiSamplePlot(datasets, "nLooseOtherMuons", 5, -0.5, 4.5, "Nr. of Loose Other Muons");

  MSPlot["ElectronEt"] = new MultiSamplePlot(datasets, "ElectronEt", 15, 0, 75, "Electron E_{T}");
  MSPlot["ElectronEta"] = new MultiSamplePlot(datasets, "ElectronEta", 13, -2.6, 2.6, "Electron #eta");
  MSPlot["ElectronRelIso"] = new MultiSamplePlot(datasets, "ElectronRelIso", 20, -0.0001, 10, "Electron RelIso");
  
  MSPlot["nLooseElectrons"] = new MultiSamplePlot(datasets, "nLooseElectrons", 5, -0.5, 4.5, "Nr. of Loose Electrons");
  
  MSPlot["allJetsPt1"] = new MultiSamplePlot(datasets, "allJetsPt1", 20, 0, 250, "Hardest Jet P_{T}");
  MSPlot["allJetsPt2"] = new MultiSamplePlot(datasets, "allJetsPt2", 20, 0, 200, "2nd Hardest Jet P_{T}");
  MSPlot["allJetsPt3"] = new MultiSamplePlot(datasets, "allJetsPt3", 25, 0, 150, "3rd Hardest Jet P_{T}");
  MSPlot["allJetsPt4"] = new MultiSamplePlot(datasets, "allJetsPt4", 25, 0, 100, "4th Hardest Jet P_{T}");
  MSPlot["allJetsEta"] = new MultiSamplePlot(datasets, "allJetsEta", 25, -2.5, 2.5, "Jet #eta");
  MSPlot["allJetsEMF"] = new MultiSamplePlot(datasets, "allJetsEMF", 20, 0, 1, "Jet EMF");
  MSPlot["allJetsn90Hits"] = new MultiSamplePlot(datasets, "allJetsn90Hits", 60, -0.5, 59.5, "Jet n90Hits");
  MSPlot["allJetsfHPD"] = new MultiSamplePlot(datasets, "allJetsfHPD", 25, 0, 1, "Jet fHPD");
  
  MSPlot["nSelectedJets"] = new MultiSamplePlot(datasets, "nSelectedJets", 10, -0.5, 9.5, "Nr. of Selected Jets");
  
  // Plots of selected events (before b-tagging)
  MSPlot["selectedEventsMuPt"] = new MultiSamplePlot(datasets, "selectedEventsMuPt", 8, 0, 160, "Selected Muon P_{T}");
  MSPlot["selectedEventsMuEta"] = new MultiSamplePlot(datasets, "selectedEventsMuEta", 11, -2.2, 2.2, "Selected Muon #eta");
  MSPlot["selectedEventsMuRelIso"] = new MultiSamplePlot(datasets, "selectedEventsMuRelIso", 10, -0.0001, 0.05, "Selected Muon RelIso");
  MSPlot["selectedEventsMuChi2"] = new MultiSamplePlot(datasets, "selectedEventsMuChi2", 20, 0, 10, "Selected Muon #Chi^{2}/ndf");
  MSPlot["selectedEventsMuNValidHits"] = new MultiSamplePlot(datasets, "selectedEventsMuNValidHits", 25, 9.5, 34.5, "Selected Muon Number of Valid Hits");
  MSPlot["selectedEventsMuMinDRMuJet"] = new MultiSamplePlot(datasets, "selectedEventsMuMinDRMuJet", 10, 0, 4, "Selected Muon Minimal #DeltaR(#mu, jet)");    
  MSPlot["selectedEventsMud0"] = new MultiSamplePlot(datasets, "selectedEventsMud0", 10, -0.02, 0.02, "Selected Muon d0");
  
  MSPlot["selectedEventsMET"] = new MultiSamplePlot(datasets, "selectedEventsMET", 50, 0, 500, "Selected Events MET");
  MSPlot["selectedEventsM3"] = new MultiSamplePlot(datasets, "selectedEventsM3", 20, 0, 1000, "Selected Events M3");

  MSPlot["selectedEventsJetsPt"] = new MultiSamplePlot(datasets, "selectedEventsJetsPt", 100, 0, 500, "Selected Jets P_{T}");
  MSPlot["selectedEventsJetsEta"] = new MultiSamplePlot(datasets, "selectedEventsJetsEta", 13, -2.6, 2.6, "Jet #eta");


  // Plots of selected events AFTER b-tagging
  MSPlot["btagged_selectedEventsMuPt"] = new MultiSamplePlot(datasets, "btagged_selectedEventsMuPt", 8, 0, 160, "Selected Muon P_{T}");
  MSPlot["btagged_selectedEventsMuEta"] = new MultiSamplePlot(datasets, "btagged_selectedEventsMuEta", 11, -2.2, 2.2, "Selected Muon #eta");
  MSPlot["btagged_selectedEventsMuRelIso"] = new MultiSamplePlot(datasets, "btagged_selectedEventsMuRelIso", 10, -0.0001, 0.05, "Selected Muon RelIso");
  MSPlot["btagged_selectedEventsMuChi2"] = new MultiSamplePlot(datasets, "btagged_selectedEventsMuChi2", 20, 0, 10, "Selected Muon #Chi^{2}/ndf");
  MSPlot["btagged_selectedEventsMuNValidHits"] = new MultiSamplePlot(datasets, "btagged_selectedEventsMuNValidHits", 25, 9.5, 34.5, "Selected Muon Number of Valid Hits");
  MSPlot["btagged_selectedEventsMuMinDRMuJet"] = new MultiSamplePlot(datasets, "btagged_selectedEventsMuMinDRMuJet", 10, 0, 4, "Selected Muon Minimal #DeltaR(#mu, jet)");    
  MSPlot["btagged_selectedEventsMud0"] = new MultiSamplePlot(datasets, "btagged_selectedEventsMud0", 10, -0.02, 0.02, "Selected Muon d0");
  
  MSPlot["btagged_selectedEventsMET"] = new MultiSamplePlot(datasets, "btagged_selectedEventsMET", 50, 0, 500, "Selected Events MET (GeV)");
  MSPlot["btagged_selectedEventsM3"] = new MultiSamplePlot(datasets, "btagged_selectedEventsM3", 20, 0, 1000, "Selected Events M3");

  MSPlot["btagged_selectedEventsJetsPt"] = new MultiSamplePlot(datasets, "btagged_selectedEventsJetsPt", 100, 0, 500, "Selected Jets P_{T}");
  MSPlot["btagged_selectedEventsJetsEta"] = new MultiSamplePlot(datasets, "btagged_selectedEventsJetsEta", 13, -2.6, 2.6, "Jet #eta");

  MSPlot["btagged_selectedEventsLeadingJetPt"] = new MultiSamplePlot(datasets, "btagged_selectedEventsLeadingJetPt", 100, 0, 500, "Leading jet P_{T}");
  MSPlot["btagged_selectedEventsLeadingJetEta"] = new MultiSamplePlot(datasets, "btagged_selectedEventsLeadingJetEta", 13, -2.6, 2.6, "Leading jet #eta");
  MSPlot["btagged_selectedEventsJet2Pt"] = new MultiSamplePlot(datasets, "btagged_selectedEventsJet2Pt", 100, 0, 500, "2nd jet P_{T}");
  MSPlot["btagged_selectedEventsJet2Eta"] = new MultiSamplePlot(datasets, "btagged_selectedEventsJet2Eta", 13, -2.6, 2.6, "2nd jet #eta");
  MSPlot["btagged_selectedEventsJet3Pt"] = new MultiSamplePlot(datasets, "btagged_selectedEventsJet3Pt", 100, 0, 500, "3rd jet P_{T}");
  MSPlot["btagged_selectedEventsJet3Eta"] = new MultiSamplePlot(datasets, "btagged_selectedEventsJet3Eta", 13, -2.6, 2.6, "3rd jet #eta");
  MSPlot["btagged_selectedEventsJet4Pt"] = new MultiSamplePlot(datasets, "btagged_selectedEventsJet4Pt", 100, 0, 500, "4th jet P_{T}");
  MSPlot["btagged_selectedEventsJet4Eta"] = new MultiSamplePlot(datasets, "btagged_selectedEventsJet4Eta", 13, -2.6, 2.6, "4th jet #eta");
  
  MSPlot["btagged_selectedEventsWHadronicJet1Pt"] = new MultiSamplePlot(datasets, "btagged_selectedEventsWHadronicJet1Pt", 100, 0, 500, "Hadronic W jet 1 (MVA) P_{T}");
  MSPlot["btagged_selectedEventsWHadronicJet1Eta"] = new MultiSamplePlot(datasets, "btagged_selectedEventsWHadronicJet1Eta", 13, -2.6, 2.6, "Hadronic W jet 1 (MVA) #eta");
  MSPlot["btagged_selectedEventsWHadronicJet2Pt"] = new MultiSamplePlot(datasets, "btagged_selectedEventsWHadronicJet2Pt", 100, 0, 500, "Hadronic W jet 2 (MVA) P_{T}");
  MSPlot["btagged_selectedEventsWHadronicJet2Eta"] = new MultiSamplePlot(datasets, "btagged_selectedEventsWHadronicJet2Eta", 13, -2.6, 2.6, "Hadronic W jet 2 (MVA) #eta");
  MSPlot["btagged_selectedEventsHadronicBJetPt"] = new MultiSamplePlot(datasets, "btagged_selectedEventsHadronicBJetPt", 100, 0, 500, "Hadronic b-jet (MVA) P_{T}");
  MSPlot["btagged_selectedEventsHadronicBJetEta"] = new MultiSamplePlot(datasets, "btagged_selectedEventsHadronicBJetEta", 100, 0, 500, "Hadronic b-jet (MVA) #eta");
  MSPlot["btagged_selectedEventsLeptonicBJetPt"] = new MultiSamplePlot(datasets, "btagged_selectedEventsLeptonicBJetPt", 100, 0, 500, "Leptonic b-jet (MVA) P_{T}");
  MSPlot["btagged_selectedEventsLeptonicBJetEta"] = new MultiSamplePlot(datasets, "btagged_selectedEventsLeptonicBJetEta", 100, 0, 500, "Leptonic b-jet (MVA) #eta");
	
	//control plots
	MSPlot["pT_bjet_aftermva"] = new MultiSamplePlot(datasets,"pT_bjet_aftermva",150,0,300,"b Jet P_{T} after MVA");
	MSPlot["pT_ljet_aftermva"] = new MultiSamplePlot(datasets,"pT_ljet_aftermva",150,0,300,"l Jets P_{T} after MVA");
	MSPlot["mW_aftermva"] = new MultiSamplePlot(datasets,"mW_aftermva",100,0,200,"m_{W} after MVA");
	MSPlot["mtop_aftermva"] = new MultiSamplePlot(datasets,"mtop_aftermva",200,0,400,"m_{top} after MVA");
	MSPlot["mvavalue_aftermva"] = new MultiSamplePlot(datasets,"mvavalue_aftermva",100,0,1,"maxMVA value");
	MSPlot["massmub_aftermva"] = new MultiSamplePlot(datasets,"massmub_aftermva",200,0,400,"m_{b#mu} after MVA");
	MSPlot["njets_aftermva"] = new MultiSamplePlot(datasets,"njets_aftermva",8,2.5,10.5,"Number of selected jets after MVA");		
  
  MSPlot["selectedEvents_btagvalueHadB_fromMVA"] = new MultiSamplePlot(datasets,"selectedEvents_btagvalueHadB_fromMVA",70,-10,40,"TCHEM of hadronic b after MVA");
  MSPlot["selectedEvents_btagvalueLepB_fromMVA"] = new MultiSamplePlot(datasets,"selectedEvents_btagvalueLepB_fromMVA",70,-10,40,"TCHEM of leptonic b after MVA");
  
	// maxMVA stuff
  histo1D["maxMVAweight_goodcomb"] = new TH1F("maxMVAweight_goodcomb","maxMVAweight_goodcomb",100,0,4);
  histo1D["maxMVAweight_badcomb"] = new TH1F("maxMVAweight_badcomb","maxMVAweight_badcomb",100,0,4);
  histo1D["maxMVAweight_AllMC"] = new TH1F("maxMVAweight_AllMC","maxMVAweight_AllMC",100,0,4);
	histo1D["maxMVAweight_Data"] = new TH1F("maxMVAweight_Data","maxMVAweight_Data",100,0,4);
  
  //pileup reweighting control plot
  MSPlot["lumiWeights_beforeselection"] = new MultiSamplePlot(datasets,"lumiWeights_beforeselection",100,0,4,"lumiWeights;lumiWeight;#events");

  //some (normalized)'key kinematic distributions'
  TH1F* hglobalMuPt_ttbar = new TH1F("globalMuPt_ttbar","Global Muon P_{T}", 50, 0, 300);
   TH1F* hglobalMuEta_ttbar = new TH1F("globalMuEta_ttbar","Global Muon #eta", 11, -2.2, 2.2);
   TH1F* hJetsPt_ttbar = new TH1F("JetsPt_ttbar","Jet P_{T}", 50, 0, 400);
   TH1F* hJetsEta_ttbar = new TH1F("JetsEta_ttbar","Jet #eta", 25, -2.5, 2.5);
   TH1F* hMET_ttbar = new TH1F("MET_ttbar","MET", 50, 0, 500);
  TH1F* hglobalMuPt_tprime350 = new TH1F("globalMuPt_tprime350","Global Muon P_{T}", 50, 0, 300);
   TH1F* hglobalMuEta_tprime350 = new TH1F("globalMuEta_tprime350","Global Muon #eta", 11, -2.2, 2.2);
   TH1F* hJetsPt_tprime350 = new TH1F("JetsPt_tprime350","Jet P_{T}", 50, 0, 400);
   TH1F* hJetsEta_tprime350 = new TH1F("JetsEta_tprime350","Jet #eta", 25, -2.5, 2.5);   
   TH1F* hMET_tprime350 = new TH1F("MET_tprime350","MET", 50, 0, 500);
  TH1F* hglobalMuPt_tprime400 = new TH1F("globalMuPt_tprime400","Global Muon P_{T}", 50, 0, 300);
   TH1F* hglobalMuEta_tprime400 = new TH1F("globalMuEta_tprime400","Global Muon #eta", 11, -2.2, 2.2);
   TH1F* hJetsPt_tprime400 = new TH1F("JetsPt_tprime400","Jet P_{T}", 50, 0, 400);
   TH1F* hJetsEta_tprime400 = new TH1F("JetsEta_tprime400","Jet #eta", 25, -2.5, 2.5);      
   TH1F* hMET_tprime400 = new TH1F("MET_tprime400","MET", 50, 0, 500);
  TH1F* hglobalMuPt_tprime450 = new TH1F("globalMuPt_tprime450","Global Muon P_{T}", 50, 0, 300);
   TH1F* hglobalMuEta_tprime450 = new TH1F("globalMuEta_tprime450","Global Muon #eta", 11, -2.2, 2.2);
   TH1F* hJetsPt_tprime450 = new TH1F("JetsPt_tprime450","Jet P_{T}", 50, 0, 400);
   TH1F* hJetsEta_tprime450 = new TH1F("JetsEta_tprime450","Jet #eta", 25, -2.5, 2.5);      
   TH1F* hMET_tprime450 = new TH1F("MET_tprime450","MET", 50, 0, 500);
  TH1F* hglobalMuPt_tprime500 = new TH1F("globalMuPt_tprime500","Global Muon P_{T}", 50, 0, 300);
   TH1F* hglobalMuEta_tprime500 = new TH1F("globalMuEta_tprime500","Global Muon #eta", 11, -2.2, 2.2);
   TH1F* hJetsPt_tprime500 = new TH1F("JetsPt_tprime500","Jet P_{T}", 50, 0, 400);
   TH1F* hJetsEta_tprime500 = new TH1F("JetsEta_tprime500","Jet #eta", 25, -2.5, 2.5);      
   TH1F* hMET_tprime500 = new TH1F("MET_tprime500","MET", 50, 0, 500);
   
   //for maxMVA comparison ttbar vs tprime...(normalized)
   TH1F* hMaxMVA_ttbar = new TH1F("MaxMVA_ttbar","Maximum MVA", 20, 0, 1);
   TH1F* hMaxMVA_tprime400 = new TH1F("MaxMVA_tprime400","Maximum MVA", 20, 0, 1);
   TH1F* hMaxMVA_tprime500 = new TH1F("MaxMVA_tprime500","Maximum MVA", 20, 0, 1);
   
   //generator level
   TH1F* hMCMuEta_ttbar = new TH1F("MCMuEta_ttbar","MC Muon #eta", 22, -2.2, 2.2);
   TH1F* hMCMuEta_tprime400 = new TH1F("MCMuEta_tprime400","MC Muon #eta", 22, -2.2, 2.2);
   TH1F* hMCMatchedQuarksEta_ttbar = new TH1F("MCMatchedQuarksEta_ttbar","MC MatchedQuarks #eta", 22, -2.2, 2.2);
   TH1F* hMCMatchedQuarksEta_tprime400 = new TH1F("MCMatchedQuarksEta_tprime400","MC MatchedQuarks #eta", 22, -2.2, 2.2);
   TH1F* hMCTopEta_ttbar = new TH1F("MCTopEta_ttbar","MC Top #eta", 22, -2.2, 2.2); // not anti-t
   TH1F* hMCTopEta_realTop_ttbar = new TH1F("MCTopEta_realTop_ttbar","MC Top #eta", 22, -2.2, 2.2); // not anti-t
   TH1F* hMCAntiTopEta_ttbar = new TH1F("MCTopEta_ttbar","MC Top #eta", 22, -2.2, 2.2); // anti-t, not t
   TH1F* hMCTprimeEta_tprime400 = new TH1F("MCTprimeEta_tprime400","MC Tprime #eta", 22, -2.2, 2.2); //not anti-t'
   TH1F* hMCAntiTprimeEta_tprime400 = new TH1F("MCTprimeEta_tprime400","MC Tprime #eta", 22, -2.2, 2.2); //anti-t', not t'
   TH1F* hMCTprimeEta_tprime500 = new TH1F("MCTprimeEta_tprime500","MC Tprime #eta", 22, -2.2, 2.2); //not anti-t'
  

   MSPlot["NumberofSelEvents_0tag"] = new MultiSamplePlot(datasets, "NselectedEvents_0tag", 1, 0, 1, "0 b-tags");
   MSPlot["NumberofSelEvents_1tag"] = new MultiSamplePlot(datasets, "NselectedEvents_1tag", 1, 0, 1, "1 b-tag");
   MSPlot["NumberofSelEvents_2tag"] = new MultiSamplePlot(datasets, "NselectedEvents_2tag", 1, 0, 1, "2 b-tags");
   MSPlot["NumberofSelEvents_tags"] = new MultiSamplePlot(datasets, "NselectedEvents_tags", 3, 0, 3, "selected events number of b-tags");

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
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-2);
  CutsSelecTable.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-1);
  CutsSelecTable.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets);
  CutsSelecTable.push_back(string(LabelNJets));
  CutsSelecTable.push_back(string("mva value"));
  if (verbose > 0)
    cout << " - CutsSelectionTable instantiated ..." << endl;
  SelectionTable selecTable(CutsSelecTable, datasets);
  selecTable.SetLuminosity(Luminosity);
  if (verbose > 0)
    cout << " - SelectionTable instantiated ..." << endl;
  

  /////////////////////////////////////
  // Initialize JetCombination Stuff //
  /////////////////////////////////////

  bool Tprime = true;
  bool TrainMVA = false; // If false, the previously trained MVA will be used to calculate stuff
  string MVAmethod = "Likelihood"; // MVAmethod to be used to get the good jet combi calculation (not for training! this is chosen in the jetcombiner class)
    
  JetCombiner* jetCombiner = new JetCombiner(TrainMVA, Luminosity, datasets, MVAmethod, Tprime);
  if (verbose > 0)
    cout << " - JetCombiner instantiated ..." << endl;

  ////////////////////////
  // NEW PileUp Reweighting //
  ////////////////////////

  // initialize LumiReWeighting stuff
  // Summer11 PU_S4, distribution obtained by averaging the number of interactions, taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupMCReweightingUtilities
   // in each beam crossing to estimate the true mean.  THIS IS THE RECOMMENDED ONE for reweighting.
/*
  Double_t probdistFlat10[25] = {
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0630151648,
    0.0526654164,
    0.0402754482,
    0.0292988928,
    0.0194384503,
    0.0122016783,
    0.007207042,
    0.004003637,
    0.0020278322,
    0.0010739954,
    0.0004595759,
    0.0002229748,
    0.0001028162,
    4.58337152809607E-05
  };
   
  Double_t MCLumi_f[25] = {
    0.104109,
    0.0703573,
    0.0698445,
    0.0698254,
    0.0697054,
    0.0697907,
    0.0696751,
    0.0694486,
    0.0680332,
    0.0651044,
    0.0598036,
    0.0527395,
    0.0439513,
    0.0352202,
    0.0266714,
    0.019411,
    0.0133974,
    0.00898536,
    0.0057516,
    0.00351493,
    0.00212087,
    0.00122891,
    0.00070592,
    0.000384744,
    0.000219377
  };
  
  Double_t TopDBDist2011Data_f[25] = {
    0.0127118660008111155,
    0.0273174253882752516,
    0.0647422373974094190,
    0.108494213975257103,
    0.140081296984992526,
    0.150411260268535935,
    0.142773479388604602,
    0.118012735306947752,
    0.0881395784021791473,
    0.0603740700218931975,
    0.0382939204454870868,
    0.0227366747939989136,
    0.0127228459417252551,
    0.00674674468025676568,
    0.00340977235841692389,
    0.00165292589154045016,
    0.000771798466244840342,
    0.000347480158040664431,
    0.000151563397272207710,
    0.0000642172483977206039,
    0.0000264962736283059724,
    0.0000106455374332742453,
    0.00000418355451211455042,
    0.00000161033109693768961,
    0.000000606815958689117662
  };
  
  Double_t TrueDist2011_f[25] = {
    0.019091,
    0.0293974,
    0.0667931,
    0.108859,
    0.139533,
    0.149342,
    0.138629,
    0.114582,
    0.0859364,
    0.059324,
    0.0381123,
    0.0229881,
    0.0131129,
    0.00711764,
    0.00369635,
    0.00184543,
    0.000889604,
    0.000415683,
    0.000188921,
    0.000146288,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0
  };
  
  vector<float> TrueDist2011, MClumi, Spring11MClumi, TopDBDist2011Data;
  for( int i=0; i<25; ++i)
  {
    TopDBDist2011Data.push_back(TopDBDist2011Data_f[i]);
    TrueDist2011.push_back(TrueDist2011_f[i]);
    MClumi.push_back(MCLumi_f[i]);
    Spring11MClumi.push_back(probdistFlat10[i]);
  }
*/

  cout << " Initializing LumiReWeighting stuff..." << endl;
//  LumiReWeighting LumiWeightsSpring11 = LumiReWeighting(Spring11MClumi, TrueDist2011);
//  LumiReWeighting LumiWeights = LumiReWeighting(MClumi, TopDBDist2011Data);
  LumiReWeighting LumiWeights = LumiReWeighting("ReweightHistos/pileup/WJets.root", "ReweightHistos/pileup/data_to_173692.root", "pileup", "pileup");
  cout << " Initialized LumiReWeighting stuff" << endl;
 
 

  ////////////////////////////////////
  //	Loop on datasets
  ////////////////////////////////////
  
  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
  for (unsigned int d = 0; d < datasets.size (); d++)
  { 
    string dataSetName = datasets[d]->Name();
    string dataSetTitle = datasets[d]->Title ();
    
    if (verbose > 1)
      cout << "   Dataset " << d << ": " << dataSetName << "/ title : " << dataSetTitle << endl;
    if (verbose > 1)
      std::cout<<"      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
		
		//open files and load
    cout<<"LoadDataset"<<endl;
    treeLoader.LoadDataset (datasets[d], anaEnv);
    
    string previousFilename = "";
    int iFile = -1;
       
    
    selecTable.Fill(d,0, datasets[d]->Xsection() * datasets[d]->EquivalentLumi() );
    
    /////////////////////////////////////
    /// Initialize JEC factors
    /////////////////////////////////////
   	
    vector<JetCorrectorParameters> vCorrParam;
    if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") // Data!
    {
      JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("JECFiles/Jec11V2_db_AK5PFchs_L2L3Residual.txt");
      vCorrParam.push_back(*ResJetCorPar);
    }
    /*else // MC!
    {
    }*/
    
////    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/START41_V0_AK5PF_Uncertainty.txt"); //newer
  
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/Jec11V2_db_AK5PFchs_Uncertainty.txt"); //newest (~7july2011)
    
    
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, false);

    vector<string> BranchesTtree;
    BranchesTtree = obs.ListOfVariables ();
    BranchesTtree.push_back("scaleFactor"); //not an observable...
    TTreeObservables TtreeObs(BranchesTtree,dataSetName); //first I used dataSetTitle
   
   
   
    ////////////////////////////////////
    //	Loop on events
    ////////////////////////////////////
    
    nEvents[d] = 0;
    int itrigger = -1, previousRun = -1;
    if (verbose > 1)
      cout << "	Loop over events " << endl;
    int start = 0;
    int end = datasets[d]->NofEvtsToRunOver();
    if(TrainMVA && datasets[d]->Name () == "TTbarJets_SemiMuon"){ 
      start = 0;
      end = int(datasets[d]->NofEvtsToRunOver()/3);
    }
    else if (!TrainMVA && datasets[d]->Name () == "TTbarJets_SemiMuon"){    
      start = int(datasets[d]->NofEvtsToRunOver()/3);
      end = datasets[d]->NofEvtsToRunOver();
    }
    else{
      start = 0;
      end = datasets[d]->NofEvtsToRunOver();
    } //commented for sync test
    
    for (unsigned int ievt = start; ievt < end; ievt++) //ievt < datasets[d]->NofEvtsToRunOver()
    {
      nEvents[d]++;
      if(ievt%1000 == 0)
        std::cout<<"Processing the "<<ievt<<"th event" <<flush<<"\r";
      
      //load event
      event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);

      // scale factor for the event
      float scaleFactor = 1.;
      
      //load pile-up weight for this event
      float lumiweight = 1;
      if(dataSetName != "Data" && dataSetName != "data" && dataSetName != "DATA"){
          float avPU = ( (float)event->nPu(-1) + (float)event->nPu(0) + (float)event->nPu(+1) ) / 3.;
          lumiweight = LumiWeights.ITweight(avPU);
	  if(verbose>2)cout<<"    pileup weight = "<<lumiweight<<endl;
      }


      //to check generator level of t' samples
      vector<TRootMCParticle*> mcParticles_copy;
      if(dataSetName.find("TTbarJets_SemiMu") == 0 || dataSetName.find("TTbarJets_Other") == 0 || dataSetName.find("Tprime") == 0)
      {
        mcParticles_copy = treeLoader.LoadMCPart(ievt);
        sort(mcParticles_copy.begin(),mcParticles_copy.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
      if(dataSetName.find("TTbarJets_SemiMu") == 0 || dataSetName.find("TTbarJets_Other") == 0){
         //cout<<endl<<"*****************"<<endl;
	 TRootMCParticle myW, myB, myTop;
         for(unsigned int i=0; i<mcParticles_copy.size(); i++)
         {  
	    //cout<<"mcParticles_copy["<<i<<"]->type() = "<<mcParticles_copy[i]->type()<<endl;
	    //cout<<"   status = "<<mcParticles_copy[i]->status()<<endl;
	    //cout<<"   mother = "<<mcParticles_copy[i]->motherType()<<endl;
	    //cout<<"   granny = "<<mcParticles_copy[i]->grannyType()<<endl;
	    
	    if(mcParticles_copy[i]->type() == 6){
	       hMCTopEta_realTop_ttbar->Fill(mcParticles_copy[i]->Eta()); //only top, not anti-top
	       //cout<<"mcParticles_copy["<<i<<"]->type() = "<<mcParticles_copy[i]->type()<<endl;
	    }
	    if(mcParticles_copy[i]->type() == 24 && mcParticles_copy[i]->motherType() == 6){
	      //cout<<"    status = "<<mcParticles_copy[i]->status()<<endl;
	      myW = (TRootMCParticle) (*(mcParticles_copy[i]));
	    }
	    if(mcParticles_copy[i]->type() == -5 && mcParticles_copy[i]->motherType() == 6){
	      myB = (TRootMCParticle) (*(mcParticles_copy[i]));
	    }
	    
	    if( mcParticles_copy[i]->status() != 3) continue;
	    if(fabs(mcParticles_copy[i]->type()) == 13 && fabs(mcParticles_copy[i]->motherType()) == 24 && fabs(mcParticles_copy[i]->grannyType()) == 6){
	    	hMCMuEta_ttbar->Fill(mcParticles_copy[i]->Eta());
	    }
	 
	 }
	 myTop = (myW + myB);
	 //cout<<" myTop Eta = "<<myTop.Eta()<<endl;
	 hMCTopEta_ttbar->Fill(myTop.Eta());
      }
      if(dataSetName.find("Tprime") == 0){
         //cout<<endl<<"*****************"<<endl;
	 TRootMCParticle myW, myB, myTprime;
         for(unsigned int i=0; i<mcParticles_copy.size(); i++)
         {
	    ////cout<<"mcParticles_copy["<<i<<"]->grannyType() = "<<mcParticles_copy[i]->grannyType()<<endl;
	    //cout<<"mcParticles_copy["<<i<<"]->type() = "<<mcParticles_copy[i]->type()<<endl;
	    //cout<<"   status = "<<mcParticles_copy[i]->status()<<endl;
	    //cout<<"   mother = "<<mcParticles_copy[i]->motherType()<<endl;
	    //cout<<"   granny = "<<mcParticles_copy[i]->grannyType()<<endl;
	    if(mcParticles_copy[i]->type() == 8){
	       if(dataSetName.find("Tprime400") == 0) if(dataSetName.find("Tprime400") == 0)hMCTprimeEta_tprime400->Fill(mcParticles_copy[i]->Eta()); //only t', not anti-t'
	       //cout<<"mcParticles_copy["<<i<<"]->type() = "<<mcParticles_copy[i]->type()<<endl;
	    }
	    if(mcParticles_copy[i]->type() == 24 && mcParticles_copy[i]->motherType() == 8){
	      myW = (TRootMCParticle) (*(mcParticles_copy[i]));
	    }
	    if(mcParticles_copy[i]->type() == -5 && mcParticles_copy[i]->motherType() == 8){
	      myB = (TRootMCParticle) (*(mcParticles_copy[i]));
	    }
	    
	    if( mcParticles_copy[i]->status() != 3) continue;
	    if(fabs(mcParticles_copy[i]->type()) == 13 && fabs(mcParticles_copy[i]->motherType()) == 24 && fabs(mcParticles_copy[i]->grannyType()) == 8){
	    	if(dataSetName.find("Tprime400") == 0) hMCMuEta_tprime400->Fill(mcParticles_copy[i]->Eta());
	    
	    }
	 
	 } 
	 myTprime = (myW + myB);
	 //cout<<" myTprime Eta = "<<myW.Eta()<<endl;
	 if(dataSetName.find("Tprime400") == 0) hMCTprimeEta_tprime400->Fill(myTprime.Eta());     
	 if(dataSetName.find("Tprime500") == 0) hMCTprimeEta_tprime500->Fill(myTprime.Eta());   
      } 


      vector<TRootGenJet*> genjets;
      if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
      {
        genjets = treeLoader.LoadGenJet(ievt);
        sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
      }  
      
      // Load the GenEvent and calculate the branching ratio correction
      if(dataSetName.find("TTbarJets") == 0)
      {
        TRootGenEvent* genEvt = treeLoader.LoadGenEvent(ievt);
        if( genEvt->isSemiLeptonic() )
		scaleFactor *= (0.108*9.)*(0.676*1.5);
	else if( genEvt->isFullHadronic() )
		scaleFactor *= (0.676*1.5)*(0.676*1.5);
	else if( genEvt->isFullLeptonic() )
  		scaleFactor *= (0.108*9.)*(0.108*9.);
      }

      
      // Clone the init_jets vector, otherwise the corrections will be removed
      for(unsigned int i=0; i<init_jets_corrected.size(); i++)
        if(init_jets_corrected[i]) delete init_jets_corrected[i];
      init_jets_corrected.clear();

      for(unsigned int i=0; i<init_jets.size(); i++)
        init_jets_corrected.push_back( (TRootJet*) init_jets[i]->Clone() );
      
      // Clone the mets vector, otherwise the corrections will be removed
      for(unsigned int i=0; i<mets_corrected.size(); i++)
        if(mets_corrected[i]) delete mets_corrected[i];
      mets_corrected.clear();

      for(unsigned int i=0; i<mets.size(); i++)
        mets_corrected.push_back( (TRootMET*) mets[i]->Clone() );
 
 
    // check which file in the dataset it is to have the HLTInfo right
    string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
    if(previousFilename != currentFilename)
    {
      	previousFilename = currentFilename;
	      iFile++;
	      cout<<"File changed!!! => iFile = "<<iFile<<endl;
    }
      
    int currentRun = event->runId();
    if(previousRun != currentRun)
    {
        previousRun = currentRun;
	if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
	{    
  				if (event->runId() >= 160431 && event->runId() <= 163261)
    				 itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v5"), currentRun, iFile);
  				else if (event->runId() >= 163270 && event->runId() <= 163869)
    				 itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v6"), currentRun, iFile);
  				else if (event->runId() >= 165088 && event->runId() <= 165633)
    				 itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v8"), currentRun, iFile);
  				else if (event->runId() >= 165970 && event->runId() <= 167043 && event->runId() != 166346)
    				 itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v9"), currentRun, iFile);
  				else if (event->runId() == 166346)
    				 itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v10"), currentRun, iFile);
  				else if (event->runId() >= 167078 && event->runId() <= 167913)
    				 itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v11"), currentRun, iFile);
				else if (event->runId() >= 170249 && event->runId() <= 172619) //Aug05ReReco: equivalent to the run range of PromptReco_v5 normally, but Aug05 replaces this. Warning: somewhere we last about 5/pb in this data?
				 itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_v8"), currentRun, iFile);
				else if (event->runId() >= 172620 && event->runId() <= 173198) //first part of PromptReco_v6, same as previous trigger
                                 itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_v8"), currentRun, iFile);
				else if (event->runId() >= 173236 && event->runId() <= 173692) //second part of PromptReco_v6
				 itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_v9"), currentRun, iFile);
				
  				if(itrigger == 9999){
    				 cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (DATA) IN RUN " << event->runId() << endl;
    				 exit(1);
  				}
	}
	else 
	{  
   				itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v5"), currentRun, iFile);//Summer11 MC has other triggers!	
    
  				if(itrigger == 9999)
				{
    				 cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (" << dataSetName << ") IN RUN " << event->runId() << endl;
    				 exit(1);
				}
  	}	
     }

    

     if( dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
     {       
       // Apply the scraping veto, not needed anymore, already done on PAT level in the toptree production
       /*bool isBeamBG = true;
       if(event->nTracks() > 10)
       {
          if( ( (float) event->nHighPurityTracks() ) / ( (float) event->nTracks() ) > 0.25 )  //this is changed: 0.2 and already applied on PAT level in the toptree production
            isBeamBG = false;
       }
       if(isBeamBG){
          cout<<"Beam scraping event!! Continuing..."<<endl;
	  continue;
       }*/
	
       // Apply Jet Corrections on-the-fly
       jetTools->correctJets(init_jets_corrected, vertex);
	 
       for(unsigned int i=0; i<init_jets_corrected.size(); i++){
         //cout<<"   init_jets["<<i<<"]->Pt() = "<<init_jets[i]->Pt()<<", init_jets_corrected["<<i<<"]->Pt() = "<<init_jets_corrected[i]->Pt()<<endl;
         mets_corrected[0]->SetPx(mets_corrected[0]->Px() + init_jets[i]->Px() - init_jets_corrected[i]->Px());
         mets_corrected[0]->SetPy(mets_corrected[0]->Py() + init_jets[i]->Py() - init_jets_corrected[i]->Py());
	 mets_corrected[0]->SetE(sqrt(pow(mets_corrected[0]->Px(),2) + pow(mets_corrected[0]->Py(),2)));
       }   	 
     
     }
     else if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
     {
        // Correct for the difference in muon efficiency (HLT and Id) between Data and MC
        //scaleFactor *= 0.965;  // Fall10 value...
        
        // Correct jets and MET for JES uncertainy systematics	
        if (doJESShift == 1)
	    jetTools->correctJetJESUnc(init_jets_corrected, mets_corrected[0], "minus"); //should I work with these 'corrected', these are actually in this stage the same as init_jets etc, but later on I use init_jets_corrected for the selection and so on... so I'll use 'corrected'
	else if (doJESShift == 2)
	    jetTools->correctJetJESUnc(init_jets_corrected, mets_corrected[0], "plus");      
	    	 
	 
	// Match RecoJets with GenJets
        vector< pair<size_t, size_t> > indexVector; //first index = RecoJet, second index = GenJet
        vector<bool> mLock(genjets.size(),false);   // when locked, genJet is already matched to a recoJet
        for(size_t i=0; i<init_jets_corrected.size(); i++)
        {
          pair<size_t, size_t> tmpIndex;
          float minDR = 9999.;
          for(size_t j=0; j<genjets.size(); j++)
          {
            if( ! mLock[j] )
            {
              if( init_jets_corrected[i]->DeltaR(*genjets[j]) < 0.4 && init_jets_corrected[i]->DeltaR(*genjets[j]) < minDR )
              {
                minDR = init_jets_corrected[i]->DeltaR(*genjets[j]);
                tmpIndex = pair<size_t, size_t>(i,j);
              }
            }
          }
          if(minDR < 999.)
          {
            mLock[tmpIndex.second] = true;
            indexVector.push_back(tmpIndex);
          }
        }

        // Apply correction for jet energy resolution on-the-fly, only for recoJets matched with a genJet
        for(size_t i=0; i<indexVector.size(); i++)
        {
          if( genjets[indexVector[i].second]->Pt() < 15 ) continue;
          //MSPlot["Pt_before_JERcorr"]->Fill(init_jets_corrected[indexVector[i].first]->Pt(), datasets[d], true, Luminosity*scaleFactor);
          mets_corrected[0]->SetPx(mets_corrected[0]->Px() + init_jets_corrected[indexVector[i].first]->Px());
	  mets_corrected[0]->SetPy(mets_corrected[0]->Py() + init_jets_corrected[indexVector[i].first]->Py());
	  float corrFactor = 0.1;
          float fabsEta = fabs(init_jets_corrected[indexVector[i].first]->Eta());
          if(doJERShift == 1)
          {
            if(fabsEta <= 1.5) corrFactor = 0.0;
            else if(fabsEta < 2.0 && fabsEta > 1.5) corrFactor = -0.05;
            else corrFactor = -0.1;
          }
          else if(doJERShift == 2)
          {
            if(fabsEta <= 1.5) corrFactor = 0.2;
            else if(fabsEta < 2.0 && fabsEta > 1.5) corrFactor = 0.25;
            else corrFactor = 0.3;
          }
          float deltapt = ( init_jets_corrected[indexVector[i].first]->Pt() - genjets[indexVector[i].second]->Pt() ) * corrFactor;
          float ptscale = max(0.0, ( init_jets_corrected[indexVector[i].first]->Pt() + deltapt) / init_jets_corrected[indexVector[i].first]->Pt() );
          if(ptscale > 0.0)
            init_jets_corrected[indexVector[i].first]->SetPxPyPzE(init_jets_corrected[indexVector[i].first]->Px()*ptscale, init_jets_corrected[indexVector[i].first]->Py()*ptscale,
              init_jets_corrected[indexVector[i].first]->Pz()*ptscale, init_jets_corrected[indexVector[i].first]->E()*ptscale);
          //MSPlot["Pt_after_JERcorr"]->Fill(init_jets_corrected[indexVector[i].first]->Pt(), datasets[d], true, Luminosity*scaleFactor);
          //if(init_jets_corrected[indexVector[i].first]->Pt() > 30)
          //  MSPlot["Pt_before_JERcorr_PtCorr30"]->Fill(init_jets_corrected[indexVector[i].first]->Pt()/ptscale, datasets[d], true, Luminosity*scaleFactor);
          mets_corrected[0]->SetPx(mets_corrected[0]->Px() - init_jets_corrected[indexVector[i].first]->Px());
	  mets_corrected[0]->SetPy(mets_corrected[0]->Py() - init_jets_corrected[indexVector[i].first]->Py());
	  mets_corrected[0]->SetE(sqrt(pow(mets_corrected[0]->Px(),2) + pow(mets_corrected[0]->Py(),2)));
	}
        
        //Scale jets with a certain factor
        //jetTools->scaleJets(init_jets_corrected, 1.);
	// apply PU Reweighting 
         //   cout << "before: " << scaleFactor << endl;
	 scaleFactor = scaleFactor*lumiweight; //pile-up reweighting
	 //  cout << "after: " << scaleFactor << endl;
     }
     
     if(doBtagScaling)
     {
       //b-tagging efficiencies...//remember: if the scalefactor is scaled in Tprime.cc, so in the trees, than don't scale in TprimeAnalysisTemplates a 2nd time!!
       if(datasets[d]->Name()=="TTbarJets_SemiMuon" || (datasets[d]->Name()).find("Tprime")<=(datasets[d]->Name()).size()){//check again!
	    //cout<<"Applying b-efficiency scale factor to TTbarJets_SemiMuon"<<endl;
	    float btageff_SF = 0.95;
	    scaleFactor = scaleFactor*btageff_SF;		 
       }
       if(datasets[d]->Name()=="TTbarJets_Other"){
	    //cout<<"Applying b-efficiency scale factor to TTbarJets_Other"<<endl;
	    float btageff_SF = 0.95;
	    scaleFactor = scaleFactor*btageff_SF;	 
       }
       if((datasets[d]->Name()).find("ST_")==0){
	    //cout<<"Applying b-efficiency scale factor to Single top"<<endl;
	    float btageff_SF = 0.95;
	    scaleFactor = scaleFactor*btageff_SF;	 
       }
     }     
     if(WJetsSmoothing) //remember: if the scalefactor is scaled in Tprime.cc, so in the trees, then don't scale in TprimeAnalysisTemplates a 2nd time!!
     {		
	if((datasets[d]->Name()).find("WJets")==0){
	    //cout<<"RESCALING WJETS (without btagging)"<<endl;
	    float btageff = 1.;
	    if(doJESShift==0) btageff= 0.14859; //(JES nominal)
	    else if(doJESShift==1) btageff= 0.14898; //JES minus
	    else if(doJESShift==2) btageff= 0.13874; //JES plus
	    
	    scaleFactor = scaleFactor*btageff; 
	}
     }   



      /////////////////////////////
      //   Selection
      /////////////////////////////
      
      //Declare selection instance    
      Selection selection(init_jets_corrected, init_muons, init_electrons, mets_corrected); //mets also corrected...
      
      
//      cout<<"itrigger = "<<itrigger<<endl;
      bool trigged = treeLoader.EventTrigged (itrigger);
      bool isGoodPV = selection.isPVSelected(vertex, anaEnv.PVertexNdofCut, anaEnv.PVertexZCut, anaEnv.PVertexRhoCut);
      //bool isGoodPV = selection.isPVSelected(vertex, 4,24,2.); //this is the RefSelV4 selection for the vertex, make sure this is put in the config file correctly

      vector<TRootJet*> selectedJets;
      vector<TRootMuon*> selectedMuons;
	   
      if (init_jets_corrected.size() > 0) {
	      if (init_jets_corrected[0]->jetType() == 1 || doPF2PAT) { // calojets
	        //cout << "Selecting for caloJets" << endl;
	        selectedJets = selection.GetSelectedJets(true);
	        selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
	      }	else {
	        //cout << "Selecting for PF/JPT jets" << endl;
	        vector<TRootMuon*> overlapMuons = selection.GetSelectedMuons(vertex[0]);
	        selection.setJetCuts(35.,2.4,0.01,1.,0.98,0.3,0.1);//selection.setJetCuts(30.,2.4,0.01,1.,0.98,0.3,0.1); // refSelV4 values
	        //selection.setJetCuts(30.,2.4,0.01,1.,0.98,0.3,0.1); //for sync test
		selectedJets = selection.GetSelectedJets(overlapMuons,true);
		selection.setMuonCuts(35.,2.1,0.125,10,0.02,0.3,1,1,1);//added this, 35 GeV instead of 20 GeV; RefSelV4: selection.setMuonCuts(20,2.1,0.05,10,0.02,0.3,1,1,1). //new (PF isolation): selection.setMuonCuts(20,2.1,0.125,10,0.02,0.3,1,1,1). 
		//selection.setMuonCuts(20.,2.1,0.125,10,0.02,0.3,1,1,1); //for sync test
		selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
	      }
      }
      vector<TRootMuon*> vetoMuons = selection.GetSelectedLooseMuons();
      vector<TRootElectron*> vetoElectrons = selection.GetSelectedLooseElectrons(false);

      MSPlot["nPrimaryVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);

      bool eventSelected = false;

      MSPlot["nEventsAfterCuts"]->Fill(0, datasets[d], true, Luminosity*scaleFactor);
      selecTable.Fill(d,1,1);
      
      if(trigged){
//      cout<<"TRIGGED!"<<endl;
        MSPlot["lumiWeights_beforeselection"]->Fill(lumiweight, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["nEventsAfterCuts"]->Fill(1, datasets[d], true, Luminosity*scaleFactor);
        selecTable.Fill(d,2,1);
        
        if(isGoodPV){
          MSPlot["nEventsAfterCuts"]->Fill(2, datasets[d], true, Luminosity*scaleFactor);
          selecTable.Fill(d,3,1);
      	  MSPlot["nGoodPrimaryVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
          				      
	  for(unsigned int i=0; i<init_jets_corrected.size(); i++){	  
		if(dataSetName == "TTbarJets_SemiMuon" || dataSetName == "TTbarJets_Other"){
			hJetsPt_ttbar->Fill(init_jets_corrected[i]->Pt(),scaleFactor);
			hJetsEta_ttbar->Fill(init_jets_corrected[i]->Eta(),scaleFactor);
		}
		if(dataSetName == "Tprime350"){
			hJetsPt_tprime350->Fill(init_jets_corrected[i]->Pt(),scaleFactor);
			hJetsEta_tprime350->Fill(init_jets_corrected[i]->Eta(),scaleFactor);
		}
		if(dataSetName == "Tprime400"){
			hJetsPt_tprime400->Fill(init_jets_corrected[i]->Pt(),scaleFactor);
			hJetsEta_tprime400->Fill(init_jets_corrected[i]->Eta(),scaleFactor);
		}
		if(dataSetName == "Tprime450"){
			hJetsPt_tprime450->Fill(init_jets_corrected[i]->Pt(),scaleFactor);
			hJetsEta_tprime450->Fill(init_jets_corrected[i]->Eta(),scaleFactor);
		}
		if(dataSetName == "Tprime500"){
			hJetsPt_tprime500->Fill(init_jets_corrected[i]->Pt(),scaleFactor);
			hJetsEta_tprime500->Fill(init_jets_corrected[i]->Eta(),scaleFactor);
		}
	  }
	  if(dataSetName == "TTbarJets_SemiMuon" || dataSetName == "TTbarJets_Other"){
		hMET_ttbar->Fill(mets_corrected[0]->Et(),scaleFactor);
	  }
	  if(dataSetName == "Tprime350"){
		hMET_tprime350->Fill(mets_corrected[0]->Et(),scaleFactor);
	  }
	  if(dataSetName == "Tprime400"){
		hMET_tprime400->Fill(mets_corrected[0]->Et(),scaleFactor);
	  }
	  if(dataSetName == "Tprime450"){
		hMET_tprime450->Fill(mets_corrected[0]->Et(),scaleFactor);
	  }
	  if(dataSetName == "Tprime500"){
		hMET_tprime500->Fill(mets_corrected[0]->Et(),scaleFactor);	  
	  }
	  
	  for(unsigned int i=0; i<init_muons.size(); i++)
          {
            if(init_muons[i]->isGlobalMuon() && init_muons[i]->Pt() > 18)
            {
              MSPlot["globalMuPt"]->Fill(init_muons[i]->Pt(), datasets[d], true, Luminosity*scaleFactor);
              MSPlot["globalMuEta"]->Fill(init_muons[i]->Eta(), datasets[d], true, Luminosity*scaleFactor);
              MSPlot["globalMuChi2"]->Fill(init_muons[i]->chi2(), datasets[d], true, Luminosity*scaleFactor);
              MSPlot["globalMuNValidHits"]->Fill(init_muons[i]->nofValidHits(), datasets[d], true, Luminosity*scaleFactor);
              MSPlot["globalMud0"]->Fill(init_muons[i]->d0(), datasets[d], true, Luminosity*scaleFactor);
              
	      	      
	      if(dataSetName == "TTbarJets_SemiMuon" || dataSetName == "TTbarJets_Other"){
	         hglobalMuPt_ttbar->Fill(init_muons[i]->Pt(), scaleFactor); //doesn't need to be normalized to the luminosity, I will normalize to 1
	         hglobalMuEta_ttbar->Fill(init_muons[i]->Eta(), scaleFactor);
	      }
	      if(dataSetName == "Tprime350"){
	         hglobalMuPt_tprime350->Fill(init_muons[i]->Pt(), scaleFactor); //doesn't need to be normalized to the luminosity, I will normalize to 1
	         hglobalMuEta_tprime350->Fill(init_muons[i]->Eta(), scaleFactor);
	      }
	      if(dataSetName == "Tprime400"){
	         hglobalMuPt_tprime400->Fill(init_muons[i]->Pt(), scaleFactor); //doesn't need to be normalized to the luminosity, I will normalize to 1
	         hglobalMuEta_tprime400->Fill(init_muons[i]->Eta(), scaleFactor);
	      }
	      if(dataSetName == "Tprime450"){
	         hglobalMuPt_tprime450->Fill(init_muons[i]->Pt(), scaleFactor); //doesn't need to be normalized to the luminosity, I will normalize to 1
	         hglobalMuEta_tprime450->Fill(init_muons[i]->Eta(), scaleFactor);
	      }
	      if(dataSetName == "Tprime500"){
	         hglobalMuPt_tprime500->Fill(init_muons[i]->Pt(), scaleFactor); //doesn't need to be normalized to the luminosity, I will normalize to 1
	         hglobalMuEta_tprime500->Fill(init_muons[i]->Eta(), scaleFactor);
	      }

//              if()
              {
                MSPlot["globalMuRelIso"]->Fill(init_muons[i]->relativeIso03(), datasets[d], true, Luminosity*scaleFactor);

                float mindRMuJet = 999.;
                for(unsigned int j=0;j<selectedJets.size();j++)
                {
                  float dRMuJet = init_muons[i]->DeltaR(*selectedJets[j]);
                  if(dRMuJet < mindRMuJet) mindRMuJet = dRMuJet;
                }
                if(mindRMuJet != 999.)
                  MSPlot["globalMuMinDRMuJet"]->Fill(mindRMuJet, datasets[d], true, Luminosity*scaleFactor);
              }
            }
          }
          MSPlot["nSelectedMuons"]->Fill(selectedMuons.size(), datasets[d], true, Luminosity*scaleFactor);
          
          if(selectedMuons.size()==1){
            MSPlot["nEventsAfterCuts"]->Fill(3, datasets[d], true, Luminosity*scaleFactor);
  	    selecTable.Fill(d,4,1);
  		      
  	    unsigned int nMuonsFound = 0, indexDuplicatedMuon = -1;
  	    for(unsigned int i=0; i<init_muons.size(); i++)
            {
              if(init_muons[i]->isGlobalMuon())
              {
                if((init_muons[i]->Pt() - selectedMuons[0]->Pt())<0.0000001 && (init_muons[i]->Eta() - selectedMuons[0]->Eta())<0.0000001 && (init_muons[i]->relativeIso03() - selectedMuons[0]->relativeIso03())<0.0000001 &&
                  init_muons[i]->isTrackerMuon() && init_muons[i]->idGlobalMuonPromptTight())
                {
                  nMuonsFound++;
                  if(nMuonsFound == 2) indexDuplicatedMuon = i;
                }
                else
                {
   	    	  MSPlot["otherGlobalMuPt"]->Fill(init_muons[i]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                  MSPlot["otherGlobalMuEta"]->Fill(init_muons[i]->Eta(), datasets[d], true, Luminosity*scaleFactor);
                  MSPlot["otherGlobalMuRelIso"]->Fill(init_muons[i]->relativeIso03(), datasets[d], true, Luminosity*scaleFactor);
                }
  	      }
  	    }
  	    if(nMuonsFound != 1)
  	    {
  		        /*cerr<<"nMuonsfound = "<<nMuonsFound<<" !!!!  This should be equal to 1!!!"<<endl;
  		        if(nMuonsFound>1)
  		        {
    		        cerr<<"selMu:  Pt: "<<selectedMuons[0]->Pt()<<"  Eta: "<<selectedMuons[0]->Eta()<<"  Phi: "<<selectedMuons[0]->Phi()
    		          <<"  RelIso: "<<selectedMuons[0]->relativeIso03()<<"  Chi2: "<<selectedMuons[0]->chi2()<<"  nHits: "<<selectedMuons[0]->nofValidHits()
    		          <<"  globalMuPrmptTight: "<<selectedMuons[0]->idGlobalMuonPromptTight()<<"  trackerMu: "<<selectedMuons[0]->isTrackerMuon()
    		          <<"  d0: "<<selectedMuons[0]->d0()<<endl;
    		        cerr<<"DuplMu:  Pt: "<<init_muons[indexDuplicatedMuon]->Pt()<<"  Eta: "<<init_muons[indexDuplicatedMuon]->Eta()<<"  Phi: "<<init_muons[indexDuplicatedMuon]->Phi()
    		          <<"  RelIso: "<<init_muons[indexDuplicatedMuon]->relativeIso03()<<"  Chi2: "<<init_muons[indexDuplicatedMuon]->chi2()<<"  nHits: "<<init_muons[indexDuplicatedMuon]->nofValidHits()
    		          <<"  globalMuPrmptTight: "<<init_muons[indexDuplicatedMuon]->idGlobalMuonPromptTight()<<"  trackerMu: "<<init_muons[indexDuplicatedMuon]->isTrackerMuon()
    		          <<"  d0: "<<init_muons[indexDuplicatedMuon]->d0()<<endl;
  		        }*/
  	     }
  		        
  	     MSPlot["nLooseOtherMuons"]->Fill(vetoMuons.size()-1, datasets[d], true, Luminosity*scaleFactor);
  		        
             if(vetoMuons.size()==1)
	     {
               MSPlot["nEventsAfterCuts"]->Fill(4, datasets[d], true, Luminosity*scaleFactor);
               selecTable.Fill(d,5,1);
              
               for(unsigned int i=0; i<init_electrons.size(); i++)
               {
                  MSPlot["ElectronEt"]->Fill(init_electrons[i]->Et(), datasets[d], true, Luminosity*scaleFactor);
                  MSPlot["ElectronEta"]->Fill(init_electrons[i]->Eta(), datasets[d], true, Luminosity*scaleFactor);
                  MSPlot["ElectronRelIso"]->Fill(init_electrons[i]->combinedIso(3,3,3)/init_electrons[i]->Et(), datasets[d], true, Luminosity*scaleFactor);
               }
              
               MSPlot["nLooseElectrons"]->Fill(vetoElectrons.size(), datasets[d], true, Luminosity*scaleFactor);
              
               if(vetoElectrons.size()==0)
               {
                  MSPlot["nEventsAfterCuts"]->Fill(5, datasets[d], true, Luminosity*scaleFactor);
                  selecTable.Fill(d,6,1);
                
                  sort(init_jets_corrected.begin(),init_jets_corrected.end(),HighestPt()); // HighestPt() is included from the Selection class
                
                  int previousJet = 0;
                
                  MSPlot["nSelectedJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
                
                  if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-3)
                  {
                    MSPlot["nEventsAfterCuts"]->Fill(6, datasets[d], true, Luminosity*scaleFactor);
                    selecTable.Fill(d,7,1);
                    if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-2)
                    {
                      MSPlot["nEventsAfterCuts"]->Fill(7, datasets[d], true, Luminosity*scaleFactor);
                      selecTable.Fill(d,8,1);
                      if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-1)
                      {
                        MSPlot["nEventsAfterCuts"]->Fill(8, datasets[d], true, Luminosity*scaleFactor);
                        selecTable.Fill(d,9,1);		 
                        if(selectedJets.size()>=(unsigned int)anaEnv.NofJets) //this has to be changed when you want exactly 4 selected jets...
                        {			
		          if(mets_corrected[0]->Et()>40){ //commented for sync test
                             MSPlot["nEventsAfterCuts"]->Fill(9, datasets[d], true, Luminosity*scaleFactor);
                             selecTable.Fill(d,10,1);
                             eventSelected = true;
                           //cout<<" EVENT SELECTED!"<<endl;
                        
			   // plot some properties of the selected events
                             MSPlot["selectedEventsMuPt"]->Fill(selectedMuons[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                             MSPlot["selectedEventsMuEta"]->Fill(selectedMuons[0]->Eta(), datasets[d], true, Luminosity*scaleFactor);
                             //MSPlot["selectedEventsMuRelIso"]->Fill(selectedMuons[0]->relativeIso03(), datasets[d], true, Luminosity*scaleFactor); //old reliso
                             MSPlot["selectedEventsMuChi2"]->Fill(selectedMuons[0]->chi2(), datasets[d], true, Luminosity*scaleFactor);
                             MSPlot["selectedEventsMuNValidHits"]->Fill(selectedMuons[0]->nofValidHits(), datasets[d], true, Luminosity*scaleFactor);
                             MSPlot["selectedEventsMud0"]->Fill(selectedMuons[0]->d0(), datasets[d], true, Luminosity*scaleFactor);
                             MSPlot["selectedEventsMET"]->Fill( mets_corrected[0]->Et(), datasets[d], true, Luminosity*scaleFactor);
			     
			     float mindRMuJet = 999.;
                             for(unsigned int j=0;j<selectedJets.size();j++)
                             {
                                float dRMuJet = selectedMuons[0]->DeltaR(*selectedJets[j]);
                                if(dRMuJet < mindRMuJet) mindRMuJet = dRMuJet;
                             }
                             MSPlot["selectedEventsMuMinDRMuJet"] ->Fill(mindRMuJet, datasets[d], true, Luminosity*scaleFactor);   
                        
                             if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
                                cout << event->runId() << ":" << event->eventId() << ":" << event->lumiBlockId()  << endl;
                        
                             float M3 = -1, maxPt = -1;
			     int previousJet = 0;
			     for(int i=0;i<selectedJets.size();i++)
			     {
				                  for(int j=0;j<i;j++)
				                  {
					                  for(int k=0;k<j;k++)
					                  {
						                  float combinedPt = (*selectedJets[i] + *selectedJets[j] + *selectedJets[k]).Pt();
						                  if(combinedPt > maxPt)
						                  {
							                  maxPt = combinedPt;
							                  M3 = (*selectedJets[i] + *selectedJets[j] + *selectedJets[k]).M();
						                  }
					                  }
				                  }
				                  
				                  //Fill jet plots inside the first jet-loop
				                  MSPlot["selectedEventsJetsPt"]->Fill(selectedJets[i]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                                                  MSPlot["selectedEventsJetsEta"]->Fill(selectedJets[i]->Eta(), datasets[d], true, Luminosity*scaleFactor);
			                         
						 /*added the following myself...*/
						  
					         if(previousJet == 0){
					           MSPlot["allJetsPt1"]->Fill(selectedJets[i]->Pt(), datasets[d], true, Luminosity*scaleFactor); 
					           previousJet = 1;
						 }
						 else if(previousJet == 1){
					           MSPlot["allJetsPt2"]->Fill(selectedJets[i]->Pt(), datasets[d], true, Luminosity*scaleFactor); 
					           previousJet = 2;
						 }
						 else if(previousJet == 2){
                      				   MSPlot["allJetsPt3"]->Fill(selectedJets[i]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                      				   previousJet = 3;
                    				 }
                    				 else if(previousJet == 3){
                      				   MSPlot["allJetsPt4"]->Fill(selectedJets[i]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                      				   previousJet = 4;
                    				 }	
                    				 MSPlot["allJetsEta"]->Fill(selectedJets[i]->Eta(), datasets[d], true, Luminosity*scaleFactor);
						 
			     }
                             MSPlot["selectedEventsM3"]->Fill(M3, datasets[d], true, Luminosity*scaleFactor);
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
       
      if( ! (trigged && isGoodPV) ) continue;
      // Don't look further at events with a bad primary vertex or not passing the trigger

      if(selectedJets.size() < 4) continue;

      if(selectedMuons.size() != 1) continue;

      vector<TRootMCParticle*> mcParticles;
      if(dataSetName.find("TTbarJets_SemiMu") == 0 || dataSetName.find("Tprime") == 0)
      {
        mcParticles = treeLoader.LoadMCPart(ievt);
        sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
      

      sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class
      jetCombiner->ProcessEvent(datasets[d], mcParticles, selectedJets, selectedMuons[0], vertex[0], eventSelected, init_electrons, init_muons, scaleFactor);


      if ( !TrainMVA && eventSelected )
      {
        //cout<<"Doing stuff with selected events..."<<endl;
        //get the MC matched jet combination, not the MVA best matched
	vector<unsigned int> goodCombi = jetCombiner->GetGoodJetCombination(); 
        
        // check if the b-jet contains a muon
        bool muonInBJet = false;
        for(int i=0; i<init_muons.size(); i++)
        {
          if(init_muons[i]->Pt() > 5 && fabs( init_muons[i]->Eta() ) < 2.5)
          {
  	        if (goodCombi[2] < selectedJets.size())
	            if (init_muons[i]->DeltaR(*selectedJets[goodCombi[2]]) < 0.5) muonInBJet = true;
          }
        }
//        if(muonInBJet) continue; // skip events with muon in b-jet!

        // check if the b-jet contains an electron
        bool electronInBJet = false;
        for(int i=0; i<init_electrons.size(); i++)
        {
          if(init_electrons[i]->Pt() > 5 && fabs( init_electrons[i]->Eta() ) < 2.5)
          {
  	        if (goodCombi[2] < selectedJets.size())
	            if (init_electrons[i]->DeltaR(*selectedJets[goodCombi[2]]) < 0.5) electronInBJet = true;
          }
        }
//        if(electronInBJet) continue; // skip events with electron in b-jet!

	pair<float, vector<unsigned int> > MVAvals = jetCombiner->getMVAValue(MVAmethod, 1); // 1 means the highest MVA value
				
	float mW_maxMVA = (*selectedJets[MVAvals.second[0]] + *selectedJets[MVAvals.second[1]]).M();
        float mTop_maxMVA = (*selectedJets[MVAvals.second[0]] + *selectedJets[MVAvals.second[1]] + *selectedJets[MVAvals.second[2]]).M();

 	      
				//expected corrections before mvacut
				//if( dataSetName.find("TTbarJets_SemiMu") == 0 && jetCombiner->isGoodJetCombination(MVAmethod, 1, false) == true )
		    //    jetCombiner->FillExpCorr(expCorrNoCut, vertex[0], init_muons, init_electrons, MVAvals.first);
      	  
				//if(dataSetName.find("TTbarJets_SemiMu") == 0 && goodCombi[2] >= 0 && goodCombi[0] >= 0 && goodCombi[1] >= 0)
				//XXXXXXXXX;

	   TLorentzVector Wmass = *selectedJets[MVAvals.second[0]]+*selectedJets[MVAvals.second[1]];
	   TLorentzVector topmass = Wmass+*selectedJets[MVAvals.second[2]];
	   TLorentzVector mub = *selectedJets[MVAvals.second[3]]+*selectedMuons[0];
        
	 if(dataSetName == "TTbarJets_SemiMuon" || dataSetName == "TTbarJets_Other") hMaxMVA_ttbar->Fill(MVAvals.first); //fill with max mva value among all jet combinations
	 if(dataSetName == "Tprime400") hMaxMVA_tprime400->Fill(MVAvals.first);
	 if(dataSetName == "Tprime500") hMaxMVA_tprime500->Fill(MVAvals.first);
	
				//cout << " W mass = "<< Wmass.M() << endl;         
	 bool performBtag = true, EventBtagged = false; // I mean, if there is at least one of the 'b-jets' b-tagged, EventBtagged should be true
	 string myChoice("AtleastSingleBtag");
	 if(dataSetName.find("WJets")==0 && WJetsSmoothing){
	    performBtag = false;
	 }
	 float btagvalueHadB = 0, btagvalueLepB = 0;
	 btagvalueHadB = selectedJets[MVAvals.second[2]]->btag_trackCountingHighEffBJetTags();
	 btagvalueLepB = selectedJets[MVAvals.second[3]]->btag_trackCountingHighEffBJetTags();
	 //cout<<"btagvalueHadB = "<<btagvalueHadB<<",  btagvalueLepB"<<btagvalueLepB<<endl;
	 float btagCuts[] = {1.7,3.3,10.2}; //trackcountinghighefficiency working points: loose, medium, tight. I will first work with medium. NOTE: if you want to change, also the scale factor of the b-tagging eff and the uncertainties should be changed!!

         MSPlot["selectedEvents_btagvalueHadB_fromMVA"]->Fill(btagvalueHadB,datasets[d], true, Luminosity*scaleFactor);
         MSPlot["selectedEvents_btagvalueLepB_fromMVA"]->Fill(btagvalueLepB,datasets[d], true, Luminosity*scaleFactor);

	 if(!(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA"))
	 {
	    //cout << "WILL DO THE TAGGING AND UNTAGGING GAME FOR THIS EVENT" << endl;
	    /*if(fabs(selectedJets[MVAvals.second[2]]->partonFlavour()) == 5 )
	    { // is a true b-jet
	      if(btagvalueHadB > btagCuts[1])
	      {
	        float random_btag = scalefactor_btageff->Uniform();
	        if(random_btag > scalefactorbtageff){ 
		    btagvalueHadB = -1000;//jet not tagged if random dice > 0.94 for nominal samples and > 0.85 for "-"	
		    //cout << "		the jet will be untagged because random dice is " << random_btag << endl;
		}
	      }	    
	    }
	    else{
	      if(btagvalueHadB > btagCuts[1])
	      {
	         scaleFactor = scaleFactor * mistagfactor;
              }	    
	    }
	    
	    if(fabs(selectedJets[MVAvals.second[3]]->partonFlavour()) == 5 )
	    { // is a true b-jet
	      if(btagvalueLepB > btagCuts[1])
	      {
	        float random_btag = scalefactor_btageff->Uniform();
	        if(random_btag > scalefactorbtageff){ 
		    btagvalueLepB = -1000;//jet not tagged if random dice > 0.94 for nominal samples and > 0.85 for "-"	
		    //cout << "		the jet will be untagged because random dice is " << random_btag << endl;
		}
	      }	    
	    }
	    else{
	      if(btagvalueLepB > btagCuts[1])
	      {
	         scaleFactor = scaleFactor * mistagfactor;
              }	    
	    }*/
	    
	    //for c-quark use alse b-tag eff scale factor, for now... https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagWeight	    
	    if(btagvalueHadB>btagCuts[1] && btagvalueLepB>btagCuts[1]) //2-tag
	    {
	       if((fabs(selectedJets[MVAvals.second[2]]->partonFlavour()) == 5 || fabs(selectedJets[MVAvals.second[2]]->partonFlavour()) == 4) && (fabs(selectedJets[MVAvals.second[3]]->partonFlavour()) == 5 || fabs(selectedJets[MVAvals.second[3]]->partonFlavour()) == 4))
	          scaleFactor = (1 - pow(1-scalefactorbtageff,2)*(1-scalefactorbtageff)/fabs(1-scalefactorbtageff));
	       else if((fabs(selectedJets[MVAvals.second[2]]->partonFlavour()) != 5 && fabs(selectedJets[MVAvals.second[2]]->partonFlavour()) != 4) && (fabs(selectedJets[MVAvals.second[3]]->partonFlavour()) != 5 && fabs(selectedJets[MVAvals.second[3]]->partonFlavour()) != 4))
	          scaleFactor = (1 - pow(1-mistagfactor,2)*(1-mistagfactor)/fabs(1-mistagfactor));
	       else
	          scaleFactor = scaleFactor * scalefactorbtageff * mistagfactor; //if 1 true b-jet and 1 non-true b-jet...
	    
	    }
	    else if(btagvalueHadB>btagCuts[1] || btagvalueLepB>btagCuts[1]) //else, exactly 1-tag
	    {	    
	       if(btagvalueHadB>btagCuts[1])
	       {
	          if(fabs(selectedJets[MVAvals.second[2]]->partonFlavour()) == 5 || fabs(selectedJets[MVAvals.second[2]]->partonFlavour()) == 4)
	             scaleFactor = scaleFactor * scalefactorbtageff;
	          else
	             scaleFactor = scaleFactor * mistagfactor;
	       }
	       else if(btagvalueLepB>btagCuts[1])
	       {
	          if(fabs(selectedJets[MVAvals.second[3]]->partonFlavour()) == 5 || fabs(selectedJets[MVAvals.second[3]]->partonFlavour()) == 4)
	             scaleFactor = scaleFactor * scalefactorbtageff;
	          else
	             scaleFactor = scaleFactor * mistagfactor;
	       }
	    }
	    
	 }
	 
	 if(myChoice=="AtleastSingleBtag") EventBtagged = (btagvalueHadB>btagCuts[1] || btagvalueLepB>btagCuts[1]);
	 if(myChoice=="DoubleBtag") EventBtagged = (btagvalueHadB>btagCuts[1] && btagvalueLepB>btagCuts[1]); //a lot of events may be rejected when you demand TWO jets to be tagged...	 
	 
	 //for the number of b-tags; WARNING: with the b-tag eff scale factors the different b-tag bins should be revisited... then final selection should be okay, however.
	 if(!(btagvalueHadB>btagCuts[1]) && !(btagvalueLepB>btagCuts[1]))
	 {
	    //0 tags
	    MSPlot["NumberofSelEvents_0tag"]->Fill(0,datasets[d], true, Luminosity*scaleFactor);
	    MSPlot["NumberofSelEvents_tags"]->Fill(0,datasets[d], true, Luminosity*scaleFactor);
	 }
	 if((btagvalueHadB>btagCuts[1] || btagvalueLepB>btagCuts[1]) && !(btagvalueHadB>btagCuts[1] && btagvalueLepB>btagCuts[1]))
	 {
	     //1 tag
	     MSPlot["NumberofSelEvents_1tag"]->Fill(0,datasets[d], true, Luminosity*scaleFactor);
	     MSPlot["NumberofSelEvents_tags"]->Fill(1,datasets[d], true, Luminosity*scaleFactor);
	 }
	 if(btagvalueHadB>btagCuts[1] && btagvalueLepB>btagCuts[1])
	 {
	     //2 tags
	     MSPlot["NumberofSelEvents_2tag"]->Fill(0,datasets[d], true, Luminosity*scaleFactor);
	     MSPlot["NumberofSelEvents_tags"]->Fill(2,datasets[d], true, Luminosity*scaleFactor);
	 }
	 
	 if(!performBtag || EventBtagged)
	 {
	      //cout<<"Event b-tagged"<<endl;
	     //filled without b-tag eff scale factor... and Wjets heavy flavour factor (1.5)... to be changed?
	      MSPlot["btagged_selectedEventsMuPt"]->Fill(selectedMuons[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
              MSPlot["btagged_selectedEventsMuEta"]->Fill(selectedMuons[0]->Eta(), datasets[d], true, Luminosity*scaleFactor);
              MSPlot["btagged_selectedEventsMET"]->Fill( mets_corrected[0]->Et(), datasets[d], true, Luminosity*scaleFactor);     	      

              for(unsigned int i=0;i<selectedJets.size();i++){
	        MSPlot["btagged_selectedEventsJetsPt"]->Fill(selectedJets[i]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["btagged_selectedEventsJetsEta"]->Fill(selectedJets[i]->Eta(), datasets[d], true, Luminosity*scaleFactor);
	      }

              MSPlot["btagged_selectedEventsLeadingJetPt"]->Fill(selectedJets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
              MSPlot["btagged_selectedEventsLeadingJetEta"]->Fill(selectedJets[0]->Eta(), datasets[d], true, Luminosity*scaleFactor);
              MSPlot["btagged_selectedEventsJet2Pt"]->Fill(selectedJets[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
              MSPlot["btagged_selectedEventsJet2Eta"]->Fill(selectedJets[1]->Eta(), datasets[d], true, Luminosity*scaleFactor);
              MSPlot["btagged_selectedEventsJet3Pt"]->Fill(selectedJets[2]->Pt(), datasets[d], true, Luminosity*scaleFactor);
              MSPlot["btagged_selectedEventsJet3Eta"]->Fill(selectedJets[2]->Eta(), datasets[d], true, Luminosity*scaleFactor);
              MSPlot["btagged_selectedEventsJet4Pt"]->Fill(selectedJets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
              MSPlot["btagged_selectedEventsJet4Eta"]->Fill(selectedJets[3]->Eta(), datasets[d], true, Luminosity*scaleFactor);
	      
	      MSPlot["btagged_selectedEventsWHadronicJet1Pt"]->Fill(selectedJets[MVAvals.second[0]]->Pt(), datasets[d], true, Luminosity*scaleFactor);
              MSPlot["btagged_selectedEventsWHadronicJet1Eta"]->Fill(selectedJets[MVAvals.second[0]]->Eta(), datasets[d], true, Luminosity*scaleFactor);
              MSPlot["btagged_selectedEventsWHadronicJet2Pt"]->Fill(selectedJets[MVAvals.second[1]]->Pt(), datasets[d], true, Luminosity*scaleFactor);
              MSPlot["btagged_selectedEventsWHadronicJet2Eta"]->Fill(selectedJets[MVAvals.second[1]]->Eta(), datasets[d], true, Luminosity*scaleFactor);
              MSPlot["btagged_selectedEventsHadronicBJetPt"]->Fill(selectedJets[MVAvals.second[2]]->Pt(), datasets[d], true, Luminosity*scaleFactor);
              MSPlot["btagged_selectedEventsHadronicBJetEta"]->Fill(selectedJets[MVAvals.second[2]]->Eta(), datasets[d], true, Luminosity*scaleFactor);
              MSPlot["btagged_selectedEventsLeptonicBJetPt"]->Fill(selectedJets[MVAvals.second[3]]->Pt(), datasets[d], true, Luminosity*scaleFactor);
              MSPlot["btagged_selectedEventsLeptonicBJetEta"]->Fill(selectedJets[MVAvals.second[3]]->Eta(), datasets[d], true, Luminosity*scaleFactor);
	    
	    //piece ~ from ObsTreeMaker.cc
	    //cout<<"LOOP OVER VARIABLES"<<endl;
	    for (unsigned int v = 0; v < lstVar.size (); v++) {	     			
		   //	if (!anaEnv.runOnObs ((int) v))
	           //	  continue;		   
	      Observables obsEvt (*(selectedMuons[0]),*selectedJets[MVAvals.second[0]],*selectedJets[MVAvals.second[1]], *selectedJets[MVAvals.second[2]], *selectedJets[MVAvals.second[3]], *(mets_corrected[0]));
	      //the following is to set the scaleFactor of the event, to be able to store this in a branch in the end... the scaleFactor is not an observable of course, so in the future, a lot of care has to be taken w.r.t. e.g. the ObservablesRanker		   		   		   
	      obsEvt.setEventscaleFactor(scaleFactor); //for data, scaleFactor should be 1	        		
	      //cout<<"... Event is going to be filled"<<endl;
	      bool fill;fill=0;
	      if (v==lstVar.size()-1)  //just to fill TtreeObs only once instead of in each iteration of the loop over variables
		 fill=1;
              TtreeObs.Fill (obsEvt,fill);
	      if(lstVar[v]=="MassHadTop" &&  (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")){
	      	if(obsEvt.Variable("MassHadTop") > 1000) cout<<" ------------------>>>>>>>>>>>> Event with MassHadTop  = "<<obsEvt.Variable("MassHadTop")<<": "<< event->runId() << ":" << event->eventId() << ":" << event->lumiBlockId() <<"  <<<<<<<<<<<<------------------"<<endl;
	      }
	      //if(lstVar[v]=="MassHadTop"){
	      //  hMassHadTop->Fill(obsEvt.Variable("MassHadTop"));
              //}                		
	    }			//loop over variables	
	 }
	 else continue;
	
      }// end !TrainMVA && eventSelected
      
        //delete selection;
    }	//loop on events
    cout<<endl;
    
    //if( anaEnv.nPseudoExp==0 || dataSetTitle=="Data")
      TtreeObs.Write (true); //write all tree OBS for each dataset
      
    //important: free memory
    //treeLoader.UnLoadDataset();
    if(jetTools) delete jetTools;
  }				//loop on datasets
  
  //Once everything is filled ...
  if (verbose > 0)
    cout << " We ran over all the data ;-)" << endl;
  
  //Selection tables
  selecTable.TableCalculator(false, true, true, true);
  string selectiontable = "SelectionTable_Tprime";
	if (argc >= 3){
		string sample=string(argv[2]);
 		selectiontable = selectiontable +"_"+sample;
 	}
	selectiontable = selectiontable +".tex"; 	
	selecTable.Write(selectiontable.c_str());

  // Do some special things with certain plots (normalize, BayesDivide, ... )
  if (verbose > 0)
    cout << "Treating the special plots." << endl;
  
  ///////////////////
  // Writing
  //////////////////
  if (verbose > 1)
  	cout << " - Writing outputs to the files ..." << endl;

	string pathPNGJetCombi = pathPNG+"JetCombination/";
	string pathPNGExpCorr = pathPNG+"ExpCorr/";

  mkdir((pathPNG+"MSPlot/").c_str(),0777);
  mkdir(pathPNGJetCombi.c_str(),0777);
  mkdir(pathPNGExpCorr.c_str(),0777);
  
  jetCombiner->Write(fout, true, pathPNGJetCombi);
  
  // 1D 
  TDirectory* th1dir = fout->mkdir("1D_histograms");
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
    MultiSamplePlot *temp = it->second;
    string name = it->first;
    temp->Draw(false, name, true, true, true, true, true,5); //ttbar stond op true, maar door iets in MS werd het toch niet samen genomen
    temp->Write(fout, name, true, pathPNG+"MSPlot/");
  }

  //Write histograms
  fout->cd();
  th1dir->cd();

  fout->cd();

	for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
	{
		TH1F *temp = it->second;
//		int N = temp->GetNbinsX();
//  	temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
//  	temp->SetBinContent(N+1,0);
//		temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
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

	//Write TGraphAsymmErrors
	for(map<string,TGraphAsymmErrors*>::const_iterator it = graphAsymmErr.begin(); it != graphAsymmErr.end(); it++)
	{
	  TGraphAsymmErrors *temp = it->second;
	  temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
		tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
	}
	
  fout->cd();
  //add configuration info
  fout->cd();
  configTree->Fill();
  configTree->Write();


  //Write TGraphErrors
  fout->cd();
  for(map<string,TGraphErrors*>::const_iterator it = graphErr.begin(); it != graphErr.end(); it++)
  {
    TGraphErrors *temp = it->second;
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
  }

  //write the normalized 'key kinematic distributions to compare ttbar with t' pair
  fout->cd();
  
  TCanvas *c1 = new TCanvas("cglobalMuPt","cglobalMuPt",500,500);
  TLegend* legend_c1 = new TLegend(0.65,0.80,0.89,0.70);
  ModifyHist(hglobalMuPt_ttbar,1);
  ModifyHist(hglobalMuPt_tprime350,kGreen);
  ModifyHist(hglobalMuPt_tprime400,kBlue);
  ModifyHist(hglobalMuPt_tprime450,kMagenta);
  ModifyHist(hglobalMuPt_tprime500,kYellow);  
  hglobalMuPt_ttbar->Draw();
  hglobalMuPt_tprime350->Draw("SAME");
  hglobalMuPt_tprime400->Draw("SAME");
  hglobalMuPt_tprime450->Draw("SAME");
  hglobalMuPt_tprime500->Draw("SAME");  
  legend_c1->SetTextFont(70);
  legend_c1->SetTextSize(0.03);
  legend_c1->AddEntry(hglobalMuPt_ttbar,"ttbar+jets","L");
  legend_c1->AddEntry(hglobalMuPt_tprime350,"t' pair, m(t') = 350 GeV","L");
  legend_c1->AddEntry(hglobalMuPt_tprime400,"t' pair, m(t') = 400 GeV","L");
  legend_c1->AddEntry(hglobalMuPt_tprime450,"t' pair, m(t') = 450 GeV","L");
  legend_c1->AddEntry(hglobalMuPt_tprime500,"t' pair, m(t') = 500 GeV","L");
  legend_c1->Draw("SAME");
  c1->Write();
   
   
  TCanvas *c2 = new TCanvas("cglobalMuEta","cglobalMuEta",500,500);
  TLegend* legend_c2 = new TLegend(0.65,0.80,0.89,0.70);
  ModifyHist(hglobalMuEta_ttbar,1);
  ModifyHist(hglobalMuEta_tprime350,kGreen);
  ModifyHist(hglobalMuEta_tprime400,kBlue);
  ModifyHist(hglobalMuEta_tprime450,kMagenta);
  ModifyHist(hglobalMuEta_tprime500,kYellow);
  hglobalMuEta_tprime500->Draw();  //first this, maximum of plot...
  hglobalMuEta_ttbar->Draw("SAME");
  hglobalMuEta_tprime350->Draw("SAME");
  hglobalMuEta_tprime400->Draw("SAME");
  hglobalMuEta_tprime450->Draw("SAME");
  legend_c2->SetTextFont(70);
  legend_c2->SetTextSize(0.03);
  legend_c2->AddEntry(hglobalMuEta_ttbar,"ttbar+jets","L");
  legend_c2->AddEntry(hglobalMuEta_tprime350,"t' pair, m(t') = 350 GeV","L");
  legend_c2->AddEntry(hglobalMuEta_tprime400,"t' pair, m(t') = 400 GeV","L");
  legend_c2->AddEntry(hglobalMuEta_tprime450,"t' pair, m(t') = 450 GeV","L");
  legend_c2->AddEntry(hglobalMuEta_tprime500,"t' pair, m(t') = 500 GeV","L");
  legend_c2->Draw("SAME");
  c2->Write();
  
  TCanvas *c3 = new TCanvas("cJetsPt","cJetsPt",500,500);
  TLegend* legend_c3 = new TLegend(0.65,0.80,0.89,0.70);
  ModifyHist(hJetsPt_ttbar,1);
  ModifyHist(hJetsPt_tprime350,kGreen);
  ModifyHist(hJetsPt_tprime400,kBlue);
  ModifyHist(hJetsPt_tprime450,kMagenta);
  ModifyHist(hJetsPt_tprime500,kYellow);  
  hJetsPt_ttbar->Draw();
  hJetsPt_tprime350->Draw("SAME");
  hJetsPt_tprime400->Draw("SAME");
  hJetsPt_tprime450->Draw("SAME");
  hJetsPt_tprime500->Draw("SAME");  
  legend_c3->SetTextFont(70);
  legend_c3->SetTextSize(0.03);
  legend_c3->AddEntry(hJetsPt_ttbar,"ttbar+jets","L");
  legend_c3->AddEntry(hJetsPt_tprime350,"t' pair, m(t') = 350 GeV","L");
  legend_c3->AddEntry(hJetsPt_tprime400,"t' pair, m(t') = 400 GeV","L");
  legend_c3->AddEntry(hJetsPt_tprime450,"t' pair, m(t') = 450 GeV","L");
  legend_c3->AddEntry(hJetsPt_tprime500,"t' pair, m(t') = 500 GeV","L");
  legend_c3->Draw("SAME");
  c3->Write();
  
  TCanvas *c3p1 = new TCanvas("cJetsEta","cJetsEta",500,500);
  TLegend* legend_c3p1 = new TLegend(0.65,0.80,0.89,0.70);
  ModifyHist(hJetsEta_ttbar,1);
  ModifyHist(hJetsEta_tprime350,kGreen);
  ModifyHist(hJetsEta_tprime400,kBlue);
  ModifyHist(hJetsEta_tprime450,kMagenta);
  ModifyHist(hJetsEta_tprime500,kYellow);  
  hJetsEta_tprime500->Draw(); //first this, maximum of plot...
  hJetsEta_ttbar->Draw("SAME");
  hJetsEta_tprime350->Draw("SAME");
  hJetsEta_tprime400->Draw("SAME");
  hJetsEta_tprime450->Draw("SAME");  
  legend_c3p1->SetTextFont(70);
  legend_c3p1->SetTextSize(0.03);
  legend_c3p1->AddEntry(hJetsEta_ttbar,"ttbar+jets","L");
  legend_c3p1->AddEntry(hJetsEta_tprime350,"t' pair, m(t') = 350 GeV","L");
  legend_c3p1->AddEntry(hJetsEta_tprime400,"t' pair, m(t') = 400 GeV","L");
  legend_c3p1->AddEntry(hJetsEta_tprime450,"t' pair, m(t') = 450 GeV","L");
  legend_c3p1->AddEntry(hJetsEta_tprime500,"t' pair, m(t') = 500 GeV","L");
  legend_c3p1->Draw("SAME");
  c3p1->Write();
    
  TCanvas *c4 = new TCanvas("cMET","cMET",500,500);
  TLegend* legend_c4 = new TLegend(0.65,0.80,0.89,0.70);
  ModifyHist(hMET_ttbar,1);
  ModifyHist(hMET_tprime350,kGreen);
  ModifyHist(hMET_tprime400,kBlue);
  ModifyHist(hMET_tprime450,kMagenta);
  ModifyHist(hMET_tprime500,kYellow);  
  hMET_ttbar->Draw();
  hMET_tprime350->Draw("SAME");
  hMET_tprime400->Draw("SAME");
  hMET_tprime450->Draw("SAME");
  hMET_tprime500->Draw("SAME");  
  legend_c4->SetTextFont(70);
  legend_c4->SetTextSize(0.03);
  legend_c4->AddEntry(hMET_ttbar,"ttbar+jets","L");
  legend_c4->AddEntry(hMET_tprime350,"t' pair, m(t') = 350 GeV","L");
  legend_c4->AddEntry(hMET_tprime400,"t' pair, m(t') = 400 GeV","L");
  legend_c4->AddEntry(hMET_tprime450,"t' pair, m(t') = 450 GeV","L");
  legend_c4->AddEntry(hMET_tprime500,"t' pair, m(t') = 500 GeV","L");
  legend_c4->Draw("SAME");
  c4->Write();  

  TCanvas *c5 = new TCanvas("cMaxMVA","cMaxMVA",500,500);
  TLegend* legend_c5 = new TLegend(0.65,0.80,0.89,0.70);
  ModifyHist(hMaxMVA_ttbar,1);
  ModifyHist(hMaxMVA_tprime400,kBlue);
  ModifyHist(hMaxMVA_tprime500,kYellow);  
  hMaxMVA_tprime500->Draw(); //first this, maximum of plot...
  hMaxMVA_ttbar->Draw("SAME");
  hMaxMVA_tprime400->Draw("SAME");  
  legend_c5->SetTextFont(70);
  legend_c5->SetTextSize(0.03);
  legend_c5->AddEntry(hMaxMVA_ttbar,"ttbar+jets","L");
  legend_c5->AddEntry(hMaxMVA_tprime400,"t' pair, m(t') = 400 GeV","L");
  legend_c5->AddEntry(hMaxMVA_tprime500,"t' pair, m(t') = 500 GeV","L");
  legend_c5->Draw("SAME");
  c5->Write();
  
  TCanvas *c6 = new TCanvas("cMCMuEta","cMCMuEta",500,500);
  TLegend* legend_c6 = new TLegend(0.65,0.80,0.89,0.70);
  ModifyHist(hMCMuEta_ttbar,1);
  ModifyHist(hMCMuEta_tprime400,kBlue); 
  hMCMuEta_tprime400->Draw(); //first this, maximum of plot...
  hMCMuEta_ttbar->Draw("SAME");
  legend_c6->SetTextFont(70);
  legend_c6->SetTextSize(0.03);
  legend_c6->AddEntry(hMCMuEta_ttbar,"ttbar+jets","L");
  legend_c6->AddEntry(hMCMuEta_tprime400,"t' pair, m(t') = 400 GeV","L");
  legend_c6->Draw("SAME");
  c6->Write();
     
  
  TCanvas *c7 = new TCanvas("cMCTopEta","cMCTopEta",500,500);
  TLegend* legend_c7 = new TLegend(0.65,0.80,0.89,0.70);
  ModifyHist(hMCTopEta_ttbar,1);
  ModifyHist(hMCTopEta_realTop_ttbar,kMagenta);
  ModifyHist(hMCTprimeEta_tprime400,kBlue);
  ModifyHist(hMCTprimeEta_tprime500,kRed);
  hMCTprimeEta_tprime400->Draw(); //first this, maximum of plot...
  hMCTprimeEta_tprime500->Draw("SAME");
  hMCTopEta_ttbar->Draw("SAME");
  hMCTopEta_realTop_ttbar->Draw("SAME");
  legend_c7->SetTextFont(70);
  legend_c7->SetTextSize(0.03);
  legend_c7->AddEntry(hMCTopEta_ttbar,"ttbar+jets ((W+b)->Eta())","L");
  legend_c7->AddEntry(hMCTopEta_realTop_ttbar,"ttbar+jets (Top->Eta())","L");
  legend_c7->AddEntry(hMCTprimeEta_tprime400,"t' pair, m(t') = 400 GeV ((W+b)->Eta())","L");
  legend_c7->AddEntry(hMCTprimeEta_tprime500,"t' pair, m(t') = 500 GeV ((W+b)->Eta())","L");
  legend_c7->Draw("SAME");
  c7->Write();  
  
  
  
  
  //
  if (verbose > 1)
    cout << " - Done with writing the module outputs in the ouput file ..." << endl;
  cout << " - Closing the output file now..." << endl;
//  fout->Write();
  fout->Close();

  //delete
  if(jetCombiner) delete jetCombiner;
  
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

void ModifyHist (TH1F* &h, Color_t lcolor)
{
	double temp_integral;

	h->SetLineColor(lcolor);
	temp_integral = h->Integral();
	//cout << temp_integral << endl;
	if(temp_integral!=0) h->Scale(pow(temp_integral,-1));
}

