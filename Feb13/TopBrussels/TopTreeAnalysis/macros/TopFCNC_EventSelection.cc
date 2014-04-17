#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include "TRandom3.h"

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "../Selection/interface/SelectionTable.h"
#include "../Content/interface/AnalysisEnvironment.h"
#include "../Content/interface/Dataset.h"
#include "../Tools/interface/JetTools.h"
#include "../Tools/interface/PlottingTools.h"
#include "../Tools/interface/MultiSamplePlot.h"
#include "../Tools/interface/TTreeLoader.h"
#include "../Tools/interface/AnalysisEnvironmentLoader.h"
#include "../Reconstruction/interface/JetCorrectorParameters.h"
#include "../Reconstruction/interface/JetCorrectionUncertainty.h"
#include "../Reconstruction/interface/MakeBinning.h"
#include "../MCInformation/interface/Lumi3DReWeighting.h"
#include "../TopFCNC/interface/TopFCNC_Evt.h"

#include "Style.C"

using namespace std;
using namespace TopTree;

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

struct HighestTCHEBtag{
    bool operator()( TRootJet* j1, TRootJet* j2 ) const{
    	return j1->btag_trackCountingHighEffBJetTags() > j2->btag_trackCountingHighEffBJetTags();
    }
};
struct HighestCVSBtag{
    bool operator()( TRootJet* j1, TRootJet* j2 ) const{
    	return j1->btag_combinedSecondaryVertexBJetTags() > j2->btag_combinedSecondaryVertexBJetTags();
    }
};

int main (int argc, char *argv[])
{
  int doJESShift = 0; // 0: off 1: minus 2: plus
  cout << "doJESShift: " << doJESShift << endl;

  int doJERShift = 0; // 0: off (except nominal scalefactor for jer) 1: minus 2: plus
  cout << "doJERShift: " << doJERShift << endl;

  int dobTagEffShift = 0; //0: off (except nominal scalefactor for btag eff) 1: minus 2: plus
  cout << "dobTagEffShift: " << dobTagEffShift << endl;

  int domisTagEffShift = 0; //0: off (except nominal scalefactor for mistag eff) 1: minus 2: plus
  cout << "domisTagEffShift: " << domisTagEffShift << endl;

  int doPUShift = 0; //0: off (except nominal PU reweighting) 1: minus 2: plus
  cout << "doPUShift: " << doPUShift << endl;

  string btagger = "TCHEM";
// b-tag scalefactor => TCHEL: data/MC scalefactor = 0.95 +- 0.10,    TCHEM: data/MC scalefactor = 0.94 +- 0.09
// mistag scalefactor => TCHEL: data/MC scalefactor = 1.11 +- 0.12,    TCHEM: data/MC scalefactor = 1.21 +- 0.17
  float scalefactorbtageff, mistagfactor;
  if(btagger == "TCHEL") //track counting high eff loose working point
  {
	  if(dobTagEffShift == 0)
		scalefactorbtageff = 0.95;
	  if(dobTagEffShift == 1)
		scalefactorbtageff = 0.85;
	  if(dobTagEffShift == 2)
		scalefactorbtageff = 1.05;
	  
	  if(domisTagEffShift == 0)
		mistagfactor = 1.11;
	  if(domisTagEffShift == 1)
		mistagfactor = 0.99;
	  if(domisTagEffShift == 2)
		mistagfactor = 1.23;
		
  }
  else if(btagger == "TCHEM") //track counting high eff medium working point
  {
  	  if(dobTagEffShift == 0)
		scalefactorbtageff = 0.94;
	  if(dobTagEffShift == 1)
		scalefactorbtageff = 0.85;
	  if(dobTagEffShift == 2)
		scalefactorbtageff = 1.03;
		
	  if(domisTagEffShift == 0)
		mistagfactor = 1.21;
	  if(domisTagEffShift == 1)
		mistagfactor = 1.04;
	  if(domisTagEffShift == 2)
		mistagfactor = 1.38;
  }
  float workingpointvalue = 9999; //{1.7,3.3,10.2}; trackcountinghighefficiency working points: loose, medium, tight
  if(btagger == "TCHEL")
     workingpointvalue = 1.7;
  else if(btagger == "TCHEM")
     workingpointvalue = 3.3;	

  clock_t start = clock();

  cout << "*************************************************************" << endl;
  cout << " Beginning of the program for the FCNC search ! " << endl;
  cout << "*************************************************************" << endl;

  //SetStyle if needed
  //setTDRStyle();
  setGregStyle();
  //setMyStyle();

  string postfix = "_EventSelection"; // to relabel the names of the output file

  if (doJESShift == 1)
    postfix= postfix+"_JESMinus";
  if (doJESShift == 2)
    postfix= postfix+"_JESPlus";
  if (doJERShift == 1)
    postfix= postfix+"_JERMinus";
  if (doJERShift == 2)
    postfix= postfix+"_JERPlus";
  if (dobTagEffShift == 1)
    postfix= postfix+"_bTagMinus";
  if (dobTagEffShift == 2)
    postfix= postfix+"_bTagPlus";
  if(domisTagEffShift == 1)
    postfix= postfix+"_misTagMinus";
  if(domisTagEffShift == 2)
    postfix= postfix+"_misTagPlus";

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////// Configuration ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string channelpostfix = "";
  string xmlFileName = "";

  bool diElectron = false; // use diElectron channel?
  bool diMuon = true; // use diMuon channel?
  if(diElectron && diMuon){
	cout << "  --> Using both diMuon and diElectron channel? Choose only one (for the moment, since this requires running on different samples/skims)!" << endl;
	exit(1);
  }

  if(diMuon){
	cout << " --> Using the diMuon channel..." << endl;
	channelpostfix = "_DiMu";
	xmlFileName = "../config/myTopFCNCconfig_Muon.xml";
  }
  else if(diElectron){
	cout << " --> Using the diElectron channel..." << endl;
	channelpostfix = "_DiEl";
	xmlFileName = "../config/myTopFCNCconfig_Electron.xml";
  }

  const char *xmlfile = xmlFileName.c_str();
  cout << "used config file: " << xmlfile << endl;    
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////// AnalysisEnvironment /////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  AnalysisEnvironment anaEnv;
  cout<<" - Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  int verbose = 2;//anaEnv.Verbose;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////// Load Datasets ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TTreeLoader treeLoader;
  vector < Dataset* > datasets;
  cout << " - Loading datasets ..." << endl;
  treeLoader.LoadDatasets (datasets, xmlfile);
  float Luminosity = 5000;
  
  for (unsigned int d = 0; d < datasets.size (); d++)
  {
    	string dataSetName = datasets[d]->Name();
	if(dataSetName.find("Data")<=0 || dataSetName.find("data")<=0 || dataSetName.find("DATA")<=0)
	{
		  Luminosity = datasets[d]->EquivalentLumi();
		  break;
	 }
  }
  cout << "Rescaled to an integrated luminosity of "<< Luminosity << endl;

  //Output ROOT file
  string rootFileName ("TopFCNC"+postfix+channelpostfix+".root");
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");

  //vector of objects
  cout << " - Variable declaration ..." << endl;
  vector < TRootVertex* >   vertex;
  vector < TRootMuon* >     init_muons;
  vector < TRootElectron* > init_electrons;
  vector < TRootJet* >      init_jets;
  vector < TRootMET* >      mets;

  //Global variable
  TRootEvent* event = 0;
  TopFCNC_Evt* MyTopFCNC_EvtCand = 0;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////// Histograms /////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////
  ////////////////// MultiSample plots  //////////////////////////////
  ////////////////////////////////////////////////////////////////////

  MSPlot["RhoCorrection"]              = new MultiSamplePlot(datasets, "RhoCorrection", 100, 0, 100, "#rho");
  MSPlot["NbOfVertices"]               = new MultiSamplePlot(datasets, "NbOfVertices", 20, 0, 20, "Nb. of vertices");

  MSPlot["1stLeadingMuonRelIsolation"] = new MultiSamplePlot(datasets, "1stLeadingMuonRelIsolation", 500, 0, 0.5, "RelIso");
  MSPlot["2ndLeadingMuonRelIsolation"] = new MultiSamplePlot(datasets, "2ndLeadingMuonRelIsolation", 500, 0, 0.5, "RelIso");
  MSPlot["3rdLeadingMuonRelIsolation"] = new MultiSamplePlot(datasets, "3rdLeadingMuonRelIsolation", 500, 0, 0.5, "RelIso");

  MSPlot["NbOfIsolatedMuons"]          = new MultiSamplePlot(datasets, "NbOfIsolatedMuons", 5, 0, 5, "Nb. of isolated muons");
  MSPlot["NbOfIsolatedElectrons"]      = new MultiSamplePlot(datasets, "NbOfIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");

  MSPlot["DiLeptonInvMass"]            = new MultiSamplePlot(datasets, "DiLeptonInvMass", 400, 50, 130, "m_{ll}");

  MSPlot["NbOfExtraIsolatedMuons"]     = new MultiSamplePlot(datasets, "NbOfExtraIsolatedMuons", 5, 0, 5, "Nb. of isolated muons");
  MSPlot["NbOfExtraIsolatedElectrons"] = new MultiSamplePlot(datasets, "NbOfExtraIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");

  MSPlot["NbOfSelectedJets_Before3rdLeptCut"]       = new MultiSamplePlot(datasets, "NbOfSelectedJets_Before3rdLeptCut", 15, 0, 15, "Nb. of jets");
  MSPlot["NbOfSelectedJets_mm_ch"]                  = new MultiSamplePlot(datasets, "NbOfSelectedJets_mm_ch", 15, 0, 15, "Nb. of jets");
  MSPlot["NbOfSelectedJets_mme_ch"]                 = new MultiSamplePlot(datasets, "NbOfSelectedJets_mme_ch", 15, 0, 15, "Nb. of jets");
  MSPlot["NbOfSelectedJets_mmm_ch"]                 = new MultiSamplePlot(datasets, "NbOfSelectedJets_mmm_ch", 15, 0, 15, "Nb. of jets");
  MSPlot["NbOfSelectedJets_ee_ch"]                  = new MultiSamplePlot(datasets, "NbOfSelectedJets_ee_ch", 15, 0, 15, "Nb. of jets");
  MSPlot["NbOfSelectedJets_eee_ch"]                 = new MultiSamplePlot(datasets, "NbOfSelectedJets_eee_ch", 15, 0, 15, "Nb. of jets");
  MSPlot["NbOfSelectedJets_eem_ch"]                 = new MultiSamplePlot(datasets, "NbOfSelectedJets_eem_ch", 15, 0, 15, "Nb. of jets");

  MSPlot["BdiscBJetCand_mm_ch_CVS"]           = new MultiSamplePlot(datasets, "BdiscBJetCand_mm_ch_CVS", 100, 0, 1, "CSV b-disc.");
  MSPlot["BdiscBJetCand_ee_ch_CVS"]           = new MultiSamplePlot(datasets, "BdiscBJetCand_ee_ch_CVS", 100, 0, 1, "CSV b-disc.");
  MSPlot["BdiscBJetCand_mm_ch_TCHE"]          = new MultiSamplePlot(datasets, "BdiscBJetCand_mm_ch_TCHE", 100, 0, 50, "TCHE b-disc.");
  MSPlot["BdiscBJetCand_ee_ch_TCHE"]          = new MultiSamplePlot(datasets, "BdiscBJetCand_ee_ch_TCHE", 100, 0, 50, "TCHE b-disc.");

  MSPlot["BdiscBJetCand_mme_ch_CVS"]          = new MultiSamplePlot(datasets, "BdiscBJetCand_mme_ch_CVS", 100, 0, 1, "CSV b-disc.");
  MSPlot["BdiscBJetCand_eee_ch_CVS"]          = new MultiSamplePlot(datasets, "BdiscBJetCand_eee_ch_CVS", 100, 0, 1, "CSV b-disc.");
  MSPlot["BdiscBJetCand_mme_ch_TCHE"]         = new MultiSamplePlot(datasets, "BdiscBJetCand_mme_ch_TCHE", 100, 0, 50, "TCHE b-disc.");
  MSPlot["BdiscBJetCand_eee_ch_TCHE"]         = new MultiSamplePlot(datasets, "BdiscBJetCand_eee_ch_TCHE", 100, 0, 50, "TCHE b-disc.");

  MSPlot["BdiscBJetCand_mmm_ch_CVS"]          = new MultiSamplePlot(datasets, "BdiscBJetCand_mmm_ch_CVS", 100, 0, 1, "CSV b-disc.");
  MSPlot["BdiscBJetCand_eem_ch_CVS"]          = new MultiSamplePlot(datasets, "BdiscBJetCand_eem_ch_CVS", 100, 0, 1, "CSV b-disc.");
  MSPlot["BdiscBJetCand_mmm_ch_TCHE"]         = new MultiSamplePlot(datasets, "BdiscBJetCand_mmm_ch_TCHE", 100, 0, 50, "TCHE b-disc.");
  MSPlot["BdiscBJetCand_eem_ch_TCHE"]         = new MultiSamplePlot(datasets, "BdiscBJetCand_eem_ch_TCHE", 100, 0, 50, "TCHE b-disc.");

  MSPlot["HighestBdisc_mm_ch_CVS"]            = new MultiSamplePlot(datasets, "HighestBdisc_mm_ch_CVS", 100, 0, 1, "CSV b-disc.");
  MSPlot["HighestBdisc_ee_ch_CVS"]            = new MultiSamplePlot(datasets, "HighestBdisc_ee_ch_CVS", 100, 0, 1, "CSV b-disc.");
  MSPlot["HighestBdisc_mm_ch_TCHE"]           = new MultiSamplePlot(datasets, "HighestBdisc_mm_ch_TCHE",100, 0, 50, "TCHE b-disc.");
  MSPlot["HighestBdisc_ee_ch_TCHE"]           = new MultiSamplePlot(datasets, "HighestBdisc_ee_ch_TCHE",100, 0, 50, "TCHE b-disc.");
  
  MSPlot["MET_mm_ch"]                         = new MultiSamplePlot(datasets, "MET_mm_ch",  100, 0, 200, "MET");
  MSPlot["MET_mme_ch"]                        = new MultiSamplePlot(datasets, "MET_mme_ch", 100, 0, 200, "MET");
  MSPlot["MET_mmm_ch"]                        = new MultiSamplePlot(datasets, "MET_mmm_ch", 100, 0, 200, "MET");
  MSPlot["MET_ee_ch"]                         = new MultiSamplePlot(datasets, "MET_ee_ch",  100, 0, 200, "MET");
  MSPlot["MET_eee_ch"]                        = new MultiSamplePlot(datasets, "MET_eee_ch", 100, 0, 200, "MET");
  MSPlot["MET_eem_ch"]                        = new MultiSamplePlot(datasets, "MET_eem_ch", 100, 0, 200, "MET");

  MSPlot["Mtt_mm_ch"]                         = new MultiSamplePlot(datasets, "Mtt_mm_ch",  100, 0, 1000, "m_{t#bar{t}}");
  MSPlot["Mtt_mme_ch"]                        = new MultiSamplePlot(datasets, "Mtt_mme_ch", 100, 0, 1000, "m_{t#bar{t}}");
  MSPlot["Mtt_mmm_ch"]                        = new MultiSamplePlot(datasets, "Mtt_mmm_ch", 100, 0, 1000, "m_{t#bar{t}}");
  MSPlot["Mtt_ee_ch"]                         = new MultiSamplePlot(datasets, "Mtt_ee_ch",  100, 0, 1000, "m_{t#bar{t}}");
  MSPlot["Mtt_eee_ch"]                        = new MultiSamplePlot(datasets, "Mtt_eee_ch", 100, 0, 1000, "m_{t#bar{t}}");
  MSPlot["Mtt_eem_ch"]                        = new MultiSamplePlot(datasets, "Mtt_eem_ch", 100, 0, 1000, "m_{t#bar{t}}");

  MSPlot["Mzq_mm_ch"]                         = new MultiSamplePlot(datasets, "Mzq_mm_ch",  100, 80, 300, "m_{Zq}");
  MSPlot["Mzq_mme_ch"]                        = new MultiSamplePlot(datasets, "Mzq_mme_ch", 100, 80, 300, "m_{Zq}");
  MSPlot["Mzq_mmm_ch"]                        = new MultiSamplePlot(datasets, "Mzq_mmm_ch", 100, 80, 300, "m_{Zq}");
  MSPlot["Mzq_ee_ch"]                         = new MultiSamplePlot(datasets, "Mzq_ee_ch",  100, 80, 300, "m_{Zq}");
  MSPlot["Mzq_eee_ch"]                        = new MultiSamplePlot(datasets, "Mzq_eee_ch", 100, 80, 300, "m_{Zq}");
  MSPlot["Mzq_eem_ch"]                        = new MultiSamplePlot(datasets, "Mzq_eem_ch", 100, 80, 300, "m_{Zq}");

  MSPlot["MassChi2_mm_ch"]                    = new MultiSamplePlot(datasets, "MassChi2_mm_ch",  20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_mmm_ch"]                   = new MultiSamplePlot(datasets, "MassChi2_mmm_ch", 20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_mme_ch"]                   = new MultiSamplePlot(datasets, "MassChi2_mme_ch", 20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_ee_ch"]                    = new MultiSamplePlot(datasets, "MassChi2_ee_ch",  20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_eem_ch"]                   = new MultiSamplePlot(datasets, "MassChi2_eem_ch", 20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_eee_ch"]                   = new MultiSamplePlot(datasets, "MassChi2_eee_ch", 20, 0, 20, "\\chi^{2}");

  MSPlot["MassChi2_SM_mm_ch"]                 = new MultiSamplePlot(datasets, "MassChi2_SM_mm_ch",  20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_SM_mmm_ch"]                = new MultiSamplePlot(datasets, "MassChi2_SM_mmm_ch", 20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_SM_mme_ch"]                = new MultiSamplePlot(datasets, "MassChi2_SM_mme_ch", 20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_SM_ee_ch"]                 = new MultiSamplePlot(datasets, "MassChi2_SM_ee_ch",  20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_SM_eem_ch"]                = new MultiSamplePlot(datasets, "MassChi2_SM_eem_ch", 20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_SM_eee_ch"]                = new MultiSamplePlot(datasets, "MassChi2_SM_eee_ch", 20, 0, 20, "\\chi^{2}");

  MSPlot["MassChi2_FCNC_mm_ch"]               = new MultiSamplePlot(datasets, "MassChi2_FCNC_mm_ch",  20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_FCNC_mmm_ch"]              = new MultiSamplePlot(datasets, "MassChi2_FCNC_mmm_ch", 20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_FCNC_mme_ch"]              = new MultiSamplePlot(datasets, "MassChi2_FCNC_mme_ch", 20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_FCNC_ee_ch"]               = new MultiSamplePlot(datasets, "MassChi2_FCNC_ee_ch",  20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_FCNC_eem_ch"]              = new MultiSamplePlot(datasets, "MassChi2_FCNC_eem_ch", 20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_FCNC_eee_ch"]              = new MultiSamplePlot(datasets, "MassChi2_FCNC_eee_ch", 20, 0, 20, "\\chi^{2}");

  MSPlot["TriLeptonInvMass_mmm_ch"]           = new MultiSamplePlot(datasets, "TriLeptonInvMass_mmm_ch", 160, 50, 130, "m_{lll}");
  MSPlot["TriLeptonInvMass_mme_ch"]           = new MultiSamplePlot(datasets, "TriLeptonInvMass_mme_ch", 160, 50, 130, "m_{lll}");
  MSPlot["TriLeptonInvMass_eem_ch"]           = new MultiSamplePlot(datasets, "TriLeptonInvMass_eem_ch", 160, 50, 130, "m_{lll}");
  MSPlot["TriLeptonInvMass_eee_ch"]           = new MultiSamplePlot(datasets, "TriLeptonInvMass_eee_ch", 160, 50, 130, "m_{lll}");
//  MSPlot["NbOfLooseMuon"]     = new MultiSamplePlot(datasets, "NbOfLooseMuon", 10, 0, 10, "Nb. of loose muons");
//  MSPlot["NbOfLooseElectron"] = new MultiSamplePlot(datasets, "NbOfLooseElectron", 10, 0, 10, "Nb. of loose electrons");

  ////////////////////////////////////////////////////////////////////
  ////////////////// 1D histograms  //////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  histo1D["lumiWeights"] = new TH1F("lumiWeights","lumiWeights;lumiWeight;#events",100,0,4);
  for (unsigned int d = 0; d < datasets.size(); d++){
	histo2D[("MET_vs_Mzq_mm_ch_"+datasets[d]->Name()).c_str()] = new TH2F(("MET_vs_Mzq_mm_ch_"+datasets[d]->Name()).c_str(),"MET:m_{zq}",100,0,200,100,80,300);
	histo2D[("MET_vs_Mzq_mmm_ch_"+datasets[d]->Name()).c_str()] = new TH2F(("MET_vs_Mzq_mmm_ch_"+datasets[d]->Name()).c_str(),"MET:m_{zq}",100,0,200,100,80,300);
	histo2D[("MET_vs_Mzq_mme_ch_"+datasets[d]->Name()).c_str()] = new TH2F(("MET_vs_Mzq_mme_ch_"+datasets[d]->Name()).c_str(),"MET:m_{zq}",100,0,200,100,80,300);
  }
  cout << " - Declared histograms ..." <<  endl;
	
  ////////////////////////////////////////////////////////////////////
  ////////////////// Plots  //////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  string pathPNG = "TopFCNC"+postfix+channelpostfix;
  pathPNG += "_MSPlots/"; 	
//  pathPNG = pathPNG +"/"; 	
  mkdir(pathPNG.c_str(),0777);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////// Selection Tables ///////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////
  ///////////////////// Channel : µµ  ////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  vector<string> CutsSelecTableDiMu;
  CutsSelecTableDiMu.push_back(string("initial"));
  CutsSelecTableDiMu.push_back(string("PU reweighting"));
  CutsSelecTableDiMu.push_back(string("Trigger"));
  CutsSelecTableDiMu.push_back(string("Good PV"));
  CutsSelecTableDiMu.push_back(string("$\\geq$ 2 isolated muons"));
  CutsSelecTableDiMu.push_back(string("$|m_{ll}-m_Z|<30$ GeV"));
  CutsSelecTableDiMu.push_back(string("Veto on 3rd iso. lept."));
  CutsSelecTableDiMu.push_back(string("$\\geq$ 1 jet"));
  CutsSelecTableDiMu.push_back(string("$\\geq$ 2 jet"));
  CutsSelecTableDiMu.push_back(string("$\\geq$ 3 jet"));
  CutsSelecTableDiMu.push_back(string("$\\geq$ 4 jet"));
  CutsSelecTableDiMu.push_back(string("$b\\texttt{-}disc \\geq 0.7$ (CSV)"));
//  CutsSelecTableDiMu.push_back(string("$MET \\leq 30$ GeV"));  

  SelectionTable selecTableDiMu(CutsSelecTableDiMu, datasets);
  selecTableDiMu.SetLuminosity(Luminosity);
  selecTableDiMu.SetPrecision(1);

  ////////////////////////////////////////////////////////////////////
  ///////////////////// Channel : µµµ ////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  vector<string> CutsSelecTableTriMu;
  CutsSelecTableTriMu.push_back(string("initial"));
  CutsSelecTableTriMu.push_back(string("PU reweighting"));
  CutsSelecTableTriMu.push_back(string("Trigger"));
  CutsSelecTableTriMu.push_back(string("Good PV"));
  CutsSelecTableTriMu.push_back(string("$\\geq$ 2 isolated muons"));
  CutsSelecTableTriMu.push_back(string("$|m_{ll}-m_Z|<30$ GeV"));
  CutsSelecTableTriMu.push_back(string("3rd iso. mu."));
  CutsSelecTableTriMu.push_back(string("$\\geq$ 1 jet"));
  CutsSelecTableTriMu.push_back(string("$\\geq$ 2 jet"));
  CutsSelecTableTriMu.push_back(string("$\\geq$ 3 jet"));
//  CutsSelecTableTriMu.push_back(string("$MET \\leq 30$ GeV"));  

  SelectionTable selecTableTriMu(CutsSelecTableTriMu, datasets);
  selecTableTriMu.SetLuminosity(Luminosity);
  selecTableTriMu.SetPrecision(1);

  ////////////////////////////////////////////////////////////////////
  ///////////////////// Channel : µµe ////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  vector<string> CutsselecTableDiMuEl;
  CutsselecTableDiMuEl.push_back(string("initial"));
  CutsselecTableDiMuEl.push_back(string("PU reweighting"));
  CutsselecTableDiMuEl.push_back(string("Trigger"));
  CutsselecTableDiMuEl.push_back(string("Good PV"));
  CutsselecTableDiMuEl.push_back(string("$\\geq$ 2 isolated muons"));
  CutsselecTableDiMuEl.push_back(string("$|m_{ll}-m_Z|<30$ GeV"));
  CutsselecTableDiMuEl.push_back(string("3rd iso. elec."));
  CutsselecTableDiMuEl.push_back(string("$\\geq$ 1 jet"));
  CutsselecTableDiMuEl.push_back(string("$\\geq$ 2 jet"));
  CutsselecTableDiMuEl.push_back(string("$\\geq$ 3 jet"));
//  CutsselecTableDiMuEl.push_back(string("$MET \\leq 30$ GeV"));  

  SelectionTable selecTableDiMuEl(CutsselecTableDiMuEl, datasets);
  selecTableDiMuEl.SetLuminosity(Luminosity);
  selecTableDiMuEl.SetPrecision(1);

  ////////////////////////////////////////////////////////////////////
  ///////////////////// Channel : ee /////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  vector<string> CutsSelecTableDiEl;
  CutsSelecTableDiEl.push_back(string("initial"));
  CutsSelecTableDiEl.push_back(string("PU reweighting"));
  CutsSelecTableDiEl.push_back(string("Trigger"));
  CutsSelecTableDiEl.push_back(string("Good PV"));
  CutsSelecTableDiEl.push_back(string("$\\geq$ 2 isolated electrons"));
  CutsSelecTableDiEl.push_back(string("$|m_{ll}-m_Z|<30$ GeV"));
  CutsSelecTableDiEl.push_back(string("Veto on 3rd iso. lept."));
  //CutsSelecTableDiEl.push_back(string("Conversion veto"));
  CutsSelecTableDiEl.push_back(string("$\\geq$ 1 jet"));
  CutsSelecTableDiEl.push_back(string("$\\geq$ 2 jet"));
  CutsSelecTableDiEl.push_back(string("$\\geq$ 3 jet"));
  CutsSelecTableDiEl.push_back(string("$\\geq$ 4 jet"));
  CutsSelecTableDiEl.push_back(string("$b\\texttt{-}disc \\geq 0.7$ (CSV)"));
  //CutsSelecTableDiEl.push_back(string("MET > 40 GeV"));


  SelectionTable selecTableDiEl(CutsSelecTableDiEl, datasets);
  selecTableDiEl.SetLuminosity(Luminosity);
  selecTableDiEl.SetPrecision(1);
  
  ////////////////////////////////////////////////////////////////////
  ///////////////////// Channel : eee ////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  vector<string> CutsSelecTableTriEl;
  CutsSelecTableTriEl.push_back(string("initial"));
  CutsSelecTableTriEl.push_back(string("PU reweighting"));
  CutsSelecTableTriEl.push_back(string("Trigger"));
  CutsSelecTableTriEl.push_back(string("Good PV"));
  CutsSelecTableTriEl.push_back(string("$\\geq$ 2 isolated electrons"));
  CutsSelecTableTriEl.push_back(string("$|m_{ll}-m_Z|<30$ GeV"));
  CutsSelecTableTriEl.push_back(string("3rd iso. elec."));
  CutsSelecTableTriEl.push_back(string("$\\geq$ 1 jet"));
  CutsSelecTableTriEl.push_back(string("$\\geq$ 2 jet"));
  CutsSelecTableTriEl.push_back(string("$\\geq$ 3 jet"));
//  CutsSelecTableTriEl.push_back(string("$MET \\leq 30$ GeV"));  

  SelectionTable selecTableTriEl(CutsSelecTableTriEl, datasets);
  selecTableTriEl.SetLuminosity(Luminosity);
  selecTableTriEl.SetPrecision(1);
  

  ////////////////////////////////////////////////////////////////////
  ///////////////////// Channel : eeµ ////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  vector<string> CutsSelecTableDiElMu;
  CutsSelecTableDiElMu.push_back(string("initial"));
  CutsSelecTableDiElMu.push_back(string("PU reweighting"));
  CutsSelecTableDiElMu.push_back(string("Trigger"));
  CutsSelecTableDiElMu.push_back(string("Good PV"));
  CutsSelecTableDiElMu.push_back(string("$\\geq$ 2 isolated electrons"));
  CutsSelecTableDiElMu.push_back(string("$|m_{ll}-m_Z|<30$ GeV"));
  CutsSelecTableDiElMu.push_back(string("3rd iso. muon"));
  CutsSelecTableDiElMu.push_back(string("$\\geq$ 1 jet"));
  CutsSelecTableDiElMu.push_back(string("$\\geq$ 2 jet"));
  CutsSelecTableDiElMu.push_back(string("$\\geq$ 3 jet"));
  //CutsSelecTableDiElMu.push_back(string("Conversion veto"));


  SelectionTable selecTableDiElMu(CutsSelecTableDiElMu, datasets);
  selecTableDiElMu.SetLuminosity(Luminosity);
  selecTableDiElMu.SetPrecision(1);

  cout << " - SelectionTable instantiated ..." << endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////// PileUp Reweighting - 3D ///////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  Lumi3DReWeighting Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Fall11.root","PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root", "pileup", "pileup");
  Lumi3DWeights.weight3D_init(1.0);

  if(doPUShift == 1)
  	Lumi3DWeights.weight3D_init(0.92);
  else if(doPUShift == 2)
  	Lumi3DWeights.weight3D_init(1.08);
  else
  	Lumi3DWeights.weight3D_init(1.0);

  cout << " - Initialized LumiReWeighting stuff" << endl;
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////// Loop on datasets ///////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;

  for (unsigned int d = 0; d < datasets.size(); d++) //d < datasets.size()
  {
    if (verbose > 1){
      cout << "   Dataset " << d << " name : " << datasets[d]->Name () << " / title : " << datasets[d]->Title () << endl;
      cout << " - Cross section = " << datasets[d]->Xsection() << endl;
      cout << " - IntLumi = " << datasets[d]->EquivalentLumi() << "  NormFactor = " << datasets[d]->NormFactor() << endl;
      cout << " - Nb of events : " << datasets[d]->NofEvtsToRunOver() << endl;
    }
    //open files and load
    cout<<"Load Dataset"<<endl;
    treeLoader.LoadDataset (datasets[d], anaEnv);
		
    string previousFilename = "";
    int iFile = -1;
    
    string dataSetName = datasets[d]->Name();	
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////// Initialize JEC factors /////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
    vector<JetCorrectorParameters> vCorrParam;

    // Create the JetCorrectorParameter objects, the order does not matter.
    // YYYY is the first part of the txt files: usually the global tag from which they are retrieved
    JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters("JECFiles/START42_V17_AK5PFchs_L1FastJet.txt");
    JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters("JECFiles/START42_V17_AK5PFchs_L2Relative.txt");
    JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters("JECFiles/START42_V17_AK5PFchs_L3Absolute.txt");

    //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
    vCorrParam.push_back(*L1JetPar);
    vCorrParam.push_back(*L2JetPar);
    vCorrParam.push_back(*L3JetPar);

    if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") // Data!
      {
	JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("JECFiles/START42_V17_AK5PFchs_L2L3Residual.txt");
	vCorrParam.push_back(*ResJetCorPar);
      }
    
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/START42_V17_AK5PFchs_Uncertainty.txt");
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true); // last boolean ('startFromRaw') = false!    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////// Loop on events //////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int itrigger = -1, previousRun = -1;
    int fourIsoLeptCounter = 0;
      
    int start = 0;
    unsigned int end = datasets[d]->NofEvtsToRunOver();
     
    if (verbose > 1) cout << " - Loop over events " << endl;      
    
    for (unsigned int ievt = start; ievt < end; ievt++)
    {        

	if(ievt%1000 == 0)
		std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC << " ("<<100*(ievt-start)/(end-start)<<"%)"<<flush<<"\r";

	//load event
	event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);

	vector<TRootGenJet*> genjets;
	if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
	{
		genjets = treeLoader.LoadGenJet(ievt);
	}

	// check which file in the dataset it is to have the HLTInfo right
	string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
	if(previousFilename != currentFilename)
	{
		previousFilename = currentFilename;
        	iFile++;
		cout<<"File changed!!! => iFile = "<<iFile<<endl;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////// trigger /////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool trigged = false;
	int currentRun = event->runId();
	if(previousRun != currentRun)
	{
      		previousRun = currentRun;
		if(diMuon)
		{
			if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
			{
				/*------------------------------------------------------------------
				Dataset : DoubleMu/Run2011A-May10ReReco-v1
				--------------------------------------------------------------------
				Trigger HLT_DoubleMu7_v1 available for runs 160431-163261
				Trigger HLT_DoubleMu7_v2 available for runs 163270-163869
				------------------------------------------------------------------*/
				if(event->runId() >= 160431 && event->runId() <= 163261) // IntLumi = 33.601(/pb)
					itrigger = treeLoader.iTrigger (string ("HLT_DoubleMu7_v1"), currentRun, iFile);
  			else if (event->runId() >= 163270 && event->runId() <= 163869) // IntLumi = 167.434(/pb)
    			itrigger = treeLoader.iTrigger (string ("HLT_DoubleMu7_v2"), currentRun, iFile);
				/*--------------------------------------------------------------------
				Sub-Total integrated luminosity = 201,035(/pb)
				    Total integrated luminosity = 201,035(/pb)
				------------------------------------------------------------------*/

				/*------------------------------------------------------------------
				Dataset : DoubleMu/Run2011A-PromptReco-v4
				--------------------------------------------------------------------
				Trigger HLT_Mu13_Mu8_v2 available for runs 165088-167043
				Trigger HLT_Mu13_Mu8_v3 available for runs 166346-166346
				Trigger HLT_Mu13_Mu8_v4 available for runs 167078-167913
				------------------------------------------------------------------*/
  		  else if (event->runId() >= 165088 && event->runId() <= 167043 && event->runId() != 166346) // IntLumi = 619,142(/pb)
          itrigger = treeLoader.iTrigger (string ("HLT_Mu13_Mu8_v2"), currentRun, iFile);
  		  else if (event->runId() == 166346) // IntLumi = 4.291(/pb)
          itrigger = treeLoader.iTrigger (string ("HLT_Mu13_Mu8_v3"), currentRun, iFile);
  		  else if (event->runId() >= 167078 && event->runId() <= 167913) // IntLumi = 226.454(/pb)
          itrigger = treeLoader.iTrigger (string ("HLT_Mu13_Mu8_v4"), currentRun, iFile);
				/*--------------------------------------------------------------------
				Sub-Total integrated luminosity =  849,887(/pb)
				    Total integrated luminosity = 1050,922(/pb)
				------------------------------------------------------------------*/

				/*------------------------------------------------------------------
				Dataset : DoubleMu/Run2011A-05Aug2011-v1
				--------------------------------------------------------------------
				Trigger HLT_Mu13_Mu8_v6 available for runs 170826-172619
				------------------------------------------------------------------*/
				else if (event->runId() >= 170249 && event->runId() <= 172619)  // IntLumi = //Aug05ReReco equivalent to PromptReco_v5
				  itrigger = treeLoader.iTrigger (string ("HLT_Mu13_Mu8_v6"), currentRun, iFile);
				/*------------------------------------------------------------------
				Sub-Total integrated luminosity =  XXXX(/pb)
				    Total integrated luminosity = XXXX(/pb)
				------------------------------------------------------------------*/

				/*------------------------------------------------------------------
				Dataset : DoubleMu/Run2011A-PromptReco-v6
				--------------------------------------------------------------------
				Trigger HLT_Mu13_Mu8_v6 available for runs 172620-173198
				Trigger HLT_Mu13_Mu8_v7 available for runs 173236-173692
				------------------------------------------------------------------*/
				else if (event->runId() >= 172620 && event->runId() <= 173198) // IntLumi = 
          itrigger = treeLoader.iTrigger (string ("HLT_Mu13_Mu8_v6"), currentRun, iFile);
				else if (event->runId() >= 173236 && event->runId() <= 173692) // IntLumi = 
				  itrigger = treeLoader.iTrigger (string ("HLT_Mu13_Mu8_v7"), currentRun, iFile);
				/*------------------------------------------------------------------
				Sub-Total integrated luminosity =  XXXX(/pb)
				    Total integrated luminosity = XXXX(/pb)
				------------------------------------------------------------------*/

				/*------------------------------------------------------------------
				Dataset : DoubleMu/Run2011B-PromptReco-v1
				--------------------------------------------------------------------
				Trigger HLT_Mu13_Mu8_v7  available for runs 175860-178380
				Trigger HLT_Mu13_Mu8_v10 available for runs 178420-179889
				Trigger HLT_Mu13_Mu8_v11 available for runs 179959-180252
				------------------------------------------------------------------*/
        else if( event->runId() >=  175860 && event->runId() <= 178380 ) // IntLumi = 
   			  itrigger = treeLoader.iTrigger (string ("HLT_Mu13_Mu8_v7"), currentRun, iFile);
        else if( event->runId() >=  178420 && event->runId() <= 179889 ) // IntLumi = 
					itrigger = treeLoader.iTrigger (string ("HLT_Mu13_Mu8_v10"),currentRun, iFile);
        else if( event->runId() >=  179959 && event->runId() <=  180252 ) // IntLumi = 
					itrigger = treeLoader.iTrigger (string ("HLT_Mu13_Mu8_v11"),currentRun, iFile); 
									   
  		  if(itrigger == 9999)
				{
    		  cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (DATA) IN RUN " << event->runId() << endl;
    		  exit(1);
 	 			}
				//trigged = treeLoader.EventTrigged (itrigger);			

	   		}
	   		else 
	   		{
				if(dataSetName != "ttbar_fcnc") itrigger = treeLoader.iTrigger (string ("HLT_Mu13_Mu8_v7"), currentRun, iFile);
				else itrigger = treeLoader.iTrigger (string ("HLT_DoubleMu7_v1"), currentRun, iFile);
    
  				if(itrigger == 9999)
				{
    			  		cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (" << dataSetName << ") IN RUN " << event->runId() << endl;
    			  		exit(1);
				}
				cout<<"Trigger bit nr : "<<itrigger<<endl;
			}
		} //end if diMuon
		else if(diElectron)
		{
			if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
			{
				/*------------------------------------------------------------------
				Dataset : DoubleElectron/Run2011A-May10ReReco-v1
				--------------------------------------------------------------------
				Trigger HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1 available for runs 160431-161016
				Trigger HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2 available for runs 162762-163261
				Trigger HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3 available for runs 163270-163869
				Trigger HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2
				->	available for runs 162762-163261
				Trigger HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v3
				->	available for runs 163270-163869
				------------------------------------------------------------------*/
				if(currentRun >= 150000 && currentRun <= 161176)      // IntLumi = 5.281(/pb)
					itrigger = treeLoader.iTrigger ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1", currentRun);
				else if(currentRun >= 161179 && currentRun <= 163261) // IntLumi = 28.321(/pb)
					itrigger = treeLoader.iTrigger ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2", currentRun);
				else if(currentRun >= 163262 && currentRun <= 164237) // IntLumi = 167.518(/pb)
					itrigger = treeLoader.iTrigger ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3", currentRun);
				/*--------------------------------------------------------------------
				Sub-Total integrated luminosity = 201,12(/pb)
				    Total integrated luminosity = 201,12(/pb)
				------------------------------------------------------------------*/

				/*--------------------------------------------------------------------
				Dataset : DoubleElectron/Run2011A-PromptReco-v4
				--------------------------------------------------------------------
				Trigger HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4 available for runs 165088-165633
				Trigger HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5 available for runs 165970-166967
				Trigger HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6 available for runs 167039-167913
				Trigger HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v4
				->	available for runs 165088-165633
				Trigger HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v5
				->	available for runs 165970-166967
				Trigger HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5
				->	available for runs 167039-167913
				--------------------------------------------------------------------*/
				else if(currentRun >= 165085 && currentRun <= 165888) // IntLumi = 139.027(/pb)
					itrigger = treeLoader.iTrigger ("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v4", currentRun);
				else if(currentRun >= 165900 && currentRun <= 167043) // IntLumi = 524.904(/pb)
					itrigger = treeLoader.iTrigger ("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v5", currentRun);
				else if(currentRun >= 167044 && currentRun <= 170053) // IntLumi = 243.081(/pb)
					itrigger = treeLoader.iTrigger ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5", currentRun);
				/*--------------------------------------------------------------------
				Sub-Total integrated luminosity =  907,012(/pb)
				    Total integrated luminosity = 1108,132(/pb)
				------------------------------------------------------------------*/

				/*--------------------------------------------------------------------
				Dataset : DoubleElectron/Run2011A-05Aug2011-v1
				--------------------------------------------------------------------
				Trigger HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v8 available for runs 170826-172619
				Trigger HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7 
				->	available for runs 170826-172619
				------------------------------------------------------------------*/
				else if(currentRun >= 170054 && currentRun <= 170759) // IntLumi = 
					itrigger = treeLoader.iTrigger ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6", currentRun);
				else if(currentRun >= 170760 && currentRun <= 172619) // IntLumi = 
					itrigger = treeLoader.iTrigger ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7", currentRun);
				/*--------------------------------------------------------------------
				Sub-Total integrated luminosity = 907,012(/pb)
				    Total integrated luminosity = 1108,132(/pb)
				------------------------------------------------------------------*/

				/*--------------------------------------------------------------------
				Dataset : DoubleElectron/Run2011A-PromptReco-v6
				--------------------------------------------------------------------
				Trigger HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v8 available for runs 172620-173198
				Trigger HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7
				->	available for runs 172620-173198
				Trigger HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8
				->	available for runs 173236-173692
				------------------------------------------------------------------*/
				else if(currentRun >= 172620 && currentRun <= 173198) // IntLumi = 410.523(/pb)
					itrigger = treeLoader.iTrigger ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7", currentRun);
				else if(currentRun >= 173199 && currentRun <= 173692) // IntLumi = 
				itrigger = treeLoader.iTrigger ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8", currentRun);
				/*--------------------------------------------------------------------
				Sub-Total integrated luminosity = 907,012(/pb)
				    Total integrated luminosity = 1108,132(/pb)
				------------------------------------------------------------------*/

				/*--------------------------------------------------------------------
				Dataset : DoubleElectron/Run2011B-PromptReco-v1
				--------------------------------------------------------------------
				Trigger HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8
				->	available for runs 175860-178380
				Trigger HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9
				->	available for runs 178420-179889
				Trigger HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10
				->	available for runs 179959-180252
				------------------------------------------------------------------*/

				else if(currentRun >= 175860 && currentRun <= 178380) // IntLumi = 
				itrigger = treeLoader.iTrigger ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8", currentRun);
				else if(currentRun >= 178420 && currentRun <= 179899) // IntLumi = 
				itrigger = treeLoader.iTrigger ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9", currentRun);
				else if(currentRun >= 179900 && currentRun <= 999999) // IntLumi = 
				itrigger = treeLoader.iTrigger ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10", currentRun);

			}
	   		else 
	   		{
				if(dataSetName != "ttbar_fcnc") itrigger = treeLoader.iTrigger (string ("HLT_"), currentRun, iFile);
				else itrigger = treeLoader.iTrigger (string ("HLT_"), currentRun, iFile);
    
  				if(itrigger == 9999)
				{
    			  		cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (" << dataSetName << ") IN RUN " << event->runId() << endl;
    			  		exit(1);
				}
				cout<<"Trigger bit nr : "<<itrigger<<endl;
			}
		} //end if diElectron
	} //end previousRun != currentRun

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////// Jet energy scale corrections     /////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Apply Jet Corrections on-the-fly
	//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JES correction on the fly:");
	if( dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
		jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho(),true); //last boolean: isData (needed for L2L3Residual...)
	else
		jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho(),false); //last boolean: isData (needed for L2L3Residual...)
	//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JES correction on the fly:");

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////// Type I MET corrections: (Only for |eta| <=4.7 ) //////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before MET type I correction:");      
	if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
		jetTools->correctMETTypeOne(init_jets,mets[0],true);
	else
		jetTools->correctMETTypeOne(init_jets,mets[0],false);
	//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After MET type I correction:");

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////// Jet energy smearing and systematic uncertainty ///////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
	{
		// JER systematic! 
		if(doJERShift == 1)
			jetTools->correctJetJER(init_jets, genjets, mets[0], "minus");
		else if(doJERShift == 2)
			jetTools->correctJetJER(init_jets, genjets, mets[0], "plus");
		else
			jetTools->correctJetJER(init_jets, genjets, mets[0], "nominal");
	  
		//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JER correction:");	       

		// JES systematic! 
		if (doJESShift == 1)
			jetTools->correctJetJESUnc(init_jets, mets[0], "minus");
		else if (doJESShift == 2)
			jetTools->correctJetJESUnc(init_jets, mets[0], "plus");
	
		//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JES correction:");
	
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////// Beam scrapping veto and PU reweighting ///////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// scale factor for the event
	float scaleFactor = 1.;

	if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
	{
		// Apply the scraping veto. (Is it still needed?)
        	bool isBeamBG = true;
        	if(event->nTracks() > 10)
        	{
			if( ( (float) event->nHighPurityTracks() ) / ( (float) event->nTracks() ) > 0.25 )
			isBeamBG = false;
		}
      		if(isBeamBG) continue;
	}
	else{
		// Apply pile-up reweighting
		double lumiWeight3D = Lumi3DWeights.weight3D(event->nPu(-1),event->nPu(0),event->nPu(+1));
	 	scaleFactor *= lumiWeight3D;
	}
	histo1D["lumiWeights"]->Fill(scaleFactor);	
	MSPlot["RhoCorrection"]->Fill(event->kt6PFJetsPF2PAT_rho(), datasets[d], true, Luminosity*scaleFactor);
			
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////// Event selection ////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	selecTableDiMu.Fill(d,0, 1.);//datasets[d]->Xsection() * datasets[d]->EquivalentLumi() );
	selecTableTriMu.Fill(d,0, 1.);
	selecTableDiMuEl.Fill(d,0, 1.);
	selecTableDiEl.Fill(d,0, 1.);
	selecTableTriEl.Fill(d,0, 1.);
	selecTableDiElMu.Fill(d,0, 1.);

	selecTableDiMu.Fill(d,1,scaleFactor);
	selecTableTriMu.Fill(d,1,scaleFactor);
	selecTableDiMuEl.Fill(d,1,scaleFactor);
	selecTableDiEl.Fill(d,1,scaleFactor);
	selecTableTriEl.Fill(d,1,scaleFactor);
	selecTableDiElMu.Fill(d,1,scaleFactor);
		
	// Apply trigger selection
	trigged = treeLoader.EventTrigged (itrigger);

	if(!trigged)		   continue;
	selecTableDiMu.Fill(d,2,scaleFactor);
	selecTableTriMu.Fill(d,2,scaleFactor);
	selecTableDiMuEl.Fill(d,2,scaleFactor);
	selecTableDiEl.Fill(d,2,scaleFactor);
	selecTableTriEl.Fill(d,2,scaleFactor);
	selecTableDiElMu.Fill(d,2,scaleFactor);

	// Declare selection instance    
	Selection selection(init_jets, init_muons, init_electrons, mets); //mets can also be corrected...
/*            
	vector<TRootMCParticle*> mcParticles;
	if(dataSetName.find("ttbar_fcnc") == 0)
	{
		treeLoader.LoadMCEvent(ievt, 0, 0, mcParticles,false);  
		sort(mcParticles.begin(),mcParticles.end(),HighestPt());
	}
*/
	// Define object selection cuts
	selection.setJetCuts(20.,2.4,0.01,1.,0.98,0.3,0.1);

	selection.setDiElectronCuts(20,2.5,0.17,0.02,1); //Et,Eta,RelIso,d0,DistVzPVz
	//selection.setLooseElectronCuts(15,2.5,0.2);

	selection.setDiMuonCuts(20,2.4,0.125,10,0.02); //Et,Eta,RelIso,NValidHits,d0
	//selection.setLooseMuonCuts(15,2.4,0.2);
	  
	//Select objects 
	vector<TRootElectron*> selectedElectrons_NoIso = selection.GetSelectedDiElectrons(20,2.4,999.,vertex[0]);
	vector<TRootElectron*> selectedElectrons       = selection.GetSelectedDiElectrons(vertex[0]);
	vector<TRootElectron*> selectedExtraElectrons;

	vector<TRootMuon*>     selectedMuons_NoIso = selection.GetSelectedDiMuons(20,2.4,999.);
	vector<TRootMuon*>     selectedMuons       = selection.GetSelectedDiMuons();
	vector<TRootMuon*>     selectedExtraMuons;

	vector<TRootJet*>      selectedJets        = selection.GetSelectedJets(true); // ApplyJetId

	//vector<TRootMuon*>     looseMuons     = selection.GetSelectedLooseMuons();
	//vector<TRootElectron*> looseElectrons = selection.GetSelectedLooseElectrons(true); // VBTF Id

	MSPlot["NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);

	// Apply primary vertex selection
	bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2);
  if(!isGoodPV) continue;
	selecTableDiMu.Fill(d,3,scaleFactor);
	selecTableTriMu.Fill(d,3,scaleFactor);
	selecTableDiMuEl.Fill(d,3,scaleFactor);
	selecTableDiEl.Fill(d,3,scaleFactor);
	selecTableTriEl.Fill(d,3,scaleFactor);
	selecTableDiElMu.Fill(d,3,scaleFactor);

	// Select events with at least two leptons
	if(diMuon){
		if(selectedMuons_NoIso.size()<2) continue;
	}
	else if(diElectron){
		if(selectedElectrons_NoIso.size()<2) continue;
	}

	if(diMuon){
		MSPlot["1stLeadingMuonRelIsolation"]->Fill(selectedMuons_NoIso[0]->relativePfIso(), datasets[d], true, Luminosity*scaleFactor);
		MSPlot["2ndLeadingMuonRelIsolation"]->Fill(selectedMuons_NoIso[1]->relativePfIso(), datasets[d], true, Luminosity*scaleFactor);
		if(selectedMuons_NoIso.size()>2)
			MSPlot["3rdLeadingMuonRelIsolation"]->Fill(selectedMuons_NoIso[2]->relativePfIso(), datasets[d], true, Luminosity*scaleFactor);
	}
	else if(diElectron){
		MSPlot["1stLeadingElectronRelIsolation"]->Fill(selectedElectrons_NoIso[0]->relativePfIso(), datasets[d], true, Luminosity*scaleFactor);
		MSPlot["2ndLeadingElectronRelIsolation"]->Fill(selectedElectrons_NoIso[1]->relativePfIso(), datasets[d], true, Luminosity*scaleFactor);
		if(selectedElectrons_NoIso.size()>2)
			MSPlot["3rdLeadingElectronRelIsolation"]->Fill(selectedElectrons_NoIso[2]->relativePfIso(), datasets[d], true, Luminosity*scaleFactor);
	}

	MSPlot["NbOfIsolatedMuons"]->Fill(selectedMuons.size(), datasets[d], true, Luminosity*scaleFactor);
	MSPlot["NbOfIsolatedElectrons"]->Fill(selectedElectrons.size(), datasets[d], true, Luminosity*scaleFactor);

	// Select events with at least two isolated muons
	if(diMuon){
		if(selectedMuons.size()<2) continue;
	}
	else if(diElectron){
		if(selectedElectrons.size()<2) continue;
	}

	selecTableDiMu.Fill(d,4,scaleFactor);
	selecTableTriMu.Fill(d,4, scaleFactor);
	selecTableDiMuEl.Fill(d,4, scaleFactor);
	selecTableDiEl.Fill(d,4,scaleFactor);
	selecTableTriEl.Fill(d,4,scaleFactor);
	selecTableDiElMu.Fill(d,4,scaleFactor);

  bool foundZ = false;
	int idx_Z_1 = -1, idx_Z_2 = -1;
	float Zmass = 91.;
	float Zwindowsize = 30.;
	// Calculate the invariant mass for each isolated lepton pairs
	// - return true if the mass is the Z boson mass window 
	// - return the indices of the lepton candidates
	if(diMuon){
    for(unsigned int i=0;i<selectedMuons.size()-1;i++)
    {
      for(unsigned int j=i+1;j<selectedMuons.size();j++)
      {
        TRootMuon* mu1 = (TRootMuon*) selectedMuons[i];
        TRootMuon* mu2 = (TRootMuon*) selectedMuons[j];
        if(mu1->charge() == mu2->charge()) continue;
        double invMass = (*mu1 + *mu2).M();
        MSPlot["DiLeptonInvMass"]->Fill(invMass, datasets[d], true, Luminosity*scaleFactor);
        if( invMass >= (Zmass-Zwindowsize) && invMass <= (Zmass+Zwindowsize) )
        {
          idx_Z_1 = i;
          idx_Z_2 = j;
          foundZ = true;
        }
      }
    }
	}
	else if(diElectron){
	  	for(unsigned int i=0;i<selectedElectrons.size()-1;i++)
	  	{
	  		for(unsigned int j=i+1;j<selectedElectrons.size();j++)
	  		{
	   			TRootElectron* el1 = (TRootElectron*) selectedElectrons[i];
	   			TRootElectron* el2 = (TRootElectron*) selectedElectrons[j];
				if(el1->charge() == el2->charge()) continue;
				double invMass = (*el1 + *el2).M();
				MSPlot["DiLeptonInvMass"]->Fill(invMass, datasets[d], true, Luminosity*scaleFactor);
				if( invMass >= (Zmass-Zwindowsize) && invMass <= (Zmass+Zwindowsize) )
				{
					idx_Z_1 = i;
					idx_Z_2 = j;
					foundZ = true;
				}
	    		}
	  	}
	}
	// Select events with at least one pair of opposite charge leptons with |mll-mz|<windowsize
	if(!foundZ) continue; 
	selecTableDiMu.Fill(d,5,scaleFactor);
	selecTableTriMu.Fill(d,5, scaleFactor);
	selecTableDiMuEl.Fill(d,5, scaleFactor);
	selecTableDiEl.Fill(d,5,scaleFactor);
	selecTableTriEl.Fill(d,5,scaleFactor);
	selecTableDiElMu.Fill(d,5,scaleFactor);

	selectedExtraMuons = selectedMuons;
	selectedExtraElectrons = selectedElectrons;
	// Erase Z boson lepton candidates
	if(diMuon){
		selectedExtraMuons.erase(selectedExtraMuons.begin()+idx_Z_2);
		selectedExtraMuons.erase(selectedExtraMuons.begin()+idx_Z_1);
	}
	else if(diElectron){
		selectedExtraElectrons.erase(selectedExtraElectrons.begin()+idx_Z_2);
		selectedExtraElectrons.erase(selectedExtraElectrons.begin()+idx_Z_1);

	}

	MSPlot["NbOfExtraIsolatedMuons"]->Fill(selectedExtraMuons.size(), datasets[d], true, Luminosity*scaleFactor);
	MSPlot["NbOfExtraIsolatedElectrons"]->Fill(selectedExtraElectrons.size(), datasets[d], true, Luminosity*scaleFactor);

	MSPlot["NbOfSelectedJets_Before3rdLeptCut"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);

	MyTopFCNC_EvtCand = 0;
	double invMass = 0;
	double highestbtagdisc = 0;
	// Select events based on the presence of *exactly one* extra isolated lepton
	if(diMuon && selectedExtraMuons.size()==0 && selectedElectrons.size()==0){
		selecTableDiMu.Fill(d,6,scaleFactor);
		MSPlot["NbOfSelectedJets_mm_ch"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
		if(selectedJets.size()>0){ //at least 1 jet
			selecTableDiMu.Fill(d,7,scaleFactor); 
			if(selectedJets.size()>1){ //at least 2 jets
				selecTableDiMu.Fill(d,8,scaleFactor);
				if(selectedJets.size()>2){ //at least 3 jets
					selecTableDiMu.Fill(d,9,scaleFactor);
					if(selectedJets.size()>3){ //at least 4 jets
						selecTableDiMu.Fill(d,10,scaleFactor);
						sort(selectedJets.begin(),selectedJets.end(),HighestCVSBtag());
						MSPlot["HighestBdisc_mm_ch_CVS"]->Fill(selectedJets[0]->btag_combinedSecondaryVertexBJetTags(),datasets[d], true, Luminosity*scaleFactor);
						highestbtagdisc = selectedJets[0]->btag_combinedSecondaryVertexBJetTags();
						sort(selectedJets.begin(),selectedJets.end(),HighestTCHEBtag());
						MSPlot["HighestBdisc_mm_ch_TCHE"]->Fill(selectedJets[0]->btag_trackCountingHighEffBJetTags(),datasets[d], true, Luminosity*scaleFactor);
						MyTopFCNC_EvtCand = new TopFCNC_Evt(TopFCNC_Evt::kMuon);
						MyTopFCNC_EvtCand->ReconstructDiLeptEvt(selectedMuons[idx_Z_1], selectedMuons[idx_Z_2], selectedJets, true);
						MSPlot["BdiscBJetCand_mm_ch_CVS"]->Fill(MyTopFCNC_EvtCand->B().btag_combinedSecondaryVertexBJetTags(),datasets[d],true, Luminosity*scaleFactor);
						MSPlot["BdiscBJetCand_mm_ch_TCHE"]->Fill(MyTopFCNC_EvtCand->B().btag_trackCountingHighEffBJetTags(),datasets[d],true, Luminosity*scaleFactor);
						if(highestbtagdisc<0.7) continue;
						selecTableDiMu.Fill(d,11,scaleFactor);
						MSPlot["MET_mm_ch"]->Fill(mets[0]->Et(),datasets[d], true, Luminosity*scaleFactor);
						MSPlot["Mtt_mm_ch"]->Fill((MyTopFCNC_EvtCand->smDecayTop()+MyTopFCNC_EvtCand->fcncDecayTop()).M(),datasets[d],true,Luminosity*scaleFactor);
						MSPlot["Mzq_mm_ch"]->Fill(MyTopFCNC_EvtCand->fcncDecayTop().M(),datasets[d], true, Luminosity*scaleFactor);
						MSPlot["MassChi2_mm_ch"]->Fill(MyTopFCNC_EvtCand->MassChi2(),datasets[d], true, Luminosity*scaleFactor);
						MSPlot["MassChi2_SM_mm_ch"]->Fill(MyTopFCNC_EvtCand->MassChi2_SM(),datasets[d], true, Luminosity*scaleFactor);
						MSPlot["MassChi2_FCNC_mm_ch"]->Fill(MyTopFCNC_EvtCand->MassChi2_FCNC(),datasets[d], true, Luminosity*scaleFactor);
						histo2D[("MET_vs_Mzq_mm_ch_"+datasets[d]->Name()).c_str()]->Fill(mets[0]->Et(),MyTopFCNC_EvtCand->fcncDecayTop().M());
					}
				}
			}
		}
	}
	else if(diMuon && selectedExtraMuons.size()==0 && selectedElectrons.size()==1){
		selecTableDiMuEl.Fill(d,6,scaleFactor);
		MSPlot["NbOfSelectedJets_mme_ch"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
		if(selectedJets.size()>0){ //at least 1 jet
			selecTableDiMuEl.Fill(d,7,scaleFactor); 
			if(selectedJets.size()>1){ //at least 2 jets
				selecTableDiMuEl.Fill(d,8,scaleFactor);
				if(selectedJets.size()>2){ //at least 3 jets
					selecTableDiMuEl.Fill(d,9,scaleFactor);
					MyTopFCNC_EvtCand = new TopFCNC_Evt(TopFCNC_Evt::kMuon,TopFCNC_Evt::kElec);
					MyTopFCNC_EvtCand->ReconstructTriLeptEvt(selectedMuons[idx_Z_1], selectedMuons[idx_Z_2], selectedElectrons[0], selectedJets, mets[0], true);
					invMass = (*selectedMuons[idx_Z_1]+*selectedMuons[idx_Z_2]+*selectedElectrons[0]).M();
					MSPlot["BdiscBJetCand_mme_ch_CVS"]->Fill(MyTopFCNC_EvtCand->B().btag_combinedSecondaryVertexBJetTags(), datasets[d], true, Luminosity*scaleFactor);
					MSPlot["BdiscBJetCand_mme_ch_TCHE"]->Fill(MyTopFCNC_EvtCand->B().btag_trackCountingHighEffBJetTags(), datasets[d], true, Luminosity*scaleFactor);
					MSPlot["MET_mme_ch"]->Fill(mets[0]->Et(),datasets[d], true, Luminosity*scaleFactor);
					MSPlot["Mtt_mme_ch"]->Fill((MyTopFCNC_EvtCand->smDecayTop()+MyTopFCNC_EvtCand->fcncDecayTop()).M(),datasets[d], true, Luminosity*scaleFactor);
					MSPlot["Mzq_mme_ch"]->Fill(MyTopFCNC_EvtCand->fcncDecayTop().M(),datasets[d], true, Luminosity*scaleFactor);
					MSPlot["MassChi2_mme_ch"]->Fill(MyTopFCNC_EvtCand->MassChi2(),datasets[d], true, Luminosity*scaleFactor);
					MSPlot["MassChi2_SM_mme_ch"]->Fill(MyTopFCNC_EvtCand->MassChi2_SM(),datasets[d],true,Luminosity*scaleFactor);
					MSPlot["MassChi2_FCNC_mme_ch"]->Fill(MyTopFCNC_EvtCand->MassChi2_FCNC(),datasets[d],true, Luminosity*scaleFactor);
					MSPlot["TriLeptonInvMass_mme_ch"]->Fill(invMass, datasets[d], true, Luminosity*scaleFactor);
					histo2D[("MET_vs_Mzq_mme_ch_"+datasets[d]->Name()).c_str()]->Fill(mets[0]->Et(),MyTopFCNC_EvtCand->fcncDecayTop().M());
				}
			}
		}
	}
	else if(diMuon && selectedExtraMuons.size()==1 && selectedElectrons.size()==0){
		selecTableTriMu.Fill(d,6,scaleFactor);
		MSPlot["NbOfSelectedJets_mmm_ch"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
		if(selectedJets.size()>0){ //at least 1 jet
			selecTableTriMu.Fill(d,7,scaleFactor); 
			if(selectedJets.size()>1){ //at least 2 jets
				selecTableTriMu.Fill(d,8,scaleFactor);
				if(selectedJets.size()>2){ //at least 3 jets
					selecTableTriMu.Fill(d,9,scaleFactor);
					MSPlot["MET_mmm_ch"]->Fill(mets[0]->Et(),datasets[d], true, Luminosity*scaleFactor);
					MyTopFCNC_EvtCand = new TopFCNC_Evt(TopFCNC_Evt::kMuon,TopFCNC_Evt::kMuon);
					MyTopFCNC_EvtCand->ReconstructTriLeptEvt(selectedMuons[idx_Z_1], selectedMuons[idx_Z_2], selectedExtraMuons[0], selectedJets, mets[0], true);
					invMass = (*selectedMuons[idx_Z_1]+*selectedMuons[idx_Z_2]+*selectedExtraMuons[0]).M();
					MSPlot["BdiscBJetCand_mmm_ch_CVS"]->Fill(MyTopFCNC_EvtCand->B().btag_combinedSecondaryVertexBJetTags(),datasets[d], true, Luminosity*scaleFactor);
					MSPlot["BdiscBJetCand_mmm_ch_TCHE"]->Fill(MyTopFCNC_EvtCand->B().btag_trackCountingHighEffBJetTags(),datasets[d], true, Luminosity*scaleFactor);
					MSPlot["Mtt_mmm_ch"]->Fill((MyTopFCNC_EvtCand->smDecayTop()+MyTopFCNC_EvtCand->fcncDecayTop()).M(),datasets[d], true, Luminosity*scaleFactor);
					MSPlot["Mzq_mmm_ch"]->Fill(MyTopFCNC_EvtCand->fcncDecayTop().M(),datasets[d], true, Luminosity*scaleFactor);
					MSPlot["MassChi2_mmm_ch"]->Fill(MyTopFCNC_EvtCand->MassChi2(),datasets[d], true, Luminosity*scaleFactor);
					MSPlot["MassChi2_SM_mmm_ch"]->Fill(MyTopFCNC_EvtCand->MassChi2_SM(),datasets[d],true,Luminosity*scaleFactor);
					MSPlot["MassChi2_FCNC_mmm_ch"]->Fill(MyTopFCNC_EvtCand->MassChi2_FCNC(),datasets[d],true, Luminosity*scaleFactor);
					MSPlot["TriLeptonInvMass_mmm_ch"]->Fill(invMass, datasets[d], true, Luminosity*scaleFactor);
					histo2D[("MET_vs_Mzq_mmm_ch_"+datasets[d]->Name()).c_str()]->Fill(mets[0]->Et(),MyTopFCNC_EvtCand->fcncDecayTop().M());
				}
			}
		}
	}
	else if(diElectron && selectedExtraElectrons.size()==0 && selectedMuons.size()==0){
		selecTableDiEl.Fill(d,6,scaleFactor);
		MSPlot["NbOfSelectedJets_ee_ch"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
		if(selectedJets.size()>0){ //at least 1 jet
			selecTableDiEl.Fill(d,7,scaleFactor); 
			if(selectedJets.size()>1){ //at least 2 jets
				selecTableDiEl.Fill(d,8,scaleFactor);
				if(selectedJets.size()>2){ //at least 3 jets
					selecTableDiEl.Fill(d,9,scaleFactor);
					if(selectedJets.size()>3){ //at least 4 jets
						selecTableDiEl.Fill(d,10,scaleFactor);
						sort(selectedJets.begin(),selectedJets.end(),HighestCVSBtag());
						MSPlot["HighestBdisc_ee_ch_CVS"]->Fill(selectedJets[0]->btag_combinedSecondaryVertexBJetTags(),datasets[d], true, Luminosity*scaleFactor);
						highestbtagdisc = selectedJets[0]->btag_combinedSecondaryVertexBJetTags();
						sort(selectedJets.begin(),selectedJets.end(),HighestTCHEBtag());
						MSPlot["HighestBdisc_ee_ch_TCHE"]->Fill(selectedJets[0]->btag_trackCountingHighEffBJetTags(),datasets[d], true, Luminosity*scaleFactor);
						MyTopFCNC_EvtCand = new TopFCNC_Evt(TopFCNC_Evt::kElec);
						MyTopFCNC_EvtCand->ReconstructDiLeptEvt(selectedElectrons[idx_Z_1], selectedElectrons[idx_Z_2], selectedJets);
						MSPlot["BdiscBJetCand_ee_ch_CVS"]->Fill(MyTopFCNC_EvtCand->B().btag_combinedSecondaryVertexBJetTags(),datasets[d], true, Luminosity*scaleFactor);
						MSPlot["BdiscBJetCand_ee_ch_TCHE"]->Fill(MyTopFCNC_EvtCand->B().btag_trackCountingHighEffBJetTags(),datasets[d],true, Luminosity*scaleFactor);
						if(highestbtagdisc<0.7) continue;
						selecTableDiEl.Fill(d,11,scaleFactor);
						MSPlot["MET_ee_ch"]->Fill(mets[0]->Et(),datasets[d], true, Luminosity*scaleFactor);
						MSPlot["Mtt_ee_ch"]->Fill((MyTopFCNC_EvtCand->smDecayTop()+MyTopFCNC_EvtCand->fcncDecayTop()).M(),datasets[d],true,Luminosity*scaleFactor);
						MSPlot["Mzq_ee_ch"]->Fill(MyTopFCNC_EvtCand->fcncDecayTop().M(),datasets[d], true, Luminosity*scaleFactor);
						MSPlot["MassChi2_ee_ch"]->Fill(MyTopFCNC_EvtCand->MassChi2(),datasets[d], true, Luminosity*scaleFactor);
						MSPlot["MassChi2_SM_ee_ch"]->Fill(MyTopFCNC_EvtCand->MassChi2_SM(),datasets[d], true, Luminosity*scaleFactor);
						MSPlot["MassChi2_FCNC_ee_ch"]->Fill(MyTopFCNC_EvtCand->MassChi2_FCNC(),datasets[d], true, Luminosity*scaleFactor);
						histo2D[("MET_vs_Mzq_ee_ch_"+datasets[d]->Name()).c_str()]->Fill(mets[0]->Et(),MyTopFCNC_EvtCand->fcncDecayTop().M());
					}
				}
			}
		}
	}
	else if(diElectron && selectedExtraElectrons.size()==0 && selectedMuons.size()==1){
		selecTableDiElMu.Fill(d,6,scaleFactor);
		MSPlot["NbOfSelectedJets_eem_ch"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
		if(selectedJets.size()>0){ //at least 1 jet
			selecTableDiElMu.Fill(d,7,scaleFactor); 
			if(selectedJets.size()>1){ //at least 2 jets
				selecTableDiElMu.Fill(d,8,scaleFactor);
				if(selectedJets.size()>2){ //at least 3 jets
					selecTableDiElMu.Fill(d,9,scaleFactor);
					MyTopFCNC_EvtCand = new TopFCNC_Evt(TopFCNC_Evt::kElec,TopFCNC_Evt::kMuon);
					MyTopFCNC_EvtCand->ReconstructTriLeptEvt(selectedElectrons[idx_Z_1], selectedElectrons[idx_Z_2], selectedMuons[0], selectedJets, mets[0]);
					invMass = (*selectedElectrons[idx_Z_1]+*selectedElectrons[idx_Z_2]+*selectedMuons[0]).M();
					MSPlot["BdiscBJetCand_eem_ch_CVS"]->Fill(MyTopFCNC_EvtCand->B().btag_combinedSecondaryVertexBJetTags(),datasets[d],true,Luminosity*scaleFactor);
					MSPlot["BdiscBJetCand_eem_ch_TCHE"]->Fill(MyTopFCNC_EvtCand->B().btag_trackCountingHighEffBJetTags(),datasets[d],true,Luminosity*scaleFactor);
					MSPlot["MET_eem_ch"]->Fill(mets[0]->Et(),datasets[d], true, Luminosity*scaleFactor);
					MSPlot["Mtt_eem_ch"]->Fill((MyTopFCNC_EvtCand->smDecayTop()+MyTopFCNC_EvtCand->fcncDecayTop()).M(),datasets[d], true, Luminosity*scaleFactor);
					MSPlot["Mzq_eem_ch"]->Fill(MyTopFCNC_EvtCand->fcncDecayTop().M(),datasets[d], true, Luminosity*scaleFactor);
					MSPlot["MassChi2_eem_ch"]->Fill(MyTopFCNC_EvtCand->MassChi2(),datasets[d], true, Luminosity*scaleFactor);
					MSPlot["MassChi2_SM_eem_ch"]->Fill(MyTopFCNC_EvtCand->MassChi2_SM(),datasets[d],true,Luminosity*scaleFactor);
					MSPlot["MassChi2_FCNC_eem_ch"]->Fill(MyTopFCNC_EvtCand->MassChi2_FCNC(),datasets[d],true,Luminosity*scaleFactor); 					
					MSPlot["TriLeptonInvMass_eem_ch"]->Fill(invMass, datasets[d], true, Luminosity*scaleFactor);
          histo2D[("MET_vs_Mzq_eem_ch_"+datasets[d]->Name()).c_str()]->Fill(mets[0]->Et(),MyTopFCNC_EvtCand->fcncDecayTop().M());
				}
			}
		}
	}
	else if(diElectron && selectedExtraElectrons.size()==1 && selectedMuons.size()==0){
		selecTableTriEl.Fill(d,6,scaleFactor);
		MSPlot["NbOfSelectedJets_eee_ch"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
		if(selectedJets.size()>0){ //at least 1 jet
			selecTableTriEl.Fill(d,7,scaleFactor); 
			if(selectedJets.size()>1){ //at least 2 jets
				selecTableTriEl.Fill(d,8,scaleFactor);
				if(selectedJets.size()>2){ //at least 3 jets
					selecTableTriEl.Fill(d,9,scaleFactor);
					MSPlot["MET_eee_ch"]->Fill(mets[0]->Et(),datasets[d], true, Luminosity*scaleFactor);
					MyTopFCNC_EvtCand = new TopFCNC_Evt(TopFCNC_Evt::kElec,TopFCNC_Evt::kElec);
					MyTopFCNC_EvtCand->ReconstructTriLeptEvt(selectedElectrons[idx_Z_1], selectedElectrons[idx_Z_2], selectedExtraElectrons[0], selectedJets, mets[0]);
					invMass = (*selectedElectrons[idx_Z_1]+*selectedElectrons[idx_Z_2]+*selectedExtraElectrons[0]).M();
					MSPlot["BdiscBJetCand_eee_ch_CVS"]->Fill(MyTopFCNC_EvtCand->B().btag_combinedSecondaryVertexBJetTags(),datasets[d],true, Luminosity*scaleFactor);
					MSPlot["BdiscBJetCand_eee_ch_TCHE"]->Fill(MyTopFCNC_EvtCand->B().btag_trackCountingHighEffBJetTags(),datasets[d],true,Luminosity*scaleFactor);
					MSPlot["Mtt_eee_ch"]->Fill((MyTopFCNC_EvtCand->smDecayTop()+MyTopFCNC_EvtCand->fcncDecayTop()).M(),datasets[d], true, Luminosity*scaleFactor);
					MSPlot["Mzq_eee_ch"]->Fill(MyTopFCNC_EvtCand->fcncDecayTop().M(),datasets[d], true, Luminosity*scaleFactor);
					MSPlot["MassChi2_eee_ch"]->Fill(MyTopFCNC_EvtCand->MassChi2(),datasets[d], true, Luminosity*scaleFactor);
					MSPlot["MassChi2_SM_eee_ch"]->Fill(MyTopFCNC_EvtCand->MassChi2_SM(),datasets[d],true,Luminosity*scaleFactor);
					MSPlot["MassChi2_FCNC_eee_ch"]->Fill(MyTopFCNC_EvtCand->MassChi2_FCNC(),datasets[d],true,Luminosity*scaleFactor); 					
					MSPlot["TriLeptonInvMass_eee_ch"]->Fill(invMass, datasets[d], true, Luminosity*scaleFactor);
          histo2D[("MET_vs_Mzq_eee_ch_"+datasets[d]->Name()).c_str()]->Fill(mets[0]->Et(),MyTopFCNC_EvtCand->fcncDecayTop().M());
				}
			}
		}
	}
	else fourIsoLeptCounter++;


// opposite charge leptons
//if(selectedMuons[0]->charge()== selectedMuons[1]->charge())
//require that there are no two electrons forming the Z mass
//if( selection.foundZCandidate(selectedMuons, selectedMuons, 20.) )
//it should not be an electron from a conversion!


	//delete selection;
	if(MyTopFCNC_EvtCand) delete MyTopFCNC_EvtCand;
    }//loop on events
    
    cout<<endl;
    cout<<"FYI ; nb of events with at least four isolated leptons = "<<fourIsoLeptCounter<<endl;
    
    //important: free memory
    treeLoader.UnLoadDataset();

    if(jetTools) delete jetTools;
    
  } //loop on datasets
    
  //Once everything is filled ...
  cout << " We ran over all the data ;-)" << endl;
  
  ///////////////////
  // Writing
  //////////////////
  cout << " - Writing outputs to the files ..." << endl;

  //Selection tables
  if(diMuon){ 
	//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
	selecTableDiMu.TableCalculator(  false, true, true, true, true);
	selecTableDiMuEl.TableCalculator(false, true, true, true, true);
	selecTableTriMu.TableCalculator( false, true, true, true, true);
  //Options : WithError (false), writeMerged (true), useBookTabs (false), addRawNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false)
	selecTableDiMu.Write(  "TopFCNC"+postfix+"Table_DiMu.tex",    true,true,true,true,false,false,true);
	selecTableDiMuEl.Write("TopFCNC"+postfix+"Table_DiMuElec.tex",true,true,true,true,false,false,true);
	selecTableTriMu.Write( "TopFCNC"+postfix+"Table_TriMu.tex",   true,true,true,true,false,false,true);
  }
  else if(diElectron){
	//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
	selecTableDiEl.TableCalculator(  false, true, true, true, true);
	selecTableDiElMu.TableCalculator(false, true, true, true, true);
	selecTableTriEl.TableCalculator( false, true, true, true, true);
  //Options : WithError (false), writeMerged (true), useBookTabs (false), addRawNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false)
	selecTableDiEl.Write(  "TopFCNC"+postfix+"Table_DiEl.tex",  true,true,true,true,false,false,true);
	selecTableDiElMu.Write("TopFCNC"+postfix+"Table_DiElMu.tex",true,true,true,true,false,false,true);
	selecTableTriEl.Write( "TopFCNC"+postfix+"Table_TriEl.tex", true,true,true,true,false,false,true);
  }
  fout->cd();
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
	MultiSamplePlot *temp = it->second;
	//temp->addText("CMS preliminary");
	string name = it->first;
	temp->Draw(false, name, true, true, true, true, true,1,false); // merge TT/QCD/W/Z/ST/
	//Draw(bool addRandomPseudoData = false, string label = string("CMSPlot"), bool mergeTT = false, bool mergeQCD = false, bool mergeW = false, bool mergeZ = false, bool mergeST = false, int scaleNPSignal = 1, bool addRatio = false, bool mergeVV = false, bool mergeTTV = false);
	temp->Write(fout, name, true, pathPNG, "pdf");
  }
  TDirectory* th1dir = fout->mkdir("Histos1D");
  th1dir->cd();
  for(map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
  {
	TH1F *temp = it->second;
	temp->Write();
	//TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
	//tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
  }
  TDirectory* th2dir = fout->mkdir("Histos2D");
  th2dir->cd();
   for(map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
  {
	TH2F *temp = it->second;
	temp->Write();
	//TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
	//tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
  }

  
  //delete  
  delete fout;

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;

  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;

  return 0;
}

