#include "TStyle.h"
#include "TF2.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include "TRandom3.h"

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
#include "../Tools/interface/TwoDimTemplateTools.h"
#include "../Tools/interface/InclFourthGenSearchTools.h"
#include "../JESMeasurement/interface/JetCombiner.h"
#include "../Reconstruction/interface/JetCorrectorParameters.h"
#include "../Reconstruction/interface/JetCorrectionUncertainty.h"
#include "../Reconstruction/interface/MakeBinning.h"
#include "../Reconstruction/interface/TTreeObservables.h"
#include "../MCInformation/interface/Lumi3DReWeighting.h"
//#include "../MCInformation/interface/LumiReWeighting.h" 
//for Kinematic Fit
#include "../MCInformation/interface/ResolutionFit.h"
#include "../KinFitter/interface/TKinFitter.h"
#include "../KinFitter/interface/TFitConstraintM.h"
#include "../KinFitter/interface/TFitParticleEtThetaPhi.h"

#include "Style.C"

using namespace std;
using namespace TopTree;
using namespace RooFit;
//using namespace reweight;

//defined in InclFourthGenSearchTools.h
//std::string IntToStr( int n ){	std::ostringstream result;	result << n;	return result.str();}

struct sort_pair_decreasing
{
    bool operator()(const std::pair<int,float> &left, const std::pair<int,float> &right)
    {
        return left.second > right.second;
    }
};

//To cout the Px, Py, Pz, E and Pt of objects
void coutObjectsFourVector(vector < TRootMuon* > init_muons, vector < TRootElectron* > init_electrons, vector < TRootJet* > init_jets, vector < TRootMET* > mets, string Comment);

float jetprob(float jetpt, float btagvalue);

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

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

  string btagger = "TCHPM";
//from BTV-11-001; WARNING: btag SFs are for jets of 20-240 GeV, averaged over eta; and mistag SFs are for jets of 50-80 GeV
// b-tag scalefactor => TCHEL: data/MC scalefactor = 0.95 +- 0.10,    TCHEM: data/MC scalefactor = 0.94 +- 0.09,	TCHPM: data/MC scalefactor =  0.91 +- 0.09
// mistag scalefactor => TCHEL: data/MC scalefactor = 1.11 +- 0.12,    TCHEM: data/MC scalefactor = 1.21 +- 0.17,	TCHPM: data/MC scalefactor =  1.27 +- 0.15
  float scalefactorbtageff = 1, mistagfactor = 1;
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
	else if(btagger == "TCHPM") //track counting high pur medium working point
  {
  	  if(dobTagEffShift == 0)
		scalefactorbtageff = 0.91;
	  if(dobTagEffShift == 1)
		scalefactorbtageff = 0.8194;
	  if(dobTagEffShift == 2)
		scalefactorbtageff = 1.0006;
		
	  if(domisTagEffShift == 0)
		mistagfactor = 1.27;
	  if(domisTagEffShift == 1)
		mistagfactor = 1.1187;
	  if(domisTagEffShift == 2)
		mistagfactor = 1.4213;
  }
	else if(btagger == "TCHPT") //track counting high pur tight working point
  {
	  cout<<"WARNING: look up SFs for TCHPT"<<endl;
  	 /* if(dobTagEffShift == 0)
		scalefactorbtageff = 0.91;
	  if(dobTagEffShift == 1)
		scalefactorbtageff = 0.8194;
	  if(dobTagEffShift == 2)
		scalefactorbtageff = 1.0006;
		
	  if(domisTagEffShift == 0)
		mistagfactor = 1.27;
	  if(domisTagEffShift == 1)
		mistagfactor = 1.1187;
	  if(domisTagEffShift == 2)
		mistagfactor = 1.4213;*/
  }
  float workingpointvalue = 9999; //{1.7,3.3,10.2}; trackcountinghighefficiency working points: loose, medium, tight
  if(btagger == "TCHEL")
     workingpointvalue = 1.7;
  else if(btagger == "TCHEM")
     workingpointvalue = 3.3;
  else if(btagger == "TCHPM")
     workingpointvalue = 1.93;
	else if(btagger == "TCHPT")
     workingpointvalue = 3.41;

  clock_t start = clock();

  cout << "*************************************************************" << endl;
  cout << " Beginning of the program for the fourth generation search ! " << endl;
  cout << "*************************************************************" << endl;

  //SetStyle if needed
  setTDRStyle();
  //setMyStyle();

  string postfix = "JetMultiplicityTEST_25Jan12"; // to relabel the names of the output file

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
  if(doPUShift == 1)
    postfix= postfix+"_PUMinus";
  if(doPUShift == 2)
    postfix= postfix+"_PUPlus";


  /////////////////////
  // Configuration
  /////////////////////
  bool useMassesAndResolutions = true;
  bool doMVAjetcombination = true; //when false, the jet combination and the top mass will not be reconstructed, and nothing will be trained
  bool TrainMVA = false; // If false, the previously trained MVA will be used to calculate stuff. Note: there is an MVA output file with the training, but also some files in the ./weights directory!!
  if (doJESShift != 0 || doJERShift != 0 || dobTagEffShift != 0 || domisTagEffShift != 0 || doPUShift != 0)
  {
		useMassesAndResolutions = true;
		doMVAjetcombination = true;
		TrainMVA = false;
  }
  
  bool TrainwithTprime = false;
  string MVAmethod = "Likelihood"; // MVAmethod to be used to get the good jet combi calculation (not for training! this is chosen in the jetcombiner class)
  string channelpostfix = "";
  bool semiElectron = false; // use semiElectron channel?
  bool semiMuon = true; // use semiMuon channel?
  if(semiElectron && semiMuon)
  {
     cout << "  --> Using both semiMuon and semiElectron channel? Choose only one (for the moment, since this requires running on different samples/skims)!" << endl;
     exit(1);
  }
  else
  {
    if(semiMuon){
       cout << " --> Using the semiMuon channel..." << endl;
       channelpostfix = "_semiMu";
    }
    else if(semiElectron){
       cout << " --> Using the semiElectron channel..." << endl;
       channelpostfix = "_semiEl";
    }
  }
  
  //xml file
  string xmlFileName = "";
  if(semiElectron) xmlFileName = "../config/myFourthGenconfig_Electron.xml";
  else if(semiMuon) xmlFileName = "../config/myFourthGenconfig.xml";
  const char *xmlfile = xmlFileName.c_str();
  cout << "used config file: " << xmlfile << endl;    
  
  //Configuration and variables for kinematic fit
  bool doKinematicFit = false;
  bool CalculateResolutions = false; // If false, the resolutions will be loaded from a previous calculation and the kinfit can be done
  ResolutionFit *resFitLightJets = 0, *resFitBJets = 0;
  if(CalculateResolutions && doKinematicFit)
  {
    resFitLightJets = new ResolutionFit("LightJet");
    resFitBJets = new ResolutionFit("BJet"); //BJet resolutions not needed, but will still be made since the input of a function requires this...
  }
  else if(!CalculateResolutions && doKinematicFit)
  {
    resFitLightJets = new ResolutionFit("LightJet");
    resFitLightJets->LoadResolutions("lightJetReso.root");// if you do not want to calculate the resolutions, copy these files from Stijn or so and put them in your macros/ directory
    resFitBJets = new ResolutionFit("BJet");
    resFitBJets->LoadResolutions("bJetReso.root");
  }
  
  ////////////////////////////////////
  /// AnalysisEnvironment  
  ////////////////////////////////////

  AnalysisEnvironment anaEnv;
  cout<<" - Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  int verbose = anaEnv.Verbose;
  float anaEnvLuminosity = anaEnv.Luminosity;	// in 1/pb 
  cout << "analysis environment luminosity for rescaling "<< anaEnvLuminosity << endl;

  /////////////////////
  // Load Datasets
  /////////////////////

  TTreeLoader treeLoader;
  vector < Dataset* > datasets;
  cout << " - Loading datasets ..." << endl;
  treeLoader.LoadDatasets (datasets, xmlfile);
  
  //is this block needed?
  float Luminosity = anaEnvLuminosity;
  for (unsigned int d = 0; d < datasets.size (); d++)
  {
    	string dataSetName = datasets[d]->Name();
	if(dataSetName.find("Data")<=0 || dataSetName.find("data")<=0 || dataSetName.find("DATA")<=0)
	{
		  Luminosity = datasets[d]->EquivalentLumi();
		  break;
	 }
  }
  if(Luminosity != anaEnvLuminosity) cout << "changed analysis environment luminosity to "<< Luminosity << endl;
  
  //Output ROOT file
  string rootFileName ("InclFourthGenSearch"+postfix+channelpostfix+".root");
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");

 
  //vector of objects
  cout << " - Variable declaration ..." << endl;
  vector < TRootVertex* > vertex;
  vector < TRootMuon* > init_muons;
  vector < TRootElectron* > init_electrons;
  vector < TRootJet* > init_jets;
  vector < TRootMET* > mets;

  //Global variable
  TRootEvent* event = 0;

  string pathPNG = "InclFourthGenSearchPlots"+postfix+channelpostfix;
  pathPNG = pathPNG +"/"; 	
  pathPNG = pathPNG +"/"; 	
  mkdir(pathPNG.c_str(),0777);

  //Most 1D and MS plots are declared inside this class
  InclFourthGenSearchTools myInclFourthGenSearchTools(semiMuon, semiElectron, datasets, Luminosity,doKinematicFit);
  if(!CalculateResolutions && doKinematicFit) myInclFourthGenSearchTools.SetResolutionFit(resFitLightJets);
  
  MSPlot["allDiJetMasses"] = new MultiSamplePlot(datasets, "allDiJetMasses", 50, 0, 1000, "m_{jj}");
  MSPlot["hadronicRecoWMass_chosenWjets"] = new MultiSamplePlot(datasets, "hadronicRecoWMass_chosenWjets", 50, 0, 500, "m_{W}"); 
  histo1D["hadronicPartonWMass"] = new TH1F("hadronicPartonWMass","Hadronic W Mass, using the Partons",100,0,200);
  histo1D["hadronicRecoWMass"] = new TH1F("hadronicRecoWMass","Hadronic W Mass, using the RecoJets",100,0,200);
  
  histo1D["lumiWeights"] = new TH1F("lumiWeights","lumiWeights;lumiWeight;#events",100,0,4);
  histo1D["LeptonPt_TTbar"] = new TH1F("leptonspt ttbar","leptonspt ttbar;pt leptons;#events",250,0,500);
  histo1D["LeptonPt_Tprime500"] = new TH1F("leptonspt tprime500","leptonspt tprime500;pt leptons;#events",250,0,500);
  histo1D["LeptonPt_Bprime500"] = new TH1F("leptonspt bprime500","leptonspt bprime500;pt leptons;#events",250,0,500);
  histo1D["LeptonPt_SBprime500"] = new TH1F("leptonspt sbprime500","leptonspt sbprime500;pt leptons;#events",250,0,500);
  
  string multileptons[2] = {"SSLeptons","TriLeptons"};
  string histoName,histo_dataset;
  for(int i = 0; i<2; i++)
  {
		histoName = "NbEvents_"+multileptons[i];
		for(unsigned int d = 0; d < datasets.size (); d++){
			histo_dataset = histoName+(datasets[d]->Name()).c_str(); 
			histo1D[histo_dataset.c_str()] = new TH1F(histo_dataset.c_str(),histo_dataset.c_str(), 1, 0.5, 1.5);
		}
  }
	
  MSPlot["MS_NbSSevents"] = new MultiSamplePlot(datasets,"# events with SS leptons", 1, 0.5, 1.5, "");
  MSPlot["MS_NbTrievents"] = new MultiSamplePlot(datasets,"# events with 3 leptons", 1, 0.5, 1.5, "");
  MSPlot["MS_MET"] = new MultiSamplePlot(datasets,"MET", 75, 0, 150, "");
  MSPlot["MS_LeptonPt"] = new MultiSamplePlot(datasets,"lepton pt", 150, 0, 300, "");
  MSPlot["MS_nPV"] = new MultiSamplePlot(datasets, "nPrimaryVertices", 21, -0.5, 20.5, "Nr. of primary vertices");
  MSPlot["MS_JetMultiplicity_SingleLepton"] = new MultiSamplePlot(datasets, "JetMultiplicity", 10, -0.5, 9.5, "Jet Multiplicity");
  MSPlot["MS_BtaggedJetMultiplicity_SingleLepton"] = new MultiSamplePlot(datasets, "BtaggedJetMultiplicity", 7, -0.5, 6.5, "b-tagged jet multiplicity");
	MSPlot["MS_JetMultiplicityAtleast1Btag_SingleLepton"] = new MultiSamplePlot(datasets, "JetMultiplicityAtleast1Btag", 10, -0.5, 9.5, "Jet multiplicity (>=1 b-tag)");

  MSPlot["MS_JetPt_all_SingleLepton"] = new MultiSamplePlot(datasets,"JetPt_all", 50, 0, 300, "Pt of all jets (GeV)");
	MSPlot["MS_JetPt_btagged_SingleLepton"] = new MultiSamplePlot(datasets,"JetPt_btagged", 50, 0, 300, "Pt of b-tagged jets (GeV)");
	MSPlot["MS_JetPt_nonbtagged_SingleLepton"] = new MultiSamplePlot(datasets,"JetPt_nonbtagged", 50, 0, 300, "Pt of non b-tagged jets (GeV)");
	
	cout << " - Declared histograms ..." <<  endl;

  float NbSSevents = 0;  
  float NbSSevents_1B_2W = 0; 
  float NbSSevents_1B_3W = 0; 
  float NbSSevents_1B_4W = 0; 
  float NbSSevents_2B_2W = 0; 
  float NbSSevents_2B_3W = 0; 
  float NbSSevents_2B_4W = 0; 
  float NP_NbSSevents_1B_2W = 0; 
  float NP_NbSSevents_1B_3W = 0; 
  float NP_NbSSevents_1B_4W = 0; 
  float NP_NbSSevents_2B_2W = 0; 
  float NP_NbSSevents_2B_3W = 0; 
  float NP_NbSSevents_2B_4W = 0; 
  float NbTrievents = 0; 
  float NbTrievents_1B_3W = 0; 
  float NbTrievents_1B_4W = 0; 
  float NbTrievents_2B_3W = 0; 
  float NbTrievents_2B_4W = 0; 
  float NP_NbTrievents_1B_3W = 0; 
  float NP_NbTrievents_1B_4W = 0; 
  float NP_NbTrievents_2B_3W = 0; 
  float NP_NbTrievents_2B_4W = 0; 
  
  /////////////////////////////////////////////////////////
  //Configuration and variables for 2D HT-Mtop distribution. There is a 2D distribution for 2 boxes seperately: 1B_2W and 2B_2W
  /////////////////////////////////////////////////////////
  string xvariable = "HT", yvariable = "MTop"; //these are the two variables for which the 2D plane is made
  int nbinsxvariable = 10, nbinsyvariable = 12; //if the binning is already created, make sure these are the same as before!
  string binningFileName_HTvsMTop_1B_2W = "Binning_InclFourthGenSearch_1B_2W_TTbarJetsFlat"+channelpostfix+".root";
  string binningFileName_HTvsMTop_2B_2W = "Binning_InclFourthGenSearch_2B_2W_TTbarJetsFlat"+channelpostfix+".root";
  TwoDimTemplateTools HTvsMTop_1B_2W("1B_2W",xvariable,nbinsxvariable,yvariable,nbinsyvariable);
  TwoDimTemplateTools HTvsMTop_2B_2W("2B_2W",xvariable,nbinsxvariable,yvariable,nbinsyvariable); 
  if(doMVAjetcombination && !TrainMVA && useMassesAndResolutions) //Note: the boolean 'useMassesAndResolutions'is not needed when a Wmass file exists
  { 
    HTvsMTop_1B_2W.SetDatasets(datasets);
    HTvsMTop_2B_2W.SetDatasets(datasets);
    if(anaEnv.MCRound!=0)
    {
       HTvsMTop_1B_2W.LoadTwoDimBinning(binningFileName_HTvsMTop_1B_2W);
       HTvsMTop_2B_2W.LoadTwoDimBinning(binningFileName_HTvsMTop_2B_2W);
    }
    else
    {
       cout<<" ... The 2D binnings will be created!"<<endl;
    }
  }	
	
  ////////////////////////////////////
  /// Selection table
  ////////////////////////////////////

  vector<string> CutsSelecTableSemiMu;
  CutsSelecTableSemiMu.push_back(string("initial")); //0
  CutsSelecTableSemiMu.push_back(string("preselected"));
  CutsSelecTableSemiMu.push_back(string("trigged"));
  CutsSelecTableSemiMu.push_back(string("Good PV"));
  CutsSelecTableSemiMu.push_back(string("$\\geq$ 1 muon"));
  CutsSelecTableSemiMu.push_back(string("$\\geq$ 1 b-tagged jet"));
  CutsSelecTableSemiMu.push_back(string("MET > 40 GeV"));  
  CutsSelecTableSemiMu.push_back(string("SS $\\geq$ 1 muon"));
  CutsSelecTableSemiMu.push_back(string("trileptons $\\geq$ 1 muon"));
  CutsSelecTableSemiMu.push_back(string("box 1B 1W"));
  CutsSelecTableSemiMu.push_back(string("box 1B 2W"));//10
  CutsSelecTableSemiMu.push_back(string("box 1B 3W"));
  CutsSelecTableSemiMu.push_back(string("box 1B $\\geq$ 4W"));
  CutsSelecTableSemiMu.push_back(string("box 2B 1W"));
  CutsSelecTableSemiMu.push_back(string("box 2B 2W"));
  CutsSelecTableSemiMu.push_back(string("box 2B 3W"));
  CutsSelecTableSemiMu.push_back(string("box 2B $\\geq$ 4W"));
  CutsSelecTableSemiMu.push_back(string("SS leptons in all boxes combined")); //17
  CutsSelecTableSemiMu.push_back(string("trileptons in all boxes combined")); //18

  vector<string> CutsSelecTableSemiEl;
  CutsSelecTableSemiEl.push_back(string("initial")); //0
  CutsSelecTableSemiEl.push_back(string("preselected"));
  CutsSelecTableSemiEl.push_back(string("trigged"));
  CutsSelecTableSemiEl.push_back(string("Good PV"));
  CutsSelecTableSemiEl.push_back(string("$\\geq$ 1 selected electron"));
  CutsSelecTableSemiEl.push_back(string("Veto muon"));
  CutsSelecTableSemiEl.push_back(string("Conversion veto"));
  CutsSelecTableSemiEl.push_back(string("$\\geq$ 1 b-tagged jet"));
  CutsSelecTableSemiEl.push_back(string("MET > 40 GeV"));
  CutsSelecTableSemiEl.push_back(string("SS 2 electrons"));
  CutsSelecTableSemiEl.push_back(string("trileptons 3 electrons")); //10
  CutsSelecTableSemiEl.push_back(string("box 1B 1W"));
  CutsSelecTableSemiEl.push_back(string("box 1B 2W"));
  CutsSelecTableSemiEl.push_back(string("box 1B 3W"));
  CutsSelecTableSemiEl.push_back(string("box 1B $\\geq$ 4W"));
  CutsSelecTableSemiEl.push_back(string("box 2B 1W"));
  CutsSelecTableSemiEl.push_back(string("box 2B 2W"));
  CutsSelecTableSemiEl.push_back(string("box 2B 3W"));
  CutsSelecTableSemiEl.push_back(string("box 2B $\\geq$ 4W"));
  CutsSelecTableSemiEl.push_back(string("SS leptons in all boxes combined")); //19
  CutsSelecTableSemiEl.push_back(string("trileptons in all boxes combined")); //20

  SelectionTable selecTableSemiMu(CutsSelecTableSemiMu, datasets);
  selecTableSemiMu.SetLuminosity(Luminosity);
  SelectionTable selecTableSemiEl(CutsSelecTableSemiEl, datasets);
  selecTableSemiEl.SetLuminosity(Luminosity);
  
  cout << " - SelectionTable instantiated ..." << endl;

  ////////////////////////////////////////////////////
  // PileUp Reweighting - 3D//
  ////////////////////////////////////////////////////
  Lumi3DReWeighting Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Flat10PlusTail.root","PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root", "pileup", "pileup");
  Lumi3DWeights.weight3D_init(1.0);
	if(doPUShift == 1)
  	Lumi3DWeights.weight3D_init(0.92);
	else if(doPUShift == 2)
  	Lumi3DWeights.weight3D_init(1.08);
	else
  	Lumi3DWeights.weight3D_init(1.0);

//  LumiReWeighting LumiWeights = LumiReWeighting("PileUpReweighting/pileup_WJets_36bins.root", "PileUpReweighting/pileup_2011Data_UpToRun180252.root", "pileup2", "pileup");
//  PoissonMeanShifter PShiftUp = PoissonMeanShifter(0.6); // PU-systematic
//  PoissonMeanShifter PShiftDown = PoissonMeanShifter(-0.6); // PU-systematic
  cout << " - Initialized LumiReWeighting stuff" << endl;
  
  
  ///////////////
  // JetCombiner
  ///////////////
  JetCombiner* jetCombiner;
  if(!doMVAjetcombination) TrainMVA = false;
  else if(doMVAjetcombination)
  {
    jetCombiner = new JetCombiner(TrainMVA, Luminosity, datasets, MVAmethod, true, "",channelpostfix); //last bool is basically to use also the W mass as constraint    
  }
  

  if(doMVAjetcombination && TrainMVA) useMassesAndResolutions = true; //just to make sure the W mass plot is not produced (is not used anyway)
  
  cout << " - JetCombiner instantiated ..." << endl;

    
  ////////////////////////////////////
  ////////////////////////////////////
  ///////// Loop on datasets
  ////////////////////////////////////
  ////////////////////////////////////
  cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;

  for (unsigned int d = 0; d < datasets.size(); d++) //d < datasets.size()
  {
    if (verbose > 1)
      cout << "   Dataset " << d << " name : " << datasets[d]->Name () << " / title : " << datasets[d]->Title () << endl;
    if (verbose > 1)
      std::cout<<"      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
		
    //open files and load
    cout<<"Load Dataset"<<endl;
    treeLoader.LoadDataset (datasets[d], anaEnv);
		
    string previousFilename = "";
    int iFile = -1;
    
    string dataSetName = datasets[d]->Name();	
    selecTableSemiMu.Fill(d,0, datasets[d]->Xsection() * datasets[d]->EquivalentLumi() );
    selecTableSemiEl.Fill(d,0, datasets[d]->Xsection() * datasets[d]->EquivalentLumi() );
    
    /////////////////////////////////////
    /// Initialize JEC factors 
    /////////////////////////////////////
		
    //L2L3 residual corrections already in data Toptrees now! (because a global tag is used where these corrections are included)
    vector<JetCorrectorParameters> vCorrParam;
    /*if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") // Data!
    {
	   JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("JECFiles/Jec11V2_db_AK5PFchs_L2L3Residual.txt");
	   vCorrParam.push_back(*ResJetCorPar);
    }*/
    
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/Jec11V2_db_AK5PFchs_Uncertainty.txt");
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, false); // last boolean ('startFromRaw') = false!

    /*if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") // Data!
    {
     	JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("JECFiles/GR_R_38X_V15_AK5PF_L2Relative.txt");
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("JECFiles/GR_R_38X_V15_AK5PF_L3Absolute.txt");
      JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("JECFiles/GR_R_38X_V15_AK5PF_L2L3Residual.txt");
      vCorrParam.push_back(*L2JetCorPar);
      vCorrParam.push_back(*L3JetCorPar);
      vCorrParam.push_back(*ResJetCorPar);
    }
    else // MC!
    {
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("JECFiles/START38_V14_AK5PF_L2Relative.txt");
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("JECFiles/START38_V14_AK5PF_L3Absolute.txt");
      vCorrParam.push_back(*L2JetCorPar);
      vCorrParam.push_back(*L3JetCorPar);
    }*/  
    

    ///////////////////////////////////////////////////
    // calculate W reconstructed mass and resolution: load the mass values and resolutions for chi2 jetcombination
    ///////////////////////////////////////////////////
    
    float Wmass = -9999;
    float SigmaWmass = -9999;
    
    if(useMassesAndResolutions)
    {
      string resfilename = "MassPlot.root";
      
      cout << "INFO: Using mass and width for reconstructed W mass from: " << resfilename << endl;
      TFile* res = new TFile(resfilename.c_str(),"READ");
      res->cd();
	
      TF1* WmassFit = (TF1*)res->Get("hadronicRecoWMass_Fitted");
      if (WmassFit)
      {
	 Wmass = WmassFit->GetParameter(1); 
	 SigmaWmass = WmassFit->GetParameter(2); 	
      }

      res->Close();
      delete res;
      //delete WmassFit;
    }

    ////////////////////////////////////
    ////////////////////////////////////
    ///////// Loop on events
    ////////////////////////////////////
    ////////////////////////////////////
    int itrigger = -1, previousRun = -1;
      
    //block of code to arrange that only on 1/3th of the ttbar semimu or semiel sample is run for the training, and on 2/3th of the sample for the evaluation. 
    //Don't forget to change the eqLumi in the config depending on training or evaluation (not safe, to be changed? Not so trivial...)
    int start = 0;
    int end = datasets[d]->NofEvtsToRunOver();
    if(TrainMVA && (datasets[d]->Name () == "TTbarJets_SemiMuon" || datasets[d]->Name () == "TTbarJets_SemiElectron"))
    { 
        start = 0;
        end = int(datasets[d]->NofEvtsToRunOver()/3);
    }
    else if (!TrainMVA && (datasets[d]->Name () == "TTbarJets_SemiMuon" || datasets[d]->Name () == "TTbarJets_SemiElectron"))
    {    
        start = int(datasets[d]->NofEvtsToRunOver()/3);
        end = datasets[d]->NofEvtsToRunOver();
    }
    else{
        start = 0;
        end = datasets[d]->NofEvtsToRunOver();
    }
     
    if (verbose > 1)
      cout << " - Loop over events " << endl;      
    
    for (int ievt = start; ievt < end; ievt++)
    {        

      if(ievt%1000 == 0)
        std::cout<<"Processing the "<<ievt<<"th event ("<<100*(ievt-start)/(end-start)<<"%)"<<flush<<"\r";
        
      //load event
      event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);
      vector<TRootGenJet*> genjets;
      if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
      {
        genjets = treeLoader.LoadGenJet(ievt,false);
        sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
      
      // scale factor for the event
      float scaleFactor = 1.;
                
      // check which file in the dataset it is to have the HLTInfo right
      string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
      if(previousFilename != currentFilename)
      {
      	previousFilename = currentFilename;
        iFile++;
	 			cout<<"File changed!!! => iFile = "<<iFile<<endl;
      }
      
      ///////////////////////////////
      // trigger
      ///////////////////////////////
			int currentRun = event->runId();
			if(previousRun != currentRun)
			{
			previousRun = currentRun;
				if(semiMuon)
				{
					if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
					{
						if (event->runId() >= 160431 && event->runId() <= 163261)//May10ReReco
							itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v5"), currentRun, iFile);
  						else if (event->runId() >= 163270 && event->runId() <= 163869)
    				  			itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v6"), currentRun, iFile);
  						else if (event->runId() >= 165088 && event->runId() <= 165633)//PromptReco_v4; splitted over 2 toptrees: 565 and 641
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
				
        			// RUN2011B (promptv1)
   					else if( event->runId() >= 175860 && event->runId() <= 177452 )// TopTree ID 722
   				  	itrigger = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v3"), currentRun, iFile);
   					else if( event->runId() >=  177718 && event->runId() <=  178380 ) // TopTree ID 804
   				  	itrigger = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v3"), currentRun, iFile);
   					else if( event->runId() >=  178420 && event->runId() <=  178479 )
   				  	itrigger = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v6"), currentRun, iFile);								
						else if( event->runId() >=  178703 && event->runId() <=  179889 ) // TopTree ID 816
							itrigger = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v6"), currentRun, iFile);
						else if( event->runId() >=  179959 && event->runId() <=  180252 )
							itrigger = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v7"), currentRun, iFile); 
									   
  					if(itrigger == 9999)
						{
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
				} //end if semiMuon
	 			else if(semiElectron)
				{
	  			if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
	   			{      		
						// /SingleElectron/Run2011A-May10ReReco-v1/AOD 
						if (event->runId() >= 160404 && event->runId() < 161217)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1"), currentRun, iFile);
						else if (event->runId() >= 161217 && event->runId() < 163270)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2"), currentRun, iFile);
     				else if (event->runId() >= 163270 && event->runId() <= 163869)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3"), currentRun, iFile);				 
						// /ElectronHad/Run2011A-PromptReco-v4/AOD
						else if (event->runId() >= 165088 && event->runId() < 165970)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_TrkIdT_CentralJet30_BTagIP_v4"), currentRun, iFile);
						else if (event->runId() >= 165970 && event->runId() < 167038)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP_v1"), currentRun, iFile);
						else if (event->runId() >= 167038 && event->runId() <= 167913)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP_v2"), currentRun, iFile);			  
						// /ElectronHad/Run2011A-05Aug2011-v1/AOD
						else if (event->runId() >= 170249 && event->runId() <= 172619)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP_v4"), currentRun, iFile);  
						// /ElectronHad/Run2011A-PromptReco-v6/AOD 
						else if (event->runId() >= 172620 && event->runId() < 173212)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP_v4"), currentRun, iFile);  
						else if (event->runId() >= 173212 && event->runId() <= 173692)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP_v5"), currentRun, iFile);  				   				   	
						// RUN2011B (promptv1)
						else if (event->runId() >= 175832 && event->runId() < 178411)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP_v5"), currentRun, iFile);  				   					else if (event->runId() >= 178411 && event->runId() < 179942)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP_v8"), currentRun, iFile);  				   					else if (event->runId() >= 179942 && event->runId() <= 180296)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP_v9"), currentRun, iFile);  				   											   
  					if(itrigger == 9999)
						{
    					cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (DATA) IN RUN " << event->runId() << endl;
    					exit(1);
  					}// semi-electron
 	   			}
	   			else 
	   			{
   					itrigger = treeLoader.iTrigger (string ("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2"), currentRun, iFile);//Summer11 MC has other triggers!	
						if(itrigger == 9999)
						{
							cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (" << dataSetName << ") IN RUN " << event->runId() << endl;	
							exit(1);
						}
					}	 
				} //end if semiElectron	
			} //end previousRun != currentRun

       ///////////////////////////
       // correct MET with jet energy corrections
       ///////////////////////////
       // CHECK IF THIS IS STILL TO BE DONE!!! -> no
       /*if( dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ){
				for(unsigned int i=0; i<init_jets.size(); i++){
					mets[0]->SetPx(mets[0]->Px() + init_jets[i]->Px());
					mets[0]->SetPy(mets[0]->Py() + init_jets[i]->Py());
					mets[0]->SetE(sqrt(pow(mets[0]->Px(),2) + pow(mets[0]->Py(),2)));
				}	
				jetTools->correctJets(init_jets, vertex);
				for(unsigned int i=0; i<init_jets.size(); i++){
					mets[0]->SetPx(mets[0]->Px() - init_jets[i]->Px());
					mets[0]->SetPy(mets[0]->Py() - init_jets[i]->Py());
					mets[0]->SetE(sqrt(pow(mets[0]->Px(),2) + pow(mets[0]->Py(),2)));
				}	
      }*/

      //ordering is relevant; most probably 1) Type I MET correction, 2) JER where jet corrections are propagated to MET, 3) JES systematics where jet corrections are propagated to MET
      //----------------------------------------------------------
      // Apply type I MET corrections:  (Only for |eta| <= 4.7 )
      //---------------------------------------------------------
      //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before MET type I correction:");      
      if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
        jetTools->correctMETTypeOne(init_jets,mets[0],true);
      else
        jetTools->correctMETTypeOne(init_jets,mets[0],false);
      //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After MET type I correction:");
      
      if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
      {	
      				//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JER correction:");
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
      }

			double lumiWeight3D = 1.0;
			if(!(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")){
				////////////////////////////
      	// apply trigger Reweighting
      	////////////////////////////
				float mceventtriggerweight = 1;
				float NEWmceventtriggerweight = 1;
				float prob=0;
      	if(semiElectron){
					std::vector<float> probabilities;
					for(size_t i=0; i<init_jets.size(); ++i){
    				if(fabs(init_jets[i]->Eta())>2.6) continue;
    				probabilities.push_back(jetprob(init_jets[i]->Pt(),init_jets[i]->btag_trackCountingHighPurBJetTags())); //highEff when this is used offline
					}
				
					//use binary code for objects to be triggered or not triggered
					for(int i=0; i<pow(2.,(double)probabilities.size());++i){
    				int ntrigobj=0;
    				for(unsigned int j=0; j<probabilities.size();++j){
							if((int)(i/pow(2.,(double)j))%2) ntrigobj++;
						}
						if(ntrigobj<1) continue;  
						float newprob=1;
						for(unsigned int j=0; j<probabilities.size();++j){
							if((int)(i/pow(2.,(double)j))%2) newprob*=probabilities[j];
							else newprob*=1-probabilities[j];
						}
						prob+=newprob;
					}
					mceventtriggerweight*=prob;

					//stupid workaround, because part (single electron path corresponding to 200/pb) of the MC needs SF of 1 for trigger and part (btag trigger) needs the procedure above
					NEWmceventtriggerweight = (218./Luminosity +(mceventtriggerweight*(1-(218./Luminosity))));

					scaleFactor = scaleFactor*NEWmceventtriggerweight;
 					//cout << "mcevent triggerweight " << mceventtriggerweight << endl;
 					//cout << "scalefactor (only triggerweight) " << scaleFactor << endl;
      	}
        
      	////////////////////////////
      	// apply PU Reweighting
      	////////////////////////////
				lumiWeight3D = Lumi3DWeights.weight3D(event->nPu(-1),event->nPu(0),event->nPu(+1));
	 			scaleFactor = scaleFactor*lumiWeight3D;
      	histo1D["lumiWeights"]->Fill(scaleFactor);	
			}
						
								
      /////////////////////////////
      // Selection
      /////////////////////////////

     //Declare selection instance    
      Selection selection(init_jets, init_muons, init_electrons, mets); //mets can also be corrected...
      Selection nonstandard_selection(init_jets, init_muons, init_electrons, mets); //mets can also be corrected... 
      
      if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
      {
        // Apply the scraping veto. Note: should be checked if still necessary, maybe already done in toptree production
        bool isBeamBG = true;
        if(event->nTracks() > 10)
        {
          if( ( (float) event->nHighPurityTracks() ) / ( (float) event->nTracks() ) > 0.25 )
            isBeamBG = false;
      	}
      	if(isBeamBG) continue;
      }	

      bool trigged, isGoodPV;
      trigged = treeLoader.EventTrigged (itrigger);			
      isGoodPV = selection.isPVSelected(vertex, 4, 24., 2); //in the past this was put in the config, but this is not very useful, since the PV cuts are quite standard
			
      bool eventSelected = false;
      bool isSingleLepton = false;
      bool isSSLepton = false;
      bool isTriLepton = false;
      bool isSingleMuon = false;
      bool isSingleElectron = false;
      bool isSSMuon = false;
      bool isSSElectron = false;
      bool isSSMuEl = false;
      bool isTriMuon = false;
      bool isTriElectron = false;
      bool isTriMu2El1 = false;
      bool isTriMu1El2 = false;
      
      vector<TRootJet*> selectedJets,selectedJetsFromW,selectedJetsFromW_DropUsedJets,selectedJetsFromW_DropUsedJets_tmp;
      vector<TRootJet*> selectedForwardJets, selectedJetsLargeEtaRange;
      vector<TRootMuon*> selectedMuons;
      vector<TRootElectron*> selectedElectrons;
      vector<TRootElectron*> selectedLooseElectronsNoVBTFid;
      vector<TRootElectron*> selectedLooseElectronsVBTFid;
      vector<TRootMuon*> selectedLooseMuons;
      vector<TRootMCParticle*> mcParticles;
      
      if(dataSetName.find("TTbarJets_SemiMu") == 0 || dataSetName.find("TTbarJets_SemiElectron") == 0 || dataSetName.find("NP_Tprime")==0 || dataSetName.find("NP_overlay_Tprime")==0)
      {
        treeLoader.LoadMCEvent(ievt, 0, 0, mcParticles,false);  
        sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
      
      float METCut = 40;
      selection.setJetCuts(30.,2.4,0.01,1.,0.98,0.3,0.1);
      nonstandard_selection.setJetCuts(30.,4.7,0.01,1.,0.98,0.3,0.1); //only difference: larger eta acceptance 
      selection.setMuonCuts(40,2.1,0.1,10,0.02,0.3,1,1,1);
      selection.setLooseMuonCuts(10,2.5,0.2);
      selection.setElectronCuts(40,2.5,0.125,0.02,1,0.3);
      selection.setLooseElectronCuts(15,2.5,0.2);	 				 
      
      if (init_jets.size() > 0)
      {
	    	selectedJets = selection.GetSelectedJets(true);				
	    	selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
	    	selectedElectrons = selection.GetSelectedElectrons(vertex[0],selectedJets);
	    	selectedJetsLargeEtaRange = nonstandard_selection.GetSelectedJets(true);	    	    
      }
      selectedLooseElectronsNoVBTFid = selection.GetSelectedLooseElectrons(false); //no vbtfid is required
      selectedLooseElectronsVBTFid = selection.GetSelectedLooseElectrons(true); //vbtfid is required 
      selectedLooseMuons = selection.GetSelectedLooseMuons(); //veto muons	
      
      for(unsigned int i = 0; i<selectedJetsLargeEtaRange.size(); i++)
      {
	    	if(selectedJetsLargeEtaRange[i]->Eta() > 2.4)
	      	selectedForwardJets.push_back(selectedJetsLargeEtaRange[i]);
      }
     
      selecTableSemiMu.Fill(d,1,scaleFactor);
      selecTableSemiEl.Fill(d,1,scaleFactor);
		
			//// EVENTS TRIGGERED BY MUON TRIGGER			
      if(trigged && semiMuon)
      { 
				selecTableSemiMu.Fill(d,2,scaleFactor);
        if(isGoodPV)
				{
					selecTableSemiMu.Fill(d,3,scaleFactor);
					if(selectedMuons.size()>=1)
					{
						selecTableSemiMu.Fill(d,4,scaleFactor);
						sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class
						
						if(selectedJets.size()>=(unsigned int)anaEnv.NofJets)
						{  //at least 1 jet!
						
						  //block for the jet multiplicity plot
							if(mets[0]->Et()> METCut)
							{							
							      if(selectedMuons.size() == 1 && selectedLooseMuons.size() == selectedMuons.size() && selectedLooseElectronsNoVBTFid.size() == 0)
										{
										  int nBtags = 0;
											for(unsigned int j=0;j<selectedJets.size();j++)
											{
											   if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
												 {
												    nBtags++;
												 }
											}
											MSPlot["MS_JetMultiplicity_SingleLepton"]->Fill(selectedJets.size(),datasets[d], true, Luminosity*scaleFactor);
											MSPlot["MS_BtaggedJetMultiplicity_SingleLepton"]->Fill(nBtags,datasets[d], true, Luminosity*scaleFactor);
										  if(nBtags>0)
											  MSPlot["MS_JetMultiplicityAtleast1Btag_SingleLepton"]->Fill(selectedJets.size(),datasets[d], true, Luminosity*scaleFactor);
										}
							}
								
							//continuing for the selection 	
							for(unsigned int j=0;j<selectedJets.size();j++)
							{
								//now require at least a b-tagged jet larger than a certain pre-defined cut
								if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
								{
									selecTableSemiMu.Fill(d,5,scaleFactor); 
									if(mets[0]->Et()> METCut)
									{
										selecTableSemiMu.Fill(d,6,scaleFactor);
										eventSelected = true; 

										//cout << "event is selected according to the baseline selection!" << endl;
										
										////for single muon require exactly 1 muon, veto for other loose muons and veto for very loose electrons
										if(selectedMuons.size() == 1 && selectedLooseMuons.size() == selectedMuons.size() && selectedLooseElectronsNoVBTFid.size() == 0)
										{
											isSingleLepton = true; // we have a single muon event
											isSingleMuon = true;
											//cout << "is single muon!" << endl;
											if(dataSetName.find("TTbar")==0) histo1D["LeptonPt_TTbar"]->Fill(selectedMuons[0]->Pt(),scaleFactor);	
											if(dataSetName.find("NP_Tprime500")==0) histo1D["LeptonPt_Tprime500"]->Fill(selectedMuons[0]->Pt(),scaleFactor);
											if(dataSetName.find("NP_Bprime500")==0) histo1D["LeptonPt_Bprime500"]->Fill(selectedMuons[0]->Pt(),scaleFactor);							
											if(dataSetName.find("NP_SBprime500")==0) histo1D["LeptonPt_SBprime500"]->Fill(selectedMuons[0]->Pt(),scaleFactor);							
			     					}
			     					
										////for same-sign muons require exactly 2 muons, veto for very loose electrons
										else if(selectedMuons.size() == 2 && selectedLooseMuons.size() == selectedMuons.size() && selectedLooseElectronsNoVBTFid.size() == 0)
										{ 
											//require that there are no two muons forming the Z mass
											if( !selection.foundZCandidate(selectedMuons, selectedLooseMuons, 10.) )
											{
												//require the same charge
												if(selectedMuons[0]->charge()== selectedMuons[1]->charge())
												{
													isSSLepton = true; // we have two same-sign muons
													isSSMuon = true;
													//cout << "is same-sign muon!" << endl;
													if(dataSetName.find("TTbar")==0){
														histo1D["LeptonPt_TTbar"]->Fill(selectedMuons[0]->Pt(),scaleFactor);
														histo1D["LeptonPt_TTbar"]->Fill(selectedMuons[1]->Pt(),scaleFactor);
													}
													if(dataSetName.find("NP_Tprime500")==0){
														histo1D["LeptonPt_Tprime500"]->Fill(selectedMuons[0]->Pt(),scaleFactor);
														histo1D["LeptonPt_Tprime500"]->Fill(selectedMuons[1]->Pt(),scaleFactor);
													}
													if(dataSetName.find("NP_Bprime500")==0){
														histo1D["LeptonPt_Bprime500"]->Fill(selectedMuons[0]->Pt(),scaleFactor);
														histo1D["LeptonPt_Bprime500"]->Fill(selectedMuons[1]->Pt(),scaleFactor);
													}	  	
													if(dataSetName.find("NP_SBprime500")==0){
														histo1D["LeptonPt_SBprime500"]->Fill(selectedMuons[0]->Pt(),scaleFactor);
														histo1D["LeptonPt_SBprime500"]->Fill(selectedMuons[1]->Pt(),scaleFactor);
													}	  	
												}
											}
										}
										
										////for same-sign muon and electron require exactly 1 muon, veto for other loose muons
										else if(selectedMuons.size() == 1 && selectedLooseMuons.size() == selectedMuons.size())
										{
											//require exactly 1 electron
											if(selectedElectrons.size() == 1 && selectedLooseElectronsVBTFid.size() == selectedElectrons.size())
											{
				    						//require that there are no two electrons forming the Z mass
												if(!selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid,10.))
												{
													//it should not be an electron from a conversion!
													if( selection.passConversionRejection(selectedElectrons[0]) )
													{
														//require the same charge for muon and electron
														if(selectedElectrons[0]->charge()== selectedMuons[0]->charge())
														{
															isSSLepton = true; // we have a same-sign electron and muon
															isSSMuEl = true;
															//cout << "is same-sign muon+electron!" << endl;
															if(dataSetName.find("TTbar")==0){
																histo1D["LeptonPt_TTbar"]->Fill(selectedMuons[0]->Pt(),scaleFactor);
																histo1D["LeptonPt_TTbar"]->Fill(selectedElectrons[0]->Pt(),scaleFactor);	
															}
															if(dataSetName.find("NP_Tprime500")==0){
																histo1D["LeptonPt_Tprime500"]->Fill(selectedMuons[0]->Pt(),scaleFactor);
																histo1D["LeptonPt_Tprime500"]->Fill(selectedElectrons[0]->Pt(),scaleFactor);
															}
															if(dataSetName.find("NP_Bprime500")==0){
																histo1D["LeptonPt_Bprime500"]->Fill(selectedMuons[0]->Pt(),scaleFactor);
																histo1D["LeptonPt_Bprime500"]->Fill(selectedElectrons[0]->Pt(),scaleFactor);
															}	  	
															if(dataSetName.find("NP_SBprime500")==0){
																histo1D["LeptonPt_SBprime500"]->Fill(selectedMuons[0]->Pt(),scaleFactor);
																histo1D["LeptonPt_SBprime500"]->Fill(selectedElectrons[0]->Pt(),scaleFactor);
															}	  	
														}
													}
												}
											}
										}
										
										////three leptons
										
										//require at least 3 muons, veto for very loose electrons
										else if(selectedMuons.size() >= 3  && selectedLooseMuons.size() == selectedMuons.size() && selectedLooseElectronsNoVBTFid.size() == 0)
										{
											//require that there are no two muons forming the Z mass
											if( !selection.foundZCandidate(selectedMuons, selectedMuons, 10.) )
											{
												isTriLepton = true; // at least three muons
												isTriMuon = true;
												//cout << "is trilepton: 3 muons!" << endl;
												if(dataSetName.find("TTbar")==0){
													histo1D["LeptonPt_TTbar"]->Fill(selectedMuons[0]->Pt(),scaleFactor);
													histo1D["LeptonPt_TTbar"]->Fill(selectedMuons[1]->Pt(),scaleFactor);
													histo1D["LeptonPt_TTbar"]->Fill(selectedMuons[2]->Pt(),scaleFactor);
												}
												if(dataSetName.find("NP_Tprime500")==0){
													histo1D["LeptonPt_Tprime500"]->Fill(selectedMuons[0]->Pt(),scaleFactor);
													histo1D["LeptonPt_Tprime500"]->Fill(selectedMuons[1]->Pt(),scaleFactor);
													histo1D["LeptonPt_Tprime500"]->Fill(selectedMuons[2]->Pt(),scaleFactor);
												}
												if(dataSetName.find("NP_Bprime500")==0){
													histo1D["LeptonPt_Bprime500"]->Fill(selectedMuons[0]->Pt(),scaleFactor);
													histo1D["LeptonPt_Bprime500"]->Fill(selectedMuons[1]->Pt(),scaleFactor);
													histo1D["LeptonPt_Bprime500"]->Fill(selectedMuons[2]->Pt(),scaleFactor);
												}	  	
												if(dataSetName.find("NP_SBprime500")==0){
													histo1D["LeptonPt_SBprime500"]->Fill(selectedMuons[0]->Pt(),scaleFactor);
													histo1D["LeptonPt_SBprime500"]->Fill(selectedMuons[1]->Pt(),scaleFactor);
													histo1D["LeptonPt_SBprime500"]->Fill(selectedMuons[2]->Pt(),scaleFactor);
												}	  	
											}
					 					}
										
										//require at least 2 muons
										else if(selectedMuons.size() >= 2  && selectedLooseMuons.size() == selectedMuons.size() )
										{
											//require that there are no two muons forming the Z mass
											if( !selection.foundZCandidate(selectedMuons, selectedMuons, 10.) )
											{
												//require exactly 1 electron
												if(selectedElectrons.size() == 1 && selectedLooseElectronsVBTFid.size() == selectedElectrons.size())
												{
													//require that there are no two electrons forming the Z mass
													if(!selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid,10.))
													{
														//it should not be an electron from a conversion!
														if( selection.passConversionRejection(selectedElectrons[0]) )
														{
															isTriLepton = true; //at least two muons and one electron
															isTriMu2El1 = true;
															//cout << "is trilepton: 2 muons + 1 electron!" << endl;
															if(dataSetName.find("TTbar")==0){
																histo1D["LeptonPt_TTbar"]->Fill(selectedMuons[0]->Pt(),scaleFactor);
																histo1D["LeptonPt_TTbar"]->Fill(selectedMuons[1]->Pt(),scaleFactor);
																histo1D["LeptonPt_TTbar"]->Fill(selectedElectrons[0]->Pt(),scaleFactor);
															}
															if(dataSetName.find("NP_Tprime500")==0){
																histo1D["LeptonPt_Tprime500"]->Fill(selectedMuons[0]->Pt(),scaleFactor);
																histo1D["LeptonPt_Tprime500"]->Fill(selectedMuons[1]->Pt(),scaleFactor);
																histo1D["LeptonPt_Tprime500"]->Fill(selectedElectrons[0]->Pt(),scaleFactor);
															}
															if(dataSetName.find("NP_Bprime500")==0){
																histo1D["LeptonPt_Bprime500"]->Fill(selectedMuons[0]->Pt(),scaleFactor);
																histo1D["LeptonPt_Bprime500"]->Fill(selectedMuons[1]->Pt(),scaleFactor);
																histo1D["LeptonPt_Bprime500"]->Fill(selectedElectrons[0]->Pt(),scaleFactor);
															}	  	
															if(dataSetName.find("NP_SBprime500")==0){
																histo1D["LeptonPt_SBprime500"]->Fill(selectedMuons[0]->Pt(),scaleFactor);
																histo1D["LeptonPt_SBprime500"]->Fill(selectedMuons[1]->Pt(),scaleFactor);
																histo1D["LeptonPt_SBprime500"]->Fill(selectedElectrons[0]->Pt(),scaleFactor);
															}	  	
														}
													}
												}
											}
										}
										
										//require exactly 1 muon
										else if(selectedMuons.size() == 1 && selectedLooseMuons.size() == selectedMuons.size() )
										{
											//require at least 2 electrons
											if(selectedElectrons.size() >= 2)
											{
												//require that there are no two electrons forming the Z mass
												if(!selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid, 10.))
												{
													if(selection.passConversionRejection(selectedElectrons[0]) && selection.passConversionRejection(selectedElectrons[1]))
													{
														isTriLepton = true; //one muon and at least two electrons
														isTriMu1El2 = true;
														//cout << "is trilepton: 1 muon + 2 electrons!" << endl;
														if(dataSetName.find("TTbar")==0){
															histo1D["LeptonPt_TTbar"]->Fill(selectedMuons[0]->Pt(),scaleFactor);
															histo1D["LeptonPt_TTbar"]->Fill(selectedElectrons[0]->Pt(),scaleFactor);
															histo1D["LeptonPt_TTbar"]->Fill(selectedElectrons[1]->Pt(),scaleFactor);
														}
														if(dataSetName.find("NP_Tprime500")==0){
															histo1D["LeptonPt_Tprime500"]->Fill(selectedMuons[0]->Pt(),scaleFactor);
															histo1D["LeptonPt_Tprime500"]->Fill(selectedElectrons[0]->Pt(),scaleFactor);
															histo1D["LeptonPt_Tprime500"]->Fill(selectedElectrons[1]->Pt(),scaleFactor);
														}
														if(dataSetName.find("NP_Bprime500")==0){
															histo1D["LeptonPt_Bprime500"]->Fill(selectedMuons[0]->Pt(),scaleFactor);
															histo1D["LeptonPt_Bprime500"]->Fill(selectedElectrons[0]->Pt(),scaleFactor);
															histo1D["LeptonPt_Bprime500"]->Fill(selectedElectrons[1]->Pt(),scaleFactor);
														}	  	
														if(dataSetName.find("NP_SBprime500")==0){
															histo1D["LeptonPt_SBprime500"]->Fill(selectedMuons[0]->Pt(),scaleFactor);
															histo1D["LeptonPt_SBprime500"]->Fill(selectedElectrons[0]->Pt(),scaleFactor);
															histo1D["LeptonPt_SBprime500"]->Fill(selectedElectrons[1]->Pt(),scaleFactor);
														}	  	
													}
												}
											}
										}
									} // end MET cut
								} // end requirement of at least a b-tagged jet larger than a certain pre-defined cut																		
							} // end 'loop' on jets
						} //end 'at least one jet'  
          } // end if selectedMuons.size()>=1
        } // end good PV
      }// end trigged & semiMuon
      
			///// EVENTS TRIGGERED BY ELECTRON TRIGGER
			else if(trigged && semiElectron)
      {
        selecTableSemiEl.Fill(d,2,scaleFactor);
    	 	
        if( isGoodPV )
        {
          selecTableSemiEl.Fill(d,3,scaleFactor);
          if( selectedElectrons.size() >= 1 )
          {
            selecTableSemiEl.Fill(d,4,scaleFactor);
            if( selectedLooseMuons.size() == 0 ) 
            {
              selecTableSemiEl.Fill(d,5,scaleFactor);
              if( selection.passConversionRejection(selectedElectrons[0]) )
              {
								selecTableSemiEl.Fill(d,6,scaleFactor);
								sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class

								if( selectedJets.size()>=(unsigned int)anaEnv.NofJets)
								{
									//block for the jet multiplicity plot
									if(mets[0]->Et()> METCut)
									{							
										if(selectedElectrons.size() == 1 && !selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid, 10.))
										{
										  int nBtags = 0;
											for(unsigned int j=0;j<selectedJets.size();j++)
											{
											   if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
												 {
												    nBtags++;
												 }
											}
											MSPlot["MS_JetMultiplicity_SingleLepton"]->Fill(selectedJets.size(),datasets[d], true, Luminosity*scaleFactor);
											MSPlot["MS_BtaggedJetMultiplicity_SingleLepton"]->Fill(nBtags,datasets[d], true, Luminosity*scaleFactor);
											if(nBtags>0)
											  MSPlot["MS_JetMultiplicityAtleast1Btag_SingleLepton"]->Fill(selectedJets.size(),datasets[d], true, Luminosity*scaleFactor);
										}
									}
									
									for(unsigned int j=0;j<selectedJets.size();j++)
									{
										//now require at least a b-tagged jet larger than a certain pre-defined cut
										if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
										{
		             			selecTableSemiEl.Fill(d,7,scaleFactor);
											if(mets[0]->Et()> METCut)
											{
			        					selecTableSemiEl.Fill(d,8,scaleFactor);
                        eventSelected = true;
		
												//cout << "event is selected according to the baseline selection!" << endl;
														
												//// single electron
												if(selectedElectrons.size() == 1 && !selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid, 10.))
												{
													isSingleLepton = true;
													isSingleElectron = true;
													//cout << "is single electron!" << endl;
													if(dataSetName.find("TTbar")==0) histo1D["LeptonPt_TTbar"]->Fill(selectedElectrons[0]->Pt(),scaleFactor);
													if(dataSetName.find("NP_Tprime500")==0) histo1D["LeptonPt_Tprime500"]->Fill(selectedElectrons[0]->Pt(),scaleFactor);	
													if(dataSetName.find("NP_Bprime500")==0) histo1D["LeptonPt_Bprime500"]->Fill(selectedElectrons[0]->Pt(),scaleFactor);
													if(dataSetName.find("NP_SBprime500")==0) histo1D["LeptonPt_SBprime500"]->Fill(selectedElectrons[0]->Pt(),scaleFactor);			
												}
												
												//// two same-sign electrons
												else if(selectedElectrons.size() == 2 && !selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid, 10.))
												{
													if(selection.passConversionRejection(selectedElectrons[1]))
													{
														if(selectedElectrons[0]->charge()== selectedElectrons[1]->charge())
														{
															isSSLepton = true;
															isSSElectron = true;
															//cout << "is same-sign electron!" << endl;
															if(dataSetName.find("TTbar")==0){
																histo1D["LeptonPt_TTbar"]->Fill(selectedElectrons[0]->Pt(),scaleFactor);
																histo1D["LeptonPt_TTbar"]->Fill(selectedElectrons[1]->Pt(),scaleFactor);
															}
															if(dataSetName.find("NP_Tprime500")==0){
																histo1D["LeptonPt_Tprime500"]->Fill(selectedElectrons[0]->Pt(),scaleFactor);
																histo1D["LeptonPt_Tprime500"]->Fill(selectedElectrons[1]->Pt(),scaleFactor);
															}
															if(dataSetName.find("NP_Bprime500")==0){
																histo1D["LeptonPt_Bprime500"]->Fill(selectedElectrons[0]->Pt(),scaleFactor);
																histo1D["LeptonPt_Bprime500"]->Fill(selectedElectrons[1]->Pt(),scaleFactor);
															}	  	
															if(dataSetName.find("NP_SBprime500")==0){
																histo1D["LeptonPt_SBprime500"]->Fill(selectedElectrons[0]->Pt(),scaleFactor);
																histo1D["LeptonPt_SBprime500"]->Fill(selectedElectrons[1]->Pt(),scaleFactor);
															}	  	
														}
													}
												}
												
												//// three electrons
												else if(selectedElectrons.size() >= 3 && !selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid, 10.))
												{
													if(selection.passConversionRejection(selectedElectrons[1]) && selection.passConversionRejection(selectedElectrons[2]))
													{
														isTriLepton = true; // at least three electrons
														isTriElectron = true;
														//cout << "is trilepton: 3 electrons!" << endl;
														if(dataSetName.find("TTbar")==0){
															histo1D["LeptonPt_TTbar"]->Fill(selectedElectrons[0]->Pt(),scaleFactor);
															histo1D["LeptonPt_TTbar"]->Fill(selectedElectrons[1]->Pt(),scaleFactor);
															histo1D["LeptonPt_TTbar"]->Fill(selectedElectrons[2]->Pt(),scaleFactor);
														}
														if(dataSetName.find("NP_Tprime500")==0){
															histo1D["LeptonPt_Tprime500"]->Fill(selectedElectrons[0]->Pt(),scaleFactor);
															histo1D["LeptonPt_Tprime500"]->Fill(selectedElectrons[1]->Pt(),scaleFactor);
															histo1D["LeptonPt_Tprime500"]->Fill(selectedElectrons[2]->Pt(),scaleFactor);
														}
														if(dataSetName.find("NP_Bprime500")==0){
															histo1D["LeptonPt_Bprime500"]->Fill(selectedElectrons[0]->Pt(),scaleFactor);
															histo1D["LeptonPt_Bprime500"]->Fill(selectedElectrons[1]->Pt(),scaleFactor);
															histo1D["LeptonPt_Bprime500"]->Fill(selectedElectrons[2]->Pt(),scaleFactor);
														}	  	
														if(dataSetName.find("NP_SBprime500")==0){
															histo1D["LeptonPt_SBprime500"]->Fill(selectedElectrons[0]->Pt(),scaleFactor);
															histo1D["LeptonPt_SBprime500"]->Fill(selectedElectrons[1]->Pt(),scaleFactor);
															histo1D["LeptonPt_SBprime500"]->Fill(selectedElectrons[2]->Pt(),scaleFactor);
														}	  	
													}
												}
											} // end MET cut
										} // end requirement of at least a b-tagged jet larger than a certain pre-defined cut
									} // end 'loop' on jets
								} // end 'at least one jet'
							} // end conversion rejection for leading electron
            } // end veto muons
          } // end if selectedElectrons.size()>=1
        } // end good PV
      } // end trigged & semiElectron


						
      if(!isSingleLepton && !isSSLepton && !isTriLepton) continue; //same as all cuts just above (baseline selection is there) 

			MSPlot["MS_nPV"]->Fill(vertex.size(),datasets[d], true, Luminosity*scaleFactor);				

			if(isSingleLepton){
				MSPlot["MS_MET"]->Fill(mets[0]->Et(),datasets[d], true, Luminosity*scaleFactor);				
				if(semiElectron) MSPlot["MS_LeptonPt"]->Fill(selectedElectrons[0]->Pt(),datasets[d], true, Luminosity*scaleFactor);				
				if(semiMuon) MSPlot["MS_LeptonPt"]->Fill(selectedMuons[0]->Pt(),datasets[d], true, Luminosity*scaleFactor);				
			
				for(unsigned int j=0;j<selectedJets.size();j++)
				{
				  MSPlot["MS_JetPt_all_SingleLepton"]->Fill(selectedJets[j]->Pt(),datasets[d], true, Luminosity*scaleFactor);
					if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
					  MSPlot["MS_JetPt_btagged_SingleLepton"]->Fill(selectedJets[j]->Pt(),datasets[d], true, Luminosity*scaleFactor);
					else
					  MSPlot["MS_JetPt_nonbtagged_SingleLepton"]->Fill(selectedJets[j]->Pt(),datasets[d], true, Luminosity*scaleFactor);
				}			
			}

			if(isSSLepton)
			{
				//cout << "IS SAME-SIGN LEPTON EVENT" << endl;
				NbSSevents = NbSSevents + datasets[d]->NormFactor()*Luminosity*scaleFactor;
				if(semiElectron) selecTableSemiEl.Fill(d,9,scaleFactor);
				if(semiMuon) selecTableSemiMu.Fill(d,7,scaleFactor);				
			}
			if(isTriLepton)
			{
				//cout << "IS TRI-LEPTON EVENT" << endl;
				NbTrievents = NbTrievents + datasets[d]->NormFactor()*Luminosity*scaleFactor;
				if(semiElectron) selecTableSemiEl.Fill(d,10,scaleFactor);
				if(semiMuon) selecTableSemiMu.Fill(d,8,scaleFactor);				
			}

			//////////////////////////////////////////////////////////////////////////
      // MVA training
      //////////////////////////////////////////////////////////////////////////
      if(doMVAjetcombination && selectedJets.size()>=4 && (TrainMVA || (doKinematicFit && CalculateResolutions)) && ((dataSetName.find("TTbarJets_SemiMu") == 0 && semiMuon) || (dataSetName.find("TTbarJets_SemiElectron") == 0 && semiElectron))) //otherwise, if the jets vector only has 2 jets selected the jetcombiner crashes... 
      {
				if(!isSingleLepton) continue;
	
        			TRootGenEvent* genEvt = treeLoader.LoadGenEvent(ievt,false);
        			sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class 
				if(semiMuon) 
					jetCombiner->ProcessEvent(datasets[d],mcParticles,selectedJets,selectedMuons[0],init_electrons,init_muons,genEvt,scaleFactor,TrainwithTprime);	        
				else if(semiElectron)
					jetCombiner->ProcessEvent(datasets[d],mcParticles,selectedJets,selectedElectrons[0],init_electrons,init_muons,genEvt,scaleFactor,TrainwithTprime);
      }

       
      if(TrainMVA) continue; //for the training, only the jetcombiner is relevant, so the following can be skipped (to the next event in the event loop)
      
			
      //relevant for the kinematic fit: calculate the resolutions
      if(doMVAjetcombination && selectedJets.size()>=4 && doKinematicFit && CalculateResolutions && ((dataSetName.find("TTbarJets_SemiMu") == 0 && semiMuon) || (dataSetName.find("TTbarJets_SemiElectron") == 0 && semiElectron)))
      {
          jetCombiner->FillResolutions(resFitLightJets, resFitBJets);
	  			continue;
      }
      

			//////////////////////////////////////////////////////////////////////////
      // jet-parton matching needed to make the Wmassplot for the W counting
      //////////////////////////////////////////////////////////////////////////
      if(selectedJets.size()>=4 && !useMassesAndResolutions && ((dataSetName.find("TTbarJets_SemiMu") == 0 && semiMuon) || (dataSetName.find("TTbarJets_SemiElectron") == 0 && semiElectron)))
			{
				sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class)     
				bool Wbosonpartonsmatched = false; // True if the Wboson ttbar semi-mu or semi-el partons are matched to their 2 jets (not necessarily the 4 highest pt jets)
				int MCPermutation[4]; for (unsigned int i=0;i<4;i++) MCPermutation[i] = -1;

				pair<unsigned int, unsigned int> leptonicBJet_, hadronicBJet_, hadronicWJet1_, hadronicWJet2_; //First index is the JET number, second one is the parton
				vector<TRootMCParticle*> mcParticlesMatching_; 
				vector<TRootJet*> selectedJets_; 

				leptonicBJet_ = hadronicBJet_ = hadronicWJet1_ = hadronicWJet2_ = pair<unsigned int,unsigned int>(9999,9999);
				mcParticlesMatching_.clear();
				selectedJets_.clear();

				if((dataSetName.find("TTbarJets_SemiMu") == 0 && semiMuon) || (dataSetName.find("TTbarJets_SemiElectron") == 0 && semiElectron))
				{
      
					vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV;
            
					bool muPlusFromTop = false, muMinusFromTop = false;
					bool elPlusFromTop = false, elMinusFromTop = false;
					int leptonPDG, muonPDG = 13, electronPDG = 11;
					if(semiMuon)
					  leptonPDG = muonPDG;
					else if(semiElectron)
					  leptonPDG = electronPDG; 
					for(unsigned int i=0; i<mcParticles.size(); i++) {
						if( mcParticles[i]->status() != 3) continue;
						
						if( mcParticles[i]->type() == leptonPDG && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -6 )
					 	{
							//if(muMinusFromTop) cerr<<"muMinusFromTop was already true"<<endl;
							if(leptonPDG==muonPDG) muMinusFromTop = true;
							else if(leptonPDG==electronPDG) elMinusFromTop = true;
						}
						if( mcParticles[i]->type() == -leptonPDG && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == 6 )
						{
							//if(muPlusFromTop) cerr<<"muPlusFromTop was already true"<<endl;
							if(leptonPDG==muonPDG) muPlusFromTop = true;
							else if(leptonPDG==electronPDG) elPlusFromTop = true;
						}
						
						if( abs(mcParticles[i]->type()) < 6 || abs(mcParticles[i]->type()) == 21 )
						{
							mcParticlesTLV.push_back(*mcParticles[i]);
							mcParticlesMatching_.push_back(mcParticles[i]);
						}
					}
	
					//if(muPlusFromTop && muMinusFromTop)
	  				//cerr<<"muPlusFromTop and muMinusFromTop are both true ?!\nCheck if you are using the right sample..."<<endl;
	
					for(unsigned int i=0; i<4; i++)
	  				  selectedJetsTLV.push_back(*selectedJets[i]);
	
					JetPartonMatching matching = JetPartonMatching(mcParticlesTLV, selectedJetsTLV, 2, true, true, 0.3);
	
					if(matching.getNumberOfAvailableCombinations() != 1)
	  				  cerr << "matching.getNumberOfAvailableCombinations() = "<<matching.getNumberOfAvailableCombinations()<<"  This should be equal to 1 !!!"<<endl;
	
					vector< pair<unsigned int, unsigned int> > JetPartonPair; // First one is jet number, second one is mcParticle number
	
					for(unsigned int i=0; i<mcParticlesTLV.size(); i++) 
					{
	  				  int matchedJetNumber = matching.getMatchForParton(i, 0);
	  				  if(matchedJetNumber != -1)
	    				    JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
					}

					for(unsigned int i=0; i<JetPartonPair.size(); i++)
					{
	  				  unsigned int j = JetPartonPair[i].second;	  
	  				  if( fabs(mcParticlesMatching_[j]->type()) < 6 )
	    			          {
	      			            if( ( ((muPlusFromTop && semiMuon) || (elPlusFromTop && semiElectron)) && mcParticlesMatching_[j]->motherType() == -24 && mcParticlesMatching_[j]->grannyType() == -6 )
		  				|| ( ((muMinusFromTop && semiMuon) || (elMinusFromTop && semiElectron)) && mcParticlesMatching_[j]->motherType() == 24 && mcParticlesMatching_[j]->grannyType() == 6 ) )
					    {
		  				if(hadronicWJet1_.first == 9999)
						{
		    					hadronicWJet1_ = JetPartonPair[i];
		    					MCPermutation[0] = JetPartonPair[i].first;
		  				}
						else if(hadronicWJet2_.first == 9999)
						{
		   	 				hadronicWJet2_ = JetPartonPair[i];
		    					MCPermutation[1] = JetPartonPair[i].first;
		  				} 
						else cerr<<"Found a third jet coming from a W boson which comes from a top quark..."<<endl;
					    }
	    			          }
	  				  if( fabs(mcParticlesMatching_[j]->type()) == 5 )
	    			          {
	      			            if( ( ((muPlusFromTop && semiMuon) || (elPlusFromTop && semiElectron)) && mcParticlesMatching_[j]->motherType() == -6 )
		    			|| ( ((muMinusFromTop && semiMuon) || (elMinusFromTop && semiElectron)) && mcParticlesMatching_[j]->motherType() == 6 ) )
					    {
								hadronicBJet_ = JetPartonPair[i];
								MCPermutation[2] = JetPartonPair[i].first;
	      				     }
					     else if( ( ((muPlusFromTop && semiMuon) || (elPlusFromTop && semiElectron)) && mcParticlesMatching_[j]->motherType() == 6 )
			 				|| ( ((muMinusFromTop && semiMuon) || (elMinusFromTop && semiElectron)) && mcParticlesMatching_[j]->motherType() == -6 ) )
					     {
			 					leptonicBJet_ = JetPartonPair[i];
			 					MCPermutation[3] = JetPartonPair[i].first;
	      				     }
	    				  }
					}					
	
					if(hadronicWJet1_.first != 9999 && hadronicWJet2_.first != 9999)
					{
	  				  histo1D["hadronicPartonWMass"]->Fill((*mcParticlesMatching_[hadronicWJet1_.second]+*mcParticlesMatching_[hadronicWJet2_.second]).M(),scaleFactor);	  
	  				  Wbosonpartonsmatched = true;	  
					}
				} //if dataset ttbar

      	if(Wbosonpartonsmatched)
	  		{
					float WMassmatched_ = (*selectedJets[hadronicWJet1_.first]+*selectedJets[hadronicWJet2_.first]).M();
      		histo1D["hadronicRecoWMass"]->Fill(WMassmatched_,scaleFactor);
	  		}
			} //end selectedJets.size()>=4 !useMassesAndResolutions && dataSetName.find("TTbarJets_SemiMu") == 0

      
			

			/////////////////////////////////////////////////////////////////////////////
			//// core of the analysis, need to check b-tagging and reconstruct W bosons!	
			////////////////////////////////////////////////////////////////////////////
      if(useMassesAndResolutions)
      {
				vector<TRootJet*> selectedJets_MVAinput; //for the MVA jet combination for the mass reconstruction in the 1B_2W and 2B_2W boxes
			
				//////////////////////////
				//b-tagging stuff here
				//////////////////////////
				vector<TRootJet*> selectedJetsForBtagging; //need the jets within the tracker acceptance
				selectedJetsForBtagging = selection.GetSelectedJets(true);
			
				vector< pair< int, float > > jetindex_btagvalue;
				vector< pair< int, bool > > jetindex_isb;
				for(unsigned int i = 0; i<selectedJetsForBtagging.size(); i++)
				{
					pair<int,float> dummy (i,selectedJetsForBtagging[i]->btag_trackCountingHighPurBJetTags());
					jetindex_btagvalue.push_back(dummy);
				}
				/*cout<<"BEFORE SORTING"<<endl;
				for(int k=0;k<jetindex_btagvalue.size();k++)
				{
					cout<<" jetindex_btagvalue["<<k<<"].first = "<<jetindex_btagvalue[k].first<<", jetindex_btagvalue["<<k<<"].second = "<<jetindex_btagvalue[k].second<<endl;
				}*/
 				
				// this is sorted according to the second value of the pair and not the first
 				std::sort(jetindex_btagvalue.begin(), jetindex_btagvalue.end(), sort_pair_decreasing());

				// check if jet is truly from a b-quark or not and make a pair with the result
				if(! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ))
				{
					for(unsigned int i = 0; i<jetindex_btagvalue.size(); i++)
					{
						if(fabs(selectedJetsForBtagging[jetindex_btagvalue[i].first]->partonFlavour()) == 5 )
						{// is a true b-jet
							pair<int,float> dummy (jetindex_btagvalue[i].first,true); //not: pair<int,float> dummy (i,true)
							jetindex_isb.push_back(dummy);
						}
						else
						{// is not a b-jet
							pair<int,float> dummy (jetindex_btagvalue[i].first,false); //not: pair<int,float> dummy (i,false)
							jetindex_isb.push_back(dummy);							
						}
					}
				}
							

				//////////////////
				// rescale events according to b-tag/mistag eff scalefactors
				//////////////////
      	double HT = 0.;
      	//cout << "begin processing new selected event with current HT: " << HT << endl;
				HT = HT + mets[0]->Et();

				int nbOfBtags = 0;
				int bJet1 = 999;
				int bJet2 = 999;
				for(unsigned int i = 0; i<jetindex_btagvalue.size(); i++)
				{
			  	//cout<<" jetindex_btagvalue["<<i<<"].first = "<<jetindex_btagvalue[i].first<<", jetindex_btagvalue["<<i<<"].second = "<<jetindex_btagvalue[i].second<<", selectedJetsForBtagging[jetindex_btagvalue["<<i<<"].first]->partonFlavour() = "<<selectedJetsForBtagging[jetindex_btagvalue[i].first]->partonFlavour()<<endl;
					if(jetindex_btagvalue[i].second > workingpointvalue)
					{
						nbOfBtags++;
						if(nbOfBtags<=2)
						{
					  	HT = HT + selectedJetsForBtagging[jetindex_btagvalue[i].first]->Pt(); //only add pt of b-tagged jet for the 2 jets with the highest b-tag value (if more than 2 b-tagged jets, only the 2 highest are taken into account)
					  	selectedJets_MVAinput.push_back((TRootJet*) selectedJetsForBtagging[jetindex_btagvalue[i].first]->Clone());
						}
						if(nbOfBtags == 1)
						{
							bJet1 = jetindex_btagvalue[i].first;
							for(unsigned int j = 0; j < jetindex_isb.size(); j++)
							{
								if(bJet1 == jetindex_isb[j].first)
								{
									if(!jetindex_isb[j].second)
									{ // jet is b-tagged but NOT a true b
										scaleFactor = scaleFactor * mistagfactor;
									}
									else
									{ // jet is b-tagged and a true b
										scaleFactor = scaleFactor * scalefactorbtageff;
									}										
								}
							}
						}
						if(nbOfBtags == 2)
						{
							bJet2 = jetindex_btagvalue[i].first;
							for(unsigned int j = 0; j < jetindex_isb.size(); j++)
							{								
								if(bJet2 == jetindex_isb[j].first)
								{       
							  	//cout<<"    jetindex_isb["<<j<<"].second = "<<jetindex_isb[j].second<<endl;
									if(!jetindex_isb[j].second)
									{ // jet is b-tagged but NOT a true b
										//cout<<"jet is b-tagged but NOT a true b"<<endl;
										scaleFactor = scaleFactor * mistagfactor;
									} 
									else
									{ // jet is b-tagged and a true b
								  	//cout<<"jet is b-tagged and a true b"<<endl;
										scaleFactor = scaleFactor * scalefactorbtageff;
									}									
								}
							}
						}
					}
				} //end loop over jets (~ b-tagging)
    
			///////////////////////
			//copy the collection of selected jets except for the one or two b-tagged jets (used in the HT calculation)
			//this collection is used for W counting				
			//////////////////////
			for(unsigned int i = 0; i<selectedJetsForBtagging.size(); i++)
			{
				if((int)i != bJet1 && (int)i != bJet2)
				{
					selectedJetsFromW.push_back((TRootJet*) selectedJetsForBtagging[i]->Clone());
				}
			}
	
			///////////////////////////////
			// start the W counting procedure here
			///////////////////////////////
			int nbOfWs = 0; 
			if(isSingleLepton) nbOfWs = 1;
			if(isSSLepton) nbOfWs = 2;
			if(isTriLepton) nbOfWs = 3;
			
			if(nbOfWs == 0){ cout << "WARNING - no W counted yet, this can NOT be!!!" << endl; break;}
			
			float recoWmass = 0.; float newrecoWmass = 0.;
			int indexWjet1 = 999;
			int indexWjet2 = 999;
			float massdifference = 100;
			int previoussize = 0;
			
			for(unsigned int i = 0; i<selectedJetsFromW.size(); i++)
			{
				selectedJetsFromW_DropUsedJets.push_back((TRootJet*) selectedJetsFromW[i]->Clone());
			}
			
			int firstWboson_jet1 = 9999; int firstWboson_jet2 = 9999;
			bool firsttime = true;
			do
			{
				previoussize = selectedJetsFromW_DropUsedJets.size();
			  	indexWjet1 = 999; indexWjet2 = 999;
				massdifference = 100;
				for(unsigned int i = 0; i<selectedJetsFromW_DropUsedJets.size(); i++)
				{
					for(unsigned int j = 0; j<i; j++)
					{
						recoWmass = (*selectedJetsFromW_DropUsedJets[i]+*selectedJetsFromW_DropUsedJets[j]).M();
						if(!(dataSetName.find("NP_")<=dataSetName.size()))
 							MSPlot["allDiJetMasses"]->Fill(recoWmass, datasets[d], true, Luminosity*scaleFactor);
						if( fabs(recoWmass-Wmass)<massdifference)
						{
							massdifference = fabs(recoWmass-Wmass);
							newrecoWmass = recoWmass;
							indexWjet1 = i;
							indexWjet2 = j;
						}
					}
				}
				
				// fill a plot with the W mass from the chosen jets (only for the best jet couple, even if there are more W's reconstructed)
				if(firsttime && indexWjet1<999 && indexWjet2<999)
				{
					firsttime = false;
					firstWboson_jet1 = indexWjet1;
					firstWboson_jet2 = indexWjet2;
					float chosenWMass = (*selectedJetsFromW_DropUsedJets[firstWboson_jet1]+*selectedJetsFromW_DropUsedJets[firstWboson_jet2]).M();
					MSPlot["hadronicRecoWMass_chosenWjets"]->Fill(chosenWMass, datasets[d], true, Luminosity*scaleFactor);
				}

				if(massdifference < SigmaWmass && indexWjet1<999 && indexWjet2<999)
				{
					nbOfWs++;

					HT = HT + (*selectedJetsFromW_DropUsedJets[indexWjet1]+*selectedJetsFromW_DropUsedJets[indexWjet2]).Pt();
					
					if(selectedJets_MVAinput.size()<=2) //you're only using it for the nbOfWs=2 boxes anyway...
					{
					        selectedJets_MVAinput.push_back((TRootJet*) selectedJetsFromW_DropUsedJets[indexWjet1]->Clone());
					        selectedJets_MVAinput.push_back((TRootJet*) selectedJetsFromW_DropUsedJets[indexWjet2]->Clone());
					}
					
					for(unsigned int i = 0; i<selectedJetsFromW_DropUsedJets.size(); i++)
					{
						if((unsigned int)indexWjet1!=i && (unsigned int)indexWjet2!=i)
						{
							selectedJetsFromW_DropUsedJets_tmp.push_back((TRootJet*) selectedJetsFromW_DropUsedJets[i]->Clone());
						}
					}
					selectedJetsFromW_DropUsedJets.clear();
					for(unsigned int i = 0; i<selectedJetsFromW_DropUsedJets_tmp.size(); i++)
					{
						selectedJetsFromW_DropUsedJets.push_back((TRootJet*) selectedJetsFromW_DropUsedJets_tmp[i]->Clone());
					}
					selectedJetsFromW_DropUsedJets_tmp.clear();
				}

			}
			while(selectedJetsFromW_DropUsedJets.size()>=2 && selectedJetsFromW_DropUsedJets.size()!= (unsigned int)previoussize);
				
			if(nbOfBtags>=2) nbOfBtags = 2; //to make sure that events with more than 2 b-jets end up in the 2 b-jet bin			
			
			////////////////////////////////////
			//// here we put the events into boxes of #b's and #W's
			////////////////////////////////////
			if(nbOfWs>=4) nbOfWs = 4; // to make sure that events with more than 4 W bosons end up in the 4W bin

			string histo1_dataset = "NbEvents_SSLeptons";
			histo1_dataset = histo1_dataset+(datasets[d]->Name()).c_str();
			string histo2_dataset = "NbEvents_TriLeptons";
			histo2_dataset = histo2_dataset+(datasets[d]->Name()).c_str();
			
			//cout<<"nbOfBtags = "<<nbOfBtags<<endl;						
			if(nbOfBtags==1)
			{				
				if(nbOfWs==1)
				{
					//Requirements: a leptonically decaying W, and a forward jet and a 'central' one
					if(selectedForwardJets.size() != 1) continue;
					if(selectedJets.size() != 1) continue;
				
					float DeltaPhi = fabs(selectedJets[0]->DeltaPhi(*selectedForwardJets[0]));
					float RelPt = fabs(selectedJets[0]->Pt()-selectedForwardJets[0]->Pt())/(selectedJets[0]->Pt()+selectedForwardJets[0]->Pt());
				
					double const Pi=4*atan(1);
					if(DeltaPhi > ((Pi/2)+RelPt*Pi)) continue;
					//cout << "in 1B 1W box" << endl;
					
					if(isSingleLepton) 
					{						
						HT = HT + selectedForwardJets[0]->Pt();

						myInclFourthGenSearchTools.FillPlots(d,nbOfBtags,nbOfWs,HT,selectedMuons,selectedElectrons,mets,selectedJets,scaleFactor);
						if(isSingleElectron) selecTableSemiEl.Fill(d,11,scaleFactor);
						if(isSingleMuon) selecTableSemiMu.Fill(d,9,scaleFactor);
					}
					//cout << "done in 1B 1W box" << endl;
				}
				else if(nbOfWs==2)
				{
					if(isSingleLepton && selectedJets.size()>=4)
					{
						//cout << "in 1B 2W box" << endl;
						HT = HT + selectedJetsFromW_DropUsedJets[0]->Pt();

						myInclFourthGenSearchTools.FillPlots(d,nbOfBtags,nbOfWs,HT,selectedMuons,selectedElectrons,mets,selectedJets,scaleFactor);
						if(isSingleElectron) selecTableSemiEl.Fill(d,12,scaleFactor);
						if(isSingleMuon) selecTableSemiMu.Fill(d,10,scaleFactor);
			
						//////////////
						// MVA stuff, only to be done when this is decided in the configuration of this macro (the doMVAjetcombination boolean)
						//////////////
						pair<float, vector<unsigned int> > MVAvals;	
						if(doMVAjetcombination)
						{	
					   //you should put the highest Pt jet in here THAT IS NOT YET IN THE VECTOR!
       			 if(selectedJets_MVAinput.size()==3)
					   {
					     selectedJets_MVAinput.push_back((TRootJet*) selectedJetsFromW_DropUsedJets[0]->Clone());
					   }
					   else
					     cout<<"Something is not right, the vector of selected jets for MVA input in this stage should be 3"<<endl;
					   
					   if(selectedJets_MVAinput.size()==4)
					   {	
					
					     sort(selectedJets_MVAinput.begin(),selectedJets_MVAinput.end(),HighestPt()); // HighestPt() is included from the Selection class
                                                
						   TRootGenEvent* genEvt = treeLoader.LoadGenEvent(ievt,false);
			   			 if(semiMuon) jetCombiner->ProcessEvent(datasets[d],mcParticles,selectedJets_MVAinput,selectedMuons[0],init_electrons,init_muons,genEvt,scaleFactor);	//OLD WAY (class has changed since then): jetCombiner->ProcessEvent(datasets[d], mcParticles, selectedJets, selectedMuons[0], vertex[0], eventSelected, init_electrons, init_muons, scaleFactor);
               else if(semiElectron) jetCombiner->ProcessEvent(datasets[d],mcParticles,selectedJets_MVAinput,selectedElectrons[0],init_electrons,init_muons,genEvt,scaleFactor);
			  
			   			 //vector<unsigned int> goodCombi = jetCombiner_1B_2W->GetGoodJetCombination(); //get the MC matched jet combination, not the MVA best matched		   	
			   			 MVAvals = jetCombiner->getMVAValue(MVAmethod, 1); // 1 means the highest MVA value
					   }
					   else
					     cout<<"WARNING: vector of selected jets for MVA input is not equal to 4 (but to "<<selectedJets_MVAinput.size()<<"); fix this!!"<<endl;					   
					
					   //reconstructing the mass of the top/t', and fill vector for binning
					   myInclFourthGenSearchTools.CalculateTopMass(selectedJets_MVAinput[MVAvals.second[0]],selectedJets_MVAinput[MVAvals.second[1]],selectedJets_MVAinput[MVAvals.second[2]]); //(light 1, light 2, hadronic b-jet);
					   myInclFourthGenSearchTools.FillMassPlots(d,nbOfBtags,nbOfWs,scaleFactor);
					   	
					   //for the 2D binning (HT vs Mtop)
					   float fillweight = datasets[d]->NormFactor () * Luminosity * scaleFactor;
					   if(anaEnv.MCRound==0)
					   {
					        HTvsMTop_1B_2W.Fill_for2DBinning(myInclFourthGenSearchTools.GetHT(),myInclFourthGenSearchTools.GetMtop(),fillweight);
					   }
					   else if(anaEnv.MCRound!=0)
					   {
					        HTvsMTop_1B_2W.Fill(myInclFourthGenSearchTools.GetHT(),myInclFourthGenSearchTools.GetMtop(),fillweight,d);
					   }
					 }
					
				    } 
				    else if(isSSLepton && selectedJets.size() >=2) 
				    {
							//cout << "IS SAME-SIGN LEPTON" << endl;
							//count number of events with SS leptons
							if(dataSetName.find("NP") != 0) NbSSevents_1B_2W = NbSSevents_1B_2W + datasets[d]->NormFactor()*Luminosity*scaleFactor;
							if(dataSetName.find("NP_Tprime500")==0 || dataSetName.find("NP_STprime500")==0 || dataSetName.find("NP_Bprime500")==0 || dataSetName.find("NP_SBprime500")==0) NP_NbSSevents_1B_2W = NP_NbSSevents_1B_2W + datasets[d]->NormFactor()*Luminosity*scaleFactor ;
							MSPlot["MS_NbSSevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);
							histo1D[histo1_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
							if(semiElectron) selecTableSemiEl.Fill(d,19,scaleFactor);
						  if(semiMuon) selecTableSemiMu.Fill(d,17,scaleFactor);	
				    }
 	
					//cout << "done in 1B 2W box" << endl;
				}
				else if(nbOfWs==3)
				{
				    if(isSingleLepton && selectedJets.size() >=6) 
				    {
				      myInclFourthGenSearchTools.FillPlots(d,nbOfBtags,nbOfWs,HT,selectedMuons,selectedElectrons,mets,selectedJets,scaleFactor);				
							if(isSingleElectron) selecTableSemiEl.Fill(d,13,scaleFactor);
							if(isSingleMuon) selecTableSemiMu.Fill(d,11,scaleFactor);
				    }
				    else if(isSSLepton && selectedJets.size() >=4) 
			      {
							//cout << "IS SAME-SIGN LEPTON" << endl;
							//count number of events with SS leptons
							if(dataSetName.find("NP") != 0) NbSSevents_1B_3W = NbSSevents_1B_3W + datasets[d]->NormFactor()*Luminosity*scaleFactor;
							if(dataSetName.find("NP_Tprime500")==0 || dataSetName.find("NP_STprime500")==0 || dataSetName.find("NP_Bprime500")==0 || dataSetName.find("NP_SBprime500")==0) NP_NbSSevents_1B_3W = NP_NbSSevents_1B_3W + datasets[d]->NormFactor()*Luminosity*scaleFactor;
							MSPlot["MS_NbSSevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);
							histo1D[histo1_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
							if(semiElectron) selecTableSemiEl.Fill(d,19,scaleFactor);
						  if(semiMuon) selecTableSemiMu.Fill(d,17,scaleFactor);
						}
				    else if(isTriLepton && selectedJets.size() >=2) 
				    {
							//cout << "IS TRI-LEPTON" << endl;
							//count number of events with three leptons
							if(dataSetName.find("NP") != 0) NbTrievents_1B_3W = NbTrievents_1B_3W + datasets[d]->NormFactor()*Luminosity*scaleFactor;
							if(dataSetName.find("NP_Tprime500")==0 || dataSetName.find("NP_STprime500")==0 || dataSetName.find("NP_Bprime500")==0 || dataSetName.find("NP_SBprime500")==0) NP_NbSSevents_1B_3W = NP_NbSSevents_1B_3W + datasets[d]->NormFactor()*Luminosity*scaleFactor;
							MSPlot["MS_NbTrievents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);
							histo1D[histo2_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
				    	if(semiElectron) selecTableSemiEl.Fill(d,20,scaleFactor);
						  if(semiMuon) selecTableSemiMu.Fill(d,18,scaleFactor);
						}
				    //cout << "done in 1B 3W box" << endl;
				}
				else
				{
				   if(isSingleLepton && selectedJets.size() >=8) 
				   {
							//cout << "in 1B 4W box" << endl;
							myInclFourthGenSearchTools.FillPlots(d,nbOfBtags,nbOfWs,HT,selectedMuons,selectedElectrons,mets,selectedJets,scaleFactor);
							if(isSingleElectron) selecTableSemiEl.Fill(d,14,scaleFactor);
							if(isSingleMuon) selecTableSemiMu.Fill(d,12,scaleFactor);
					
				   }
				   else if(isSSLepton && selectedJets.size() >=6) 
				   {
							//cout << "IS SAME-SIGN LEPTON" << endl;
							//count number of events with SS leptons
							if(dataSetName.find("NP") != 0) NbSSevents_1B_4W = NbSSevents_1B_4W + datasets[d]->NormFactor()*Luminosity*scaleFactor;
							if(dataSetName.find("NP_Tprime500")==0 || dataSetName.find("NP_STprime500")==0 || dataSetName.find("NP_Bprime500")==0 || dataSetName.find("NP_SBprime500")==0) NP_NbSSevents_1B_4W = NP_NbSSevents_1B_4W + datasets[d]->NormFactor()*Luminosity*scaleFactor;
							MSPlot["MS_NbSSevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);
							histo1D[histo1_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
			     		if(semiElectron) selecTableSemiEl.Fill(d,19,scaleFactor);
						  if(semiMuon) selecTableSemiMu.Fill(d,17,scaleFactor);
					 }
				   else if(isTriLepton && selectedJets.size() >=4) 
				   {
							//cout << "IS TRI-LEPTON" << endl;
							//count number of events with three leptons
							if(dataSetName.find("NP") != 0) NbTrievents_1B_4W = NbTrievents_1B_4W + datasets[d]->NormFactor()*Luminosity*scaleFactor;
							if(dataSetName.find("NP_Tprime500")==0 || dataSetName.find("NP_STprime500")==0 || dataSetName.find("NP_Bprime500")==0 || dataSetName.find("NP_SBprime500")==0) NP_NbTrievents_1B_4W = NP_NbTrievents_1B_4W + datasets[d]->NormFactor()*Luminosity*scaleFactor;
							MSPlot["MS_NbTrievents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);
							histo1D[histo2_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
				      if(semiElectron) selecTableSemiEl.Fill(d,20,scaleFactor);
						  if(semiMuon) selecTableSemiMu.Fill(d,18,scaleFactor);
					  }
				   //cout << "done in 1B 4W box" << endl;
				}
				
			} //end number of btags == 1
			else if(nbOfBtags == 2)
			{
				if(nbOfWs>=4) nbOfWs = 4;// to make sure that events with more than 4 W bosons end up in the 4W bin
				
				if(nbOfWs==1)
				{
				   if(selectedForwardJets.size() != 0) continue;
				   if(selectedJets.size() != 2) continue;

				   float DeltaPhi = fabs(selectedJets[0]->DeltaPhi(*selectedJets[1]));
				   float RelPt = fabs(selectedJets[0]->Pt()-selectedJets[1]->Pt())/(selectedJets[0]->Pt()+selectedJets[1]->Pt());

				   double const Pi=4*atan(1);
				   if(DeltaPhi > ((Pi/2)+RelPt*Pi)) continue;
				   
				   //cout << "in 2B 1W box" << endl;

				   if(isSingleLepton) 
				   {
				        myInclFourthGenSearchTools.FillPlots(d,nbOfBtags,nbOfWs,HT,selectedMuons,selectedElectrons,mets,selectedJets,scaleFactor);
							if(isSingleElectron) selecTableSemiEl.Fill(d,15,scaleFactor);
							if(isSingleMuon) selecTableSemiMu.Fill(d,13,scaleFactor);
				   }
					//cout << "done in 2B 1W box" << endl;
				}
				else if(nbOfWs==2)
				{
				   if(isSingleLepton && selectedJets.size() >=4) 
				   {
					//cout << "in 2B 2W box" << endl;
					
					myInclFourthGenSearchTools.FillPlots(d,nbOfBtags,nbOfWs,HT,selectedMuons,selectedElectrons,mets,selectedJets,scaleFactor);
							if(isSingleElectron) selecTableSemiEl.Fill(d,16,scaleFactor);
							if(isSingleMuon) selecTableSemiMu.Fill(d,14,scaleFactor);
					
					//////////////
					// MVA stuff, only to be done when this is decided in the configuration of this macro (the doMVAjetcombination boolean)
					//////////////
					pair<float, vector<unsigned int> > MVAvals;
					bool TprimeEvaluation = false;
					if(dataSetName.find("NP_Tprime")==0 || dataSetName.find("NP_overlay_Tprime")==0)
					  TprimeEvaluation = true;
					if(doMVAjetcombination)
					{	
					   if(selectedJets_MVAinput.size()==4)
					   {
       			   			sort(selectedJets_MVAinput.begin(),selectedJets_MVAinput.end(),HighestPt()); // HighestPt() is included from the Selection class
                           			
						TRootGenEvent* genEvt = treeLoader.LoadGenEvent(ievt,false);
			   			if(semiMuon) jetCombiner->ProcessEvent(datasets[d],mcParticles,selectedJets_MVAinput,selectedMuons[0],init_electrons,init_muons,genEvt,scaleFactor,TprimeEvaluation);	//OLD WAY (class has changed since then): jetCombiner->ProcessEvent(datasets[d], mcParticles, selectedJets, selectedMuons[0], vertex[0], eventSelected, init_electrons, init_muons, scaleFactor);
              else if(semiElectron) jetCombiner->ProcessEvent(datasets[d],mcParticles,selectedJets_MVAinput,selectedElectrons[0],init_electrons,init_muons,genEvt,scaleFactor,TprimeEvaluation);
			  
			   			//vector<unsigned int> goodCombi = jetCombiner_2B_2W->GetGoodJetCombination(); //get the MC matched jet combination, not the MVA best matched		   	
			   			MVAvals = jetCombiner->getMVAValue(MVAmethod, 1); // 1 means the highest MVA value

	   					//coutObjectsFourVector(init_muons,init_electrons,selectedJets_MVAinput,mets,"*** Before kinematic fit ***");
						
						////to test 'purity' of good jet combinations with a 'simple' method (not MVA)... for the moment only testable in 2B_2W box
					        //myInclFourthGenSearchTools.TestPurityGoodCombinations(d,nbOfBtags,nbOfWs,genEvt,mcParticles,selectedJets_MVAinput,TprimeEvaluation,scaleFactor);
					        
						
					        //reconstructing the mass of the top/t', and fill vector for binning
					   	myInclFourthGenSearchTools.CalculateTopMass(selectedJets_MVAinput[MVAvals.second[0]],selectedJets_MVAinput[MVAvals.second[1]],selectedJets_MVAinput[MVAvals.second[2]]); //(light 1, light 2, hadronic b-jet);
					        myInclFourthGenSearchTools.FillMassPlots(d,nbOfBtags,nbOfWs,scaleFactor);
					   	
						//for the 2D binning (HT vs Mtop)
					   	float fillweight = datasets[d]->NormFactor () * Luminosity * scaleFactor;
					   	if(anaEnv.MCRound==0)
					   	{
					   	  HTvsMTop_2B_2W.Fill_for2DBinning(myInclFourthGenSearchTools.GetHT(),myInclFourthGenSearchTools.GetMtop(),fillweight);
					   	}
					   	else if(anaEnv.MCRound!=0)
					   	{
					   	  HTvsMTop_2B_2W.Fill(myInclFourthGenSearchTools.GetHT(),myInclFourthGenSearchTools.GetMtop(),fillweight,d);
					   	}
					   }
					   else
					     cout<<"WARNING: vector of selected jets for MVA input is not equal to 4; fix this!!"<<endl;	
					}
			             } 
				     else if(isSSLepton && selectedJets.size() >=2) 
				     {
							//cout << "IS SAME-SIGN LEPTON" << endl;
							//count number of events with SS leptons
							if(dataSetName.find("NP") != 0) NbSSevents_2B_2W = NbSSevents_2B_2W + datasets[d]->NormFactor()*Luminosity*scaleFactor;
							if(dataSetName.find("NP_Tprime500")==0 || dataSetName.find("NP_STprime500")==0 || dataSetName.find("NP_Bprime500")==0 || dataSetName.find("NP_SBprime500")==0) NP_NbSSevents_2B_2W = NP_NbSSevents_2B_2W + datasets[d]->NormFactor()*Luminosity*scaleFactor;
							MSPlot["MS_NbSSevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);
							histo1D[histo1_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
				      if(semiElectron) selecTableSemiEl.Fill(d,19,scaleFactor);
						  if(semiMuon) selecTableSemiMu.Fill(d,17,scaleFactor);
					  } 
				}
				else if(nbOfWs==3)
				{    
				     //cout << "in 2B 3W box" << endl;
				     if(isSingleLepton && selectedJets.size() >=6) 
				     {
					myInclFourthGenSearchTools.FillPlots(d,nbOfBtags,nbOfWs,HT,selectedMuons,selectedElectrons,mets,selectedJets,scaleFactor);
							if(isSingleElectron) selecTableSemiEl.Fill(d,17,scaleFactor);
							if(isSingleMuon) selecTableSemiMu.Fill(d,15,scaleFactor);
				     }
				     else if(isSSLepton && selectedJets.size() >=4) 
				     {
							//cout << "IS SAME-SIGN LEPTON" << endl;
							//count number of events with SS leptons
							if(dataSetName.find("NP") != 0) NbSSevents_2B_3W = NbSSevents_2B_3W + datasets[d]->NormFactor()*Luminosity*scaleFactor;
							if(dataSetName.find("NP_Tprime500")==0 || dataSetName.find("NP_STprime500")==0 || dataSetName.find("NP_Bprime500")==0 || dataSetName.find("NP_SBprime500")==0) NP_NbSSevents_2B_3W = NP_NbSSevents_2B_3W + datasets[d]->NormFactor()*Luminosity*scaleFactor;
							MSPlot["MS_NbSSevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);
							histo1D[histo1_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
				      if(semiElectron) selecTableSemiEl.Fill(d,19,scaleFactor);
						  if(semiMuon) selecTableSemiMu.Fill(d,17,scaleFactor);
						 }
				     else if(isTriLepton && selectedJets.size() >=2) 
				     {
							//cout << "IS TRI-LEPTON" << endl;
							//count number of events with three leptons
							if(dataSetName.find("NP") != 0) NbTrievents_2B_3W = NbTrievents_2B_3W + datasets[d]->NormFactor()*Luminosity*scaleFactor;
							if(dataSetName.find("NP_Tprime500")==0 || dataSetName.find("NP_STprime500")==0 || dataSetName.find("NP_Bprime500")==0 || dataSetName.find("NP_SBprime500")==0) NP_NbTrievents_2B_3W = NP_NbTrievents_2B_3W + datasets[d]->NormFactor()*Luminosity*scaleFactor;
							MSPlot["MS_NbTrievents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);
							histo1D[histo2_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
				      if(semiElectron) selecTableSemiEl.Fill(d,20,scaleFactor);
						  if(semiMuon) selecTableSemiMu.Fill(d,18,scaleFactor);
						 }
					//cout << "done in 2B 3W box" << endl;
				}
				else
				{
				     //cout << "in 2B 4W box" << endl;
				     if(isSingleLepton && selectedJets.size() >=8 )
				     {
					myInclFourthGenSearchTools.FillPlots(d,nbOfBtags,nbOfWs,HT,selectedMuons,selectedElectrons,mets,selectedJets,scaleFactor);
							if(isSingleElectron) selecTableSemiEl.Fill(d,18,scaleFactor);
							if(isSingleMuon) selecTableSemiMu.Fill(d,16,scaleFactor);
				     }
				     else if(isSSLepton && selectedJets.size() >=6) 
				     {
							//cout << "IS SAME-SIGN LEPTON" << endl;
							//count number of events with SS leptons
							if(dataSetName.find("NP") != 0) NbSSevents_2B_4W = NbSSevents_2B_4W + datasets[d]->NormFactor()*Luminosity*scaleFactor;
							if(dataSetName.find("NP_Tprime500")==0 || dataSetName.find("NP_STprime500")==0 || dataSetName.find("NP_Bprime500")==0 || dataSetName.find("NP_SBprime500")==0) NP_NbSSevents_2B_4W = NP_NbSSevents_2B_4W + datasets[d]->NormFactor()*Luminosity*scaleFactor;
							MSPlot["MS_NbSSevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);
							histo1D[histo1_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
				      if(semiElectron) selecTableSemiEl.Fill(d,19,scaleFactor);
						  if(semiMuon) selecTableSemiMu.Fill(d,17,scaleFactor);
						 }
				     else if(isTriLepton && selectedJets.size() >=4) 
				     {
						//cout << "IS TRI-LEPTON" << endl;
								//count number of events with three leptons
							if(dataSetName.find("NP") != 0) NbTrievents_2B_4W = NbTrievents_2B_4W + datasets[d]->NormFactor()*Luminosity*scaleFactor;
							if(dataSetName.find("NP_Tprime500")==0 || dataSetName.find("NP_STprime500")==0 || dataSetName.find("NP_Bprime500")==0 || dataSetName.find("NP_SBprime500")==0) NP_NbTrievents_2B_4W = NP_NbTrievents_2B_4W + datasets[d]->NormFactor()*Luminosity*scaleFactor;
							MSPlot["MS_NbTrievents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);
							histo1D[histo2_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
				      if(semiElectron) selecTableSemiEl.Fill(d,20,scaleFactor);
						  if(semiMuon) selecTableSemiMu.Fill(d,18,scaleFactor);
						}
				     //cout << "done in 2B 4W box" << endl;
				}
			} //end number of btags == 2
		
	}//useMassesAndResolutions = true
	//delete selection;
    }//loop on events

    ////to test 'purity' of good jet combinations with a 'simple' method (not MVA)
    //myInclFourthGenSearchTools.PrintPurityGoodCombinations();
    
    cout<<endl;
    
    //important: free memory
    treeLoader.UnLoadDataset();

    if(jetTools) delete jetTools;
    
    if (!useMassesAndResolutions && ((dataSetName.find("TTbarJets_SemiMu") == 0 && semiMuon) ||  (dataSetName.find("TTbarJets_SemiElectron") == 0 && semiElectron)))
    {
	
      string resfilename = "MassPlot.root";
     
      cout << "INFO: Creating output file to store the mass of the W boson: " << resfilename << endl;    

      TFile* outfile = new TFile(resfilename.c_str(),"RECREATE");
      outfile->cd();
      
      histo1D["hadronicRecoWMass"]->Write();

      // fit the distributions
      TH1F* histo = histo1D["hadronicRecoWMass"];
      TF1 *fitfunc;
      string func_title = string(histo->GetName())+"_Fitted";

      double rms = histo->GetRMS();
      double maxbin =  histo->GetBinCenter(histo->GetMaximumBin());

      fitfunc = new TF1(func_title.c_str(),"gaus");
      fitfunc->SetRange(maxbin-rms,maxbin+rms);
      histo->Fit(fitfunc,"RQ");

      fitfunc->Write();
      delete fitfunc;
      
      outfile->Close();
     //delete outfile; 
    }
  } //loop on datasets
  
  cout << "NbSSevents_1B_2W " << NbSSevents_1B_2W << endl;
  cout << "NbSSevents_1B_3W " << NbSSevents_1B_3W << endl;
  cout << "NbSSevents_1B_4W " << NbSSevents_1B_4W << endl;
  cout << "NbSSevents_2B_2W " << NbSSevents_2B_2W << endl;
  cout << "NbSSevents_2B_3W " << NbSSevents_2B_3W << endl;
  cout << "NbSSevents_2B_4W " << NbSSevents_2B_4W << endl;
  cout << "signal (500) NbSSevents_1B_2W " << NP_NbSSevents_1B_2W << endl;
  cout << "signal (500) NbSSevents_1B_3W " << NP_NbSSevents_1B_3W << endl;
  cout << "signal (500) NbSSevents_1B_4W " << NP_NbSSevents_1B_4W << endl;
  cout << "signal (500) NbSSevents_2B_2W " << NP_NbSSevents_2B_2W << endl;
  cout << "signal (500) NbSSevents_2B_3W " << NP_NbSSevents_2B_3W << endl;
  cout << "signal (500) NbSSevents_2B_4W " << NP_NbSSevents_2B_4W << endl;
  cout << "NbTrievents_1B_3W " << NbTrievents_1B_3W << endl;
  cout << "NbTrievents_1B_4W " << NbTrievents_1B_4W << endl;
  cout << "NbTrievents_2B_3W " << NbTrievents_2B_3W << endl;
  cout << "NbTrievents_2B_4W " << NbTrievents_2B_4W << endl;
  cout << "signal (500) NbTrievents_1B_3W " << NP_NbTrievents_1B_3W << endl;
  cout << "signal (500) NbTrievents_1B_4W " << NP_NbTrievents_1B_4W << endl;
  cout << "signal (500) NbTrievents_2B_3W " << NP_NbTrievents_2B_3W << endl;
  cout << "signal (500) NbTrievents_2B_4W " << NP_NbTrievents_2B_4W << endl;
  
  //Once everything is filled ...
  cout << " We ran over all the data ;-)" << endl;
  
  ///////////////////
  // Writing
  //////////////////
  cout << " - Writing outputs to the files ..." << endl;

  if(doMVAjetcombination)
  {
    string pathPNGJetCombi = pathPNG+"JetCombination/";
    mkdir(pathPNGJetCombi.c_str(),0777);
    jetCombiner->Write(fout, true, pathPNGJetCombi);
  }
  
  // Fill the resolution histograms and calculate the resolutions
  if(CalculateResolutions)
  {
      fout->cd();
      cout<<" Write resolutions"<<endl; 
      mkdir((pathPNG+"resFit_LightJet/").c_str(),0777);
      mkdir((pathPNG+"resFit_BJet/").c_str(),0777);
    
      resFitLightJets->WritePlots(fout, true, pathPNG+"resFit_LightJet/");
      resFitLightJets->WriteResolutions("lightJetReso.root");
      //resFitBJets->WritePlots(fout, true, pathPNG+"resFit_BJet/");
      resFitBJets->WriteResolutions("bJetReso.root");
  }

  
  if(!TrainMVA)
  {
    fout->cd();
    //Write histograms: MSPlots
    mkdir((pathPNG+"MSPlot/").c_str(),0777);
    //cout << "mkdir " << (pathPNG+"MSPlot/").c_str()<< endl;
    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
 				cout << "run over all MS plots 1" << endl;
        MultiSamplePlot *temp = it->second;
 				cout << "run over all MS plots 2" << endl;
        string name = it->first;
 				cout << "run over all MS plots 3 " << name << endl;
        temp->Draw(false, name, true, true, true, true, true,5);//(bool addRandomPseudoData, string label, bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST,int scaleNPsignal)
 				cout << "run over all MS plots 4 " << endl;
        temp->Write(fout, name, true, pathPNG+"MSPlot/");//bool savePNG
 				cout << "run over all MS plots - writing done" << endl;
    }
    cout << "MultiSamplePlots written" << endl;
	
    //Write histograms: 1D 
    TDirectory* th1dir = fout->mkdir("1D_histograms");
    fout->cd();
    th1dir->cd();
    for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
    {
		TH1F *temp = it->second;
 //		int N = temp->GetNbinsX();
 //  	temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
 //  	temp->SetBinContent(N+1,0);
 //		temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
		temp->Write();
		TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
		tempCanvas->SaveAs( (pathPNG+it->first+".pdf").c_str() );
    }    
    cout << "1D plots written" << endl;
    
    myInclFourthGenSearchTools.WritePlots(fout,th1dir,pathPNG);
    
    //delete th1dir;
    // 2D
    TDirectory* th2dir = fout->mkdir("2D_histograms_graphs");
    th2dir->cd();
    for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
    {
		TH2F *temp = it->second;
		temp->Write();
		TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
		tempCanvas->SaveAs( (pathPNG+it->first+".pdf").c_str() );
    }
    //delete th2dir;
    cout << "2D plots written" << endl;
    fout->cd();
    
    //Selection tables
    selecTableSemiMu.TableCalculator(false, true, true, true, true);//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
    string selectiontableSemiMu = "InclFourthGenSearch_SelectionTable_SemiMu"+postfix;
    if (argc >= 3)
    {
		string sample=string(argv[2]);
 		selectiontableSemiMu = selectiontableSemiMu +"_"+sample;
    }
    selectiontableSemiMu = selectiontableSemiMu +".tex"; 	
    if(semiMuon) selecTableSemiMu.Write(selectiontableSemiMu.c_str());
	
    selecTableSemiEl.TableCalculator(false, true, true, true, true);//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
    string selectiontableSemiEl = "InclFourthGenSearch_SelectionTable_SemiEl"+postfix;
    if (argc >= 3)
    {
		string sample=string(argv[2]);
 		selectiontableSemiEl = selectiontableSemiEl +"_"+sample;
    }
    selectiontableSemiEl = selectiontableSemiEl +".tex"; 	
    if(semiElectron) selecTableSemiEl.Write(selectiontableSemiEl.c_str());
	
     //regarding binning, which is only relevant when you do the MVA jet combination to reconstruct the top mass   
     if (doMVAjetcombination && anaEnv.MCRound==0)
     {	
        HTvsMTop_1B_2W.Write_for2DBinning(binningFileName_HTvsMTop_1B_2W);
				HTvsMTop_2B_2W.Write_for2DBinning(binningFileName_HTvsMTop_2B_2W);
     }    
     else if(doMVAjetcombination && anaEnv.MCRound==1)
     {  
				HTvsMTop_1B_2W.Convert2Dto1D(postfix);
        HTvsMTop_2B_2W.Convert2Dto1D(postfix);
        HTvsMTop_1B_2W.Write(fout,th1dir);
				HTvsMTop_2B_2W.Write(fout,th1dir);
     }
     
    fout->cd();
    
    cout << " - Closing the output file now..." << endl;
    fout->Close();
  } //end !trainMVA
  
  //delete
  if(jetCombiner) delete jetCombiner; //IMPORTANT!! file for training otherwise not filled... (?) //crashes when calculating resolutions for kinfit
  
  delete fout;

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;

  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;

  return 0;
}

//To cout the Px, Py, Pz, E and Pt of objects
void coutObjectsFourVector(vector < TRootMuon* > init_muons, vector < TRootElectron* > init_electrons, vector < TRootJet* > init_jets, vector < TRootMET* > mets, string Comment)
{
     cout<<Comment<<endl;
     
     for(unsigned int k=0; k<init_muons.size(); k++)
     {
	   cout<<" init_muons["<<k<<"] -> Px() = "<<init_muons[k]->Px()<<endl;
	   cout<<"              -> Py() = "<<init_muons[k]->Py()<<endl;
	   cout<<"              -> Pz() = "<<init_muons[k]->Pz()<<endl;
	   cout<<"                -> Pt() = "<<init_muons[k]->Pt()<<endl;
	   cout<<"              -> E() = "<<init_muons[k]->E()<<endl;   
     }
     for(unsigned int k=0; k<init_electrons.size(); k++)
     {
	   cout<<" init_electrons["<<k<<"] -> Px() = "<<init_electrons[k]->Px()<<endl;
	   cout<<"              -> Py() = "<<init_electrons[k]->Py()<<endl;
	   cout<<"              -> Pz() = "<<init_electrons[k]->Pz()<<endl;
	   cout<<"                -> Pt() = "<<init_electrons[k]->Pt()<<endl;
	   cout<<"              -> E() = "<<init_electrons[k]->E()<<endl;   
     }         
     for(unsigned int k=0; k<init_jets.size(); k++)
     {
	   cout<<" init_jets["<<k<<"] -> Px() = "<<init_jets[k]->Px()<<endl;
	   cout<<"              -> Py() = "<<init_jets[k]->Py()<<endl;
	   cout<<"              -> Pz() = "<<init_jets[k]->Pz()<<endl;
	   cout<<"                -> Pt() = "<<init_jets[k]->Pt()<<endl;
	   cout<<"              -> E() = "<<init_jets[k]->E()<<endl;	   
     }
     for(unsigned int k=0; k<mets.size(); k++)
     {
           cout<<" mets["<<k<<"] -> Px() = "<<mets[k]->Px()<<endl;
           cout<<"         ->  Py() = "<<mets[k]->Py()<<endl;
	   cout<<"         ->  Pz() = "<<mets[k]->Pz()<<endl;
	   cout<<"              -> Pt() = "<<mets[k]->Pt()<<endl;
	   cout<<"         ->  E() = "<<mets[k]->E()<<endl;
	   cout<<"              -> Et() = "<<mets[k]->Et()<<endl;
     }
};

//https://twiki.cern.ch/twiki/bin/viewauth/CMS/SingleTopTurnOnCurves
float jetprob(float jetpt, float btagvalue){
	float prob=0.982*exp(-30.6*exp(-0.151*jetpt));
  prob*=0.844*exp((-6.72*exp(-0.720*btagvalue))); //"for the offline TCHP tagger"
	//prob*=0.736*exp((-8.01*exp(-0.540*btagvalue))); //"for the offline TCHE tagger"
	return prob;
};
