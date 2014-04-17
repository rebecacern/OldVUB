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
#include "../JESMeasurement/interface/JetCombiner.h"
#include "../Reconstruction/interface/JetCorrectorParameters.h"
#include "../Reconstruction/interface/JetCorrectionUncertainty.h"
#include "../Reconstruction/interface/MakeBinning.h"
#include "../Reconstruction/interface/TTreeObservables.h"
#include "../MCInformation/interface/Lumi3DReWeighting.h"
#include "../InclFourthGenSearch/interface/InclFourthGenTree.h"
#include "../InclFourthGenSearch/interface/InclFourthGenSearchTools.h"
#include "../InclFourthGenSearch/interface/TwoDimTemplateTools.h"
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

//To cout the Px, Py, Pz, E and Pt of objects
void coutObjectsFourVector(vector < TRootMuon* > init_muons, vector < TRootElectron* > init_electrons, vector < TRootJet* > init_jets, vector < TRootMET* > mets, string Comment);

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

int main (int argc, char *argv[])
{ 

  //which systematic to run?
  string systematic = "Nominal";
  if (argc >= 2)
		systematic = string(argv[1]);
  cout << "Systematic to be used: " << systematic << endl;
  if( ! (systematic == "Nominal"  || systematic == "JESPlus" || systematic == "JESMinus" || systematic == "JERPlus" || systematic == "JERMinus") )
  {
    cout << "Unknown systematic!!!" << endl;
    cout << "Possible options are: Nominal JESPlus JESMinus JERPlus JERMinus" << endl;
    exit(-1);
  }
	
	//which channel to run?
	string channelpostfix = "";
  bool semiElectron = false; // use semiElectron channel?
  bool semiMuon = true; // use semiMuon channel?
	if (argc >= 3)
	{	
	  semiMuon = atoi(argv[2]);
		semiElectron = !semiMuon;
	}
	
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

	//which additional option to run?
  string option = "";
  bool datadriven = false;
	if (argc >= 4)
	{
		option = string(argv[3]);
  	cout << "Option to be used: " << option << endl;
  	datadriven = true;
	}
	if( ! (option == "" || option == "ChargeMisId" || option == "FakeLepton") )
  {
    cout << "Unknown option!!!" << endl;
    cout << "Possible options are: ChargeMisId , FakeLepton" << endl;
    exit(-1);
  }
	
	
	float CM_EB = 0.00140; //charge misid probability barrel electron
	float CM_EB_UNC = 0.00015;
	float CM_EE = 0.01393; //charge misid probability endcap electron
	float CM_EE_UNC = 0.00176;

	float Eff_TL_e = 0.0821;	
	float Eff_TL_e_unc = 0.0071;
	float Eff_TL_m = 0.0331;	
	float Eff_TL_m_unc = 0.0019;

	//btagger to use and corresponding workingpoints
  string btagger = "TCHPM";
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

	int NbOSMuons_data = 0;
	int NbOSElectrons_data = 0; int NbOSElectrons_EBEB_data = 0; int NbOSElectrons_EBEE_data = 0; int NbOSElectrons_EEEE_data = 0;
	int NbOSElMu_data = 0; int NbOSElMu_EB_data = 0; int NbOSElMu_EE_data = 0;
	int NbSSLooseMuonTightMuon_data = 0; int NbSSLooseElectronTightMuon_data = 0;
	int NbSSLooseElectronTightElectron_data = 0; int NbSSLooseMuonTightElectron_data = 0;

  //SetStyle if needed
  setTDRStyle();
  //setMyStyle();

  string postfix = ""; // to relabel the names of the output file  
	postfix= postfix+"_"+systematic;

  string Treespath = "InclFourthGenTrees_Fall11_dataonly";
  Treespath = Treespath +"/";
  if(!datadriven) mkdir(Treespath.c_str(),0777);
	bool savePNG = false;
	
  /////////////////////
  // Configuration
  /////////////////////
  
  //xml file
  string xmlFileName = "";
	if(semiElectron) xmlFileName = "../config/myFourthGenconfig_Electron_Fall11.xml";
  else if(semiMuon) xmlFileName = "../config/myFourthGenconfig_Muon_Fall11.xml";
  const char *xmlfile = xmlFileName.c_str();
  cout << "used config file: " << xmlfile << endl;    
  
  
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
  string rootFileName = (Treespath+"InclFourthGenSearch_TreeCreator"+postfix+channelpostfix+".root");
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

  string pathPNG = Treespath+"InclFourthGenSearchPlots_TreeCreator"+postfix+channelpostfix;
  pathPNG = pathPNG +"/"; 	
  pathPNG = pathPNG +"/"; 	
  if(savePNG) mkdir(pathPNG.c_str(),0777);


  MSPlot["allDiJetMasses"] = new MultiSamplePlot(datasets, "allDiJetMasses", 50, 0, 1000, "m_{jj}");
  MSPlot["hadronicRecoWMass_chosenWjets"] = new MultiSamplePlot(datasets, "hadronicRecoWMass_chosenWjets", 50, 0, 500, "m_{W}"); 
  histo1D["hadronicPartonWMass"] = new TH1F("hadronicPartonWMass","Hadronic W Mass, using the Partons",100,0,200);
  histo1D["hadronicRecoWMass"] = new TH1F("hadronicRecoWMass","Hadronic W Mass, using the RecoJets",100,0,200);
  
  histo1D["lumiWeights"] = new TH1F("lumiWeights","lumiWeights;lumiWeight;#events",100,0,4);
  histo1D["LeptonPt_TTbar"] = new TH1F("leptonspt ttbar","leptonspt ttbar;pt leptons;#events",250,0,500);
  histo1D["LeptonPt_Tprime500"] = new TH1F("leptonspt tprime500","leptonspt tprime500;pt leptons;#events",250,0,500);
  histo1D["LeptonPt_Bprime500"] = new TH1F("leptonspt bprime500","leptonspt bprime500;pt leptons;#events",250,0,500);
  histo1D["LeptonPt_SBprime500"] = new TH1F("leptonspt sbprime500","leptonspt sbprime500;pt leptons;#events",250,0,500);
	
	MSPlot["Reliso_Lepton"] = new MultiSamplePlot(datasets, "Lepton reliso", 50, 0, 10, "Lepton reliso");
	
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
  MSPlot["MS_MET"] = new MultiSamplePlot(datasets,"MET", 75, 0, 150, "Missing transverse energy (GeV)");
  MSPlot["MS_LeptonPt"] = new MultiSamplePlot(datasets,"lepton pt", 150, 0, 300, "Lepton Pt (GeV)");
  MSPlot["MS_nPV"] = new MultiSamplePlot(datasets, "nPrimaryVertices", 21, -0.5, 20.5, "Nr. of primary vertices");
  MSPlot["MS_JetMultiplicity_SingleLepton"] = new MultiSamplePlot(datasets, "JetMultiplicity", 10, -0.5, 9.5, "Jet Multiplicity");
  MSPlot["MS_BtaggedJetMultiplicity_SingleLepton"] = new MultiSamplePlot(datasets, "BtaggedJetMultiplicity", 7, -0.5, 6.5, "b-tagged jet multiplicity");
	MSPlot["MS_JetMultiplicityAtleast1Btag_SingleLepton"] = new MultiSamplePlot(datasets, "JetMultiplicityAtleast1Btag", 10, -0.5, 9.5, "Jet multiplicity (>=1 b-tag)");

  MSPlot["MS_JetPt_all_SingleLepton"] = new MultiSamplePlot(datasets,"JetPt_all", 50, 0, 300, "Pt of all jets (GeV)");
	MSPlot["MS_JetPt_btagged_SingleLepton"] = new MultiSamplePlot(datasets,"JetPt_btagged", 50, 0, 300, "Pt of b-tagged jets (GeV)");
	MSPlot["MS_JetPt_nonbtagged_SingleLepton"] = new MultiSamplePlot(datasets,"JetPt_nonbtagged", 50, 0, 300, "Pt of non b-tagged jets (GeV)");
	
	MSPlot["MS_LeptonRelIso"] = new MultiSamplePlot(datasets, "LeptonRelIso", 50, 0, 1, "Lepton RelIso");
	
	cout << " - Declared histograms ..." <<  endl;

  float NbSSevents = 0;   
  float NbTrievents = 0;  
	
  ////////////////////////////////////
  /// Selection table
  ////////////////////////////////////

  vector<string> CutsSelecTableSemiLep;
  CutsSelecTableSemiLep.push_back(string("initial")); //0
  CutsSelecTableSemiLep.push_back(string("preselected"));
  CutsSelecTableSemiLep.push_back(string("trigged"));
  CutsSelecTableSemiLep.push_back(string("Good PV"));
  CutsSelecTableSemiLep.push_back(string("$\\geq$ 1 muon/electron"));
  CutsSelecTableSemiLep.push_back(string("$\\geq$ 1 b-tagged jet"));
  CutsSelecTableSemiLep.push_back(string("MET $>$ 40 GeV"));  
  CutsSelecTableSemiLep.push_back(string("single muon/electron"));

  vector<string> CutsSelecTableMultiLepton;
  CutsSelecTableMultiLepton.push_back(string("SS leptons"));
  CutsSelecTableMultiLepton.push_back(string("trileptons"));

  vector<string> CutsSelecTableChargeMisId_2El;
  CutsSelecTableChargeMisId_2El.push_back(string("2 electrons"));
  CutsSelecTableChargeMisId_2El.push_back(string("2 electrons EB"));
  CutsSelecTableChargeMisId_2El.push_back(string("2 electrons EE"));
  CutsSelecTableChargeMisId_2El.push_back(string("2 electrons EB+EE"));
  CutsSelecTableChargeMisId_2El.push_back(string("SS el EB"));
  CutsSelecTableChargeMisId_2El.push_back(string("SS el EE"));
  CutsSelecTableChargeMisId_2El.push_back(string("SS el EB+EE"));
  CutsSelecTableChargeMisId_2El.push_back(string("OS el EB"));
  CutsSelecTableChargeMisId_2El.push_back(string("OS el EE"));
  CutsSelecTableChargeMisId_2El.push_back(string("OS el EB+EE"));
	
  vector<string> CutsSelecTableChargeMisId_ElMu;
  CutsSelecTableChargeMisId_ElMu.push_back(string("electron+muon"));
  CutsSelecTableChargeMisId_ElMu.push_back(string("electron+muon EB"));
  CutsSelecTableChargeMisId_ElMu.push_back(string("electron+muon EE"));
  CutsSelecTableChargeMisId_ElMu.push_back(string("SS el EB"));
  CutsSelecTableChargeMisId_ElMu.push_back(string("SS el EE"));
  CutsSelecTableChargeMisId_ElMu.push_back(string("OS el EB"));
  CutsSelecTableChargeMisId_ElMu.push_back(string("OS el EE"));
	
  SelectionTable selecTableSemiLep(CutsSelecTableSemiLep, datasets);
 	selecTableSemiLep.SetLuminosity(Luminosity);
  selecTableSemiLep.SetPrecision(1);
  SelectionTable selecTableMultiLepton(CutsSelecTableMultiLepton, datasets);
  selecTableMultiLepton.SetLuminosity(Luminosity);
  selecTableMultiLepton.SetPrecision(1);
  SelectionTable selecTableChargeMisId_2El(CutsSelecTableChargeMisId_2El, datasets);
  selecTableChargeMisId_2El.SetLuminosity(Luminosity);
  selecTableChargeMisId_2El.SetPrecision(2);
  SelectionTable selecTableChargeMisId_ElMu(CutsSelecTableChargeMisId_ElMu, datasets);
  selecTableChargeMisId_ElMu.SetLuminosity(Luminosity);
  selecTableChargeMisId_ElMu.SetPrecision(2);
  
  cout << " - SelectionTable instantiated ..." << endl;

  ////////////////////////////////////////////////////
  // PileUp Reweighting - 3D//
  ////////////////////////////////////////////////////
//  Lumi3DReWeighting Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Flat10PlusTail.root","PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root", "pileup", "pileup");
  Lumi3DReWeighting Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Fall11.root","PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root", "pileup", "pileup");
  Lumi3DWeights.weight3D_init(1.0);

  cout << " - Initialized LumiReWeighting stuff" << endl;  
    
	
		
  ////////////////////////////////////
  ////////////////////////////////////
  ///////// Loop on datasets
  ////////////////////////////////////
  ////////////////////////////////////
  cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
  ofstream myfileSS, myfileLLL;
	string mySSFile = Treespath+"InterestingEvents_SS"+channelpostfix+".txt";
	if(systematic=="Nominal") myfileSS.open(mySSFile.c_str());
	string myLLLFile = Treespath+"InterestingEvents_lll"+channelpostfix+".txt";
	if(systematic=="Nominal") myfileLLL.open(myLLLFile.c_str());

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
    selecTableSemiLep.Fill(d,0, datasets[d]->Xsection() * datasets[d]->EquivalentLumi() );
    
		
		
    /////////////////////////////////////
    /// Initialize JEC factors 
    /////////////////////////////////////
		
		//"OLD way, no JEC corrections on the fly"
    /*//L2L3 residual corrections already in data Toptrees now! (because a global tag is used where these corrections are included)
    vector<JetCorrectorParameters> vCorrParam;    
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/Jec11V2_db_AK5PFchs_Uncertainty.txt");
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, false);
    */
   	    
    vector<JetCorrectorParameters> vCorrParam;

    // Create the JetCorrectorParameter objects, the order does not matter.
   // YYYY is the first part of the txt files: usually the global tag from which they are retrieved
    JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters("JECFiles/START42_V17_AK5PFchs_L3Absolute.txt");
    JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters("JECFiles/START42_V17_AK5PFchs_L2Relative.txt");
    JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters("JECFiles/START42_V17_AK5PFchs_L1FastJet.txt");

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
    
    // true means redo also the L1
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);// last boolean ('startFromRaw') = true!
		
		

    string TreeFileName;
		TreeFileName = Treespath+"InclFourthGenTree_"+dataSetName+postfix+channelpostfix+".root";
    cout << "INFO: creating InclFourthGenTree file "+TreeFileName << endl;        
    TFile* treeFile;
		TTree* myInclFourthGenTree;
    InclFourthGenTree* myBranch_selectedEvents = 0;		
		if(!datadriven)
		{
			treeFile = new TFile(TreeFileName.c_str(),"RECREATE");
			myInclFourthGenTree = new TTree("myInclFourthGenTree","Tree containing the FourthGen information");    
			myInclFourthGenTree->Branch("InclFourthGenBranch_selectedEvents","InclFourthGenTree",&myBranch_selectedEvents);
		}


    ////////////////////////////////////
    ////////////////////////////////////
    ///////// Loop on events
    ////////////////////////////////////
    ////////////////////////////////////
    int itrigger = -1, previousRun = -1;
     
    if (verbose > 1)
      cout << " - Loop over events " << endl;      
    
    for (int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
//    for (int ievt = 0; ievt < 200; ievt++)
    {        

      if(ievt%1000 == 0)
        std::cout<<"Processing the "<<ievt<<"th event ("<<100*ievt/datasets[d]->NofEvtsToRunOver()<<"%)"<<flush<<"\r";
      
			//load event
      event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);

			//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"First cout after loading event:");
			
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
			//cout<<" currentFilename = "<<currentFilename<<", previousFilename = "<<previousFilename<<endl;
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
			//cout<<" currentRun = "<<currentRun<<", previousRun = "<<previousRun<<endl;
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
   					if(dataSetName == "ttW" || dataSetName == "ttZ" || dataSetName == "samesignWWjj" || dataSetName == "TTbarJets_scaleup" || dataSetName == "TTbarJets_scaledown" || dataSetName == "TTbarJets_matchingup" || dataSetName == "TTbarJets_matchingdown" || dataSetName.find("QCD")<=0 )
						  itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v5"), currentRun, iFile);//Summer11 MC! also the TTJets systematic samples...!
						else
						  itrigger = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v3"), currentRun, iFile);//Fall11 MC!
						
    
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
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP_v5"), currentRun, iFile);  				   					
						else if (event->runId() >= 178411 && event->runId() < 179942)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP_v8"), currentRun, iFile);  				   					
						else if (event->runId() >= 179942 && event->runId() <= 180296)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP_v9"), currentRun, iFile);  				   											   
  					if(itrigger == 9999)
						{
    					cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (DATA) IN RUN " << event->runId() << endl;
    					exit(1);
  					}// semi-electron
 	   			}
	   			else 
	   			{
					  //Problem: a trigger reweighting procedure for MC should be done when using the summer11 electron trigger...
					  if(dataSetName == "ttW" || dataSetName == "ttZ" || dataSetName == "samesignWWjj" || dataSetName == "TTbarJets_scaleup" || dataSetName == "TTbarJets_scaledown" || dataSetName == "TTbarJets_matchingup" || dataSetName == "TTbarJets_matchingdown")
   						itrigger = treeLoader.iTrigger (string ("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2"), currentRun, iFile);//Summer11 MC has other triggers!	
						else
						  itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP_v5"), currentRun, iFile);//Fall11 MC!
						
						if(itrigger == 9999)
						{
							cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (" << dataSetName << ") IN RUN " << event->runId() << endl;	
							exit(1);
						}
					}	 
				} //end if semiElectron	
			} //end previousRun != currentRun

			// JES CORRECTION   
      // Apply Jet Corrections on-the-fly: not if already in toptrees! (our first Fall11 round)
			//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JES correction on the fly:");
//			if( dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
//				jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho(),true); //last boolean: isData (needed for L2L3Residual...)
//			else
//				jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho(),false); //last boolean: isData (needed for L2L3Residual...)
		  //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JES correction on the fly:");

      //ordering is relevant; most probably 1) Type I MET correction, 2) JER where jet corrections are propagated to MET, 3) JES systematics where jet corrections are propagated to MET
      //----------------------------------------------------------
      // Apply type I MET corrections:  (Only for |eta| <= 4.7 )
      //---------------------------------------------------------
      
			//not if already in toptrees! (our first Fall11 round)
			//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before MET type I correction:");      
//      if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
//        jetTools->correctMETTypeOne(init_jets,mets[0],true);
//      else
//        jetTools->correctMETTypeOne(init_jets,mets[0],false);
      //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After MET type I correction:");
     	 
		  
      if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
      {	
	
      	//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JER correction:");
				if(systematic == "JERMinus")
					jetTools->correctJetJER(init_jets, genjets, mets[0], "minus",false); //false means don't use old numbers but newer ones...
				else if(systematic == "JERPlus")
					jetTools->correctJetJER(init_jets, genjets, mets[0], "plus",false);
				else
					jetTools->correctJetJER(init_jets, genjets, mets[0], "nominal",false);
				//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JER correction:");	       
		

				// JES systematic! 
				if (systematic == "JESMinus")
					jetTools->correctJetJESUnc(init_jets, mets[0], "minus");
				else if (systematic == "JESPlus")
					jetTools->correctJetJESUnc(init_jets, mets[0], "plus");	       
      }



			double lumiWeight3D = 1.0;
			if(!(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA"))
			{				
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
      Selection selectionFakeLepton(init_jets, init_muons, init_electrons, mets); //mets can also be corrected...

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
			bool isDiLepton = false;
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
			bool isEB = false;
			bool isEE = false;
			bool isEBEB = false;
			bool isEBEE = false;
			bool isEEEE = false;
			
      
      vector<TRootJet*> selectedJets, selectedForwardJets, selectedJetsLargeEtaRange;
      vector<TRootMuon*> selectedMuons, selectedLooseMuons, selectedMuons_FL, selectedLooseMuons_FL, selectedOnlyLooseMuons_FL;
      vector<TRootElectron*> selectedElectrons, selectedLooseElectronsNoVBTFid, selectedLooseElectronsVBTFid, selectedElectrons_FL, selectedLooseElectrons_FL, selectedOnlyLooseElectrons_FL;
      vector<TRootMCParticle*> mcParticles;


      //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before treeLoader.LoadMCEvent:");      
      //if(dataSetName.find("TTbarJets_SemiMu") == 0 || dataSetName.find("TTbarJets_SemiElectron") == 0 || dataSetName.find("NP_Tprime")==0 || dataSetName.find("NP_overlay_Tprime")==0)
			if(dataSetName != "data" && dataSetName != "Data" && dataSetName != "Data")
      {
        treeLoader.LoadMCEvent(ievt, 0, 0, mcParticles,false);  
        sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
			//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After treeLoader.LoadMCEvent:");

      
      float METCut = 40;
      selection.setJetCuts(30.,2.4,0.01,1.,0.98,0.3,0.1);
      nonstandard_selection.setJetCuts(30.,4.7,0.01,1.,0.98,0.3,0.1); //only difference: larger eta acceptance 
      selection.setMuonCuts(20,2.1,0.125,10,0.02,0.3,1,1,1); //selection.setMuonCuts(20,2.1,0.125,10,0.02,0.3,1,1,1);
      selection.setLooseMuonCuts(10,2.5,0.2);
      selection.setElectronCuts(20,2.5,0.1,0.02,1,0.3);
      selection.setLooseElectronCuts(15,2.5,0.2);	 				 

      if (init_jets.size() > 0)
      {
	    	selectedJets = selection.GetSelectedJets(true);				
	    	selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
	    	selectedElectrons = selection.GetSelectedElectrons(vertex[0],selectedJets);
	    	selectedJetsLargeEtaRange = nonstandard_selection.GetSelectedJets(true);	    	    
      	for(unsigned int i = 0; i<selectedJetsLargeEtaRange.size(); i++)
      	{
	    		if(selectedJetsLargeEtaRange[i]->Eta() > 2.4)
	      		selectedForwardJets.push_back(selectedJetsLargeEtaRange[i]);
      	}
      }
      selectedLooseElectronsNoVBTFid = selection.GetSelectedLooseElectrons(false); //no vbtfid is required
      selectedLooseElectronsVBTFid = selection.GetSelectedLooseElectrons(true); //vbtfid is required 
      selectedLooseMuons = selection.GetSelectedLooseMuons(); //veto muons	

      if(datadriven)
			{
      	selectionFakeLepton.setMuonCuts(20,2.1,0.125,10,0.02,0.3,1,1,1); //selection.setMuonCuts(20,2.1,0.125,10,0.02,0.3,1,1,1);
      	selectionFakeLepton.setElectronCuts(20,2.5,0.1,0.02,1,0.3);
				selectedMuons_FL = selectionFakeLepton.GetSelectedMuons(vertex[0],selectedJets);
	    	selectedElectrons_FL = selectionFakeLepton.GetSelectedElectrons(vertex[0],selectedJets);
      	selectionFakeLepton.setLooseMuonCuts(10,2.5,9999.);
      	selectionFakeLepton.setLooseElectronCuts(15,2.5,9999.);	 				
      	selectedLooseElectrons_FL = selectionFakeLepton.GetSelectedLooseElectrons(false); //for data-driven background estimation
      	selectedLooseMuons_FL = selectionFakeLepton.GetSelectedLooseMuons(); //for data-driven background estimation
				for(unsigned int i = 0; i<selectedLooseMuons_FL.size(); i++)
      	{
      		int IsTight = 0;
      		for(unsigned int j = 0; j<selectedMuons.size(); j++)
      		{
						if(fabs(selectedLooseMuons_FL[i]->Pt()-selectedMuons[j]->Pt()) < 0.00001) // the loose muon is tight
							IsTight = IsTight+1;
					}
					if(IsTight==0) selectedOnlyLooseMuons_FL.push_back(selectedLooseMuons_FL[i]); // we found no tight muon, so this loose muon is not tight!
				}
				for(unsigned int i = 0; i<selectedLooseElectrons_FL.size(); i++)
      	{
      		int IsTight = 0;
      		for(unsigned int j = 0; j<selectedElectrons.size(); j++)
      		{
						if(fabs(selectedLooseElectrons_FL[i]->Pt()-selectedElectrons[j]->Pt()) < 0.00001)
							IsTight = IsTight+1;
					}
					if(IsTight==0) selectedOnlyLooseElectrons_FL.push_back(selectedLooseElectrons_FL[i]);
				}
			}

     
      selecTableSemiLep.Fill(d,1,scaleFactor);		
			
			
			//// EVENTS TRIGGERED BY MUON TRIGGER			
      if(trigged && semiMuon)
      { 
				selecTableSemiLep.Fill(d,2,scaleFactor);
        if(isGoodPV)
				{
					selecTableSemiLep.Fill(d,3,scaleFactor);
					if(selectedMuons.size()>=1 && selectedMuons[0]->Pt()>40)
					{
						selecTableSemiLep.Fill(d,4,scaleFactor);
						sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class
						
						if(selectedJets.size()>=(unsigned int)anaEnv.NofJets)
						{  //at least 1 jet!
						
						  //block for the jet multiplicity plot
							if(mets[0]->Et()> METCut)
							{							
							      if(selectedMuons.size() == 1 && selectedLooseMuons.size() == selectedMuons.size() && selectedLooseElectronsVBTFid.size() == 0)
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
							//for(unsigned int j=0;j<selectedJets.size();j++)
							//{
								//now require at least a b-tagged jet larger than a certain pre-defined cut //update: not requiring offline b-tag in treecreator
								//if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue && !eventSelected)
								//{
									selecTableSemiLep.Fill(d,5,scaleFactor); 
									if(mets[0]->Et()> METCut)
									{
										selecTableSemiLep.Fill(d,6,scaleFactor);
										eventSelected = true; 

										//cout << "event is selected according to the baseline selection!" << endl;
										
										////for single muon require exactly 1 muon, veto for other loose muons and veto for very loose electrons
										if(selectedMuons.size() == 1 && selectedLooseMuons.size() == selectedMuons.size() && selectedLooseElectronsVBTFid.size() == 0)
										{
											isSingleLepton = true; // we have a single muon event
											isSingleMuon = true;
											//std::cout<<"Processing the "<<ievt<<"th event" << endl;
											//cout << "is single muon!" << endl;
											//cout << "-> muon pt: " << selectedMuons[0]->Pt() << endl;
			     					}
			     					
										
										////for same-sign muons require exactly 2 muons, veto for very loose electrons
										else if(selectedMuons.size() == 2 && selectedLooseMuons.size() == selectedMuons.size() && selectedLooseElectronsVBTFid.size() == 0)
										{ 
											//require that there are the two muons do not form the Z mass
											if( !selection.foundZCandidate(selectedMuons, selectedMuons, 10.) )
											{
											  isDiLepton = true;
												
												if(option=="ChargeMisId")
												{
													if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
													{
														bool eventBtag = false;
														for(unsigned int j = 0; j < selectedJets.size(); j++)
															if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
																eventBtag = true;
														if(selectedMuons[0]->charge() != selectedMuons[1]->charge() && eventBtag) 
															NbOSMuons_data++;
													}
												}
												//require the same charge
												if(selectedMuons[0]->charge() == selectedMuons[1]->charge())
												{
													isSSLepton = true; // we have two same-sign muons
													isSSMuon = true;
													bool eventBtag = false;
													for(unsigned int j = 0; j < selectedJets.size(); j++)
														if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
															eventBtag = true;
													if(eventBtag && dataSetName=="Data" && systematic=="Nominal")
														myfileSS << "Run: " << event->runId() << " Evt: " << event->eventId() << " Lumi: " << event->lumiBlockId() << " mm" << "\n";
													//std::cout<<"Processing the "<<ievt<<"th event" << endl;
													//cout << "is same-sign muon!" << endl;
													//cout << "-> muon 1 pt: " << selectedMuons[0]->Pt() << endl;
													//cout << "-> muon 2 pt: " << selectedMuons[1]->Pt() << endl;
												}
											}
										}
										
										//// for data-driven part: same-sign muons with 1 loose and 1 tight muon 
										else if(datadriven && selectedMuons.size() == 1 && selectedOnlyLooseMuons_FL.size() == 1  && selectedLooseElectronsVBTFid.size() == 0)
										{
											if(option =="FakeLepton")
											{
												//require that there are the two muons do not form the Z mass
												if( !selectionFakeLepton.foundZCandidate(selectedMuons_FL, selectedOnlyLooseMuons_FL, 10.) )
												{
													if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
													{
														bool eventBtag = false;
														for(unsigned int j = 0; j < selectedJets.size(); j++)
															if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
																eventBtag = true;
														if(selectedMuons[0]->charge() == selectedOnlyLooseMuons_FL[0]->charge() && eventBtag) 
															NbSSLooseMuonTightMuon_data++;
													}
										 		
												}
											}
										}
										
										////for same-sign muon and electron require exactly 1 muon, veto for other loose muons
										else if(selectedMuons.size() == 1 && selectedLooseMuons.size() == selectedMuons.size() && selectedElectrons.size() == 1 && selectedLooseElectronsVBTFid.size() == selectedElectrons.size())
										{
				    						//require that there are no two electrons forming the Z mass
												if(!selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid,10.))
												{
													//it should not be an electron from a conversion!
													if( selection.passConversionRejection(selectedElectrons[0]) )
													{
														selecTableChargeMisId_ElMu.Fill(d,0,scaleFactor);
														if(fabs(selectedElectrons[0]->superClusterEta())<1.4442){
															selecTableChargeMisId_ElMu.Fill(d,1,scaleFactor);
														}else if(fabs(selectedElectrons[0]->superClusterEta())>1.5660){
															selecTableChargeMisId_ElMu.Fill(d,2,scaleFactor);
														}
														isDiLepton = true;																											
														
														if(option=="ChargeMisId")
														{
															if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
															{
																bool eventBtag = false;
																for(unsigned int j = 0; j < selectedJets.size(); j++)
																	if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
																		eventBtag = true;
																if(selectedElectrons[0]->charge() != selectedMuons[0]->charge() && eventBtag)
																{
																	NbOSElMu_data++;
																	if(fabs(selectedElectrons[0]->superClusterEta())<1.4442) NbOSElMu_EB_data++;
																	else if(fabs(selectedElectrons[0]->superClusterEta())>1.5660) NbOSElMu_EE_data++;
																}
															}
														}
														//require the same charge for muon and electron
														if(selectedElectrons[0]->charge()== selectedMuons[0]->charge())
														{
															isSSLepton = true; // we have a same-sign electron and muon
															isSSMuEl = true;
															//std::cout<<"Processing the "<<ievt<<"th event" << endl;
															//cout << "is same-sign muon+electron!" << endl;
															//cout << "-> muon pt: " << selectedMuons[0]->Pt() << endl;
															//cout << "-> electron pt: " << selectedElectrons[0]->Pt() << endl;
															if(fabs(selectedElectrons[0]->superClusterEta())<1.4442){
																selecTableChargeMisId_ElMu.Fill(d,3,scaleFactor);
																isEB = true;
															}else if(fabs(selectedElectrons[0]->superClusterEta())>1.5660){
																selecTableChargeMisId_ElMu.Fill(d,4,scaleFactor);
																isEE = true;
															}
															
															bool eventBtag = false;
															for(unsigned int j = 0; j < selectedJets.size(); j++)
																if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
																	eventBtag = true;
															if(eventBtag && dataSetName=="Data" && systematic=="Nominal")
																myfileSS << "Run: " << event->runId() << " Evt: " << event->eventId() << " Lumi: " << event->lumiBlockId() << " me" << "\n";
														}else{ //opposite charge!!!
															if(fabs(selectedElectrons[0]->superClusterEta())<1.4442){
																selecTableChargeMisId_ElMu.Fill(d,5,scaleFactor);
															}else if(fabs(selectedElectrons[0]->superClusterEta())>1.5660){
																selecTableChargeMisId_ElMu.Fill(d,6,scaleFactor);
															}
														}
													}
												}
										}
										
										//// for data-driven part: same-sign muon and electron with 1 tight muon and 1 loose electron 
										else if(datadriven && selectedMuons.size() == 1 && selectedOnlyLooseElectrons_FL.size() == 1 && selectedElectrons.size()==0  && selectedOnlyLooseMuons_FL.size() == 0)
										{
											if(option =="FakeLepton")
											{
												//it should not be an electron from a conversion!
												if( selectionFakeLepton.passConversionRejection(selectedOnlyLooseElectrons_FL[0]) )
												{
													if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
													{
														bool eventBtag = false;
														for(unsigned int j = 0; j < selectedJets.size(); j++)
															if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
																eventBtag = true;
														if(selectedMuons[0]->charge() == selectedOnlyLooseElectrons_FL[0]->charge() && eventBtag) 
															NbSSLooseElectronTightMuon_data++;
													}
												}
											}
										}
										
										////three leptons
										
										//require at least 3 muons, veto for very loose electrons
										else if(selectedMuons.size() == 3  && selectedLooseMuons.size() == selectedMuons.size() && selectedLooseElectronsVBTFid.size() == 0)
										{
											//require that there are no two muons forming the Z mass
											if( !selection.foundZCandidate(selectedMuons, selectedMuons, 10.) )
											{
												isTriLepton = true; // at least three muons
												isTriMuon = true;
												//std::cout<<"Processing the "<<ievt<<"th event" << endl;
												bool eventBtag = false;
												for(unsigned int j = 0; j < selectedJets.size(); j++)
													if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
														eventBtag = true;
												if(eventBtag && dataSetName=="Data" && systematic=="Nominal")
													myfileLLL << "Run: " << event->runId() << " Evt: " << event->eventId() << " Lumi: " << event->lumiBlockId() << " mmm" << "\n";
												//cout << "is trilepton: 3 muons!" << endl;
												//cout << "-> muon 1 pt: " << selectedMuons[0]->Pt() << endl;
												//cout << "-> muon 2 pt: " << selectedMuons[1]->Pt() << endl;
												//cout << "-> muon 3 pt: " << selectedMuons[2]->Pt() << endl;
											}
					 					}
										
										//require 2 muons and an electron
										else if(selectedMuons.size() == 2  && selectedLooseMuons.size() == selectedMuons.size() && selectedElectrons.size() == 1 && selectedLooseElectronsVBTFid.size() == selectedElectrons.size())
										{
											//require that there are no two muons forming the Z mass
											if( !selection.foundZCandidate(selectedMuons, selectedMuons, 10.) )
											{
													//require that there are no two electrons forming the Z mass
													if(!selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid,10.))
													{
														//it should not be an electron from a conversion!
														if( selection.passConversionRejection(selectedElectrons[0]) )
														{
															isTriLepton = true; //at least two muons and one electron
															isTriMu2El1 = true;
															bool eventBtag = false;
															for(unsigned int j = 0; j < selectedJets.size(); j++)
																if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
																	eventBtag = true;
															if(eventBtag && dataSetName=="Data" && systematic=="Nominal")
																myfileLLL << "Run: " << event->runId() << " Evt: " << event->eventId() << " Lumi: " << event->lumiBlockId() << " mme" << "\n";
															//std::cout<<"Processing the "<<ievt<<"th event" << endl;
															//cout << "is trilepton: 2 muons + 1 electron!" << endl;
															//cout << "-> muon 1 pt: " << selectedMuons[0]->Pt() << endl;
															//cout << "-> muon 2 pt: " << selectedMuons[1]->Pt() << endl;
															//cout << "-> electron pt: " << selectedElectrons[0]->Pt() << endl;
														}
													}
											}
										}
										
										//require exactly 1 muon and 2 electrons
										else if(selectedMuons.size() == 1 && selectedLooseMuons.size() == selectedMuons.size() && selectedElectrons.size() == 2 && selectedLooseElectronsVBTFid.size() == selectedElectrons.size())
										{
												//require that there are no two electrons forming the Z mass
												if(!selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid, 10.))
												{
													if(selection.passConversionRejection(selectedElectrons[0]) && selection.passConversionRejection(selectedElectrons[1]))
													{
														isTriLepton = true; //one muon and at least two electrons
														isTriMu1El2 = true;
														bool eventBtag = false;
														for(unsigned int j = 0; j < selectedJets.size(); j++)
															if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
																eventBtag = true;
														if(eventBtag && dataSetName=="Data" && systematic=="Nominal")
															myfileLLL << "Run: " << event->runId() << " Evt: " << event->eventId() << " Lumi: " << event->lumiBlockId() << " mee" << "\n";
														//std::cout<<"Processing the "<<ievt<<"th event" << endl;
														//cout << "is trilepton: 1 muon + 2 electrons!" << endl;
														//cout << "-> muon pt: " << selectedMuons[0]->Pt() << endl;
														//cout << "-> electron 1 pt: " << selectedElectrons[0]->Pt() << endl;
														//cout << "-> electron 2 pt: " << selectedElectrons[1]->Pt() << endl;
													}
												}
										}
									} // end MET cut
								//} // end requirement of at least a b-tagged jet larger than a certain pre-defined cut																		
							//} // end 'loop' on jets
						} //end 'at least one jet'  
          } // end if selectedMuons.size()>=1
        } // end good PV
      }// end trigged & semiMuon
      
			///// EVENTS TRIGGERED BY ELECTRON TRIGGER
			else if(trigged && semiElectron)
      {
        selecTableSemiLep.Fill(d,2,scaleFactor);
    	 	
        if( isGoodPV )
        {
          selecTableSemiLep.Fill(d,3,scaleFactor);
          if( selectedElectrons.size() >= 1 && selectedElectrons[0]->Pt()>40)
          {
              if( selection.passConversionRejection(selectedElectrons[0]) )
              {
            		selecTableSemiLep.Fill(d,4,scaleFactor);
								sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class
							
								if( selectedJets.size()>=(unsigned int)anaEnv.NofJets)
								{
									//block for the jet multiplicity plot
									if(mets[0]->Et()> METCut)
									{							
										if(selectedElectrons.size() == 1 && !selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid, 10.) && selectedLooseMuons.size() == 0 && selectedLooseElectronsVBTFid.size() == selectedElectrons.size())
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
									
									//for(unsigned int j=0;j<selectedJets.size();j++)
									//{
										//now require at least a b-tagged jet larger than a certain pre-defined cut//update: not required offline in treecreator
										//if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue && !eventSelected)
										//{
		             			selecTableSemiLep.Fill(d,5,scaleFactor);
											if(mets[0]->Et()> METCut)
											{
			        					selecTableSemiLep.Fill(d,6,scaleFactor);
                        eventSelected = true;
		
												//cout << "event is selected according to the baseline selection!" << endl;
														
												//// single electron
												if(selectedElectrons.size() == 1 && !selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid, 10.) && selectedLooseMuons.size() == 0 && selectedLooseElectronsVBTFid.size() == selectedElectrons.size())
												{
													isSingleLepton = true;
													isSingleElectron = true;
													//std::cout<<"Processing the "<<ievt<<"th event" << endl;
													//cout << "is single electron!" << endl;
													//cout << "-> electron pt: " << selectedElectrons[0]->Pt() << endl;
												}
												
												//// two same-sign electrons
												else if(selectedElectrons.size() == 2 && selectedLooseElectronsVBTFid.size() == selectedElectrons.size() && selectedLooseMuons.size() == 0)
												{
													if(selection.passConversionRejection(selectedElectrons[1]))
													{
														if(!selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid, 10.))
														{ 
															selecTableChargeMisId_2El.Fill(d,0,scaleFactor);
															
															if(fabs(selectedElectrons[0]->superClusterEta())<1.4442 && fabs(selectedElectrons[1]->superClusterEta())<1.4442){
																selecTableChargeMisId_2El.Fill(d,1,scaleFactor);
															}else if(fabs(selectedElectrons[0]->superClusterEta())>1.5660 && fabs(selectedElectrons[1]->superClusterEta())>1.5660){
																selecTableChargeMisId_2El.Fill(d,2,scaleFactor);
															}else if((fabs(selectedElectrons[0]->superClusterEta())>1.5660 && fabs(selectedElectrons[1]->superClusterEta())<1.4442) || (fabs(selectedElectrons[1]->superClusterEta())>1.5660 && fabs(selectedElectrons[0]->superClusterEta())<1.4442)){
																selecTableChargeMisId_2El.Fill(d,3,scaleFactor);
															}																
															isDiLepton = true;
															
															if(option=="ChargeMisId")
															{
																if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
																{
																	bool eventBtag = false;
																	for(unsigned int j = 0; j < selectedJets.size(); j++)
																		if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
																			eventBtag = true;
																	if(selectedElectrons[0]->charge() != selectedElectrons[1]->charge() && eventBtag)
																	{ 
																		NbOSElectrons_data++;
																		if(fabs(selectedElectrons[0]->superClusterEta())<1.4442 && fabs(selectedElectrons[1]->superClusterEta())<1.4442) NbOSElectrons_EBEB_data++;
																		else if(fabs(selectedElectrons[0]->superClusterEta())>1.5660 && fabs(selectedElectrons[1]->superClusterEta())>1.5660) NbOSElectrons_EEEE_data++;
																		else if((fabs(selectedElectrons[0]->superClusterEta())>1.5660 && fabs(selectedElectrons[1]->superClusterEta())<1.4442) || (fabs(selectedElectrons[1]->superClusterEta())>1.5660 && fabs(selectedElectrons[0]->superClusterEta())<1.4442)) NbOSElectrons_EBEE_data++;
																	
																	}
																}
															}
															
															if(selectedElectrons[0]->charge()== selectedElectrons[1]->charge())
															{ 
																isSSLepton = true;
																isSSElectron = true;
																//std::cout<<"Processing the "<<ievt<<"th event" << endl;
																//cout << "is same-sign electron!" << endl;
																//cout << "-> electron 1 pt: " << selectedElectrons[0]->Pt() << endl;
																//cout << "-> electron 2 pt: " << selectedElectrons[1]->Pt() << endl;
																if(fabs(selectedElectrons[0]->superClusterEta())<1.4442 && fabs(selectedElectrons[1]->superClusterEta())<1.4442){
																	selecTableChargeMisId_2El.Fill(d,4,scaleFactor);
																	isEBEB = true;
																}else if(fabs(selectedElectrons[0]->superClusterEta())>1.5660 && fabs(selectedElectrons[1]->superClusterEta())>1.5660){
																	selecTableChargeMisId_2El.Fill(d,5,scaleFactor);
																	isEEEE = true;
																}else if((fabs(selectedElectrons[0]->superClusterEta())>1.5660 && fabs(selectedElectrons[1]->superClusterEta())<1.4442) || (fabs(selectedElectrons[1]->superClusterEta())>1.5660 && fabs(selectedElectrons[0]->superClusterEta())<1.4442)){
																	selecTableChargeMisId_2El.Fill(d,6,scaleFactor);
																	isEBEE = true;
																}																
																
																bool eventBtag = false;
																for(unsigned int j = 0; j < selectedJets.size(); j++)
																	if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
																	eventBtag = true;
																if(eventBtag && dataSetName=="Data" && systematic=="Nominal")
																	myfileSS << "Run: " << event->runId() << " Evt: " << event->eventId() << " Lumi: " << event->lumiBlockId() << " ee" << "\n";
															}else{ //opposite charge!!!
																if(fabs(selectedElectrons[0]->superClusterEta())<1.4442 && fabs(selectedElectrons[1]->superClusterEta())<1.4442){
																	selecTableChargeMisId_2El.Fill(d,7,scaleFactor);
																}else if(fabs(selectedElectrons[0]->superClusterEta())>1.5660 && fabs(selectedElectrons[1]->superClusterEta())>1.5660){
																	selecTableChargeMisId_2El.Fill(d,8,scaleFactor);
																}else if((fabs(selectedElectrons[0]->superClusterEta())>1.5660 && fabs(selectedElectrons[1]->superClusterEta())<1.4442) || (fabs(selectedElectrons[1]->Eta())>1.5660 && fabs(selectedElectrons[0]->Eta())<1.4442)){
																	selecTableChargeMisId_2El.Fill(d,9,scaleFactor);
																}																
															}															
														}
													}
												}
										
												//// for data-driven part: same-sign electrons with 1 loose and 1 tight electron 
												else if(datadriven && selectedElectrons.size() == 1 && selectedOnlyLooseElectrons_FL.size() == 1  && selectedLooseMuons_FL.size() == 0)
												{
													if(option =="FakeLepton")
													{
														//require that there are the two electrons do not form the Z mass
														if( !selectionFakeLepton.foundZCandidate(selectedElectrons_FL, selectedOnlyLooseElectrons_FL, 10.) )
														{
															if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
															{
																bool eventBtag = false;
																for(unsigned int j = 0; j < selectedJets.size(); j++)
																	if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
																		eventBtag = true;
																if(selectedElectrons[0]->charge() == selectedOnlyLooseElectrons_FL[0]->charge() && eventBtag) 
																	NbSSLooseElectronTightElectron_data++;
															}
										 		
														}
													}
												}
												
												//// a same-sign electron and muon
												else if(selectedElectrons.size() == 1 && !selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid, 10.) && selectedLooseElectronsVBTFid.size() == selectedElectrons.size() && selectedMuons.size() == 1)
												{
													if(selectedMuons[0]->Pt()<40 && selectedLooseMuons.size() == selectedMuons.size())
													{
														selecTableChargeMisId_ElMu.Fill(d,0,scaleFactor);
														if(fabs(selectedElectrons[0]->superClusterEta())<1.4442){
															selecTableChargeMisId_ElMu.Fill(d,1,scaleFactor);
														}else if(fabs(selectedElectrons[0]->superClusterEta())>1.5660){
															selecTableChargeMisId_ElMu.Fill(d,2,scaleFactor);
														}
														isDiLepton = true;																										
														
														if(option=="ChargeMisId")
														{
															if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
															{
																bool eventBtag = false;
																for(unsigned int j = 0; j < selectedJets.size(); j++)
																	if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
																		eventBtag = true;
																if(selectedElectrons[0]->charge() != selectedMuons[0]->charge() && eventBtag)
																{ 
																	NbOSElMu_data++;
																	if(fabs(selectedElectrons[0]->superClusterEta())<1.4442) NbOSElMu_EB_data++;
																	else if(fabs(selectedElectrons[0]->superClusterEta())>1.5660) NbOSElMu_EE_data++;
																}
															}
														}	
														if(selectedElectrons[0]->charge()== selectedMuons[0]->charge())
														{
															isSSLepton = true;
															isSSMuEl = true;
															//std::cout<<"Processing the "<<ievt<<"th event" << endl;
															//cout << "is same-sign electron + muon!" << endl;
															//cout << "-> electron pt: " << selectedElectrons[0]->Pt() << endl;
															//cout << "-> muon pt: " << selectedMuons[0]->Pt() << endl;
															if(fabs(selectedElectrons[0]->superClusterEta())<1.4442){
																selecTableChargeMisId_ElMu.Fill(d,3,scaleFactor);
																isEB = true;
															}else if(fabs(selectedElectrons[0]->superClusterEta())>1.5660){
																selecTableChargeMisId_ElMu.Fill(d,4,scaleFactor);
																isEE = true;
															}																
															bool eventBtag = false;
															for(unsigned int j = 0; j < selectedJets.size(); j++)
																if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
																eventBtag = true;
															if(eventBtag && dataSetName=="Data" && systematic=="Nominal")
																myfileSS << "Run: " << event->runId() << " Evt: " << event->eventId() << " Lumi: " << event->lumiBlockId() << " me" << "\n";
														}else{ //opposite charge!!!
															if(fabs(selectedElectrons[0]->superClusterEta())<1.4442){
																selecTableChargeMisId_ElMu.Fill(d,5,scaleFactor);
															}else if(fabs(selectedElectrons[0]->superClusterEta())>1.5660){
																selecTableChargeMisId_ElMu.Fill(d,6,scaleFactor);
															}
														}
													}
												}
												
												//// for data-driven part: same-sign muon and electron with 1 tight electron and 1 loose muon 
												else if(datadriven && selectedElectrons.size() == 1 && selectedOnlyLooseMuons_FL.size() == 1 && selectedMuons.size()==0  && selectedOnlyLooseElectrons_FL.size() == 0)
												{
													if(option =="FakeLepton")
													{
															if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
															{
																bool eventBtag = false;
																for(unsigned int j = 0; j < selectedJets.size(); j++)
																	if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
																		eventBtag = true;
																if(selectedElectrons[0]->charge() == selectedOnlyLooseMuons_FL[0]->charge() && eventBtag) 
																	NbSSLooseMuonTightElectron_data++;
															}
													}
												}
											
												//// 1 electron and 2 muons
												else if(selectedElectrons.size() == 1 && !selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid, 10.) && selectedLooseElectronsVBTFid.size() == selectedElectrons.size())
												{
													if(selectedMuons.size() == 2 && selectedMuons[1]->Pt()<40 && selectedLooseMuons.size() == selectedMuons.size())
													{
														isTriLepton = true;
														isTriMu2El1 = true;
														bool eventBtag = false;
														for(unsigned int j = 0; j < selectedJets.size(); j++)
															if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
															eventBtag = true;
														if(eventBtag && dataSetName=="Data" && systematic=="Nominal")
															myfileLLL << "Run: " << event->runId() << " Evt: " << event->eventId() << " Lumi: " << event->lumiBlockId() << " mme" << "\n";
														//std::cout<<"Processing the "<<ievt<<"th event" << endl;
														//cout << "is tri-lepton: 1 electron and 2 muons!" << endl;
														//cout << "-> electron pt: " << selectedElectrons[0]->Pt() << endl;
														//cout << "-> muon 1 pt: " << selectedMuons[0]->Pt() << endl;
														//cout << "-> muon 2 pt: " << selectedMuons[1]->Pt() << endl;
													}
												}
												
												//// 2 electrons and 1 muon
												else if(selectedElectrons.size() == 2 && !selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid, 10.) && selectedLooseElectronsVBTFid.size() == selectedElectrons.size())
												{
													if(selection.passConversionRejection(selectedElectrons[1])&& selectedElectrons[1]->Pt()<40)
													{
														if(selectedMuons.size() == 1 && selectedMuons[0]->Pt()<40 && selectedLooseMuons.size() == selectedMuons.size())
														{
															isTriLepton = true;
															isTriMu1El2 = true;
															bool eventBtag = false;
															for(unsigned int j = 0; j < selectedJets.size(); j++)
																if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
																eventBtag = true;
															if(eventBtag && dataSetName=="Data" && systematic=="Nominal")
																myfileLLL << "Run: " << event->runId() << " Evt: " << event->eventId() << " Lumi: " << event->lumiBlockId() << " mee" << "\n";
															//std::cout<<"Processing the "<<ievt<<"th event" << endl;
															//cout << "is tri-lepton: 2 electrons and 1 muon!" << endl;
															//cout << "-> electron 1 pt: " << selectedElectrons[0]->Pt() << endl;
															//cout << "-> electron 2 pt: " << selectedElectrons[1]->Pt() << endl;
															//cout << "-> muon 1 pt: " << selectedMuons[0]->Pt() << endl;
														}
													}
												}
												
												//// three electrons
												else if(selectedElectrons.size() == 3 && !selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid, 10.) && selectedLooseElectronsVBTFid.size() == selectedElectrons.size())
												{
													if(selection.passConversionRejection(selectedElectrons[1]) && selection.passConversionRejection(selectedElectrons[2]))
													{
														isTriLepton = true; // at least three electrons
														isTriElectron = true;
														bool eventBtag = false;
														for(unsigned int j = 0; j < selectedJets.size(); j++)
															if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
															eventBtag = true;
														if(eventBtag && dataSetName=="Data" && systematic=="Nominal")
															myfileLLL << "Run: " << event->runId() << " Evt: " << event->eventId() << " Lumi: " << event->lumiBlockId() << " eee" << "\n";
														//std::cout<<"Processing the "<<ievt<<"th event" << endl;
														//cout << "is trilepton: 3 electrons!" << endl;
														//cout << "-> electron 1 pt: " << selectedElectrons[0]->Pt() << endl;
														//cout << "-> electron 2 pt: " << selectedElectrons[1]->Pt() << endl;
														//cout << "-> electron 3 pt: " << selectedElectrons[2]->Pt() << endl;
													}
												}
											} // end MET cut
										//} // end requirement of at least a b-tagged jet larger than a certain pre-defined cut
									//} // end 'loop' on jets
								} // end 'at least one jet'
							} // end conversion rejection for leading electron
          } // end if selectedElectrons.size()>=1
        } // end good PV
      } // end trigged & semiElectron
						
      //if(!isSingleLepton && !isSSLepton && !isTriLepton) continue; //same as all cuts just above (baseline selection is there) 
			if(!isSingleLepton && !isDiLepton && !isTriLepton) continue; //all dilepton events (SS and OS) will be stored in the trees
			if(option=="ChargeMisId" || option=="ElectronFake" || option=="MuonFake") continue;

			MSPlot["MS_nPV"]->Fill(vertex.size(),datasets[d], true, Luminosity*scaleFactor);				

			if(isSingleLepton){
				MSPlot["MS_MET"]->Fill(mets[0]->Et(),datasets[d], true, Luminosity*scaleFactor);				
				if(semiElectron) MSPlot["MS_LeptonPt"]->Fill(selectedElectrons[0]->Pt(),datasets[d], true, Luminosity*scaleFactor);				
				if(semiMuon) MSPlot["MS_LeptonPt"]->Fill(selectedMuons[0]->Pt(),datasets[d], true, Luminosity*scaleFactor);				
				
				float relIso;
				if(semiElectron) relIso = (selectedElectrons[0]->chargedHadronIso()+selectedElectrons[0]->neutralHadronIso()+selectedElectrons[0]->photonIso())/selectedElectrons[0]->Pt();
				if(semiMuon) relIso = (selectedMuons[0]->chargedHadronIso()+selectedMuons[0]->neutralHadronIso()+selectedMuons[0]->photonIso())/selectedMuons[0]->Pt();      		  				
				MSPlot["MS_LeptonRelIso"]->Fill(relIso,datasets[d], true, Luminosity*scaleFactor);
				
				for(unsigned int j=0;j<selectedJets.size();j++)
				{
				  MSPlot["MS_JetPt_all_SingleLepton"]->Fill(selectedJets[j]->Pt(),datasets[d], true, Luminosity*scaleFactor);
					if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
					  MSPlot["MS_JetPt_btagged_SingleLepton"]->Fill(selectedJets[j]->Pt(),datasets[d], true, Luminosity*scaleFactor);
					else
					  MSPlot["MS_JetPt_nonbtagged_SingleLepton"]->Fill(selectedJets[j]->Pt(),datasets[d], true, Luminosity*scaleFactor);
				}			
				selecTableSemiLep.Fill(d,7,scaleFactor);
			}

			if(isSSLepton)
			{
				//cout << "IS SAME-SIGN LEPTON EVENT" << endl;
				NbSSevents = NbSSevents + datasets[d]->NormFactor()*Luminosity*scaleFactor;
				selecTableMultiLepton.Fill(d,0,scaleFactor);
			}
			if(isTriLepton)
			{
				//cout << "IS TRI-LEPTON EVENT" << endl;
				NbTrievents = NbTrievents + datasets[d]->NormFactor()*Luminosity*scaleFactor;
				selecTableMultiLepton.Fill(d,1,scaleFactor);				
			}
			
			
			//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before treeLoader.LoadGenEvent:");					
			TRootGenEvent* genEvt = 0;
			if((dataSetName.find("TTbarJets_SemiMu") == 0 && semiMuon) || (dataSetName.find("TTbarJets_SemiElectron") == 0 && semiElectron))
			  genEvt = treeLoader.LoadGenEvent(ievt,false);
			//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After treeLoader.LoadGenEvent:");		
			
			
			bool Wbosonpartonsmatched = false; // True if the Wboson ttbar semi-mu or semi-el partons are matched to their 2 jets (not necessarily the 4 highest pt jets)
			float WMassmatched_ = -9999;
			//////////////////////////////////////////////////////////////////////////
      // jet-parton matching needed to make the Wmassplot for the W counting
      //////////////////////////////////////////////////////////////////////////
      sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class)   
			
			if(selectedJets.size()>=4 && ((dataSetName.find("TTbarJets_SemiMu") == 0 && semiMuon) || (dataSetName.find("TTbarJets_SemiElectron") == 0 && semiElectron)))
			{
			
			  pair<unsigned int, unsigned int> leptonicBJet_, hadronicBJet_, hadronicWJet1_, hadronicWJet2_; //First index is the JET number, second one is the parton
			  leptonicBJet_ = hadronicBJet_ = hadronicWJet1_ = hadronicWJet2_ = pair<unsigned int,unsigned int>(9999,9999);
			
				//sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class)     
				int MCPermutation[4]; for (unsigned int i=0;i<4;i++) MCPermutation[i] = -1;

				vector<TRootMCParticle*> mcParticlesMatching_; 
				vector<TRootJet*> selectedJets_; 

				mcParticlesMatching_.clear();
				selectedJets_.clear();
      
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

      	if(Wbosonpartonsmatched)
	  		{
					WMassmatched_ = (*selectedJets[hadronicWJet1_.first]+*selectedJets[hadronicWJet2_.first]).M();
      		histo1D["hadronicRecoWMass"]->Fill(WMassmatched_,scaleFactor);
	  		}	
			
			} //end selectedJets.size()>=4 !useMassesAndResolutions && dataSet semiMu or semiEl

     
		  //now do some stuff for later use in the MVA (obtain the quarks that should be used for matching later on)
			//bool all4JetsMatched_MCdef = false; // True if the 4 highest pt jets are matched to the 4 ttbar semi-lep partons
  		//bool hadronictopJetsMatched_MCdef = false;
			vector<TLorentzVector> mcQuarksForMatching; //ordering: hadronicWQuark1,hadronicWQuark2,hadronicbQuark,leptonicbQuark
			int pdgID_top = 6; //top quark
			bool TprimePairSample = false;
			if(dataSetName.find("NP_Tprime") <= 0 || dataSetName.find("NP_overlay_Tprime") <= 0)
			  TprimePairSample = true; //tprime pair
  		if(TprimePairSample)
    		pdgID_top = 8; //4th generation t' quark
			if(selectedJets.size()>=4 && ((dataSetName.find("TTbarJets_SemiMu") == 0 && semiMuon) || (dataSetName.find("TTbarJets_SemiElectron") == 0 && semiElectron) || TprimePairSample))
			{
    		TLorentzVector hadronicWQuark1,hadronicWQuark2,hadronicbQuark,leptonicbQuark;
    		bool muPlusFromTop = false, muMinusFromTop = false, elPlusFromTop = false, elMinusFromTop = false;
    
    		for(unsigned int i=0; i<mcParticles.size(); i++)
    		{
    		  //cout << i << ":  status: " << mcParticles[i]->status() << "  pdgId: " << mcParticles[i]->type()
    		  //  << "  motherPdgId: " << mcParticles[i]->motherType() << "  grannyPdgId: " << mcParticles[i]->grannyType() << endl;
    		  if( mcParticles[i]->status() != 3) continue;
          
      		if( mcParticles[i]->type() == 13 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -pdgID_top )
      		  muMinusFromTop = true;
      		if( mcParticles[i]->type() == -13 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == pdgID_top )
      		  muPlusFromTop = true;
      		if( mcParticles[i]->type() == 11 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -pdgID_top )
      		  elMinusFromTop = true;
      		if( mcParticles[i]->type() == -11 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == pdgID_top )
      		  elPlusFromTop = true;
				}
				
				bool WQuark1Found = false,WQuark2Found = false,hadbQuarkFound=false,lepbQuarkFound=false;
        for(unsigned int i=0; i<mcParticles.size(); i++)
    		{
			  	if( mcParticles[i]->status() != 3) continue;
					
      		if( abs(mcParticles[i]->type()) < 6) //light/b quarks, 6 should stay hardcoded (not to be changed when considering t' pair)
      		{
					  if( ((( muPlusFromTop || elPlusFromTop ) && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -pdgID_top) 
						  || (( muMinusFromTop || elMinusFromTop ) && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == pdgID_top)) 
						&& !WQuark1Found)
        		{							  											
							hadronicWQuark1 = (TLorentzVector) *mcParticles[i];
							//cout<<"hadronicWQuark1.Pt() = "<<hadronicWQuark1.Pt()<<endl;
							WQuark1Found = true;
							//cout<<"  WQuark1Found"<<endl;
						}
						else if ( (( muPlusFromTop || elPlusFromTop ) && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -pdgID_top) 
						  || (( muMinusFromTop || elMinusFromTop ) && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == pdgID_top ))
					  {
							hadronicWQuark2 = (TLorentzVector) *mcParticles[i];
							//cout<<"hadronicWQuark2.Pt() = "<<hadronicWQuark2.Pt()<<endl;
							WQuark2Found = true;
							//cout<<"  WQuark2Found"<<endl;
						}
						else if( ( ( muPlusFromTop || elPlusFromTop ) && mcParticles[i]->motherType() == -pdgID_top )
          		|| ( ( muMinusFromTop || elMinusFromTop ) && mcParticles[i]->motherType() == pdgID_top ) )
						{
						  if( abs(mcParticles[i]->type()) == 5)
							{
							  hadronicbQuark = (TLorentzVector) *mcParticles[i];
								//cout<<"hadronicbQuark.Pt() = "<<hadronicbQuark.Pt()<<endl;
								hadbQuarkFound = true;
								//cout<<"  hadbQuarkFound"<<endl;
							}
						}
						else if( ( ( muPlusFromTop || elPlusFromTop ) && mcParticles[i]->motherType() == pdgID_top )
          		|| ( ( muMinusFromTop || elMinusFromTop ) && mcParticles[i]->motherType() == -pdgID_top ) )
						{
							if( abs(mcParticles[i]->type()) == 5)
							{
								leptonicbQuark = (TLorentzVector) *mcParticles[i];
								//cout<<"leptonicbQuark.Pt() = "<<leptonicbQuark.Pt()<<endl;
								lepbQuarkFound = true;
								//cout<<"  lepbQuarkFound"<<endl;
							}
						}			
					}
				}
				
				if(WQuark1Found && WQuark2Found && hadbQuarkFound && lepbQuarkFound) //For ttjets semi-lep this will always be found (right?), and for t' pair by this you make sure that only the quarks for the semilep case is stored
				{
					//respect this ordering
					mcQuarksForMatching.push_back(hadronicWQuark1);
					mcQuarksForMatching.push_back(hadronicWQuark2);
					mcQuarksForMatching.push_back(hadronicbQuark);
					mcQuarksForMatching.push_back(leptonicbQuark);
				}
			} //end selectedJets.size()>=4 && (dataSetName.find("TTbarJets_SemiMu") == 0 || dataSetName.find("TTbarJets_SemiElectron") == 0 || TprimePairSample)



			// find quarks from W bosons 
      int nW=0;
			vector < TRootMCParticle *> mcWbosons; 
			for(unsigned int i=0; i<mcParticles.size(); i++)
			{
				if(mcParticles[i]->status() != 3) continue;
				
				if(abs(mcParticles[i]->type()) == 24) 
				{
					//cout << "THIS W BOSON HAS " << mcParticles[i]->nDau() << " DAUGHTERS" << endl;
					if(mcParticles[i]->nDau() <= 3){ //crucial to avoid W bosons "inside" diagrams instead of the actual decay product
						nW++;
						//cout << "first daughter " << mcParticles[i]->dauOneId() << endl;
						//cout << "second daughter " << mcParticles[i]->dauTwoId() << endl;
						//cout << "third daughter " << mcParticles[i]->dauThreeId() << endl;
						//cout << "fourth daughter " << mcParticles[i]->dauFourId() << endl;
						mcWbosons.push_back(mcParticles[i]);
					}
				}
			}	
			//cout << "EVENT HAS " << nW << " or " << mcWbosons.size() << " W BOSONS ON GENERATOR LEVEL "<< endl;
			
			vector < TLorentzVector > quarksFromW;
			TLorentzVector qW_tmp; 
			int Wboson = 0;
			if(dataSetName != "data" && dataSetName != "Data" && dataSetName != "DATA")
			{
		  	//cout << "mcParticles.size " << mcParticles.size() << endl;
        for(unsigned int i=0; i<mcParticles.size(); i++)
    		{
					if( mcParticles[i]->status() != 3) continue;
        	for(unsigned int j=0; j < i; j++)
    			{
			  		if( mcParticles[j]->status() != 3) continue;
					
      			if( abs(mcParticles[i]->type()) < 6 && abs(mcParticles[j]->type()) < 6) //light/b quarks
      			{
							if(mcParticles[i]->motherType() == -24 || mcParticles[i]->motherType() == 24)
							{
								TLorentzVector W = (TLorentzVector) *mcParticles[i]+(TLorentzVector) *mcParticles[j];
								if( mcParticles[i]->motherType() == mcParticles[j]->motherType() && mcParticles[i]->grannyType() == mcParticles[j]->grannyType())
								{// same W boson as mother and same granny!
									//cout << "W boson mass from quarks "<< W.M() << endl;
									int thew = 999;
									for(int w = 0; w < mcWbosons.size(); w++){
										//cout << "W boson mass from generator W " << mcWbosons[w]->M() << endl;
										if(fabs(W.M()-mcWbosons[w]->M())<0.0001) thew = w; //important to avoid that a single quark makes multiple W bosons
									}
									if(thew<999){
										//cout << "grannyType: "<< mcParticles[i]->grannyType() << endl;
										//cout << "quark type 1: " << mcParticles[i]->type() << " quark type 2: " << mcParticles[j]->type() << endl;
										Wboson++;
										qW_tmp = (TLorentzVector) *mcParticles[i];
										quarksFromW.push_back(qW_tmp);
										qW_tmp = (TLorentzVector) *mcParticles[j];
										quarksFromW.push_back(qW_tmp);
										//if(mcParticles[i]->motherType() == -24) cout << "found W boson (-) number " << Wboson << " with mass "<<W.M() <<" decaying to quarks " << i << " and " << j << endl;
										//if(mcParticles[i]->motherType() == 24) cout << "found W boson (+) number " << Wboson << " with mass "<<W.M() <<" decaying to quarks " << i << " and " << j << endl;
									}
								}
							}
						}
					}
				}				
			}
			
			//cout << "EVENT SUMMARY for dataset: " << dataSetName << endl;
			//cout << "We have "<< (float) quarksFromW.size()/2 << " hadronically decaying W bosons" << endl;
			//for(unsigned int q=0; q<quarksFromW.size(); q++)
			//{
				//cout << "quark " << q << " with pt = "<< (quarksFromW[q].first).Pt() <<" is the decay product of W " << quarksFromW[q].second << endl;
			//}
			//cout << endl;

      vector<float> bTagTCHE, bTagTCHP, InitJetsbTagTCHE, InitJetsbTagTCHP;
			vector<int> partonFlavourJet;
      vector<TLorentzVector> SelectedJetsTLV, SelectedForwardJetsTLV, SelectedMuonsTLV, SelectedElectronsTLV, InitJets;
			vector<float> SelectedMuonsRelIso, SelectedElectronsRelIso;
      for(unsigned int iJet=0; iJet<selectedJets.size(); iJet++)
      {
            SelectedJetsTLV.push_back( *selectedJets[iJet] );
            bTagTCHE.push_back(selectedJets[iJet]->btag_trackCountingHighEffBJetTags());
            bTagTCHP.push_back(selectedJets[iJet]->btag_trackCountingHighPurBJetTags());
						partonFlavourJet.push_back(selectedJets[iJet]->partonFlavour());
      }
			for(unsigned int iJet=0; iJet<selectedForwardJets.size(); iJet++)
      {
			      SelectedForwardJetsTLV.push_back( *selectedForwardJets[iJet] );			  	
			}
			for(unsigned int iJet=0; iJet<init_jets.size(); iJet++)
			{
						InitJets.push_back( *init_jets[iJet] ); //needed for MC trigger reweighting in e-channel...
						InitJetsbTagTCHE.push_back(init_jets[iJet]->btag_trackCountingHighEffBJetTags());
            InitJetsbTagTCHP.push_back(init_jets[iJet]->btag_trackCountingHighPurBJetTags());			
			}
			for(unsigned int iMuon=0; iMuon<selectedMuons.size(); iMuon++)
      {
            SelectedMuonsTLV.push_back( *selectedMuons[iMuon] );
						float relIso;
						relIso = (selectedMuons[iMuon]->chargedHadronIso()+selectedMuons[iMuon]->neutralHadronIso()+selectedMuons[iMuon]->photonIso())/selectedMuons[iMuon]->Pt();      		  
						SelectedMuonsRelIso.push_back(relIso);
			}
			for(unsigned int iElectron=0; iElectron<selectedElectrons.size(); iElectron++)
      {
            SelectedElectronsTLV.push_back( *selectedElectrons[iElectron] );
						float relIso;
						relIso = (selectedElectrons[iElectron]->chargedHadronIso()+selectedElectrons[iElectron]->neutralHadronIso()+selectedElectrons[iElectron]->photonIso())/selectedElectrons[iElectron]->Pt();      		  
						SelectedElectronsRelIso.push_back(relIso);			
      }


			if(!datadriven)
			{
				myBranch_selectedEvents = new InclFourthGenTree();
      	myBranch_selectedEvents->setEventID( event->eventId() );
      	myBranch_selectedEvents->setRunID( event->runId() );
      	myBranch_selectedEvents->setLumiBlockID( event->lumiBlockId() );
      	myBranch_selectedEvents->setNPV( vertex.size() );
      	myBranch_selectedEvents->setNPUBXm1( event->nPu(-1) );
      	myBranch_selectedEvents->setNPU( event->nPu(0) );
      	myBranch_selectedEvents->setNPUBXp1( event->nPu(1) );
				myBranch_selectedEvents->setFlavorHistoryPath( event->flavorHistoryPath() );
			
      	myBranch_selectedEvents->setSelectedSingleLepton( isSingleLepton );
				myBranch_selectedEvents->setSelectedSingleMu( isSingleMuon );
				myBranch_selectedEvents->setSelectedSingleEl( isSingleElectron );
				myBranch_selectedEvents->setSelectedSSLepton( isSSLepton );
				myBranch_selectedEvents->setSelectedSSMu( isSSMuon );
				myBranch_selectedEvents->setSelectedSSEl( isSSElectron );
				myBranch_selectedEvents->setSelectedSSMuEl( isSSMuEl );
				myBranch_selectedEvents->setSelectedTriLepton( isTriLepton );
				myBranch_selectedEvents->setSelectedMuMuMu( isTriMuon );
				myBranch_selectedEvents->setSelectedMuMuEl( isTriMu2El1 );
				myBranch_selectedEvents->setSelectedMuElEl( isTriMu1El2 );
				myBranch_selectedEvents->setSelectedElElEl( isTriElectron );
				myBranch_selectedEvents->setSelectedEE( isEE );
				myBranch_selectedEvents->setSelectedEB( isEB );
				myBranch_selectedEvents->setSelectedEBEB( isEBEB );
				myBranch_selectedEvents->setSelectedEBEE( isEBEE );
				myBranch_selectedEvents->setSelectedEEEE( isEEEE );
			
      	if(selectedJets.size()>=4 && (dataSetName.find("TTbarJets_Semi") == 0 || TprimePairSample))
      	{
			  	if(mcQuarksForMatching.size()==4) //is always the case in ttbarjets semi-lep, but for inclusive t' pair samples this is only the case when indeed all mc quarks for matching are found (ie in the semi-lep case, actually)
					{	
				  	//actually the if statement is not needed, this will be put in the branch for any event anyway (empty or not...)  
						myBranch_selectedEvents->setmcQuarksForMatching( mcQuarksForMatching );
					}
					//myBranch_selectedEvents->setAll4JetsMCMatched( all4JetsMatched_MCdef );
        	//myBranch_selectedEvents->setAllHadronicJetsMCMatched( hadronictopJetsMatched_MCdef );
					//vector<unsigned int> MatchedJetsIndices;
					//MatchedJetsIndices.push_back(hadronicWJet1_.first);
					//MatchedJetsIndices.push_back(hadronicWJet2_.first);
					//MatchedJetsIndices.push_back(hadronicBJet_.first);
					//MatchedJetsIndices.push_back(leptonicBJet_.first);
					//myBranch_selectedEvents->setMatchedJetsIndices( MatchedJetsIndices );
				}
			
				if(dataSetName.find("TTbarJets_Semi") == 0)
      	{
      		myBranch_selectedEvents->setSemiMuDecay(genEvt->isSemiLeptonic( TRootGenEvent::kMuon ));
      		myBranch_selectedEvents->setSemiElDecay(genEvt->isSemiLeptonic( TRootGenEvent::kElec ));
					myBranch_selectedEvents->setWbosonpartonsmatched(Wbosonpartonsmatched);
					if(Wbosonpartonsmatched)
				  	myBranch_selectedEvents->setWMassmatched(WMassmatched_);
      	}
				myBranch_selectedEvents->setEventWeight( scaleFactor );
				myBranch_selectedEvents->setMET( *mets[0] );
     		myBranch_selectedEvents->setSelectedJets( SelectedJetsTLV );
				myBranch_selectedEvents->setBTagTCHE( bTagTCHE );
      	myBranch_selectedEvents->setBTagTCHP( bTagTCHP );
				myBranch_selectedEvents->setpartonFlavourJet( partonFlavourJet );
				myBranch_selectedEvents->setSelectedForwardJets( SelectedForwardJetsTLV );
				myBranch_selectedEvents->setInitJets( InitJets );
				myBranch_selectedEvents->setInitJetsBTagTCHE( InitJetsbTagTCHE );
      	myBranch_selectedEvents->setInitJetsBTagTCHP( InitJetsbTagTCHP );				     						
      	myBranch_selectedEvents->setMuons( SelectedMuonsTLV );
      	myBranch_selectedEvents->setElectrons( SelectedElectronsTLV );
				myBranch_selectedEvents->setMuonsRelIso( SelectedMuonsRelIso );
      	myBranch_selectedEvents->setElectronsRelIso( SelectedElectronsRelIso );
			
				//myBranch_selectedEvents->setTopDecayedLept( topDecayedLept );
      	//myBranch_selectedEvents->setAll4JetsMCMatched( jetCombiner->All4JetsMatched_MCdef() );
      	//myBranch_selectedEvents->setAllHadronicJetsMCMatched( jetCombiner->HadronicTopJetsMatched_MCdef() );
      	//myBranch_selectedEvents->setMvaVals(mvaValsVector);
      	//myBranch_selectedEvents->setMvaResults(mvaResultsVector);
      	//myBranch_selectedEvents->setHadrBJet( hadrBJetIndex );
      	//myBranch_selectedEvents->setHadrLJet1( lightJet1Index );
      	//myBranch_selectedEvents->setHadrLJet2( lightJet2Index );
      	//myBranch_selectedEvents->setLeptBJet( leptBJetIndex );      
      	//myBranch_selectedEvents->setHadrBQuark( jetCombiner->GetHadrBQuark() );
      	//myBranch_selectedEvents->setHadrLQuark1( jetCombiner->GetLightQuark1() );
      	//myBranch_selectedEvents->setHadrLQuark2( jetCombiner->GetLightQuark2() );
      	//myBranch_selectedEvents->setLeptBQuark( jetCombiner->GetLeptBQuark() );
      	
				//cout << "EVENT SUMMARY for dataset: " << dataSetName << endl;
				//cout << "We have "<< (float) quarksFromW.size()/2 << " hadronically decaying W bosons" << endl;
				//for(unsigned int q=0; q<quarksFromW.size(); q++)
				//{
				//	cout << "quark " << q << " with pt = "<< quarksFromW[q].Pt() <<" is the decay product of W " << (int)(q/2)+1 << endl;
				//}
				//cout << endl;
				if(dataSetName != "Data" && dataSetName != "data" && dataSetName != "DATA") myBranch_selectedEvents->setQuarksFromW( quarksFromW );
			
     		myInclFourthGenTree->Fill();
	    }
			delete myBranch_selectedEvents;

	
    }//loop on events

    
    cout<<endl;
		
		
    if(!datadriven) treeFile->cd();
		
		TTree *configTreeFile = new TTree("configTreeFile","configuration Tree in tree File");
    TClonesArray* tcdatasettreefile = new TClonesArray("Dataset",1);
    configTreeFile->Branch("Dataset","TClonesArray",&tcdatasettreefile);
    TClonesArray* tcAnaEnvtreeFile = new TClonesArray("AnalysisEnvironment",1);
    configTreeFile->Branch("AnaEnv","TClonesArray",&tcAnaEnvtreeFile);
    new ((*tcAnaEnvtreeFile)[0]) AnalysisEnvironment(anaEnv);
    new ((*tcdatasettreefile)[0]) Dataset(*datasets[d]);
    
    configTreeFile->Fill();
    configTreeFile->Write();

		if(!datadriven) myInclFourthGenTree->Write();
		if(!datadriven) treeFile->Close();
		delete treeFile;
		
    if(jetTools) delete jetTools;
		
		//important: free memory
    treeLoader.UnLoadDataset();
	
    
  } //loop on datasets
	if(systematic=="Nominal") myfileSS.close();
	if(systematic=="Nominal") myfileLLL.close();

	
	//Once everything is filled ...
  cout << " We ran over all the data ;-)" << endl;
  
  ///////////////////
  // Writing
  //////////////////
  cout << " - Writing outputs to the files ..." << endl;	
	
	
    fout->cd();
    //Write histograms: MSPlots
    if(savePNG) mkdir((pathPNG+"MSPlot/").c_str(),0777);
    //cout << "mkdir " << (pathPNG+"MSPlot/").c_str()<< endl;
		cout << "Running over all MS plots" << endl;
    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
        MultiSamplePlot *temp = it->second;
        string name = it->first;
        temp->Draw(false, name, true, true, true, true, true,5,false, true, true);//(bool addRandomPseudoData, string label, bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST,int scaleNPsignal, bool addRatio, bool mergeVV, bool mergeTTV)
        temp->Write(fout, name, savePNG, pathPNG+"MSPlot/");//bool savePNG
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
			if(savePNG) tempCanvas->SaveAs( (pathPNG+it->first+".pdf").c_str() ); //well, is actually not png but pdf...
    }    
    cout << "1D plots written" << endl;
    
    //delete th1dir;
    // 2D
    TDirectory* th2dir = fout->mkdir("2D_histograms_graphs");
    th2dir->cd();
    for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
    {
			TH2F *temp = it->second;
			temp->Write();
			TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
			if(savePNG) tempCanvas->SaveAs( (pathPNG+it->first+".pdf").c_str() );
    }
    //delete th2dir;
    cout << "2D plots written" << endl;
    fout->cd();
    
    //Selection tables
    selecTableSemiLep.TableCalculator(true, true, true, true, true, true, true, true);//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST, bool mergeVV, bool mergettV, bool NP_mass)
    string selectiontableSemiLep = Treespath+"InclFourthGenSearch_SelectionTable_"+postfix+channelpostfix;
    selectiontableSemiLep = selectiontableSemiLep +".tex"; 	
    selecTableSemiLep.Write(selectiontableSemiLep.c_str(),false, true, false, false, false, false, false); //(filename, error, merged, lines, unscaled, eff, totaleff, landscape)
	    
    selecTableMultiLepton.TableCalculator(true, true, true, true, true, true, true, true);//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST, bool mergeVV, bool mergettV, bool NP_mass)
    string selectiontableMultiLepton = Treespath+"InclFourthGenSearch_SelectionTable_MultiLepton"+postfix+channelpostfix;
    selectiontableMultiLepton = selectiontableMultiLepton +".tex"; 	
    selecTableMultiLepton.Write(selectiontableMultiLepton.c_str(),false, true, false, false, false, false, false);

    selecTableChargeMisId_ElMu.TableCalculator(true, true, true, true, true, false, true, true);//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST, bool mergeVV, bool mergettV, bool NP_mass)
    string selectiontableChargeMisId_ElMu = Treespath+"InclFourthGenSearch_SelectionTable_ChargeMisIdElMu"+postfix+channelpostfix;
    selectiontableChargeMisId_ElMu = selectiontableChargeMisId_ElMu +".tex"; 	
    selecTableChargeMisId_ElMu.Write(selectiontableChargeMisId_ElMu.c_str(),false, true, false, false, false, false, false);

    selecTableChargeMisId_2El.TableCalculator(true, true, true, true, true, false, true, true);//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST, bool mergeVV, bool mergettV, bool NP_mass)
    string selectiontableChargeMisId_2El = Treespath+"InclFourthGenSearch_SelectionTable_ChargeMisId2El"+postfix;
    selectiontableChargeMisId_2El = selectiontableChargeMisId_2El +".tex"; 	
		if(semiElectron) selecTableChargeMisId_2El.Write(selectiontableChargeMisId_2El.c_str(),false, true, false, false, false, false, false);


    cout << " - Closing the output file now..." << endl;
    fout->Close();
 
 
  delete fout;

  ofstream myfile1;
	if(option=="ChargeMisId" && systematic=="Nominal")
	{
		string myRockingFile1 = "ChargeMisId_OSEvents"+channelpostfix+".txt";
		myfile1.open(myRockingFile1.c_str());
		cout << endl;

		if(semiMuon) myfile1 << "THIS IS FOR THE MUON TRIGGER PART OF THE DATA"  << "\n"; 
		else  myfile1 << "THIS IS FOR THE ELECTRON TRIGGER PART OF THE DATA" << "\n";
		myfile1 << "# OS el+mu events: " << "\n";
		myfile1 << " barrel: " << NbOSElMu_EB_data << " +- " << sqrt(NbOSElMu_EB_data) << "\n"; 	
		myfile1 << " endcap: " << NbOSElMu_EE_data << " +- " << sqrt(NbOSElMu_EE_data) << "\n"; 	 	
		myfile1 << "\n";
		
		myfile1 << "# predicted SS el+mu events: " << "\n";
		
		float Unc_EB_events = sqrt(pow(CM_EB*sqrt(NbOSElMu_EB_data),2)+pow(NbOSElMu_EB_data*CM_EB_UNC,2));
		myfile1 << " 1 barrel el: " << NbOSElMu_EB_data*CM_EB  << " +- " << Unc_EB_events << "\n"; 	
		float Unc_EE_events = sqrt(pow(CM_EE*sqrt(NbOSElMu_EE_data),2)+pow(NbOSElMu_EE_data*CM_EE_UNC,2));
		myfile1 << " 1 endcap el: " << NbOSElMu_EE_data*CM_EE << " +- " << Unc_EE_events << "\n"; 	
		myfile1 << "\n";
		
		if(!semiMuon)
		{
			myfile1 << "# OS el+el events: " << "\n";
			myfile1 << " 2 barrel: " << NbOSElectrons_EBEB_data << " +- " << sqrt(NbOSElectrons_EBEB_data) << "\n"; 	
			myfile1 << " 2 endcap: " << NbOSElectrons_EEEE_data << " +- " << sqrt(NbOSElectrons_EEEE_data) << "\n"; 	 	
			myfile1 << " barrel+endcap: " << NbOSElectrons_EBEE_data << " +- " << sqrt(NbOSElectrons_EBEE_data) << "\n"; 	 	
			myfile1 << "\n";
			
			myfile1 << "# predicted SS el+el events: " << "\n";
			
			float Unc_EBEB_events = sqrt(pow(2*CM_EB*sqrt(NbOSElectrons_EBEB_data),2)+ pow(NbOSElectrons_EBEB_data*2*CM_EB_UNC,2));
			myfile1 << " 2 barrel: " << NbOSElectrons_EBEB_data*2*CM_EB << " +- " << Unc_EBEB_events << "\n"; 	
			
			float Unc_EEEE_events = sqrt(pow(2*CM_EE*sqrt(NbOSElectrons_EEEE_data),2)+ pow(NbOSElectrons_EEEE_data*2*CM_EE_UNC,2));
			myfile1 << " 2 endcap: " << NbOSElectrons_EEEE_data*2*CM_EE << " +- " << Unc_EEEE_events << "\n"; 	 	
			
			float Unc_EBEE_events = sqrt(pow((CM_EB+CM_EE)*sqrt(NbOSElectrons_EBEE_data),2)+ pow(NbOSElectrons_EBEE_data*CM_EB_UNC,2)+pow(NbOSElectrons_EBEE_data*CM_EE_UNC,2));
			myfile1 << " barrel+endcap: " << NbOSElectrons_EBEE_data*(CM_EB+CM_EE) << " +- " << Unc_EBEE_events << "\n"; 	 	
			
			myfile1 << "\n";
			myfile1 << " total: " <<  NbOSElectrons_EBEB_data*2*CM_EB+NbOSElectrons_EEEE_data*2*CM_EE+NbOSElectrons_EBEE_data*(CM_EB+CM_EE) << " +- " << sqrt(pow(Unc_EBEB_events,2)+pow(Unc_EEEE_events,2)+pow(Unc_EBEE_events,2)) << "\n"; 	 	
			myfile1 << "\n";
			myfile1.close();
		}
	}
	if(option=="FakeLepton" && systematic=="Nominal")
	{
		
		string myRockingFile1 = "FakeLepton_Events"+channelpostfix+".txt";
		myfile1.open(myRockingFile1.c_str());
		myfile1 << "\n";

		if(semiMuon)
		{
			myfile1 << "THIS IS FOR THE MUON TRIGGER PART OF THE DATA" << "\n";
			myfile1 << "# mu+mu events with fake muon: " << NbSSLooseMuonTightMuon_data << " +- " << sqrt(NbSSLooseMuonTightMuon_data)  << "\n";
			myfile1 << "# el+mu events with fake electron: " << NbSSLooseElectronTightMuon_data << " +- " << sqrt(NbSSLooseElectronTightMuon_data) << "\n";
			myfile1 << "\n";
			
			float Unc_mm = pow(Eff_TL_m*(1-Eff_TL_m)*sqrt(NbSSLooseMuonTightMuon_data),2)+pow(NbSSLooseMuonTightMuon_data*(1-2*Eff_TL_m)*Eff_TL_m_unc,2);
			myfile1 << "# predicted SS mu+muL events: " << NbSSLooseMuonTightMuon_data*Eff_TL_m*(1-Eff_TL_m) << " +- " << sqrt(Unc_mm) << "\n";
			float Unc_me = pow(Eff_TL_e*(1-Eff_TL_e)*sqrt(NbSSLooseElectronTightMuon_data),2)+pow(NbSSLooseElectronTightMuon_data*(1-2*Eff_TL_e)*Eff_TL_e_unc,2);
			myfile1 << "# predicted SS mu+eL events: " << NbSSLooseElectronTightMuon_data*Eff_TL_e*(1-Eff_TL_e) << " +- " << sqrt(Unc_me) << "\n";
		}

		if(!semiMuon)
		{
			myfile1 << "THIS IS FOR THE ELECTRON TRIGGER PART OF THE DATA" << "\n";
			myfile1 << "# el+el events with fake electron: " << NbSSLooseElectronTightElectron_data << " +- " << sqrt(NbSSLooseElectronTightElectron_data)<< "\n";
			myfile1 << "# el+mu events with fake muon: " << NbSSLooseMuonTightElectron_data << " +- " << sqrt(NbSSLooseMuonTightElectron_data) << "\n";
			
			float Unc_ee = pow(Eff_TL_e*(1-Eff_TL_e)*sqrt(NbSSLooseElectronTightElectron_data),2)+pow(NbSSLooseElectronTightElectron_data*(1-2*Eff_TL_e)*Eff_TL_e_unc,2);
			myfile1 << "# predicted SS e+eL events: " << NbSSLooseElectronTightElectron_data*Eff_TL_e*(1-Eff_TL_e) << " +- " << sqrt(Unc_ee) << "\n";
			float Unc_em = pow(Eff_TL_m*(1-Eff_TL_m)*sqrt(NbSSLooseMuonTightElectron_data),2)+pow(NbSSLooseMuonTightElectron_data*(1-2*Eff_TL_m)*Eff_TL_m_unc,2);
			myfile1 << "# predicted SS e+muL events: " << NbSSLooseMuonTightElectron_data*Eff_TL_m*(1-Eff_TL_m) << " +- " << sqrt(Unc_em) << "\n";
		}	
		myfile1 << "\n";
		myfile1.close();
	}

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
     for(unsigned int k=0; k<init_jets.size(); k++) //init_jets.size()
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

