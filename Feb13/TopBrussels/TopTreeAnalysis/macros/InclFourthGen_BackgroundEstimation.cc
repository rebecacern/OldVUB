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

  //which systematic to run?
  string option = "ChargeMisid";
  if (argc >= 2)
		option = string(argv[1]);
  cout << "Option to be used: " << option << endl;
  if( ! (option == "ChargeMisId" || option == "MuonFake" || option == "ElectronFake") )
  {
    cout << "Unknown systematic!!!" << endl;
    cout << "Possible options are: " << endl;
    exit(-1);
  }


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

  //SetStyle if needed
  setTDRStyle();
  //setMyStyle();

  string postfix = "_16Apr2012_CHARGEMISID"; // to relabel the names of the output file  
	postfix= postfix+"_"+option;

  /////////////////////
  // Configuration
  /////////////////////
	bool verbosity = true;
	string channelpostfix = "";
  bool semiElectron = false; // use semiElectron channel?
  bool semiMuon = false; // use semiMuon channel?
  if(option=="MuonFake") semiMuon = true;
	else semiElectron = true;
	
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

  //Output ROOT file
  string rootFileName ("InclFourthGenSearch_BackgroundEstimation"+postfix+channelpostfix+".root");
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");
  
  //xml file
  string xmlFileName = "";
	if(semiElectron) xmlFileName = "../config/myFourthGenconfig_Electron_Fall11_BackgroundEstimationNEW.xml";
  else if(semiMuon) xmlFileName = "../config/myFourthGenconfig_Muon_Fall11_BackgroundEstimationNEW.xml";
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
  
 
  //vector of objects
  cout << " - Variable declaration ..." << endl;
  vector < TRootVertex* > vertex;
  vector < TRootMuon* > init_muons;
  vector < TRootElectron* > init_electrons;
  vector < TRootJet* > init_jets;
  vector < TRootMET* > mets;

  //Global variable
  TRootEvent* event = 0;
	
  histo1D["LeptonPt_loose"] = new TH1F("leptonpt loose","leptonpt;pt leptons;#events",250,0,500);
  histo1D["LeptonPt_tight"] = new TH1F("leptonpt tight","leptonpt;pt leptons;#events",250,0,500);
	histo1D["LeptonReliso_loose"] = new TH1F("Lepton reliso loose", "Lepton reliso; reliso; # events", 50, 0, 1);
	histo1D["LeptonReliso_tight"] = new TH1F("Lepton reliso tight", "Lepton reliso; reliso; # events", 50, 0, 1);

  MSPlot["MS_Zmass"] = new MultiSamplePlot(datasets,"Z mass", 150, 0, 150, "Z mass (GeV)");


	cout << " - Declared histograms ..." <<  endl;

  float NbSSevents = 0;   
  float NbTrievents = 0;  
	float Nb_Zpeak_EB_SS_MC = 0; float Nb_Zpeak_EE_SS_MC = 0; float Nb_Zpeak_EB_OS_MC = 0; float Nb_Zpeak_EE_OS_MC = 0;
	float Nb_Zpeak_EB_SS_data = 0; float Nb_Zpeak_EE_SS_data = 0; float Nb_Zpeak_EB_OS_data = 0; float Nb_Zpeak_EE_OS_data = 0;
	int NbOfLooseMuons = 0; int NbOfTightMuons = 0; int NbOfLooseElectrons = 0; int NbOfTightElectrons = 0; 
	

  ////////////////////////////////////
  /// Selection table
  ////////////////////////////////////
  vector<string> CutsSelecTableChargeMisId_2El;
  CutsSelecTableChargeMisId_2El.push_back(string("Zpeak: SS el EB"));
  CutsSelecTableChargeMisId_2El.push_back(string("Zpeak: SS el EE"));
  CutsSelecTableChargeMisId_2El.push_back(string("Zpeak: OS el EB"));
  CutsSelecTableChargeMisId_2El.push_back(string("Zpeak: OS el EE"));

  SelectionTable selecTableChargeMisId_2El(CutsSelecTableChargeMisId_2El, datasets);
  selecTableChargeMisId_2El.SetLuminosity(Luminosity);
  selecTableChargeMisId_2El.SetPrecision(2);

  vector<string> CutsSelecTableFakeLepton;
  CutsSelecTableFakeLepton.push_back(string("basic cuts"));
  CutsSelecTableFakeLepton.push_back(string("lepton pt > 35"));
  SelectionTable selecTableFakeLepton(CutsSelecTableFakeLepton, datasets);
	selecTableFakeLepton.SetLuminosity(Luminosity);
  selecTableFakeLepton.SetPrecision(1);


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
		


    ////////////////////////////////////
    ////////////////////////////////////
    ///////// Loop on events
    ////////////////////////////////////
    ////////////////////////////////////
    int itrigger = -1, previousRun = -1;
     
    if (verbose > 1)
      cout << " - Loop over events " << endl;      
    
    for (int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
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
						//Run2011A
						if (event->runId() >= 160431 && event->runId() <= 163261)//May10ReReco
							itrigger = treeLoader.iTrigger (string ("HLT_Mu8_v1"), currentRun, iFile);	
  					else if (event->runId() >= 163270 && event->runId() <= 163869)
    				  itrigger = treeLoader.iTrigger (string ("HLT_Mu8_v2"), currentRun, iFile);
  					else if (event->runId() >= 165088 && event->runId() <= 165633)//PromptReco_v4; splitted over 2 toptrees: 565 and 641
    					itrigger = treeLoader.iTrigger (string ("HLT_Mu8_v3"), currentRun, iFile);
  					else if (event->runId() >= 165970 && event->runId() <= 167043 && event->runId() != 166346)
    					itrigger = treeLoader.iTrigger (string ("HLT_Mu8_v3"), currentRun, iFile);
  					else if (event->runId() == 166346)
    				  continue;
  					else if (event->runId() >= 167078 && event->runId() <= 167913)
    				  itrigger = treeLoader.iTrigger (string ("HLT_Mu8_v5"), currentRun, iFile);
						else if (event->runId() >= 170249 && event->runId() <= 172619) //Aug05ReReco: 
				  		itrigger = treeLoader.iTrigger (string ("HLT_Mu8_v7"), currentRun, iFile);
						else if (event->runId() >= 172620 && event->runId() <= 173198) //first part of PromptReco_v6, same as previous trigger
            	itrigger = treeLoader.iTrigger (string ("HLT_Mu8_v7"), currentRun, iFile);
						else if (event->runId() >= 173236 && event->runId() <= 173692) //second part of PromptReco_v6
				  		itrigger = treeLoader.iTrigger (string ("HLT_Mu8_v8"), currentRun, iFile);
						
						//no Run2011B	   
  					
						if(itrigger == 9999)
						{
    				  cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (DATA) IN RUN " << currentRun << endl;
    				  exit(1);
  					}
	   			}			
				} //end if semiMuon
	 			else if(semiElectron)
				{
	  			if(option == "ElectronFake"){
						if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
	   				{      		
							// /SingleElectron/Run2011A-May10ReReco-v1/AOD 	
							if (event->runId() >= 160404 && event->runId() < 161217)
    						itrigger = treeLoader.iTrigger (string ("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v1"), currentRun, iFile);
							else if (event->runId() >= 161217 && event->runId() < 163270)
    						itrigger = treeLoader.iTrigger (string ("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2"), currentRun, iFile);
     					else if (event->runId() >= 163270 && event->runId() <= 163869)
    						itrigger = treeLoader.iTrigger (string ("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v3"), currentRun, iFile);				 
							// /ElectronHad/Run2011A-PromptReco-v4/AOD
							else if (event->runId() >= 165088 && event->runId() < 165970)
    						itrigger = treeLoader.iTrigger (string ("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v4"), currentRun, iFile);
							else if (event->runId() >= 165970 && event->runId() < 167038)
    						itrigger = treeLoader.iTrigger (string ("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v5"), currentRun, iFile);
							else if (event->runId() >= 167038 && event->runId() <= 167913)
    						itrigger = treeLoader.iTrigger (string ("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v6"), currentRun, iFile);			  
							// /ElectronHad/Run2011A-05Aug2011-v1/AOD
							else if (event->runId() >= 170249 && event->runId() <= 172619)
    						itrigger = treeLoader.iTrigger (string ("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v7"), currentRun, iFile);  
							// /ElectronHad/Run2011A-PromptReco-v6/AOD 
							else if (event->runId() >= 172620 && event->runId() <= 173198)
    						itrigger = treeLoader.iTrigger (string ("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v7"), currentRun, iFile);  
							else if (event->runId() >= 173199 && event->runId() <= 178380)
    						itrigger = treeLoader.iTrigger (string ("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v8"), currentRun, iFile);  				   				   	
							else if(currentRun >= 178381 && currentRun <= 179889)
              	itrigger = treeLoader.iTrigger (string ("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v11"), currentRun, iFile);
            	else if(currentRun >= 179959 && currentRun <= 999999)
              	itrigger = treeLoader.iTrigger (string ("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v12"), currentRun, iFile);
							
							
  						if(itrigger == 9999)
							{
    						cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (DATA) IN RUN " << event->runId() << endl;
    						exit(1);
  						}// semi-electron
 	   				}
					}
					else if(option == "ChargeMisId")
					{
						if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
	   				{      		
						// DoubleElectron triggers for chargeMisId
 							if(currentRun >= 150000 && currentRun <= 161176){
            	  itrigger = treeLoader.iTrigger (string ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1"), currentRun, iFile);
            	}else if(currentRun >= 161179 && currentRun <= 163261){
            	  itrigger = treeLoader.iTrigger (string ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2"), currentRun, iFile);
            	}else if(currentRun >= 163262 && currentRun <= 164237){
            	  itrigger = treeLoader.iTrigger (string ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3"), currentRun, iFile);
            	}else if(currentRun >= 165085 && currentRun <= 165888){
            	  itrigger = treeLoader.iTrigger (string ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4"), currentRun, iFile);
            	}else if(currentRun >= 165900 && currentRun <= 167037){
            	  itrigger = treeLoader.iTrigger (string ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5"), currentRun, iFile);
            	}else if(currentRun >= 167038 && currentRun <= 170053){
            	  itrigger = treeLoader.iTrigger (string ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6"), currentRun, iFile);
            	}else if(currentRun >= 170054 && currentRun <= 170759){
             	 itrigger = treeLoader.iTrigger (string ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6"), currentRun, iFile);
            	}else if(currentRun >= 170760 && currentRun <= 173198){
								itrigger = treeLoader.iTrigger (string ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7"), currentRun, iFile);
            	}else if(currentRun >= 173199 && currentRun <= 178380){
            		itrigger = treeLoader.iTrigger (string ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8"), currentRun, iFile);
            	}else if(currentRun >= 178381 && currentRun <= 179958){
              	itrigger = treeLoader.iTrigger (string ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9"), currentRun, iFile);
            	}else if(currentRun >= 179959 && currentRun <= 999999){
              	itrigger = treeLoader.iTrigger (string ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10"), currentRun, iFile);
							}
								
							if(itrigger == 9999)
							{
    						cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (DATA) IN RUN " << currentRun << endl;
    						exit(1);
  						}// semi-electron
						} //no trigger on MC
					}
				} //end if semiElectron	
			} //end previousRun != currentRun
		  
      if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
      {	
	
      	//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JER correction:");
				jetTools->correctJetJER(init_jets, genjets, mets[0], "nominal",false);
				//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JER correction:");	       
		
      }

			double lumiWeight3D = 1.0;
			if(!(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA"))
			{				
      	////////////////////////////
      	// apply PU Reweighting
      	////////////////////////////
				lumiWeight3D = Lumi3DWeights.weight3D(event->nPu(-1),event->nPu(0),event->nPu(+1));
	 			scaleFactor = scaleFactor*lumiWeight3D;
      	//histo1D["lumiWeights"]->Fill(scaleFactor);	
			}
			cout << "scaleFactor " << scaleFactor << endl;
			
					
      /////////////////////////////
      // Selection
      /////////////////////////////
      bool eventSelected = false;

      vector<TRootMCParticle*> mcParticles;
      //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before treeLoader.LoadMCEvent:");      
      if(dataSetName.find("TTbarJets_SemiMu") == 0 || dataSetName.find("TTbarJets_SemiElectron") == 0 || dataSetName.find("NP_Tprime")==0 || dataSetName.find("NP_overlay_Tprime")==0)
      {
        treeLoader.LoadMCEvent(ievt, 0, 0, mcParticles,false);  
        sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      }

			//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After treeLoader.LoadMCEvent:");
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
      trigged = false;
			if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
				trigged = treeLoader.EventTrigged (itrigger);		
			else trigged = true; // for MC (no jet trigger available)

     //Declare selection instances    
      Selection selectionFakeMuon(init_jets, init_muons, init_electrons, mets); //mets can also be corrected...
      Selection selectionFakeElectron(init_jets, init_muons, init_electrons, mets); //mets can also be corrected...
      Selection selectionChargeMisId(init_jets, init_muons, init_electrons, mets); //mets can also be corrected...

			//need at least one loose muon with pT>10 (trigger = single muon with pt 8)
      selectionFakeMuon.setJetCuts(20.,2.4,0.01,1.,0.98,0.3,0.1);
      selectionFakeMuon.setMuonCuts(20,2.1,0.125,10,0.02,0.3,1,1,1);
      selectionFakeMuon.setLooseMuonCuts(10,2.5,9999.);
      selectionFakeMuon.setLooseElectronCuts(15,2.5,9999.);	 				

			//need at least one loosely isolated electron with pT>10 and a jet with pt>40 (trigger = loose caloiso single electron with pt 8 + jet with pt 40)
      selectionFakeElectron.setJetCuts(20.,2.4,0.01,1.,0.98,0.3,0.1);
      selectionFakeElectron.setLooseMuonCuts(10,2.5,9999.);
      selectionFakeElectron.setElectronCuts(20,2.5,0.1,0.02,1,0.3);
      selectionFakeElectron.setLooseElectronCuts(15,2.5,9999.);	 				

			//need at least 2 tight isolated electrons with pT>20 (trigger = 2 isolated electrons with pt 17 and pt 8)
      selectionChargeMisId.setJetCuts(30.,2.4,0.01,1.,0.98,0.3,0.1);
      selectionChargeMisId.setElectronCuts(20,2.5,0.1,0.02,1,0.3);
      selectionChargeMisId.setLooseElectronCuts(15,2.5,0.2);	 				
      			      
 
      if(trigged && semiMuon)
      { 
      	//if(verbosity) cout <<  " event is triggered, will do fake muon probability estimation"  << endl;
      	vector<TRootJet*> selectedJets_FM;
      	vector<TRootMuon*> selectedMuons_FM;
      	vector<TRootMuon*> selectedFakeMuons_FM;
      	vector<TRootElectron*> selectedFakeElectrons_FM;
      	
				if (init_jets.size() > 0)
      	{
	    		selectedJets_FM = selectionFakeMuon.GetSelectedJets(true);				
	    		selectedMuons_FM = selectionFakeMuon.GetSelectedMuons(vertex[0],selectedJets_FM);
      	}
      	selectedFakeMuons_FM = selectionFakeMuon.GetSelectedLooseMuons();
       	selectedFakeElectrons_FM = selectionFakeMuon.GetSelectedLooseElectrons(false);

	      isGoodPV = selectionFakeMuon.isPVSelected(vertex, 4, 24., 2); //in the past this was put in the config, but this is not very useful, since the PV cuts are quite standard
					
        if(isGoodPV)
				{
					if(selectedFakeMuons_FM.size()>=1 && mets[0]->Et()<20)
					{
						sort(selectedJets_FM.begin(),selectedJets_FM.end(),HighestPt()); // HighestPt() is included from the Selection class
						
						float MT = sqrt(2*selectedFakeMuons_FM[0]->Pt()*mets[0]->Pt()*(1-cos(selectedFakeMuons_FM[0]->DeltaPhi(*mets[0]))));
						if(selectedJets_FM.size()>=(unsigned int)anaEnv.NofJets && MT<25)
						{  //at least 1 jet!
//							if(dataSetName == "Data")
//							{
  							bool foundZ = false;
  							if(selectedFakeMuons_FM.size()>=1){
									for(unsigned int j=0;j<selectedFakeMuons_FM.size();j++)
  								{
  									for(unsigned int i=0;i<selectedFakeMuons_FM.size();i++)
								 		{
   										TRootMuon* mu1 = (TRootMuon*) selectedFakeMuons_FM[j];
   										TRootMuon* mu2 = (TRootMuon*) selectedFakeMuons_FM[i];
    									if( fabs(mu2->Pt() - mu1->Pt()) > 0.001 && fabs(mu2->Eta() - mu1->Eta()) > 0.001 && mu1->charge() != mu2->charge())
    									{	
      									if( (*mu1 + *mu2).M() >= (91.-20) && (*mu1 + *mu2).M() <= (91.+20) )
        									foundZ = true;
											}
								    }
								  }
								}
 								if(selectedFakeElectrons_FM.size()>=1)
								{
  								for(unsigned int j=0;j<selectedFakeElectrons_FM.size();j++)
  								{
  									for(unsigned int i=0;i<selectedFakeElectrons_FM.size();i++)
								 		{
   										TRootElectron* el1 = (TRootElectron*) selectedFakeElectrons_FM[j];
   										TRootElectron* el2 = (TRootElectron*) selectedFakeElectrons_FM[i];
    									if( fabs(el2->Pt() - el1->Pt()) > 0.001 && fabs(el2->Eta() - el1->Eta()) > 0.001 && el1->charge() != el2->charge())
    									{	
      									if( (*el1 + *el2).M() >= (91.-20) && (*el1 + *el2).M() <= (91.+20) )
        									foundZ = true;
											}
								    }
								  }
								}
								if(!foundZ) // reject events with a Z boson
								{											
									selecTableFakeLepton.Fill(d,0,scaleFactor);
									for(unsigned int j=0;j<selectedFakeMuons_FM.size();j++)
									{
										if(selectedFakeMuons_FM[j]->Pt()<35)	
											if(dataSetName == "Data") NbOfLooseMuons++;
									}
									for(unsigned int j=0;j<selectedMuons_FM.size();j++)
									{
										if(selectedMuons_FM[j]->Pt()<35) 
											if(dataSetName == "Data") NbOfTightMuons++;
									}
									if(selectedFakeMuons_FM.size() > 0 && selectedFakeMuons_FM[0]->Pt()>35) selecTableFakeLepton.Fill(d,1,scaleFactor);
								}
							//}
							////////////////////// FAKE LEPTON RATE ESTIMATION ///////////////////
						} 
          } 
        } // end good PV
      }// end trigged & semiMuon
      
			///// EVENTS TRIGGERED BY ELECTRON TRIGGER
			else if(trigged && semiElectron)
      {
      	
				if(option=="ChargeMisId")////////////////////// CHARGE MIS-ID PROBABILITY ESTIMATION //////////////////
				{
      		//if(verbosity) cout <<  " event is triggered, will do charge misid probability estimation for electrons"  << endl;

      		vector<TRootJet*> selectedJets_CM;
      		vector<TRootElectron*> selectedElectrons_CM;
      		vector<TRootElectron*> selectedLooseElectronsVBTFid_CM;
					
					if (init_jets.size() > 0)
      		{
	    			selectedJets_CM = selectionChargeMisId.GetSelectedJets(true);				
	    			selectedElectrons_CM = selectionChargeMisId.GetSelectedElectrons(vertex[0],selectedJets_CM);
      		}
      		selectedLooseElectronsVBTFid_CM = selectionChargeMisId.GetSelectedLooseElectrons(false); //loose vbtfid is required 

      		isGoodPV = selectionChargeMisId.isPVSelected(vertex, 4, 24., 2); //in the past this was put in the config, but this is not very useful, since the PV cuts are quite standard
        	if( isGoodPV )
        	{
          	if( selectedElectrons_CM.size() == 2 && selectedLooseElectronsVBTFid_CM.size() == selectedElectrons_CM.size() && mets[0]->Et()<20)
          	{
              if( selectionChargeMisId.passConversionRejection(selectedElectrons_CM[0])  && selectionChargeMisId.passConversionRejection(selectedElectrons_CM[1]))
              {
								sort(selectedJets_CM.begin(),selectedJets_CM.end(),HighestPt()); // HighestPt() is included from the Selection class

								float MT = sqrt(2*selectedElectrons_CM[0]->Pt()*mets[0]->Pt()*(1-cos(selectedElectrons_CM[0]->DeltaPhi(*mets[0]))));
								float MT2 = sqrt(2*selectedElectrons_CM[1]->Pt()*mets[0]->Pt()*(1-cos(selectedElectrons_CM[1]->DeltaPhi(*mets[0]))));
								if(MT<25 && MT2<25)
								{
									float Zmass = ((TLorentzVector) *selectedElectrons_CM[0]+ (TLorentzVector) *selectedElectrons_CM[1]).M();
									MSPlot["MS_Zmass"]->Fill(Zmass,datasets[d],true,Luminosity*scaleFactor);
									if(selectionChargeMisId.foundZCandidate(selectedElectrons_CM, selectedElectrons_CM, 10.))
									{								
										if(selectedElectrons_CM[0]->charge() == selectedElectrons_CM[1]->charge())
										{ 
											if(fabs(selectedElectrons_CM[0]->superClusterEta())<1.4442 && fabs(selectedElectrons_CM[1]->superClusterEta())<1.4442){
												selecTableChargeMisId_2El.Fill(d,0,scaleFactor);
												if(dataSetName.find("Data") == 0) Nb_Zpeak_EB_SS_data+=scaleFactor;
												else if(dataSetName.find("ZJets") == 0) Nb_Zpeak_EB_SS_MC+=1;
											}else if(fabs(selectedElectrons_CM[0]->superClusterEta())>1.5660 && fabs(selectedElectrons_CM[1]->superClusterEta())>1.5660){
												selecTableChargeMisId_2El.Fill(d,1,scaleFactor);
												if(dataSetName.find("Data") == 0) Nb_Zpeak_EE_SS_data+=scaleFactor;
												else if(dataSetName.find("ZJets") == 0) Nb_Zpeak_EE_SS_MC+=1;
											}
										}else{
											if(fabs(selectedElectrons_CM[0]->superClusterEta())<1.4442 && fabs(selectedElectrons_CM[1]->superClusterEta())<1.4442){
												selecTableChargeMisId_2El.Fill(d,2,scaleFactor);
												if(dataSetName.find("Data") == 0) Nb_Zpeak_EB_OS_data+=scaleFactor;
												else if(dataSetName.find("ZJets") == 0) Nb_Zpeak_EB_OS_MC+=1;
											}else if(fabs(selectedElectrons_CM[0]->superClusterEta())>1.5660 && fabs(selectedElectrons_CM[1]->superClusterEta())>1.5660){
												selecTableChargeMisId_2El.Fill(d,3,scaleFactor);
												if(dataSetName.find("Data") == 0) Nb_Zpeak_EE_OS_data+=scaleFactor;
												else if(dataSetName.find("ZJets") == 0) Nb_Zpeak_EE_OS_MC+=1;
											}															
										}
									}	
								}
							}
						}
					}
				}
				else if (option=="ElectronFake")////////////////////// FAKE LEPTON PROBABILITY ESTIMATION ///////////////////
				{
      			      
      		//if(verbosity) cout <<  " event is triggered, will do fake electron probability estimation"  << endl;

      		vector<TRootJet*> selectedJets_EF;
    		  vector<TRootElectron*> selectedElectrons_EF;
     			vector<TRootMuon*> selectedFakeMuons_EF;
		      vector<TRootElectron*> selectedFakeElectrons_EF;

					if (init_jets.size() > 0)
      		{
	    			selectedJets_EF = selectionFakeElectron.GetSelectedJets(true);				
	    			selectedElectrons_EF = selectionFakeElectron.GetSelectedElectrons(vertex[0],selectedJets_EF);
      		}
      		selectedFakeMuons_EF = selectionFakeElectron.GetSelectedLooseMuons();
      		selectedFakeElectrons_EF = selectionFakeElectron.GetSelectedLooseElectrons(false); 
    	 	
      		isGoodPV = selectionFakeElectron.isPVSelected(vertex, 4, 24., 2); //in the past this was put in the config, but this is not very useful, since the PV cuts are quite standard
        	
					if( isGoodPV )
        	{
          	if( selectedFakeElectrons_EF.size() >= 1 && mets[0]->Et()<20)
          	{
              if( selectionFakeElectron.passConversionRejection(selectedFakeElectrons_EF[0]) )
              {
								sort(selectedJets_EF.begin(),selectedJets_EF.end(),HighestPt()); // HighestPt() is included from the Selection class

								float MT = sqrt(2*selectedFakeElectrons_EF[0]->Pt()*mets[0]->Pt()*(1-cos(selectedFakeElectrons_EF[0]->DeltaPhi(*mets[0]))));
								if( selectedJets_EF.size()>=(unsigned int)anaEnv.NofJets && selectedJets_EF[0]->Pt()>50 && MT<25)
								{
  								bool foundZ = false;
  								if(selectedFakeMuons_EF.size()>=1){
										for(unsigned int j=0;j<selectedFakeMuons_EF.size();j++)
  									{
  										for(unsigned int i=0;i<selectedFakeMuons_EF.size();i++)
									 		{
   											TRootMuon* mu1 = (TRootMuon*) selectedFakeMuons_EF[j];
   											TRootMuon* mu2 = (TRootMuon*) selectedFakeMuons_EF[i];
    										if( fabs(mu2->Pt() - mu1->Pt()) > 0.001 && fabs(mu2->Eta() - mu1->Eta()) > 0.001 && mu1->charge() != mu2->charge())
    										{	
      										if( (*mu1 + *mu2).M() >= (91.-20) && (*mu1 + *mu2).M() <= (91.+20) )
        										foundZ = true;
												}
									    }
									  }
									}
 									if(selectedFakeElectrons_EF.size()>=1)
									{
  									for(unsigned int j=0;j<selectedFakeElectrons_EF.size();j++)
  									{
  										for(unsigned int i=0;i<selectedFakeElectrons_EF.size();i++)
									 		{
   											TRootElectron* el1 = (TRootElectron*) selectedFakeElectrons_EF[j];
   											TRootElectron* el2 = (TRootElectron*) selectedFakeElectrons_EF[i];
    										if( fabs(el2->Pt() - el1->Pt()) > 0.001 && fabs(el2->Eta() - el1->Eta()) > 0.001 && el1->charge() != el2->charge())
    										{	
      										if( (*el1 + *el2).M() >= (91.-20) && (*el1 + *el2).M() <= (91.+20) )
        										foundZ = true;
												}
									    }
									  }
									}
									
									if(!foundZ) // reject events with a Z boson
									{											
										selecTableFakeLepton.Fill(d,0,scaleFactor);
										for(unsigned int j=0;j<selectedFakeElectrons_EF.size();j++){
											if(selectedFakeElectrons_EF[j]->Pt()<35)	
												if(dataSetName == "Data") NbOfLooseElectrons++;      
										}
										histo1D["LeptonPt_loose"]->Fill(selectedFakeElectrons_EF[0]->Pt(),scaleFactor);	
										float relIso = (selectedFakeElectrons_EF[0]->chargedHadronIso()+selectedFakeElectrons_EF[0]->neutralHadronIso()+selectedFakeElectrons_EF[0]->photonIso())/selectedFakeElectrons_EF[0]->Pt();
										histo1D["LeptonReliso_loose"]->Fill(relIso,scaleFactor);	
										for(unsigned int j=0;j<selectedElectrons_EF.size();j++){
											if(selectedElectrons_EF[j]->Pt()<35)	
												if(dataSetName == "Data") NbOfTightElectrons++;
										}
										if(selectedFakeElectrons_EF.size()>0 && selectedFakeElectrons_EF[0]->Pt()>35) selecTableFakeLepton.Fill(d,1,scaleFactor);

										if(selectedElectrons_EF.size()>=1){
											histo1D["LeptonPt_tight"]->Fill(selectedElectrons_EF[0]->Pt(),scaleFactor);	
											relIso = (selectedElectrons_EF[0]->chargedHadronIso()+selectedElectrons_EF[0]->neutralHadronIso()+selectedElectrons_EF[0]->photonIso())/selectedElectrons_EF[0]->Pt();
											histo1D["LeptonReliso_tight"]->Fill(relIso,scaleFactor);	
										}
									}
								} // end 'at least one jet'
							} // end conversion rejection for leading electron
          	} // end if selectedElectrons.size()>=1
        	} // end good PV
				}
      } // end trigged & semiElectron
						
    
	
    }//loop on events
    
    cout<<endl;
		
    if(jetTools) delete jetTools;
		
		//important: free memory
    treeLoader.UnLoadDataset();
	
    
  } //loop on datasets

	if(semiElectron && option=="ChargeMisId"){
		float chargeMisId_Barrel_MC = (float)Nb_Zpeak_EB_SS_MC/(2*(float)Nb_Zpeak_EB_OS_MC);
		float chargeMisId_Endcap_MC = (float)Nb_Zpeak_EE_SS_MC/(2*(float)Nb_Zpeak_EE_OS_MC);
		float chargeMisId_Barrel_MC_unc = sqrt(pow((sqrt((float)Nb_Zpeak_EB_SS_MC)/(2*(float)Nb_Zpeak_EB_OS_MC)),2)+pow(((float)Nb_Zpeak_EB_SS_MC*sqrt((float)Nb_Zpeak_EB_OS_MC)/pow((float)Nb_Zpeak_EB_OS_MC,2)),2));	
		float chargeMisId_Endcap_MC_unc = sqrt(pow((sqrt((float)Nb_Zpeak_EE_SS_MC)/(2*(float)Nb_Zpeak_EE_OS_MC)),2)+pow(((float)Nb_Zpeak_EE_SS_MC*sqrt((float)Nb_Zpeak_EE_OS_MC)/pow((float)Nb_Zpeak_EE_OS_MC,2)),2));	
		float chargeMisId_Barrel_data = (float)Nb_Zpeak_EB_SS_data/(2*(float)Nb_Zpeak_EB_OS_data);
		float chargeMisId_Endcap_data = (float)Nb_Zpeak_EE_SS_data/(2*(float)Nb_Zpeak_EE_OS_data);
		float chargeMisId_Barrel_data_unc = sqrt(pow((sqrt((float)Nb_Zpeak_EB_SS_data)/(2*(float)Nb_Zpeak_EB_OS_data)),2)+pow(((float)Nb_Zpeak_EB_SS_data*sqrt((float)Nb_Zpeak_EB_OS_data)/pow((float)Nb_Zpeak_EB_OS_data,2)),2));	
		float chargeMisId_Endcap_data_unc = sqrt(pow((sqrt((float)Nb_Zpeak_EE_SS_data)/(2*(float)Nb_Zpeak_EE_OS_data)),2)+pow(((float)Nb_Zpeak_EE_SS_data*sqrt((float)Nb_Zpeak_EE_OS_data)/pow((float)Nb_Zpeak_EE_OS_data,2)),2));	
		
		cout << endl;
		cout << "charge mis-reconstruction probability in MC: " << endl;
		cout << " # SS events barrel: " <<  Nb_Zpeak_EB_SS_MC << "  and endcap: " <<  Nb_Zpeak_EE_SS_MC << endl;
		cout << " # OS events barrel: " <<  Nb_Zpeak_EB_OS_MC << "  and endcap: " <<  Nb_Zpeak_EE_OS_MC << endl;
		cout << " chargeMisId_Barrel_MC " << chargeMisId_Barrel_MC << " +- " << chargeMisId_Barrel_MC_unc  << endl;
		cout << " chargeMisId_Endcap_MC " << chargeMisId_Endcap_MC << " +- " << chargeMisId_Endcap_MC_unc  << endl;
		cout << endl;
		cout << "charge mis-reconstruction probability in DATA: " << endl;
		cout << " # SS events barrel: " << Nb_Zpeak_EB_SS_data  << "  and endcap: " << Nb_Zpeak_EE_SS_data  << endl;
		cout << " # OS events barrel: " << Nb_Zpeak_EB_OS_data  << "  and endcap: " << Nb_Zpeak_EE_OS_data  << endl;
		cout << " chargeMisId_Barrel_data " << chargeMisId_Barrel_data << " +- " << chargeMisId_Barrel_data_unc << endl;
		cout << " chargeMisId_Endcap_data " << chargeMisId_Endcap_data << " +- " << chargeMisId_Endcap_data_unc << endl;
	}
	
	float MuonFakeRate = (float)NbOfTightMuons/(float)(NbOfLooseMuons);
	float ElectronFakeRate = (float)NbOfTightElectrons/(float)(NbOfLooseElectrons);

	cout << endl;
	cout << "Number of tight muons " << NbOfTightMuons << " +- " << sqrt(NbOfTightMuons) << endl;
	cout << "Number of loose muons " << NbOfLooseMuons << " +- " << sqrt(NbOfLooseMuons) << endl;
	float MuonFakeRate_unc = sqrt(pow(sqrt(NbOfTightMuons)/NbOfLooseMuons,2)+pow((NbOfTightMuons*sqrt(NbOfLooseMuons))/pow((float)NbOfLooseMuons,2),2));
	cout << "Fake rate for muons " << MuonFakeRate <<  " +- " << MuonFakeRate_unc << endl;
	cout << endl;
	cout << "Number of tight electrons " << NbOfTightElectrons << " +- " << sqrt(NbOfTightElectrons) <<  endl;
	cout << "Number of loose electrons " << NbOfLooseElectrons << " +- " << sqrt(NbOfLooseElectrons) <<  endl;
	float ElectronFakeRate_unc = sqrt(pow(sqrt(NbOfTightElectrons)/NbOfLooseElectrons,2)+pow((NbOfTightElectrons*sqrt(NbOfLooseElectrons))/pow((float)NbOfLooseElectrons,2),2));
	cout << "Fake rate for electrons " << ElectronFakeRate <<  " +- " << ElectronFakeRate_unc << endl;
	cout << endl;

	
  ///////////////////
  // Writing
  //////////////////
  cout << " - Writing outputs to the files ..." << endl;	
  selecTableChargeMisId_2El.TableCalculator(true, true, true, true, true, false, true, true);//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST, bool mergeVV, bool mergettV, bool NP_mass)
  string selectiontableChargeMisId_2El = "InclFourthGenSearch_SelectionTable_ChargeMisId2El"+postfix;
  selectiontableChargeMisId_2El = selectiontableChargeMisId_2El +".tex"; 	
	if(semiElectron) selecTableChargeMisId_2El.Write(selectiontableChargeMisId_2El.c_str(),false, true, false, false, false, false, false);

  selecTableFakeLepton.TableCalculator(true, true, true, true, true, false, true, true);//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST, bool mergeVV, bool mergettV, bool NP_mass)
  string selectiontablefakelepton = "InclFourthGenSearch_SelectionTable_FakeLepton"+postfix;
  selectiontablefakelepton = selectiontablefakelepton +".tex"; 	
	selecTableFakeLepton.Write(selectiontablefakelepton.c_str(),false, true, false, false, false, false, false);
 
    fout->cd();
    //Write histograms: MSPlots
    //if(savePNG) mkdir((pathPNG+"MSPlot/").c_str(),0777);
    //cout << "mkdir " << (pathPNG+"MSPlot/").c_str()<< endl;
		cout << "Running over all MS plots" << endl;
    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
        MultiSamplePlot *temp = it->second;
        string name = it->first;
        temp->Draw(false, name, true, true, true, true, true,5,false, true, true);//(bool addRandomPseudoData, string label, bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST,int scaleNPsignal, bool addRatio, bool mergeVV, bool mergeTTV)
        temp->Write(fout, name, false, "MSPlot/");//bool savePNG
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
			tempCanvas->SaveAs( (it->first+".pdf").c_str() ); //well, is actually not png but pdf...
    }    
    cout << "1D plots written" << endl;

    fout->Close();
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

//https://twiki.cern.ch/twiki/bin/viewauth/CMS/SingleTopTurnOnCurves
float jetprob(float jetpt, float btagvalue){
	float prob=0.982*exp(-30.6*exp(-0.151*jetpt));
  prob*=0.844*exp((-6.72*exp(-0.720*btagvalue))); //"for the offline TCHP tagger"
	//prob*=0.736*exp((-8.01*exp(-0.540*btagvalue))); //"for the offline TCHE tagger"
	return prob;
};
