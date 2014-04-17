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
#include "../StopSearchesBG/interface/StopPair_Evt.h"

#include "Style.C"

using namespace std;
using namespace TopTree;

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

  int doPUShift = 0; //0: off (except nominal PU reweighting) 1: minus 2: plus
  cout << "doPUShift: " << doPUShift << endl;

  clock_t start = clock();

  cout << "*************************************************************" << endl;
  cout << " Beginning of the program for stop pair searches ! " << endl;
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
  if (doPUShift == 1)
    postfix= postfix+"_PUMinus";
  if (doPUShift == 2)
    postfix= postfix+"_PUPlus";
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////// Configuration ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string channelpostfix = "";
  string xmlFileName = "";

  bool semiElectron = false; // use semiElectron channel?
  bool semiMuon = true; // use semiMuon channel?
  if(semiElectron && semiMuon){
	cout << "  --> Using both semiMuon and semiElectron channel? Choose only one (for the moment, since this requires running on different samples/skims)!" << endl;
	exit(1);
  }

  if(semiMuon){
	cout << " --> Using the semiMuon channel..." << endl;
	channelpostfix = "_diMu";
	xmlFileName = "../config/myStopSearchConfig_Muon.xml";
  }
  else if(semiElectron){
	cout << " --> Using the semiElectron channel..." << endl;
	channelpostfix = "_diEl";
	xmlFileName = "../config/myStopSearchConfig_Electron.xml";
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

  float Luminosity = 5000; /*Comment ;  Will use this value if MC only. */
  
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
  string rootFileName ("StopSearchesBG"+postfix+channelpostfix+".root");
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
  StopPair_Evt* MyStopPair_EvtCand = 0;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////// Histograms /////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////
  ////////////////// MultiSample plots  //////////////////////////////
  ////////////////////////////////////////////////////////////////////

  MSPlot["NbOfVertices"]               = new MultiSamplePlot(datasets, "NbOfVertices", 20, 0, 20, "Nb. of vertices");

  MSPlot["1stLeadingMuonRelIsolation"] = new MultiSamplePlot(datasets, "1stLeadingMuonRelIsolation", 500, 0, 0.5, "RelIso");
  MSPlot["2ndLeadingMuonRelIsolation"] = new MultiSamplePlot(datasets, "2ndLeadingMuonRelIsolation", 500, 0, 0.5, "RelIso");

  MSPlot["NbOfIsolatedMuons"]          = new MultiSamplePlot(datasets, "NbOfIsolatedMuons", 5, 0, 5, "Nb. of isolated muons");
  MSPlot["NbOfIsolatedElectrons"]      = new MultiSamplePlot(datasets, "NbOfIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");

  MSPlot["NbOfLooseIsolatedMuons"]     = new MultiSamplePlot(datasets, "NbOfLooseIsolatedMuons", 5, 0, 5, "Nb. of isolated muons");
  MSPlot["NbOfLooseIsolatedElectrons"] = new MultiSamplePlot(datasets, "NbOfLooseIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");

  MSPlot["NbOfSelectedJets"]           = new MultiSamplePlot(datasets, "NbOfSelectedJets", 15, 0, 15, "Nb. of jets");

  MSPlot["BdiscBJetCand_CVSMVA"]       = new MultiSamplePlot(datasets, "BdiscBJetCand_CVSMVA", 100, 0, 1, "CSV(MVA) b-disc.");
  MSPlot["BdiscBJetCand_TCHE"]         = new MultiSamplePlot(datasets, "BdiscBJetCand_TCHE", 100, 0, 50, "TCHE b-disc.");

  MSPlot["MET"]                        = new MultiSamplePlot(datasets, "MET",  100, 0, 200, "MET");

  /*Comment ; Add your own histos here. */

  ////////////////////////////////////////////////////////////////////
  ////////////////// 1D histograms  //////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  histo1D["lumiWeights"] = new TH1F("lumiWeights","lumiWeights;lumiWeight;#events",100,0,4);

  cout << " - Declared histograms ..." <<  endl;
	
  ////////////////////////////////////////////////////////////////////
  ////////////////// Plots  //////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  string pathPNG = "StopSearchesBG"+postfix+channelpostfix;
  pathPNG += "_MSPlots/"; 	
//  pathPNG = pathPNG +"/"; 	
  mkdir(pathPNG.c_str(),0777);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////// Selection Tables ///////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////
  ///////////////////// Channel : µ  /////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  vector<string> CutsSelecTableSemiMu;
  CutsSelecTableSemiMu.push_back(string("initial"));
  CutsSelecTableSemiMu.push_back(string("PU reweighting"));
  CutsSelecTableSemiMu.push_back(string("Trigger"));
  CutsSelecTableSemiMu.push_back(string("Good PV"));
  CutsSelecTableSemiMu.push_back(string("Exact. 1 isolated muon"));
  CutsSelecTableSemiMu.push_back(string("Veto on 2nd iso. lept."));
  CutsSelecTableSemiMu.push_back(string("$\\geq$ 1 jet"));
  CutsSelecTableSemiMu.push_back(string("$\\geq$ 2 jet"));
  CutsSelecTableSemiMu.push_back(string("$\\geq$ 3 jet"));
  CutsSelecTableSemiMu.push_back(string("$\\geq$ 4 jet"));
//  CutsSelecTableSemiMu.push_back(string("$MET \\leq 30$ GeV"));  
/*Comment : Add new cuts here */

  SelectionTable selecTableSemiMu(CutsSelecTableSemiMu, datasets);
  selecTableSemiMu.SetLuminosity(Luminosity);
  selecTableSemiMu.SetPrecision(1); // Nb of digits after the comma

  ////////////////////////////////////////////////////////////////////
  ///////////////////// Channel : e //////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  vector<string> CutsSelecTableSemiElec;
  CutsSelecTableSemiElec.push_back(string("initial"));
  CutsSelecTableSemiElec.push_back(string("PU reweighting"));
  CutsSelecTableSemiElec.push_back(string("Trigger"));
  CutsSelecTableSemiElec.push_back(string("Good PV"));
  CutsSelecTableSemiElec.push_back(string("Exact. 1 isolated electron"));
  CutsSelecTableSemiElec.push_back(string("Veto on 2nd iso. lept."));
  CutsSelecTableSemiElec.push_back(string("$\\geq$ 1 jet"));
  CutsSelecTableSemiElec.push_back(string("$\\geq$ 2 jet"));
  CutsSelecTableSemiElec.push_back(string("$\\geq$ 3 jet"));

  SelectionTable selecTableSemiElec(CutsSelecTableSemiElec, datasets);
  selecTableSemiElec.SetLuminosity(Luminosity);
  selecTableSemiElec.SetPrecision(1); // Nb of digits after the comma


  cout << " - SelectionTable instantiated ..." << endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////// PileUp Reweighting - 3D ///////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  Lumi3DReWeighting Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Flat10PlusTail.root","PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root", "pileup", "pileup");
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
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////// Initialize JEC factors //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
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
      
    int start = 0;
    unsigned int end = datasets[d]->NofEvtsToRunOver();
     
    if (verbose > 1) cout << " - Loop over events " << endl;      
    
    for (unsigned int ievt = start; ievt < end; ievt++)
    {        

	if(ievt%1000 == 0)
		std::cout<<"Processing the "<<ievt<<"th event ("<<100*(ievt-start)/(end-start)<<"%)"<<flush<<"\r";

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
		if(semiMuon)
		{
			if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
			{
				/*------------------------------------------------------------------
				Dataset : XXXX
				--------------------------------------------------------------------
				Trigger HLT_XXXXXXXXX_v1 available for runs 160431-163261
				Trigger HLT_XXXXXXXxX_v2 available for runs 163270-163869
				------------------------------------------------------------------*/
				if      (event->runId() >= 160431 && event->runId() <= 163261)
					itrigger = treeLoader.iTrigger (string ("HLT_XXXXXXXXX_v1"), currentRun, iFile);
  				else if (event->runId() >= 163270 && event->runId() <= 163869)
    					itrigger = treeLoader.iTrigger (string ("HLT_XXXXXXXXX_v2"), currentRun, iFile);
				/*------------------------------------------------------------------
				Dataset : XXXX
				--------------------------------------------------------------------
				Trigger HLT_XXXXXXXxX_v3 available for runs 165088-167043
				Trigger HLT_XXXXXXXxX_v4 available for runs 166346-166346
				Trigger HLT_XXXXXXXxX_v5 available for runs 167078-167913
				------------------------------------------------------------------*/
  				else if (event->runId() >= 165088 && event->runId() <= 167043 && event->runId() != 166346)
    					itrigger = treeLoader.iTrigger (string ("HLT_XXXXXXXxX_v3"), currentRun, iFile);
  				else if (event->runId() == 166346)
    					itrigger = treeLoader.iTrigger (string ("HLT_XXXXXXXxX_v4"), currentRun, iFile);
  				else if (event->runId() >= 167078 && event->runId() <= 167913)
    					itrigger = treeLoader.iTrigger (string ("HLT_XXXXXXXxX_v5"), currentRun, iFile);
				/*------------------------------------------------------------------
				Dataset : XXXXXXXxX
				--------------------------------------------------------------------
				Trigger HLT_XXXXXXXxX_v7 available for runs 170826-172619
				------------------------------------------------------------------*/
				else if (event->runId() >= 170249 && event->runId() <= 172619) //Aug05ReReco equivalent to PromptReco_v5 (about 5/pb lost)
				  	itrigger = treeLoader.iTrigger (string ("HLT_XXXXXXXxX_v7"), currentRun, iFile);
				/*------------------------------------------------------------------
				Dataset : XXXXXXXxX
				--------------------------------------------------------------------
				Trigger HLT_XXXXXXXxX_v7 available for runs 172620-173198
				Trigger HLT_XXXXXXXxX_v8 available for runs 173236-173692
				------------------------------------------------------------------*/
				else if (event->runId() >= 172620 && event->runId() <= 173198)
            				itrigger = treeLoader.iTrigger (string ("HLT_XXXXXXXxX_v7"), currentRun, iFile);
				else if (event->runId() >= 173236 && event->runId() <= 173692)
				  	itrigger = treeLoader.iTrigger (string ("HLT_XXXXXXXxX_v8"), currentRun, iFile);
				/*------------------------------------------------------------------
				Dataset : XXXXXXXxX
				--------------------------------------------------------------------
				Trigger HLT_XXXXXXXxX_v8  available for runs 175860-178380
				Trigger HLT_XXXXXXXxX_v11 available for runs 178420-179889
				Trigger HLT_XXXXXXXxX_v12 available for runs 179959-180252
				------------------------------------------------------------------*/
   				else if( event->runId() >=  175860 && event->runId() <= 178380 )
   				  	itrigger = treeLoader.iTrigger (string ("HLT_XXXXXXXxX_v8"), currentRun, iFile);
   				else if( event->runId() >=  178420 && event->runId() <= 179889 )
					itrigger = treeLoader.iTrigger (string ("HLT_XXXXXXXxX_v11"),currentRun, iFile);
				else if( event->runId() >=  179959 && event->runId() <=  180252 )
					itrigger = treeLoader.iTrigger (string ("HLT_XXXXXXXxX_v12"),currentRun, iFile); 
									   
  				if(itrigger == 9999)
				{
    				  cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (DATA) IN RUN " << event->runId() << endl;
    				  exit(1);
 	 			}
				//trigged = treeLoader.EventTrigged (itrigger);			

	   		}
	   		else 
	   		{
				itrigger = treeLoader.iTrigger (string ("HLT_XXXXXXXxX_v1"), currentRun, iFile);
    
  				if(itrigger == 9999)
				{
    			  		cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (" << dataSetName << ") IN RUN " << event->runId() << endl;
    			  		exit(1);
				}
				cout<<"Trigger bit nr : "<<itrigger<<endl;
				//trigged = treeLoader.EventTrigged (itrigger);
				//cout<<"Triggered ? : "<<trigged<<endl;
			}
		} //end if semiMuon
		else if(semiElectron)
		{
			cerr << "To BE IMPLEMENTED" << endl;
    			exit(1);
		} //end if semiElectron
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
	///////////////////////////////////////// Type I MET corrections: (Only for |eta| <=4.7 ) /////////////////////////////////////////////////
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
	
		//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JER correction:");
	
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
			
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////// Event selection ////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	selecTableSemiMu.Fill(d,0, 1.);//datasets[d]->Xsection() * datasets[d]->EquivalentLumi() );
	selecTableSemiElec.Fill(d,0, 1.);
	//selecTableDiEl.Fill(d,0, datasets[d]->Xsection() * datasets[d]->EquivalentLumi() );

	selecTableSemiMu.Fill(d,1,scaleFactor);
	selecTableSemiElec.Fill(d,1, scaleFactor);
	//selecTableSemiEl.Fill(d,1,scaleFactor);
		
	// Apply trigger selection
	trigged = treeLoader.EventTrigged (itrigger);
	if(!trigged)		   continue;
	selecTableSemiMu.Fill(d,2,scaleFactor);
	selecTableSemiElec.Fill(d,2, scaleFactor);

	// Declare selection instance    
	Selection selection(init_jets, init_muons, init_electrons, mets); // Cf. Selection class for more infos.
            
	// Define object selection cuts
	selection.setJetCuts(20.,2.4,0.01,1.,0.98,0.3,0.1);

	selection.setElectronCuts(20,2.5,0.125,0.02,1,0.3); // Et,Eta,RelIso,d0,DistVzPVz,DRJets. Apply Elec. ID VBTFWP80. Need to be updated
	selection.setLooseElectronCuts(15,2.5,0.2);

        selection.setMuonCuts(20.,2.4,0.125,10,0.02,0.3,1,1,1); // Pt, Eta, RelIso, NValidHits, d0, DRJets, NMatches, DistVzPVz, NPixelLayersWithMeas
	selection.setLooseMuonCuts(10,2.4,0.2);
	  
	//Select objects 
	vector<TRootJet*>      selectedJets        = selection.GetSelectedJets(true); // ApplyJetId

	vector<TRootElectron*> selectedElectrons   = selection.GetSelectedElectrons(vertex[0]);// Apply electron-vertex association

	vector<TRootMuon*>     selectedMuons_NoIso = selection.GetSelectedMuons(20,2.4,999.);
//	vector<TRootMuon*>     selectedMuons       = selection.GetSelectedMuons(vertex[0]);// Apply muon-vertex association
	vector<TRootMuon*>     selectedMuons       = selection.GetSelectedMuons(vertex[0],selectedJets); // Apply muon-vertex association and return only muons with DR(jets)>0.3

	vector<TRootMuon*>     looseMuons     = selection.GetSelectedLooseMuons();
	vector<TRootElectron*> vetoElectrons = selection.GetSelectedLooseElectrons(true); // VBTF Id
	//vector<TRootMCParticle*> mcParticles;

	MSPlot["NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);

	// Apply primary vertex selection
	bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2);
        if(!isGoodPV)   	   continue;
	selecTableSemiMu.Fill(d,3,scaleFactor);
	selecTableSemiElec.Fill(d,3, scaleFactor);

	MSPlot["1stLeadingMuonRelIsolation"]->Fill(selectedMuons_NoIso[0]->relativePfIso03(), datasets[d], true, Luminosity*scaleFactor);
	if(selectedMuons_NoIso.size()>1) MSPlot["2ndLeadingMuonRelIsolation"]->Fill(selectedMuons_NoIso[1]->relativePfIso03(), datasets[d], true, Luminosity*scaleFactor);

	MSPlot["NbOfIsolatedMuons"]->Fill(selectedMuons.size(), datasets[d], true, Luminosity*scaleFactor);
	MSPlot["NbOfIsolatedElectrons"]->Fill(selectedElectrons.size(), datasets[d], true, Luminosity*scaleFactor);

	// Select events with exactly one isolated lepton
	if(semiMuon){
		if(selectedMuons.size()!=1) continue;
		selecTableSemiMu.Fill(d,4,scaleFactor);
	}
	else if(semiElectron){
		if(selectedElectrons.size()!=1) continue;
		selecTableSemiElec.Fill(d,4,scaleFactor);
	}
	

	MSPlot["NbOfLooseIsolatedMuons"]    ->Fill(looseMuons.size(), datasets[d], true, Luminosity*scaleFactor);
	MSPlot["NbOfLooseIsolatedElectrons"]->Fill(looseMuons.size(), datasets[d], true, Luminosity*scaleFactor);

	// Veto events with another isolated lepton (loose isolation)
	if(semiMuon){
		if(looseMuons.size()!=1) continue; // Veto on another loose isolated muon
		selecTableSemiMu.Fill(d,5,scaleFactor);
	}
	else if(semiElectron){
	        if(vetoElectrons.size() != 0) continue; // Veto on a loose isolated electron
		selecTableSemiElec.Fill(d,5,scaleFactor);
	}

	MSPlot["NbOfSelectedJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);

	if(selectedJets.size()<1) continue; //at least 1 jet
	selecTableSemiMu.Fill(d,6,scaleFactor); 
	selecTableSemiElec.Fill(d,6, scaleFactor);
	if(selectedJets.size()<2) continue; //at least 2 jets
	selecTableSemiMu.Fill(d,7,scaleFactor);
	selecTableSemiElec.Fill(d,7, scaleFactor);
	if(selectedJets.size()<3) continue; //at least 3 jets
	selecTableSemiMu.Fill(d,8,scaleFactor);
	selecTableSemiElec.Fill(d,8, scaleFactor);
	if(selectedJets.size()<4) continue; //at least 4 jets
	selecTableSemiMu.Fill(d,9,scaleFactor);
	selecTableSemiElec.Fill(d,9, scaleFactor);

	MyStopPair_EvtCand = 0;
	// Do whatever you want with it...
	if(MyStopPair_EvtCand) delete MyStopPair_EvtCand;

    }//loop on events
    
    cout<<endl;
    
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

  //Selection tables (latex formatted output table containing the cut-flow)
  selecTableSemiMu.TableCalculator(false, true, true, true, true); //(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
  selecTableSemiElec.TableCalculator(false, true, true, true, true); //(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)

  // Options : WithError (false), writeMerged (true), useBookTabs (false), addRawNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false) (For more infos, cf. SelectableTable class)
if(semiMuon)
  selecTableSemiMu.Write("StopSearchesBG_EventSelectionTable_DiMu.tex",true, true,true,true,false,false,true);
else if(semiElectron)
  selecTableSemiElec.Write("StopSearchesBG_EventSelectionTable_DiElec.tex",true, true,true,true,false,false,true);

  fout->cd();
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
	MultiSamplePlot *temp = it->second;
	string name = it->first;
	// Options : bool addRandomPseudoData = false, string label = string("CMSPlot"), bool mergeTT = false, bool mergeQCD = false, bool mergeW = false, bool mergeZ = false, bool mergeST = false, int scaleNPSignal = 1, bool addRatio = false);
	temp->Draw(false, name, true, true, true, true, true,1,false); // merge TT/QCD/W/Z/ST/
	
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

