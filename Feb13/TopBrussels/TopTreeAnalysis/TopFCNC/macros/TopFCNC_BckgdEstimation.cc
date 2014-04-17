#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
//#include "TKey.h"
//#include "TRandom3.h"

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "../../Selection/interface/SelectionTable.h"
#include "../../Content/interface/AnalysisEnvironment.h"
#include "../../Content/interface/Dataset.h"
#include "../../Tools/interface/JetTools.h"
//#include "../../Tools/interface/PlottingTools.h"
#include "../../Tools/interface/MultiSamplePlot.h"
#include "../../Tools/interface/TTreeLoader.h"
#include "../../Tools/interface/AnalysisEnvironmentLoader.h"
#include "../../Reconstruction/interface/JetCorrectorParameters.h"
//#include "../../Reconstruction/interface/JetCorrectionUncertainty.h"
//#include "../../Reconstruction/interface/MakeBinning.h"
#include "../../MCInformation/interface/LumiReWeighting.h"
#include "../interface/TopFCNC_Evt.h"

#include "../../macros/Style.C"

using namespace std;
using namespace reweight;
using namespace TopTree;

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

struct HighestJPBtag{
  bool operator()( TRootJet* j1, TRootJet* j2 ) const{
    return j1->btag_jetProbabilityBJetTags() > j2->btag_jetProbabilityBJetTags();
  }
};
struct HighestCVSBtag{
  bool operator()( TRootJet* j1, TRootJet* j2 ) const{
    return j1->btag_combinedSecondaryVertexBJetTags() > j2->btag_combinedSecondaryVertexBJetTags();
  }
};

double MuEffSF_Id_Run2012(double eta, double pt);
double MuEffSF_Iso04_Run2012(double eta, double pt);
double MuEffSF_TrgMu8_Run2012(double eta, double pt);
double MuEffSF_TrgMu17_Run2012(double eta, double pt);

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
  
  /*https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
   Tagger name  	            WP name WP Discr cut
   TrackCountingHighPur 	    TCHPT 	3.41
   JetProbability 	          JPL 	  0.275
   JetProbability 	          JPM 	  0.545
   JetProbability 	          JPT 	  0.790
   CombinedSecondaryVertex 	CSVL 	  0.244
   CombinedSecondaryVertex 	CSVM 	  0.679
   CombinedSecondaryVertex 	CSVT 	  0.898
   */
  float btagcut     = 0.679;
	float Zmass       = 91.2;
	float Zwindowsize = 30.;
	bool applyAsymmJetPtCut = true;
	float JetPtCuts[4]={50.,30.,20.,20.};
  
  clock_t start = clock();
  
  cout << "*************************************************************" << endl;
  cout << " Beginning of the program for the FCNC search ! " << endl;
  cout << "*************************************************************" << endl;
  
  //SetStyle if needed
  //setTDRStyle();
  setGregStyle();
  //setMyStyle();
  
  string postfix = "_BckgdEstimation"; // to relabel the names of the output file
  
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
  
  string channelpostfix = "_MuElTrigger";
  string comments = "_Run2012A";
  string xmlFileName = "../config/myTopFCNCconfig_MuEl.xml";

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
  string rootFileName ("TopFCNC"+postfix+channelpostfix+comments+"_Plots.root");
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
  Float_t rho = 0;
//  Float_t     LeptIdSF = 0;
//  Float_t     LeptIsoSF= 0;
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////// Histograms /////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////////
  ////////////////// MultiSample plots  //////////////////////////////
  ////////////////////////////////////////////////////////////////////
  
  MSPlot["RhoCorrection"]                     = new MultiSamplePlot(datasets, "RhoCorrection", 100, 0, 100, "#rho");
  MSPlot["NbOfVertices"]                      = new MultiSamplePlot(datasets, "NbOfVertices", 30, 0, 30, "Nb. of vertices");
  
  MSPlot["MuonDzero"]                         = new MultiSamplePlot(datasets, "MuonDzero", 500, -0.02, 0.02, "d_{0} [cm]");
  MSPlot["ElectronDzero"]                     = new MultiSamplePlot(datasets, "ElectronDzero", 500, -0.02, 0.02, "d_{0} [cm]");
  
  MSPlot["1stLeadingMuonPt"]                  = new MultiSamplePlot(datasets, "1stLeadingMuonPt", 300, 0, 150, "p_{T} [GeV/c]");
  MSPlot["1stLeadingElectronPt"]              = new MultiSamplePlot(datasets, "1stLeadingElectronPt", 300, 0, 150, "p_{T} [GeV/c]");
  
  MSPlot["1stLeadingMuonRelIso"]              = new MultiSamplePlot(datasets, "1stLeadingMuonRelIso", 100, 0, 0.5, "RelIso");
  MSPlot["1stLeadingElectronRelIso"]          = new MultiSamplePlot(datasets, "1stLeadingElectronRelIso", 100, 0, 0.5, "RelIso");
  
  MSPlot["NbOfIsolatedMuons"]                 = new MultiSamplePlot(datasets, "NbOfIsolatedMuons", 5, 0, 5, "Nb. of isolated muons");
  MSPlot["NbOfIsolatedElectrons"]             = new MultiSamplePlot(datasets, "NbOfIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");
  
  MSPlot["DiLeptonInvMass_me_ch"]            = new MultiSamplePlot(datasets, "DiLeptonInvMass_me_ch", 400, 50, 130, "m_{ll}");
  MSPlot["DiLeptonSystPt_me_ch"]             = new MultiSamplePlot(datasets, "DiLeptonSystPt_me_ch", 400, 0, 400, "p_{T}^{ll} [GeV/c]");
  MSPlot["DiLeptonDR_me_ch"]                 = new MultiSamplePlot(datasets, "DiLeptonDR_me_ch", 400, 0, 10, "#Delta R(l^{+}l^{-})");
  
  MSPlot["DiLeptonSystPt_AtLeastFourJets_me_ch"] = new MultiSamplePlot(datasets, "DiLeptonSystPt_AtLeastFourJets_me_ch", 100, 0, 400, "p_{T}^{ll} [GeV/c]");
  MSPlot["DiLeptonDR_AtLeastFourJets_me_ch"]     = new MultiSamplePlot(datasets, "DiLeptonDR_AtLeastFourJets_me_ch", 100, 0, 10, "#Delta R(l^{+}l^{-})");
  
  MSPlot["NbOfExtraIsolatedMuons"]            = new MultiSamplePlot(datasets, "NbOfExtraIsolatedMuons", 5, 0, 5, "Nb. of isolated muons");
  MSPlot["NbOfExtraIsolatedElectrons"]        = new MultiSamplePlot(datasets, "NbOfExtraIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");
  
  MSPlot["NbOfSelectedJets_me_ch"]            = new MultiSamplePlot(datasets, "NbOfSelectedJets_me_ch",  15, 0, 15, "Nb. of jets");
  
  MSPlot["FirstLeadingJetPt_me_ch"]           = new MultiSamplePlot(datasets, "FirstLeadingJetPt_me_ch", 100, 0, 500, "Jet p_{T} [GeV/c]");
  MSPlot["SecondLeadingJetPt_me_ch"]          = new MultiSamplePlot(datasets, "SecondLeadingJetPt_me_ch", 80, 0, 400, "Jet p_{T} [GeV/c]");
  MSPlot["ThirdLeadingJetPt_me_ch"]           = new MultiSamplePlot(datasets, "ThirdLeadingJetPt_me_ch",  60, 0, 300, "Jet p_{T} [GeV/c]");
  MSPlot["FourthLeadingJetPt_me_ch"]          = new MultiSamplePlot(datasets, "FourthLeadingJetPt_me_ch", 50, 0, 250, "Jet p_{T} [GeV/c]");
  
  MSPlot["NbOfVertices_AtLeastFourJets_me_ch"]= new MultiSamplePlot(datasets, "NbOfVertices_AtLeastFourJets_me_ch", 30, 0, 30, "Nb. of vertices");
  
  MSPlot["HighestBdisc_me_ch_CVS"]            = new MultiSamplePlot(datasets, "HighestBdisc_me_ch_CVS", 50, 0, 1, "CSV b-disc.");
  
  MSPlot["MET_me_ch"]                         = new MultiSamplePlot(datasets, "MET_me_ch",  50, 0, 200, "\\slashE_{T} [GeV]");
  
  ////////////////////////////////////////////////////////////////////
  ////////////////// 1D histograms  //////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  
  histo1D["lumiWeights"] = new TH1F("lumiWeights","lumiWeights;lumiWeight;#events",50,0,4);
  for (unsigned int d = 0; d < datasets.size(); d++){
  	histo2D[("d0_vs_phi_1stleadingmuon_"+datasets[d]->Name()).c_str()] = new TH2F(("d0_vs_phi_1stleadingmuon_"+datasets[d]->Name()).c_str(),"d_{0}:#phi",500,-0.02,0.02,500,0,4);
  	histo2D[("d0_vs_phi_1stleadingelec_"+datasets[d]->Name()).c_str()] = new TH2F(("d0_vs_phi_1stleadingelec_"+datasets[d]->Name()).c_str(),"d_{0}:#phi",500,-0.02,0.02,500,0,4);
  }
  cout << " - Declared histograms ..." <<  endl;
	
  ////////////////////////////////////////////////////////////////////
  ////////////////// Plots  //////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  
  string pathPNG = "TopFCNC"+postfix+channelpostfix+comments+"_MSPlots/";
  //  pathPNG = pathPNG +"/";
  mkdir(pathPNG.c_str(),0777);
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////// Selection Tables ///////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  char zmasscutname[100];
  sprintf(zmasscutname,"$|m_{ll}-m_Z|<%f$ GeV",Zwindowsize);
  char btagcutname[100];
  sprintf(btagcutname,"$b\\texttt{-}disc \\geq %f$ (CSV)",btagcut);
  
  ////////////////////////////////////////////////////////////////////
  ///////////////////// Channel : µe  ////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  vector<string> CutsSelecTableMuEl;
  CutsSelecTableMuEl.push_back(string("initial"));
  CutsSelecTableMuEl.push_back(string("PU reweighting"));
  CutsSelecTableMuEl.push_back(string("Trigger"));
  CutsSelecTableMuEl.push_back(string("Good PV"));
  CutsSelecTableMuEl.push_back(string("$\\geq$ 2 isolated leptons"));
  CutsSelecTableMuEl.push_back(string(zmasscutname));//"$|m_{ll}-m_Z|<30$ GeV"));
  CutsSelecTableMuEl.push_back(string("Veto on 3rd iso. lept."));
  CutsSelecTableMuEl.push_back(string("$\\geq$ 1 jet"));
  CutsSelecTableMuEl.push_back(string("$\\geq$ 2 jet"));
  CutsSelecTableMuEl.push_back(string("$\\geq$ 3 jet"));
  CutsSelecTableMuEl.push_back(string("$\\geq$ 4 jet"));
  CutsSelecTableMuEl.push_back(string(btagcutname));//"$b\\texttt{-}disc \\geq 0.7$ (CSV)"));
  //  CutsSelecTableMuEl.push_back(string("$MET \\leq 30$ GeV"));
  
  SelectionTable selecTableMuEl(CutsSelecTableMuEl, datasets);
  selecTableMuEl.SetLuminosity(Luminosity);
  selecTableMuEl.SetPrecision(1);
    
  cout << " - SelectionTable instantiated ..." << endl;
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////// PileUp Reweighting ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  LumiReWeighting LumiWeights = LumiReWeighting("../pileup/pileup_MC_S10.root", "../pileup/MuEG_Run2012A_13Jul_TopTreeID_2009_PileupHistogram.root", "pileup", "pileup");
  
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
    cout<<"Load Dataset"<<endl;
		
    string previousFilename = "";
    int iFile = -1;
    
    string dataSetName = datasets[d]->Name();
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////// Initialize JEC factors /////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
    vector<JetCorrectorParameters> vCorrParam;
    
    // Create the JetCorrectorParameter objects, the order does not matter.
    // YYYY is the first part of the txt files: usually the global tag from which they are retrieved
    JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters("../../macros/JECFiles/Jec11V2_db_AK5PFchs_L1FastJet.txt");
    JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters("../../macros/JECFiles/Jec11V2_db_AK5PFchs_L2Relative.txt");
    JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters("../../macros/JECFiles/Jec11V2_db_AK5PFchs_L3Absolute.txt");
    
    //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!!
    vCorrParam.push_back(*L1JetPar);
    vCorrParam.push_back(*L2JetPar);
    vCorrParam.push_back(*L3JetPar);
    if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") // Data!
    {
      JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("../../macros/JECFiles/Jec11V2_db_AK5PFchs_L2L3Residual.txt");
      vCorrParam.push_back(*ResJetCorPar);
    }
    
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("../../macros/JECFiles/Jec11V2_db_AK5PFchs_Uncertainty.txt");
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true); // last boolean ('startFromRaw') = false!
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////// Create TTree ///////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    string TTreeFileName ("TopFCNC"+postfix+channelpostfix+comments+"_TTree_"+dataSetName+".root");
    
    cout << "INFO: creating file : "+TTreeFileName << endl;
    
    TFile* TTreeFile = new TFile(TTreeFileName.c_str(),"RECREATE");
    
    TTree* Tree = new TTree("Tree","Tree containing the TopFCNC event candidate");
    
    TopFCNC_Evt* MyTopFCNC_EvtCand = 0;
    
    Tree->Branch("TheTopFCNC_Evt","TopFCNC_Evt",&MyTopFCNC_EvtCand);
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////// Loop on events //////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    int itrigger1 = -1, itrigger2 = -1, previousRun = -1;
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
      if(! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) )
	    {
		    genjets = treeLoader.LoadGenJet(ievt,false);
	    }
      //cout << "run: " << event->runId() << "  lumi: " << event->lumiBlockId() << "  event: " << event->eventId() << endl;
      
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
        if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
        {
          /*------------------------------------------------------------------
           Dataset : MuEG/Run2012A-13Jul2012-v1
           --------------------------------------------------------------------
           Will search for runs with trigger HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v* available
           Trigger HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4 available for runs 190645-190738
           Trigger HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5 available for runs 191043-191411
           Trigger HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6 available for runs 191695-193621
           ------------------------------------------------------------------
           Will search for runs with trigger HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v* available
           Trigger HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4 available for runs 190645-190738
           Trigger HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5 available for runs 191043-191411
           Trigger HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6 available for runs 191695-193621
           ------------------------------------------------------------------*/
          if(currentRun >= 190645 && currentRun <= 190738){
            itrigger1 = treeLoader.iTrigger (string ("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4"), currentRun, iFile);
					  itrigger2 = treeLoader.iTrigger (string ("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4"), currentRun, iFile);
            // int. lumi = XXX/pb
					}
          else if (currentRun >= 191043 && currentRun <= 191411){
            itrigger1 = treeLoader.iTrigger (string ("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5"), currentRun, iFile);
            itrigger2 = treeLoader.iTrigger (string ("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5"), currentRun, iFile);
            // int. lumi = XXX/pb
          }
          else if (currentRun >= 191695 && currentRun <= 193621){
            itrigger1 = treeLoader.iTrigger (string ("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6"), currentRun, iFile);
            itrigger2 = treeLoader.iTrigger (string ("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6"), currentRun, iFile);
            // int. lumi = XXX/pb
          }
				  /*--------------------------------------------------------------------
           Sub-Total integrated luminosity = 778,2(/pb)
           Total integrated luminosity = 778,2(/pb)
           ------------------------------------------------------------------*/
          
				  /*------------------------------------------------------------------
           Dataset : MuEG/Run2012B-13Jul2012-v1
           --------------------------------------------------------------------
           Will search for runs with trigger HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v* available
           Trigger HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7 available for runs 193834-196531
           ------------------------------------------------------------------
           Will search for runs with trigger HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v* available
           Trigger HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7 available for runs 193834-196531
           ------------------------------------------------------------------*/
          else if (currentRun >= 193806 && currentRun <= 196531){
            itrigger1 = treeLoader.iTrigger (string ("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7"), currentRun, iFile);
            itrigger2 = treeLoader.iTrigger (string ("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7"), currentRun, iFile);
            // int. lumi = XXX/pb
          }
          /*--------------------------------------------------------------------
            Sub-Total integrated luminosity = XXX(/pb)
            Total integrated luminosity = XXX(/pb)
            ------------------------------------------------------------------*/
            
          if(itrigger1 == 9999 && itrigger2 == 9999)
          {
            cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (DATA) IN RUN " << event->runId() << endl;
            exit(1);
          }
          //trigged = treeLoader.EventTrigged (itrigger);
        }
        else{
          itrigger1 = -1;
          itrigger2 = -1;
        }
	    } //end previousRun != currentRun
      
	    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    //////////////////////////////////////////// Jet energy scale corrections     /////////////////////////////////////////////////////////////
	    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
	    // Apply Jet Corrections on-the-fly
	    //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JES correction on the fly:");
      /*
       if( dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
       jetTools->correctJets(init_jets,event->kt6PFJets_rho(),true); //last boolean: isData (needed for L2L3Residual...)
       else
       jetTools->correctJets(init_jets,event->kt6PFJets_rho(),false); //last boolean: isData (needed for L2L3Residual...)
       //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JES correction on the fly:");
       
       */
	    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    //////////////////////////////////////////// Type I MET corrections: (Only for |eta| <=4.7 ) //////////////////////////////////////////////
	    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /*
       //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before MET type I correction:");
       if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
       jetTools->correctMETTypeOne(init_jets,mets[0],true);
       else
       jetTools->correctMETTypeOne(init_jets,mets[0],false);
       //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After MET type I correction:");
       */
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////// Jet energy smearing and systematic uncertainty ///////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
      //if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
      if(! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) )
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
      /*
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
       */
      if(dataSetName != "Data" && dataSetName != "data" && dataSetName != "DATA")
      {
        // Apply pile-up reweighting
        double lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );
        scaleFactor *= lumiWeight;
      }
      
      histo1D["lumiWeights"]->Fill(scaleFactor);
      
      rho = event->kt6PFJets_rho();
      MSPlot["RhoCorrection"]->Fill(rho, datasets[d], true, Luminosity*scaleFactor);
			
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////// Event selection ////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
      selecTableMuEl.Fill(d,0, 1.);
      selecTableMuEl.Fill(d,1,scaleFactor);
      
      // Apply trigger selection
      if (itrigger1 == -1 && itrigger2 == -1) {
        trigged = true;
      }
      else {
        trigged = (treeLoader.EventTrigged (itrigger1) || treeLoader.EventTrigged (itrigger2));
      }
      
      if(!trigged) continue;
      selecTableMuEl.Fill(d,2,scaleFactor);
      
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
      selection.setJetCuts(20.,2.5,0.01,1.,0.98,0.3,0.1);
      
      selection.setDiElectronCuts(20,2.5,0.15,0.04,0.,1,0.3,1); //Et,Eta,RelIso,d0,MVAId,DistVzPVz,DRJets,MaxMissHits
      //selection.setLooseElectronCuts(15,2.5,0.2);
      
      selection.setDiMuonCuts(20,2.4,0.20,999.); //Et,Eta,RelIso,d0
      //selection.setLooseMuonCuts(15,2.4,0.2);
      
      //Select objects
      vector<TRootElectron*> selectedElectrons_NoIso = selection.GetSelectedDiElectrons(20,2.5,999.); //Et,Eta,RelIso
      vector<TRootElectron*> selectedElectrons       = selection.GetSelectedDiElectrons(20,2.5,0.15);//,rho);
      vector<TRootElectron*> selectedExtraElectrons;
      
      vector<TRootMuon*>     selectedMuons_NoIso = selection.GetSelectedDiMuons(20,2.4,999.); //Et,Eta,RelIso
      vector<TRootMuon*>     selectedMuons       = selection.GetSelectedDiMuons();
      vector<TRootMuon*>     selectedExtraMuons;
      
      vector<TRootJet*>      selectedJets        = selection.GetSelectedJets(true); // ApplyJetId
      
      // Apply primary vertex selection
      bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2);
      if(!isGoodPV) continue;
      selecTableMuEl.Fill(d,3,scaleFactor);
      
      // Select events with at exactly one muon and one electron and apply muon Id/Iso SF
      // TO_DO : - appply electron Id/Iso SF
      //         - apply cross-trigger MuEG SF

      if(selectedMuons_NoIso.size()<1 && selectedElectrons_NoIso.size()<1) continue;
      scaleFactor *= MuEffSF_Id_Run2012(selectedMuons_NoIso[0]->Eta(), selectedMuons_NoIso[0]->Pt());// Id SF
      
      MSPlot["NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
      
      histo2D[("d0_vs_phi_1stleadingmuon_"+datasets[d]->Name()).c_str()]->Fill(selectedMuons_NoIso[0]->d0(),selectedMuons_NoIso[0]->Phi());
      histo2D[("d0_vs_phi_1stleadingelec_"+datasets[d]->Name()).c_str()]->Fill(selectedElectrons_NoIso[1]->d0(),selectedElectrons_NoIso[1]->Phi());
      for(unsigned int i=0;i<selectedMuons_NoIso.size();i++)
        MSPlot["MuonDzero"]->Fill(selectedMuons_NoIso[i]->d0(), datasets[d], true, Luminosity*scaleFactor);
      for(unsigned int i=0;i<selectedElectrons_NoIso.size();i++)
        MSPlot["ElecDzero"]->Fill(selectedElectrons_NoIso[i]->d0(), datasets[d], true, Luminosity*scaleFactor);


      MSPlot["1stLeadingMuonPt"]->Fill(selectedMuons_NoIso[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["1stLeadingMuonRelIso"]->Fill(selectedMuons_NoIso[0]->relativePfIso(), datasets[d], true, Luminosity*scaleFactor);
      
      MSPlot["1stLeadingElecPt"]->Fill(selectedElectrons_NoIso[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["1stLeadingElecRelIso"]->Fill(selectedElectrons_NoIso[0]->relativePfIso(), datasets[d], true, Luminosity*scaleFactor);
      
      MSPlot["NbOfIsolatedMuons"]->Fill(selectedMuons.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["NbOfIsolatedElectrons"]->Fill(selectedElectrons.size(), datasets[d], true, Luminosity*scaleFactor);
      
      // Select events with at least two isolated muons
      if(selectedMuons.size()<1 && selectedElectrons.size()<1) continue;
      scaleFactor *= MuEffSF_Iso04_Run2012(selectedMuons[0]->Eta(), selectedMuons[0]->Pt());// Iso SF
      
      selecTableMuEl.Fill(d,4,scaleFactor);
      
      bool foundZ = false;
      int idx_Z_mu = -1, idx_Z_el = -1;
      // Calculate the invariant mass for each isolated lepton pairs
      // - return true if the mass is the Z boson mass window
      // - return the indices of the lepton candidates
      for(unsigned int i=0;i<selectedMuons.size();i++)
      {
        for(unsigned int j=0;j<selectedElectrons.size();j++)
        {
          TRootMuon*     mu = (TRootMuon*)    selectedMuons[i];
          TRootElectron* el = (TRootElectron*) selectedElectrons[j];
          if(mu->charge() == el->charge()) continue;
          double invMass = (*mu + *el).M();
          MSPlot["DiLeptonInvMass_me_ch"]->Fill(invMass, datasets[d], true, Luminosity*scaleFactor);
          if( invMass >= (Zmass-Zwindowsize) && invMass <= (Zmass+Zwindowsize) )
          {
            idx_Z_mu = i;
            idx_Z_el = j;
            foundZ = true;
          }
        }
      }
      // Select events with at least one pair of opposite charge leptons with |mll-mz|<windowsize
      if(!foundZ)
        continue;
      
      float DiLeptonSystPt = (*selectedMuons[idx_Z_mu] + *selectedElectrons[idx_Z_el]).Pt();
      float DiLeptonDR     =   selectedMuons[idx_Z_mu]->DeltaR(*selectedElectrons[idx_Z_el]);
      MSPlot["DiLeptonSystPt_me_ch"]->Fill(DiLeptonSystPt, datasets[d], true, Luminosity*scaleFactor);
      MSPlot["DiLeptonDR_me_ch"]    ->Fill(DiLeptonDR,     datasets[d], true, Luminosity*scaleFactor);
      
      selecTableMuEl.Fill(d,5,scaleFactor);
      
      selectedExtraMuons = selectedMuons;
      selectedExtraElectrons = selectedElectrons;
      // Erase Z boson lepton candidates
      selectedExtraMuons.erase(selectedExtraMuons.begin()+idx_Z_mu);
      selectedExtraElectrons.erase(selectedExtraElectrons.begin()+idx_Z_el);
      
      MSPlot["NbOfExtraIsolatedMuons"]->Fill(selectedExtraMuons.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["NbOfExtraIsolatedElectrons"]->Fill(selectedExtraElectrons.size(), datasets[d], true, Luminosity*scaleFactor);
      
      MyTopFCNC_EvtCand = 0;
      double invMass = 0;
      double highestbtagdisc = 0;
      // Select events based on the presence of *exactly one* extra isolated lepton
      if(selectedExtraMuons.size()==0 && selectedElectrons.size()==0){
        selecTableMuEl.Fill(d,6,scaleFactor);
        MSPlot["NbOfSelectedJets_me_ch"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
        if(selectedJets.size()>0){ //at least 1 jet
          MSPlot["FirstLeadingJetPt_me_ch"]->Fill(selectedJets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
          if(applyAsymmJetPtCut && selectedJets[0]->Pt()<JetPtCuts[0]) continue;
          selecTableMuEl.Fill(d,7,scaleFactor);
          if(selectedJets.size()>1){ //at least 2 jets
            MSPlot["SecondLeadingJetPt_me_ch"]->Fill(selectedJets[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
            if(applyAsymmJetPtCut && selectedJets[1]->Pt()<JetPtCuts[1]) continue;
            selecTableMuEl.Fill(d,8,scaleFactor);
            if(selectedJets.size()>2){ //at least 3 jets
              MSPlot["ThirdLeadingJetPt_me_ch"]->Fill(selectedJets[2]->Pt(), datasets[d], true, Luminosity*scaleFactor);
              selecTableMuEl.Fill(d,9,scaleFactor);
              if(selectedJets.size()>3){ //at least 4 jets
                selecTableMuEl.Fill(d,10,scaleFactor);
                
                MSPlot["FourthLeadingJetPt_me_ch"]->Fill(selectedJets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["NbOfVertices_AtLeastFourJets_me_ch"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["DiLeptonDR_AtLeastFourJets_me_ch"]->Fill(DiLeptonDR, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["DiLeptonSystPt_AtLeastFourJets_me_ch"]->Fill(DiLeptonSystPt, datasets[d], true, Luminosity*scaleFactor);
                
                sort(selectedJets.begin(),selectedJets.end(),HighestCVSBtag());
                MSPlot["HighestBdisc_me_ch_CVS"]->Fill(selectedJets[0]->btag_combinedSecondaryVertexBJetTags(),datasets[d], true, Luminosity*scaleFactor);
                highestbtagdisc = selectedJets[0]->btag_combinedSecondaryVertexBJetTags();
                sort(selectedJets.begin(),selectedJets.end(),HighestJPBtag());
                MSPlot["HighestBdisc_me_ch_JP"]->Fill(selectedJets[0]->btag_jetProbabilityBJetTags(),datasets[d], true, Luminosity*scaleFactor);
                
                // Create TopFCNC_Evt object
                MyTopFCNC_EvtCand = new TopFCNC_Evt(TopFCNC_Evt::kMuon);
                MyTopFCNC_EvtCand->SetLepton1FromZ(*selectedMuons[idx_Z_mu]);
                MyTopFCNC_EvtCand->SetLepton2FromZ(*selectedElectrons[idx_Z_el]);
                MyTopFCNC_EvtCand->SetSelectedJets(selectedJets);
                MyTopFCNC_EvtCand->SetMET(*mets[0]);
                
                if(highestbtagdisc>btagcut){
                  selecTableMuEl.Fill(d,11,scaleFactor);
                  MSPlot["MET_me_ch"]->Fill(mets[0]->Et(),datasets[d], true, Luminosity*scaleFactor);
                }
              }
            }
          }
        }
      }

      if(MyTopFCNC_EvtCand){
        MyTopFCNC_EvtCand->SetEventID( event->eventId() );
        MyTopFCNC_EvtCand->SetRunID( event->runId() );
        MyTopFCNC_EvtCand->SetLumiBlockID( event->lumiBlockId() );
        MyTopFCNC_EvtCand->SetIdParton1( event->idParton1() );
        MyTopFCNC_EvtCand->SetIdParton2( event->idParton2() );
        MyTopFCNC_EvtCand->SetxParton1( event->xParton1() );
        MyTopFCNC_EvtCand->SetxParton2( event->xParton2() );
        MyTopFCNC_EvtCand->SetFactorizationScale( event->factorizationScale() );
        MyTopFCNC_EvtCand->SetnPV(vertex.size());
        MyTopFCNC_EvtCand->SetnTruePU(event->nTruePU());
        MyTopFCNC_EvtCand->SetEventWeight(scaleFactor);
      }
      Tree->Fill();
      
      //delete selection;
      if(MyTopFCNC_EvtCand) delete MyTopFCNC_EvtCand;
    }//loop on events
    
    cout<<endl;
    cout<<"FYI ; nb of events with at least four isolated leptons = "<<fourIsoLeptCounter<<endl;
    
    TTreeFile->cd();
    
    TTree *configTree = new TTree("configTree","configuration Tree");
    TClonesArray* tcdataset = new TClonesArray("Dataset",1);
    configTree->Branch("Dataset","TClonesArray",&tcdataset);
    new ((*tcdataset)[0]) Dataset(*datasets[d]);
    
    configTree->Fill();
    configTree->Write();
    Tree->Write();
    TTreeFile->Close();
    delete TTreeFile;
    
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

  //(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)	
  selecTableMuEl.TableCalculator(  false, true, true, true, true);
  //Options : WithError (false), writeMerged (true), useBookTabs (false), addRawNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false)
  selecTableMuEl.Write(  "TopFCNC"+postfix+channelpostfix+comments+"_SelectionTable_DiMu.tex",    true,true,true,true,false,false,true);

  fout->cd();
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
	  MultiSamplePlot *temp = it->second;
	  //temp->addText("CMS preliminary");
	  string name = it->first;
	  name += comments;
	  temp->Draw(false, name, true, true, true, true, true,1,true); // merge TT/QCD/W/Z/ST/
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
double MuEffSF_Id_Run2012(double eta, double pt){
  if( fabs(eta) >= 0 && fabs(eta) < 0.9 ){
    if( pt >= 10 && pt < 20)
      return 0.981239;
    else if( pt >= 20 && pt < 25)
      return 0.98043;
    else if( pt >= 25 && pt < 30)
      return 0.991182;
    else if( pt >= 30 && pt < 35)
      return 0.988984;
    else if( pt >= 35 && pt < 40)
      return 0.990422;
    else if( pt >= 40 && pt < 50)
      return 0.988821;
    else if( pt >= 50 && pt < 60)
      return 0.989954;
    else if( pt >= 60 && pt < 90)
      return 0.989478;
    else if( pt >= 90 && pt < 140)
      return 1.00023;
    else if( pt >= 140 && pt < 500)
      return 0.985242;
    else
      return 1;
  }
  if( fabs(eta) >= 0.9 && fabs(eta) < 1.2 ){
    if( pt >= 10 && pt < 20)
      return 0.936552;
    else if( pt >= 20 && pt < 25)
      return 0.995381;
    else if( pt >= 25 && pt < 30)
      return 0.991899;
    else if( pt >= 30 && pt < 35)
      return 0.990177;
    else if( pt >= 35 && pt < 40)
      return 0.990865;
    else if( pt >= 40 && pt < 50)
      return 0.986703;
    else if( pt >= 50 && pt < 60)
      return 0.990343;
    else if( pt >= 60 && pt < 90)
      return 0.985462;
    else if( pt >= 90 && pt < 140)
      return 0.999933;
    else if( pt >= 140 && pt < 500)
      return 0.922726;
    else
      return 1;
  }
  if( fabs(eta) >= 1.2 && fabs(eta) < 2.1 ){
    if( pt >= 10 && pt < 20)
      return 1.00304;
    else if( pt >= 20 && pt < 25)
      return 0.998942;
    else if( pt >= 25 && pt < 30)
      return 0.996901;
    else if( pt >= 30 && pt < 35)
      return 0.997486;
    else if( pt >= 35 && pt < 40)
      return 0.994566;
    else if( pt >= 40 && pt < 50)
      return 0.996159;
    else if( pt >= 50 && pt < 60)
      return 0.997451;
    else if( pt >= 60 && pt < 90)
      return 0.996516;
    else if( pt >= 90 && pt < 140)
      return 1.03286;
    else if( pt >= 140 && pt < 500)
      return 1.05323;
    else
      return 1;
  }
}
double MuEffSF_Iso04_Run2012(double eta, double pt){
  if( fabs(eta) >= 0 && fabs(eta) < 0.9 ){
    if( pt >= 10 && pt < 20)
      return 0.922624;
    else if( pt >= 20 && pt < 25)
      return 0.976962;
    else if( pt >= 25 && pt < 30)
      return 0.997654;
    else if( pt >= 30 && pt < 35)
      return 0.997849;
    else if( pt >= 35 && pt < 40)
      return 0.998674;
    else if( pt >= 40 && pt < 50)
      return 0.998288;
    else if( pt >= 50 && pt < 60)
      return 0.998246;
    else if( pt >= 60 && pt < 90)
      return 0.99948;
    else if( pt >= 90 && pt < 140)
      return 1.00003;
    else if( pt >= 140 && pt < 500)
      return 0.996181;
    else
      return 1;
  }
  if( fabs(eta) >= 0.9 && fabs(eta) < 1.2 ){
    if( pt >= 10 && pt < 20)
      return 0.948552;
    else if( pt >= 20 && pt < 25)
      return 0.981943;
    else if( pt >= 25 && pt < 30)
      return 0.996887;
    else if( pt >= 30 && pt < 35)
      return 0.999591;
    else if( pt >= 35 && pt < 40)
      return 1.00033;
    else if( pt >= 40 && pt < 50)
      return 0.999218;
    else if( pt >= 50 && pt < 60)
      return 0.998999;
    else if( pt >= 60 && pt < 90)
      return 0.99905;
    else if( pt >= 90 && pt < 140)
      return 0.997391;
    else if( pt >= 140 && pt < 500)
      return 1.00422;
    else
      return 1;
  }
  if( fabs(eta) >= 1.2 && fabs(eta) < 2.1 ){
    if( pt >= 10 && pt < 20)
      return 0.970175;
    else if( pt >= 20 && pt < 25)
      return 0.989697;
    else if( pt >= 25 && pt < 30)
      return 1.0003;
    else if( pt >= 30 && pt < 35)
      return 1.00058;
    else if( pt >= 35 && pt < 40)
      return 1.00088;
    else if( pt >= 40 && pt < 50)
      return 0.999595;
    else if( pt >= 50 && pt < 60)
      return 0.999906;
    else if( pt >= 60 && pt < 90)
      return 0.999467;
    else if( pt >= 90 && pt < 140)
      return 0.997148;
    else if( pt >= 140 && pt < 500)
      return 0.997978;
    else
      return 1;
  }
}
double MuEffSF_TrgMu8_Run2012(double eta, double pt){
  if( fabs(eta) >= 0 && fabs(eta) < 0.9 ){
    if( pt >= 10 && pt < 20)
      return 0.991061;
    else if( pt >= 20 && pt < 25)
      return 0.988522;
    else if( pt >= 25 && pt < 30)
      return 0.98938;
    else if( pt >= 30 && pt < 35)
      return 0.987832;
    else if( pt >= 35 && pt < 40)
      return 0.989023;
    else if( pt >= 40 && pt < 50)
      return 0.988155;
    else if( pt >= 50 && pt < 60)
      return 0.987275;
    else if( pt >= 60 && pt < 90)
      return 0.989316;
    else if( pt >= 90 && pt < 140)
      return 0.990073;
    else if( pt >= 140 && pt < 500)
      return 0.982128;
    else
      return 1;
  }
  if( fabs(eta) >= 0.9 && fabs(eta) < 1.2 ){
    if( pt >= 10 && pt < 20)
      return 1.00247;
    else if( pt >= 20 && pt < 25)
      return 0.98477;
    else if( pt >= 25 && pt < 30)
      return 0.985676;
    else if( pt >= 30 && pt < 35)
      return 0.983014;
    else if( pt >= 35 && pt < 40)
      return 0.983788;
    else if( pt >= 40 && pt < 50)
      return 0.983716;
    else if( pt >= 50 && pt < 60)
      return 0.985706;
    else if( pt >= 60 && pt < 90)
      return 0.982735;
    else if( pt >= 90 && pt < 140)
      return 0.982356;
    else if( pt >= 140 && pt < 500)
      return 0.963695;
    else
      return 1;
  }
  if( fabs(eta) >= 1.2 && fabs(eta) < 2.1 ){
    if( pt >= 10 && pt < 20)
      return 1.00883;
    else if( pt >= 20 && pt < 25)
      return 1.00035;
    else if( pt >= 25 && pt < 30)
      return 0.993731;
    else if( pt >= 30 && pt < 35)
      return 0.990587;
    else if( pt >= 35 && pt < 40)
      return 0.987497;
    else if( pt >= 40 && pt < 50)
      return 0.985698;
    else if( pt >= 50 && pt < 60)
      return 0.98527;
    else if( pt >= 60 && pt < 90)
      return 0.983774;
    else if( pt >= 90 && pt < 140)
      return 0.971552;
    else if( pt >= 140 && pt < 500)
      return 1.00464;
    else
      return 1;
  }
}
double MuEffSF_TrgMu17_Run2012(double eta, double pt){
  if( fabs(eta) >= 0 && fabs(eta) < 0.9 ){
    if( pt >= 10 && pt < 20)
      return 0.991061;
    else if( pt >= 20 && pt < 25)
      return 0.988522;
    else if( pt >= 25 && pt < 30)
      return 0.98938;
    else if( pt >= 30 && pt < 35)
      return 0.987832;
    else if( pt >= 35 && pt < 40)
      return 0.989023;
    else if( pt >= 40 && pt < 50)
      return 0.988155;
    else if( pt >= 50 && pt < 60)
      return 0.987275;
    else if( pt >= 60 && pt < 90)
      return 0.989316;
    else if( pt >= 90 && pt < 140)
      return 0.990073;
    else if( pt >= 140 && pt < 500)
      return 0.982128;
    else
      return 1;
  }
  if( fabs(eta) >= 0.9 && fabs(eta) < 1.2 ){
    if( pt >= 10 && pt < 20)
      return 1.00247;
    else if( pt >= 20 && pt < 25)
      return 0.98477;
    else if( pt >= 25 && pt < 30)
      return 0.985676;
    else if( pt >= 30 && pt < 35)
      return 0.983014;
    else if( pt >= 35 && pt < 40)
      return 0.983788;
    else if( pt >= 40 && pt < 50)
      return 0.983716;
    else if( pt >= 50 && pt < 60)
      return 0.985706;
    else if( pt >= 60 && pt < 90)
      return 0.982735;
    else if( pt >= 90 && pt < 140)
      return 0.982356;
    else if( pt >= 140 && pt < 500)
      return 0.963695;
    else
      return 1;
  }
  if( fabs(eta) >= 1.2 && fabs(eta) < 2.1 ){
    if( pt >= 10 && pt < 20)
      return 1.00883;
    else if( pt >= 20 && pt < 25)
      return 1.00035;
    else if( pt >= 25 && pt < 30)
      return 0.993731;
    else if( pt >= 30 && pt < 35)
      return 0.990587;
    else if( pt >= 35 && pt < 40)
      return 0.987497;
    else if( pt >= 40 && pt < 50)
      return 0.985698;
    else if( pt >= 50 && pt < 60)
      return 0.98527;
    else if( pt >= 60 && pt < 90)
      return 0.983774;
    else if( pt >= 90 && pt < 140)
      return 0.971552;
    else if( pt >= 140 && pt < 500)
      return 1.00464;
    else
      return 1;
  }
}
