///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//////   Skelton macro for T1tttt analysis in the multi lepton final state(s). ////////
//////    J. Keaveney, K.Deroover                                               ///////
//////                                                                          ///////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////



#include "TStyle.h"
#include "TPaveText.h"
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
#include "../MCInformation/interface/LumiReWeighting.h"
#include "../Reconstruction/interface/MEzCalculator.h"
#include "../Reconstruction/interface/TTreeObservables.h"
#include "../InclFourthGenSearch/interface/InclFourthGenTree.h"
#include "../InclFourthGenSearch/interface/InclFourthGenSearchTools.h"
#include "../macros/Style.C"

using namespace std;
using namespace TopTree;
using namespace reweight;


//Setting this bool to true will split the ttbar + jets
// dataset into tt + ll, tt + cc and tt + bb datasets.
bool split_ttbar = false;

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

/// MultiPadPlot
map<string,MultiSamplePlot*> MultiPadPlot;

//structs which return jet with largest value 
// of b-tag discriminant.

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

//Main body of analysis code starts here...
int main (int argc, char *argv[])
{

  //these variables control the application of a shift 
  // to the Jet Energy Scale (JES)
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

    
 //Choose which b-tag algorithm will be used.
  string btagger = "CSVM";
// b-tag scalefactor => TCHEL: data/MC scalefactor = 0.95 +- 0.10,    TCHEM: data/MC scalefactor = 0.94 +- 0.09
// mistag scalefactor => TCHEL: data/MC scalefactor = 1.11 +- 0.12,    TCHEM: data/MC scalefactor = 1.21 +- 0.17
  float scalefactorbtageff, mistagfactor;
  if(btagger == "TCHPM"  || btagger == "TCHET"  ||  btagger == "SSV" ){
    cout<<"This tagger ("<< btagger <<")is not commisioned in 2012, please use CSV, TCHP or JetProb"<<endl;
    exit(1);
}
  else if(btagger == "TCHEM") //redundant for now (these values need updating), but will use as skeleton for CSVM
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
    
  float workingpointvalue = 9999; //working points updated to 2012 BTV-POG recommendations.
 
  if(btagger == "TCHPM"  || btagger == "TCHET"  ||  btagger == "SSV" ){
    cout<<"This tagger ("<< btagger <<")is not commisioned in 2012, please use CSV, TCHP or JetProb"<<endl;
    exit(1); 
  }
   else if(btagger == "CSVL")
     workingpointvalue = .244;	
  else if(btagger == "CSVM")
    workingpointvalue = .679;
  else if(btagger == "CSVT")
    workingpointvalue = .898;

  clock_t start = clock();

 cout << "*************************************************************" << endl;
  cout << " Beginning of the program for the FourTop search ! "           << endl;
  cout << "*************************************************************" << endl;

  //SetStyle if needed (for plots)
  setTDRStyle();
  //setGregStyle();
  //setMyStyle();

  string postfix = "_EventSelection"; // to relabel the names of the output file

//flags which control the application of scale to the JES, BTag eff, Mistag eff
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

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////// Configuration ///////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 

  string channelpostfix = "";
  string xmlFileName = "";

  bool Electron = false; // use Electron channel?
  bool Muon = true; // use Muon channel?
  if(Electron && Muon){
	cout << "  --> Using both Muon and Electron channel? Choose only one ( since different samples/skims are required)!" << endl;
	exit(1);
  }

    
  //specify the .xml (containing lists of skimmed toptrees you want to use)
  xmlFileName = "config/test_fullsamples.xml";

  const char *xmlfile = xmlFileName.c_str();
  cout << "used config file: " << xmlfile << endl;    


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////// AnalysisEnvironment /////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //This class reads the .xml file and makes the information available
    //here in the macro
  AnalysisEnvironment anaEnv;
  cout<<" - Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  int verbose = 2;//anaEnv.Verbose;


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////// Load Datasets ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TTreeLoader treeLoader;
  vector < Dataset* > datasets;
  cout << " - Loading datasets ..." << endl;
  treeLoader.LoadDatasets (datasets, xmlfile);

//This is a default integrated lumi which should be overwritten by the lumi
//read in from the xml file.
  float Luminosity = 5343.64; //pb^-1??
  
  for (unsigned int d = 0; d < datasets.size (); d++)
  {
     cout <<"found sample with equivalent lumi "<<  datasets[d]->EquivalentLumi() <<endl;
     string dataSetName = datasets[d]->Name();
   if(dataSetName.find("Data")<=0 || dataSetName.find("data")<=0 || dataSetName.find("DATA")<=0)
	{
		  Luminosity = datasets[d]->EquivalentLumi();//reading in integrated lumi from xml
		  cout <<"found DATA sample with equivalent lumi "<<  datasets[d]->EquivalentLumi() <<endl;
		  break;
	 }
  }
    
    
  cout << "Rescaling to an integrated luminosity of "<< Luminosity <<" pb^-1" << endl;
    //Check that this figure is correct...

    //N.B. for splitting the ttbar sample, it is essential to have the ttjets sample as the last dataset loaded
    if (split_ttbar){
        //creating three new datsets to contain the tt+ll, tt+cc and tt+ bb compenents.
        int ndatasets = datasets.size() - 1 ;
        cout << " - splitting TTBar dataset ..." << ndatasets   << endl;
        vector<string> ttbar_filenames = datasets[ndatasets]->Filenames();
        cout <<"ttbar filenames =  "<< ttbar_filenames[0] <<endl;
        
        Dataset* ttbar_ll = new Dataset("TTJets_ll","tt + ll" , true, 633, 2, 2, 1, 213.4,ttbar_filenames );
        Dataset* ttbar_cc = new Dataset("TTJets_cc","tt + cc" , true, 633, 2, 2, 1, 6.9, ttbar_filenames );
        Dataset* ttbar_bb = new Dataset("TTJets_bb","tt + bb" , true, 633, 2, 2, 1, 4.8, ttbar_filenames );
        
        ttbar_ll->SetEquivalentLuminosity(30483.61);
        ttbar_cc->SetEquivalentLuminosity(30483.61);
        ttbar_bb->SetEquivalentLuminosity(30483.61);
        
        ttbar_ll->SetColor(kRed);
        ttbar_cc->SetColor(kRed-3);
        ttbar_bb->SetColor(kRed+2);
        
        datasets.pop_back();
        datasets.push_back(ttbar_bb);
        datasets.push_back(ttbar_cc);
        datasets.push_back(ttbar_ll);     
    }     
    
  //Output ROOT file, specify any names you like.
  string rootFileName ("FourTop_MultiLepton"+postfix+channelpostfix+".root");
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");

  //vectors of objects that we want to study
  cout << " - Variable declaration ..." << endl;
  vector < TRootVertex* >   vertex;
  vector < TRootMuon* >     init_muons;
  vector < TRootElectron* > init_electrons;
  vector < TRootJet* >      init_jets;
  vector < TRootMET* >      mets;

  //Global variable
  TRootEvent* event = 0;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////// Histograms /////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////// MultiSample plots: convenient class which combines multiple MC and DATA histograms into single plots.  //////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    MSPlot["RhoCorrection"]              = new MultiSamplePlot(datasets, "RhoCorrection", 100, 0, 100, "#rho");
    MSPlot["NbOfVertices"]               = new MultiSamplePlot(datasets, "NbOfVertices", 40, 0, 40, "Nb. of vertices");
    
    //Muons
    MSPlot["NbOfIsolatedMuons"]          = new MultiSamplePlot(datasets, "NbOfIsolatedMuons", 5, 0, 5, "Nb. of isolated muons");
    MSPlot["NbOfIsolatedElectrons"]      = new MultiSamplePlot(datasets, "NbOfIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");
    MSPlot["NbOfExtraIsolatedMuons"]     = new MultiSamplePlot(datasets, "NbOfExtraIsolatedMuons", 5, 0, 5, "Nb. of isolated muons");
    MSPlot["NbOfExtraIsolatedElectrons"] = new MultiSamplePlot(datasets, "NbOfExtraIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");
    MSPlot["MuonRelIsolation"] = new MultiSamplePlot(datasets, "MuonRelIsolation", 50, 0, .15, "RelIso");
    MSPlot["MuonPt"]              = new MultiSamplePlot(datasets, "MuonPt", 75, 0, 300, "PT_{#mu}");
    MSPlot["MuonEta"]              = new MultiSamplePlot(datasets, "MuonEta", 25, -2.4, 2.4, "#eta_{#mu}");
    MSPlot["MuonPhi"]              = new MultiSamplePlot(datasets, "MuonPhi", 50, -4, 4, "#phi_{#mu}");
    MSPlot["MuonNValidHits"]              = new MultiSamplePlot(datasets, "MuonNValidHits", 30, 0, 30, "NValidHits_{#mu}");
    MSPlot["Muond0"]              = new MultiSamplePlot(datasets, "Muond0", 50, 0, .05, "d0_{#mu}");
    MSPlot["MuondRJets"]              = new MultiSamplePlot(datasets, "MuondRJets", 50, 0, 10, "dRJets_{#mu}");
    MSPlot["MuonNMatchedStations"]              = new MultiSamplePlot(datasets, "MuonNMatchedStations", 10, 0, 10, "NMatchedStations_{#mu}");
    MSPlot["MuonDistVzPVz"]              = new MultiSamplePlot(datasets, "MuonDistVzPVz", 75, -.3,.3, "DistVzPVz_{#mu}");
    MSPlot["MuonTrackerLayersWithMeasurement"]    = new MultiSamplePlot(datasets, "MuonTrackerLayersWithMeasurement", 25, 0, 25, "nLayers");
    MSPlot["DiMuon_InvMass"]     = new MultiSamplePlot(datasets, "DiMuon_InvMass", 60, 0, 120, "DiMuon_InvMass");
    MSPlot["NbOfLooseMuon"]     = new MultiSamplePlot(datasets, "NbOfLooseMuon", 10, 0, 10, "Nb. of loose muons");
    
    //B-tagging discriminators
    MSPlot["BdiscBJetCand_CSV"]  = new MultiSamplePlot(datasets, "HighestBdisc_CSV", 75, 0, 1, "CSV b-disc.");
    MSPlot["HighestBdisc_m_ch_CSV"]            = new MultiSamplePlot(datasets, "HighestBdisc_mm_ch_CVS", 100, 0, 1, "CSV b-disc.");
    MSPlot["HighestBdisc_e_ch_CSV"]            = new MultiSamplePlot(datasets, "HighestBdisc_ee_ch_CVS", 100, 0, 1, "CSV b-disc.");
    
    //Jets
    MSPlot["NbOfSelectedJets"]                  = new MultiSamplePlot(datasets, "NbOfSelectedJets", 15, 0, 15, "Nb. of jets");
    MSPlot["NbOfSelectedLightJets"]                  = new MultiSamplePlot(datasets, "NbOfSelectedLightJets", 10, 0, 10, "Nb. of jets");
    MSPlot["NbOfSelectedBJets"]                  = new MultiSamplePlot(datasets, "NbOfSelectedBJets", 8, 0, 8, "Nb. of jets");
    MSPlot["JetEta"]                  = new MultiSamplePlot(datasets, "JetEta", 30,-3, 3, "Jet #eta");
    MSPlot["JetPhi"]                  = new MultiSamplePlot(datasets, "JetPhi", 50, -4,4 , "Jet #phi");
    
    MSPlot["FirstTriJetMass_1BJet_g"] = new MultiSamplePlot(datasets, "FirstTriJetMass_1BJet_g", 15, 0, 1000, "m_{bjj}");
    MSPlot["FirstTriJetMass_1BJet_All"] = new MultiSamplePlot(datasets, "FirstTriJetMass_1BJet_All", 50, 0, 1000, "m_{bjj}");
    MSPlot["FirstTriJetPt_1BJet_g"] = new MultiSamplePlot(datasets, "FirstTriJetPt_1BJet_g", 50, 0, 1000, "Pt_{bjj}");
    MSPlot["FirstDiJetMass"] = new MultiSamplePlot(datasets, "FirstDiJetMass", 50, 0, 200, "m_{jj}");
    MSPlot["FirstDiJetMass_g"] = new MultiSamplePlot(datasets, "FirstDiJetMass_g", 50, 0, 200, "m_{jj}");
      
    
    MSPlot["SelectedJetPt"] = new MultiSamplePlot(datasets, "JetPt", 50, 0, 1000, "PT_{jet}");
    MSPlot["4thJetPt"] = new MultiSamplePlot(datasets, "4thJetPt", 60, 0, 400, "PT_{jet}");
    MSPlot["5thJetPt"] = new MultiSamplePlot(datasets, "5thJetPt", 60, 0, 400, "PT_{jet}");
    MSPlot["6thJetPt"] = new MultiSamplePlot(datasets, "6thJetPt", 60, 0, 400, "PT_{jet}");
    MSPlot["7thJetPt"] = new MultiSamplePlot(datasets, "7thJetPt", 60, 0, 400, "PT_{jet}");
    MSPlot["SelectedJetPt_light"] = new MultiSamplePlot(datasets, "JetPt_light", 50, 0, 1000, "PT_{lightjet}");
    MSPlot["SelectedJetPt_b"] = new MultiSamplePlot(datasets, "JetPt_b", 50, 0, 1000, "PT_{bjet}");
    MSPlot["HT_SelectedJets"] = new MultiSamplePlot(datasets, "HT_SelectedJets", 50, 0, 1500, "HT");
    MSPlot["MHT_SelectedJets"] = new MultiSamplePlot(datasets, "MHT_SelectedJets", 75, 0, 1000, "MHT");
    MSPlot["MHTSig_SelectedJets"] = new MultiSamplePlot(datasets, "MHTSig_SelectedJets", 75, 0, 30, "MHTSig");
    MSPlot["MET"] = new MultiSamplePlot(datasets, "MET", 40, 0, 700, "MET");
    MSPlot["MET_MHT"]= new MultiSamplePlot(datasets, "MET_MHT", 75, 0, 200, "MET_MHT");
    MSPlot["STLep"] = new MultiSamplePlot(datasets, "STLep", 20, 0, 1000, "STLep");
    MSPlot["STJet"] = new MultiSamplePlot(datasets, "STJet", 20, 0, 1000, "STJet");

    //Electrons
    MSPlot["ElectronPt"]              = new MultiSamplePlot(datasets, "ElectronPt", 50, 0, 100, "PT_{e}");
    MSPlot["NbOfLooseElectron"] = new MultiSamplePlot(datasets, "NbOfLooseElectron", 10, 0, 10, "Nb. of loose electrons");

  //Declare arrays of MSPlots: We may want to look at certain variables in different event categories, i.e.
  // events with X jets and Y b tags, so we declare arrays of MultiSamplePlots, with one 
  //MultiSamplePlot corresponding to each category.
    
  Int_t minNJets=3, maxNJets=5, minNBJets=0, maxNBJets=3;

  for (Int_t q = minNJets; q <= maxNJets; q++){
    for (Int_t p = minNBJets; p<= maxNBJets; p++){

    string NJets_str = static_cast<ostringstream*>( &(ostringstream() << q) )->str();
    string NBJets_str = static_cast<ostringstream*>( &(ostringstream() << p) )->str();
   // string HT_Name = "HT_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string MET_Name = "MET_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string InclMET_Name = "InclMET_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    //string HTX_Name = "HTX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    //string EventMass_Name = "EventMass_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    //string EventMassX_Name = "EventMassX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    //string FirstTriJetMass_1BJet_g_Name = "FirstTriJetMass_1BJet_g"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    //string FirstDiJetMass_1BJet_g_Name = "FirstDiJetMass_1BJet_g"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    //string SecondDiJetMass_1BJet_g_Name = "SecondDiJetMass_1BJet_g"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    //string SecondTriJetMass_1BJet_g_Name = "SecondTriJetMass_1BJet_g"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    //string SecondTriJetMass_1BJet_g_chi2_Name = "SecondTriJetMass_1BJet_g_chi2_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    //string MuMetBMasses_g_Name = "MuMetBMasses_g_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    //string MuMetMasses_g_Name = "MuMetMasses_g_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;

        
    //MSPlot[MuMetMasses_g_Name.c_str() ] = new MultiSamplePlot(datasets, MuMetMasses_g_Name.c_str() , 50, 0, 700, "M_{muMET}");
    //MSPlot[MuMetBMasses_g_Name.c_str() ] = new MultiSamplePlot(datasets, MuMetBMasses_g_Name.c_str() , 50, 0, 700, "M_{muMETb}");
    //MSPlot[HT_Name.c_str() ] = new MultiSamplePlot(datasets, HT_Name.c_str() , 10, 0, 1700, "HT");
    //MSPlot[HTX_Name.c_str() ] = new MultiSamplePlot(datasets, HTX_Name.c_str() , 10, 0, 1200, "HTX");
    //MSPlot[EventMass_Name.c_str() ] = new MultiSamplePlot(datasets, EventMass_Name.c_str() , 10, 0, 2500, "EventMass");
    //MSPlot[EventMassX_Name.c_str() ] = new MultiSamplePlot(datasets, EventMassX_Name.c_str() , 12, 0, 1700, "EventMassX");
    //MSPlot[FirstTriJetMass_1BJet_g_Name.c_str() ] = new MultiSamplePlot(datasets, FirstTriJetMass_1BJet_g_Name.c_str() , 20, 50 ,300 , "M_{bjj}");
    //MSPlot[FirstDiJetMass_1BJet_g_Name.c_str() ] = new MultiSamplePlot(datasets, FirstDiJetMass_1BJet_g_Name.c_str() , 20, 50 ,300 , "M_{jj}");
    //MSPlot[SecondDiJetMass_1BJet_g_Name.c_str() ] = new MultiSamplePlot(datasets, SecondDiJetMass_1BJet_g_Name.c_str() , 20, 50 ,300 , "M_{jj}");
    //MSPlot[SecondTriJetMass_1BJet_g_Name.c_str() ] = new MultiSamplePlot(datasets, SecondTriJetMass_1BJet_g_Name.c_str() , 20, 50 ,300 , "M_{bjj}");
    //MSPlot[SecondTriJetMass_1BJet_g_chi2_Name.c_str() ] = new MultiSamplePlot(datasets, SecondTriJetMass_1BJet_g_chi2_Name.c_str() , 15, 0 ,250 , "#chi^{2}");
    MSPlot[MET_Name.c_str() ] = new MultiSamplePlot(datasets, MET_Name.c_str() , 40, 0, 700, "MET");
    MSPlot[InclMET_Name.c_str() ] = new MultiSamplePlot(datasets, InclMET_Name.c_str() , 40, 0, 700, "MET");
        
}
}

  ////////////////////////////////////////////////////////////////////
  ////////////////// 1D histograms  //////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  histo1D["lumiWeights"] = new TH1F("lumiWeights","lumiWeights;lumiWeight;#events",100,0,4);
  for (unsigned int d = 0; d < datasets.size(); d++){
    histo2D[("RelIso_vs_MET_"+datasets[d]->Name()).c_str()] = new TH2F(("RelIso_vs_MET_"+datasets[d]->Name()).c_str(),"RelIso:MET",100,0,1000, 100, 0,1);
    histo2D[("ThirdTopMass_vs_HT_"+datasets[d]->Name()).c_str()] = new TH2F(("ThirdTopMass_vs_HT_"+datasets[d]->Name()).c_str(),"ThirdTopMass_vs_HT",30,0,1000, 30, 0,1500);
    histo2D[("ThirdTopMass_vs_SecondTopMass"+datasets[d]->Name()).c_str()] = new TH2F(("ThirdTopMass_vs_SecondTopMass"+datasets[d]->Name()).c_str(),"ThirdTopMass_vs_SecondTopMass",30,0,1000, 30, 0,1000);
      
    histo2D[("MassChi2_vs_HT"+datasets[d]->Name()).c_str()] = new TH2F(("MassChi2_vs_HT"+datasets[d]->Name()).c_str(),"MassChi2_vs_HT",15,0,400, 15, 0,1400);
    histo2D[("EventMassX_vs_HT"+datasets[d]->Name()).c_str()] = new TH2F(("EventMassX_vs_HT"+datasets[d]->Name()).c_str(),"EventMassX_vs_HT",20,0,1500, 15, 0,1200);
    histo2D[("EventMassX_vs_HTX"+datasets[d]->Name()).c_str()] = new TH2F(("EventMassX_vs_HTX"+datasets[d]->Name()).c_str(),"EventMassX_vs_HTX",20,0,1500, 15, 0,1200);

  }
  //  cout << " - Declared histograms ..." <<  endl;
	
  ////////////////////////////////////////////////////////////////////
  ////////////////// Plots  //////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  //Defining a directory in which .png files of all the plots created will be stored.
  string pathPNG = "FourTop"+postfix+channelpostfix;
  pathPNG += "_MSPlots_DiLepton/"; 	
  //  pathPNG = pathPNG +"/"; 	
  mkdir(pathPNG.c_str(),0777);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////// Selection Table/////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //This class is very convenient for monitoring the effect of a series of cuts (a cut-flow)
    // on your data and MC. It writes out a latex (.tex) file with the number of events passing
    //each cut. The cuts indicated here apply to the single lepton + jets channel and should be
    //changed if you want to make a cut flow table.
    
  ////////////////////////////////////////////////////////////////////
  ///////////////////// Channel : µ + jets  //////////////////////////
  ////////////////////////////////////////////////////////////////////
  vector<string> CutsSelecTableMu;
  CutsSelecTableMu.push_back(string("initial"));
  //  CutsSelecTableMu.push_back(string("PU reweighting"));
  CutsSelecTableMu.push_back(string("Event cleaning and Trigger"));
  CutsSelecTableMu.push_back(string("Exactly 1 isolated muon"));
  CutsSelecTableMu.push_back(string("Loose muon veto"));
  CutsSelecTableMu.push_back(string("Electron veto"));
  CutsSelecTableMu.push_back(string("$\\geq$ 1 jet"));
  CutsSelecTableMu.push_back(string("$\\geq$ 2 jet"));
  CutsSelecTableMu.push_back(string("$\\geq$ 3 jet"));
  CutsSelecTableMu.push_back(string("$\\geq$ 4 jet"));
  CutsSelecTableMu.push_back(string("$b\\texttt{-}disc \\geq 0.679$ (CSVM)"));

  SelectionTable selecTableMu(CutsSelecTableMu, datasets);
  selecTableMu.SetLuminosity(Luminosity);
  selecTableMu.SetPrecision(1);


  ////////////////////////////////////////////////////////////////////
  ///////////////////// Channel : e + jets    ////////////////////////
  ////////////////////////////////////////////////////////////////////

  vector<string> CutsSelecTableEl;
  CutsSelecTableEl.push_back(string("initial"));
  //  CutsSelecTableEl.push_back(string("PU reweighting"));
  CutsSelecTableEl.push_back(string("Event cleaning and Trigger"));
  CutsSelecTableEl.push_back(string("Exactly 1 isolated electron"));
  CutsSelecTableEl.push_back(string("Loose muon veto"));
  CutsSelecTableEl.push_back(string("Loose electron veto"));
  CutsSelecTableEl.push_back(string("Conversion veto"));
  CutsSelecTableEl.push_back(string("$\\geq$ 1 jet"));
  CutsSelecTableEl.push_back(string("$\\geq$ 2 jet"));
  CutsSelecTableEl.push_back(string("$\\geq$ 3 jet"));
  CutsSelecTableEl.push_back(string("$\\geq$ 4 jet"));
  CutsSelecTableEl.push_back(string("$b\\texttt{-}disc \\geq 0.679$ (CSVM)"));
  //CutsSelecTableEl.push_back(string("MET > 40 GeV"));


  SelectionTable selecTableEl(CutsSelecTableEl, datasets);
  selecTableEl.SetLuminosity(Luminosity);
  selecTableEl.SetPrecision(1);
  
  //cout << " - SelectionTable instantiated ..." << endl;

       
    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// PileUp Reweighting - N True interactions, recommended method for 2012 //////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
   //PU Reweighting is an event-by-event weighting applied to MC to ensure it has roughly
    //the same  No. of vertices distribution as in data. 
    LumiReWeighting LumiWeights;
    LumiWeights = LumiReWeighting("PUReweighting/pileup_MC_S10.root", "PUReweighting/PURewtoptree_id_2014_1350565897_PileupHistogram.root", "pileup", "pileup");
    
    
  //Flags to shift PU reweighting up and down
  //reweight::PoissonMeanShifter PShiftDown_ = reweight::PoissonMeanShifter(-0.6);
  //reweight::PoissonMeanShifter PShiftUp_ = reweight::PoissonMeanShifter(0.6);

  cout << " Initialized LumiReWeighting stuff" << endl;
  
  //exit(1);
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////// Loop on datasets
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //We now start to loop over our datasets...
    
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
		
      
    // It is possible to apply the Jet Energy Correction here in the macro, however this is not neccessary at the moment
    // as the correction has already been applied in the TopTree, I leave the code here in case we need to apply a new
    // correction in the future.
    vector<JetCorrectorParameters> vCorrParam;

    // Create the JetCorrectorParameter objects, the order does not matter.
    // YYYY is the first part of the txt files: usually the global tag from which they are retrieved
    JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters("../macros/JECFiles/START42_V17_AK5PFchs_L1FastJet.txt");
    JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters("../macros/JECFiles/START42_V17_AK5PFchs_L2Relative.txt");
    JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters("../macros/JECFiles/START42_V17_AK5PFchs_L3Absolute.txt");

    //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
    vCorrParam.push_back(*L1JetPar);
    vCorrParam.push_back(*L2JetPar);
    vCorrParam.push_back(*L3JetPar);

    if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") // Data!
      {
	JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("../macros/JECFiles/START42_V17_AK5PFchs_L2L3Residual.txt");
	vCorrParam.push_back(*ResJetCorPar);
      }
    
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("../macros/JECFiles/START42_V17_AK5PFchs_Uncertainty.txt");
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true); // last boolean ('startFromRaw') = false!    


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////                      Loop on events                                                    ///////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    int itrigger = -1, previousRun = -1;
   
    int start = 0;
    unsigned int end = datasets[d]->NofEvtsToRunOver();

    cout <<"Number of events = "<<  end  <<endl;

      //flag which controls the output of debugging statements
    bool debug = false;

    if (verbose > 1) cout << " - Loop over events " << endl;  
      
      int nBBBar, nCCBar, nLLBar;
      nBBBar=  nCCBar = nLLBar = 0;
      int nSignal, nBackground =0;
      
    for (unsigned int ievt = start; ievt < end; ievt++)
    {  
	

	if(ievt%1000 == 0)
		std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC << " ("<<100*(ievt-start)/(end-start)<<"%)"<<flush<<"\r";

     selecTableMu.Fill(d,0,1.);
     selecTableEl.Fill(d,0,1.);

	//load event
	event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);

        
    //load mcparticles to check jet flavour for ttjets events
    vector<TRootMCParticle*> mcParticles_flav;
    Int_t ttbar_flav = -1;

    double nExB,nExC,nExL;
    nExB = nExC = nExL = 0.;
        
    TRootGenEvent* genEvt_flav = 0;
    genEvt_flav = treeLoader.LoadGenEvent(ievt,false);
        
    //load the MC particles of the TopTree into a vector    
    treeLoader.LoadMCEvent(ievt, genEvt_flav, 0, mcParticles_flav,false); 
        
       // cout <<" mc parts "<< mcParticles_flav.size()  <<endl;

	/*double nLept = 0.;

        if(  (dataSetName != "Data" || dataSetName != "data" || dataSetName != "DATA" ) )
        {
	     //looping over the vector of mc particles and make cuts on the number of leptons coming from W-bosons on truth-level
        for(unsigned int p=0; p<mcParticles_flav.size(); p++) {
          
            if(mcParticles_flav[p]->status()==3 && abs(mcParticles_flav[p]->type())==13 || abs(mcParticles_flav[p]->type())==11 && abs(mcParticles_flav[p]->motherType())==24) {
                nLept++;
	    }

        }
                        if (nLept !=2.)  continue;
       
        }*/
	

        if(  (dataSetName == "TTJets_ll" || dataSetName == "TTJets_cc" || dataSetName == "TTJets_bb" ) )
        {
            //looping over the vector of mc particles and count the number of light quarks (DUSG), C quarks and B quarks
            // which are produced in the decay chain of the top quark.
        for(unsigned int p=0; p<mcParticles_flav.size(); p++) {
          
            if(mcParticles_flav[p]->status()==3 && abs(mcParticles_flav[p]->type())==5 && abs(mcParticles_flav[p]->motherType())!=6) {
               // ttbar_flav=2;
                nExB++;  
            }
            
            else if (mcParticles_flav[p]->status()==3 && abs(mcParticles_flav[p]->type())==4 && abs(mcParticles_flav[p]->motherType())!=6
                && abs(mcParticles_flav[p]->motherType())!=5 && abs(mcParticles_flav[p]->motherType())!=24
                ){
           // ttbar_flav=1;
                 nExC++; 
            }
            
            else if (mcParticles_flav[p]->status()==3 && abs(mcParticles_flav[p]->type())<4 && abs(mcParticles_flav[p]->motherType())!=6){
                // ttbar_flav=1;
                nExL++; 
            }

        }
            
     //   cout <<"TTBar flav composition : " << nExL  <<"  light, " << nExC <<"  C, " << nExB<< "  B" <<  endl;
            
     //   if (ttbar_flav != 1 && ttbar_flav != 2 ) ttbar_flav = 0;
       
            // there are 2 or more b-quarks, the event is ttbar+bb, if there are two or more c quarks, it is ttbar + cc,
            //otherwise, it is ttbar + light quarks.
            
            if (nExB >= 2.){
            ttbar_flav =2; 
            nBBBar++ ; //  bbbar
            }
            else if ( nExC >=2.) {
            ttbar_flav =1; 
            nCCBar++ ; //ccbar
            }
            else{
            ttbar_flav =0.; 
                nLLBar++;  //llbar   
            }
            
            //If the dataset does not match the type of ttbar event, we skip the event
            if (ttbar_flav ==0 && (dataSetName == "TTJets_cc"  || dataSetName == "TTJets_bb"))  continue;
            if (ttbar_flav ==1 && (dataSetName == "TTJets_ll"  || dataSetName == "TTJets_bb" ))  continue;
            if (ttbar_flav ==2 && (dataSetName == "TTJets_ll"  || dataSetName == "TTJets_cc" ))  continue;
        
        }
        
	vector<TRootGenJet*> genjets;
	if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
	{
	  // loading GenJets as I need them for JER
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

    //This is where the triggers are specified, for now Im using the Single lepton (muon) plus jets triggers
    // we may need to use different triggers...
        
    bool trigged = false;
	int currentRun = event->runId();
	if(previousRun != currentRun)
	{
	  cout <<"What run? "<< currentRun<<endl;
      		previousRun = currentRun;
		if(Muon)
		{
			if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
			{
		      
			  if(event->runId() >= 190456 && event->runId() <= 190761)
				 itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2"), currentRun, iFile);
			  else if ( event->runId() >= 190762 && event->runId() <= 191511 )
			         itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v3"), currentRun, iFile);
			  else if ( event->runId() >= 191512  && event->runId() <= 193805  )
			         itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v4"), currentRun, iFile);
			  //for  top triggers change PD at this point.
			  else if ( event->runId() >= 193806  && event->runId() <= 194269  )
			         itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_v1"), currentRun, iFile);
                //extendidng the 52X 1 fb
			  else if ( event->runId() >= 194269 && event->runId() <= 196531  )
			    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1"), currentRun, iFile);
                
              else if ( event->runId() >= 196532 && event->runId() <= 198522  )
                  itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v2"), currentRun, iFile);

              else if ( event->runId() >= 198523 && event->runId() <= 199608  )
                  itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v2"), currentRun, iFile);
                
              else if ( event->runId() >= 199609 && event->runId() <= 201678  )
                  itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet45_35_25_v1"), currentRun, iFile);
                
              else if ( event->runId() >= 201679 && event->runId() <= 202016  )
                  itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet45_35_25_v1"), currentRun, iFile);
			
  		  if(itrigger == 9999)
				{
    		  cout << "NO VALID TRIGGER FOUND FOR THIS EVENT (DATA) IN RUN " << event->runId() << endl;
		    exit(1);
 	 			}
	   		}
	   		else 
	   		{
			  //for refsel, HLT_IsoMu17_eta2p1_TriCentralPFJet30_v2 recommended not HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2
			 // if(dataSetName == "TTJets" || dataSetName == "TTJets_ll" || dataSetName == "TTJets_cc"  || dataSetName == "TTJets_bb" ) itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2"), currentRun, iFile);
                
                  if(dataSetName == "TTJets" || dataSetName == "TTJets_ll" || dataSetName == "TTJets_cc"  || dataSetName == "TTJets_bb" ) itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1"), currentRun, iFile);
                
                else  if(dataSetName == "TTTT") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1"),currentRun, iFile);
                
                 else  if(dataSetName == "T1TTTT") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFJet30_v2"),currentRun, iFile);
                
			  // 	    if(dataSetName == "TTJets") itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v8"), currentRun, iFile);
      			else if (dataSetName == "WJets") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1"), currentRun, iFile);

                else if (dataSetName == "ZJets") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1"), currentRun, iFile);

			    else if (dataSetName == "SingleTop_t") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1"), currentRun, iFile);

       			else if (dataSetName == "SingleTop_tW_T") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1"), currentRun, iFile);
       			else if (dataSetName == "SingleTop_tW_TBar") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1"), currentRun, iFile);
			
			
			////////My added triggers for the new datasets
			else  if(dataSetName == "WWJets") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1"),currentRun, iFile);
			else  if(dataSetName == "ZZJets") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1"),currentRun, iFile);
			else  if(dataSetName == "TTZ") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1"),currentRun, iFile);
			else  if(dataSetName == "WZJets") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1"),currentRun, iFile);
			else  if(dataSetName == "TTW") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1"),currentRun, iFile);

			///////

				else if (dataSetName == "MultiJet") {
                                                                   
	                 	  itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1"), currentRun, iFile); 
				
			
				  }
		    
				       
  				if(itrigger == 9999)
				{
    			  		cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (" << dataSetName << ") IN RUN " << event->runId() << endl;
					//		exit(1);
				}
			}

            if (debug) cout <<"end if muon..."<<endl;
		} //end if Muon
		else if(Electron)
		{
			if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
			{
			
			  if(event->runId() >= 190456 && event->runId() <= 190738)
				 itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v8"), currentRun, iFile);
			  else if ( event->runId() >= 190762 && event->runId() <= 191511 )
			         itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v9"), currentRun, iFile);
		  else if ( event->runId() >= 191512  )
			         itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v10"), currentRun, iFile);
	  if(itrigger == 9999)
				{
    		  cout << "NO VALID TRIGGER FOUND FOR THIS EVENT (DATA) IN RUN " << event->runId() << endl;
		    exit(1);
 	 			}		
			}
	   		else 
			  {
				if(dataSetName == "TTJets") itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v8"), currentRun, iFile);
				else if (dataSetName == "WJets") itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v8"), currentRun, iFile);
    
  				if(itrigger == 9999)
				{
    			  		cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (" << dataSetName << ") IN RUN " << event->runId() << endl;
						exit(1);
				}
			}
		} //end if Electron
	} //end previousRun != currentRun

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////// Jet energy scale corrections     /////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Apply Jet Corrections on-the-fly, as mentioned in an earlier comment, there is no need to apply these corrections
    // at the moment.
        
	//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JES correction on the fly:");
    //	if( dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
    //		jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho(),true); //last boolean: isData (needed for L2L3Residual...)
    //	else
    //		jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho(),false); //last boolean: isData (needed for L2L3Residual...)
	//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JES correction on the fly:");

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////// Type I MET corrections: (Only for |eta| <=4.7 ) //////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before MET type I correction:");      
	//if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
	//	jetTools->correctMETTypeOne(init_jets,mets[0],true);
	//else
	//	jetTools->correctMETTypeOne(init_jets,mets[0],false);
	//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After MET type I correction:");

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////// Jet energy smearing and systematic uncertainty ///////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //As the jet energy resolution is observed to be worse in data than MC, the energy of
        //MC jets is smeared by 10% in MC, this is disabled for now.
        
	if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
	{
	  //JER 
	  doJERShift == 0;
	//	if(doJERShift == 1)
	//		jetTools->correctJetJER(init_jets, genjets, mets[0], "minus");
	//	else if(doJERShift == 2)
	//		jetTools->correctJetJER(init_jets, genjets, mets[0], "plus");
	//	else
	//		jetTools->correctJetJER(init_jets, genjets, mets[0], "nominal");
	  
		//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JER correction:");	       

		// JES systematic! 
	//	if (doJESShift == 1)
	//		jetTools->correctJetJESUnc(init_jets, mets[0], "minus");
	//	else if (doJESShift == 2)
	//jetTools->correctJetJESUnc(init_jets, mets[0], "plus");
	
		//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JES correction:");
	
	}


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////// Beam scrapping veto and PU reweighting ///////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// scale factor for the event
	float scaleFactor = 1.;

        
        //A simple filter to remove events which arise from the interaction LHC beam with
        // the beampipe known as "beam scraping events".
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
	double lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );
	double lumiWeightOLD=lumiWeight;
	if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
	  lumiWeight=1;
	scaleFactor = scaleFactor*lumiWeight;
	
	}

	histo1D["lumiWeights"]->Fill(scaleFactor);	
			
//////////////////////////////////// Event select/////////////////////////////////////////////////////////////
// Here an event selection is applied!
	       
	// Apply trigger selection
	trigged = treeLoader.EventTrigged (itrigger);
	if (debug)cout<<"triggered? Y/N?  "<< trigged  <<endl;
	//if(!trigged)		   continue;  // This checks if trigger is passed or not

	// Declare selection instance, this class is very convenient for applying generic cuts to 
    // objects (jets, muons, electrons, met)
	Selection selection(init_jets, init_muons, init_electrons, mets);

	// Define object selection cuts
	selection.setJetCuts(40.,2.5,0.01,1.,0.98,0,0);//Pt, Eta, EMF, n90Hits, fHPD, dRJetElectron, DRJets

	selection.setElectronCuts(20.,2.5,0.1,0.02,0.,999.,0,1); //Pt,Eta,RelIso,d0,MVAId,DistVzPVz, DRJets, MaxMissingHits
	selection.setLooseElectronCuts(20,2.5,0.2,0.);
        
    selection.setMuonCuts(20.,2.4,.2,0,0.2,0,1,0.5,5,1 ); //Pt,Eta,RelIso,NValidMuHits,d0, dRJets, NMatchedStations,DistVzPVz,
    selection.setLooseMuonCuts(10,2.5,0.2);
	  
	//Select objects 
    //Now we use the selection class to fill vectors of our selected objects.
    vector<TRootElectron*> selectedElectrons       = selection.GetSelectedElectrons(vertex[0]);
	vector<TRootElectron*> selectedExtraElectrons;
	vector<TRootMuon*>     selectedMuons_NoIso = selection.GetSelectedMuons(26,2.4,999.); 
	vector<TRootMuon*>     selectedMuons       = selection.GetSelectedMuons(vertex[0]);
	vector<TRootMuon*>     selectedExtraMuons;
	vector<TRootJet*>      selectedJets        = selection.GetSelectedJets(true); // ApplyJetId
    vector<TRootJet*>      selectedSoftJets        = selection.GetSelectedJets(20.,2.5, selectedMuons, 0., true); // ApplyJetId
	vector<TRootMuon*>     selectedLooseMuons     = selection.GetSelectedLooseMuons();
    vector<TRootElectron*> selectedLooseElectrons = selection.GetSelectedLooseElectrons(); // VBTF ID
    vector<TRootJet*>      selectedBJets; // B-Jets
    vector<TRootJet*>      selectedLightJets; // light-Jets

	//order jets wrt to Pt, then set bool corresponding to RefSel cuts.                                                                                 
    sort(selectedJets.begin(),selectedJets.end(),HighestPt()); //order muons wrt Pt.                                                                    

    int JetCut =0;
    int nMu = selectedMuons.size();
    int nEl = selectedElectrons.size();


	// Apply primary vertex selection (requiring a good PV in the event.)
	bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2);
        if(!isGoodPV) continue;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////// Sync'ing cutflow//////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////     
                
        //Writing out the cut-flow mentioned earlier, need to edit the cuts here if we want to 
        //use this...
        
	selecTableMu.Fill(d,1,1.) ;
	selecTableEl.Fill(d,1,1.) ;
	int nTags = 0;
        
        
	  if (nMu ==1){selecTableMu.Fill(d,2,1.) ;
	    if (selectedLooseMuons.size()==1){  selecTableMu.Fill(d,3,1.);
	         if (selectedLooseElectrons.size()==0)  {selecTableMu.Fill(d,4,1.);     
		   if (selectedJets.size() >= 1 && selectedJets[0]->Pt() >45.) {  selecTableMu.Fill(d,5,1.) ;
		     if (selectedJets.size() >= 2 && selectedJets[1]->Pt() >45.) { selecTableMu.Fill(d,6,1.) ;
		       if (selectedJets.size() >= 3 && selectedJets[2]->Pt() >45.){  selecTableMu.Fill(d,7,1.) ;
			 if (selectedJets.size() >= 4 && selectedJets[3]->Pt() >20.) { selecTableMu.Fill(d,8,1.) ;

			   for (int testjet=0; testjet<selectedJets.size(); testjet++) {
			     if (selectedJets[testjet]->btag_combinedSecondaryVertexBJetTags() > 0.679)  nTags++;
			   }
        if (nTags >= 1) selecTableMu.Fill(d,9,1.);
	}
		}
		}
	      }
	    }
	  }
	  }

	  //Sync'ing cutflow -electrons
	   nTags = 0;

	   if (nEl ==1){selecTableEl.Fill(d,2,1.) ;//one isolated electron
	     if (selectedLooseMuons.size()==0){  selecTableEl.Fill(d,3,1.);// loose muon veto
	       if (selectedLooseElectrons.size()==1)  {selecTableEl.Fill(d,4,1.);    // di-lepton veto 
		if( selectedElectrons[0]->passConversion() == true  ){selecTableEl.Fill(d,5,1.); //put conversion rejection cut here
		  if (selectedJets.size() >= 1 && selectedJets[0]->Pt() >45.) {  selecTableEl.Fill(d,6,1.) ; //1 jet 45
		    if (selectedJets.size() >= 2 && selectedJets[1]->Pt() >45.) { selecTableEl.Fill(d,7,1.) ; //2 jet 45
		      if (selectedJets.size() >= 3 && selectedJets[2]->Pt() >45.){  selecTableEl.Fill(d,8,1.) ; //3 jet 45
			if (selectedJets.size() >= 4 && selectedJets[3]->Pt() >20.) { selecTableEl.Fill(d,9,1.) ; //4 jet 20
		      for (int testjet=0; testjet<selectedJets.size(); testjet++) {
			if (selectedJets[testjet]->btag_combinedSecondaryVertexBJetTags() > 0.679)  nTags++;
        }
                      if (nTags >= 1) selecTableEl.Fill(d,10,1.);  //1 CSVM Btag
	}
	}
	}
		}
	      }
	    }
	  }
	  }
        
        
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Applying baseline offline event selection: finally make requirents and skip events which do not fulfill the requirements. ////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (debug)	cout <<" applying baseline event selection..."<<endl;
                
	//if (debug) cout<<" jet1 pt =  "<<selectedJets[0]->Pt() << "   "<< " jet2 pt =  "<<selectedJets[1]->Pt() << "   "<< " jet2 pt =  "<<selectedJets[2]->Pt()   <<endl;
 
        //filling vector of b-jets
	for (Int_t seljet =0; seljet < selectedJets.size(); seljet++ ){
	  if( selectedJets[seljet]->btag_combinedSecondaryVertexBJetTags() > workingpointvalue) selectedBJets.push_back(selectedJets[seljet]);
	  else selectedLightJets.push_back(selectedJets[seljet]);
	}

	if(Muon){

    if (debug) cout <<"Number of Muons, Jets, BJets, JetCut  ===>  "<< selectedMuons.size() <<"  "  << selectedJets.size()   <<"  " <<  selectedBJets.size()   <<endl;

	  //Apply the selection, requring 3 jets, 1 btag, 1 muon, 1 electon. Leptons are required to be same-sign.
	  if  (  !( selectedJets.size() >= 4 &&    selectedBJets.size() >= 1 &&  nMu == 1 &&  nEl == 1 && (selectedMuons[0]->charge() != selectedElectrons[0]->charge()  )  )) continue; 

	  if (debug) cout<< "Event passed..."<<endl;
        
        
////////////////////////////////////////////////////////////////////////////////////
//// Getting Gen Event
////////////////////////////////////////////////////////////////////////////////////
         if(dataSetName != "data" && dataSetName != "Data" && dataSetName != "Data"){
        vector<TRootMCParticle*> mcParticles;
        vector<TRootMCParticle*> mcTops;
        
        TRootGenEvent* genEvt = 0;
        genEvt = treeLoader.LoadGenEvent(ievt,false);
        sort(selectedJets.begin(),selectedJets.end(),HighestPt()); 
        treeLoader.LoadMCEvent(ievt, genEvt, 0, mcParticles,false);  
        
        int leptonPDG = 13;
        if (debug) cout <<"size   "<< mcParticles.size()<<endl;
    }
    
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////Filling histograms / plotting
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Order of plotting
	  // 0. Vertices
	  // 1. MET
	  // 2. Muons
	  // 3. Jets: per jet plots, event level variables, jet combinations,discriminants.


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////Vertices
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	MSPlot["NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////MET
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	MSPlot["MET"]->Fill(mets[0]->Et(), datasets[d], true, Luminosity*scaleFactor);	

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////Muons
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (debug) cout <<"  nMUons  "<< selectedMuons.size()<< "  n loose muons " << selectedLooseMuons.size()    <<endl;
    //loop to get inv. masses of muon-loosemuon combinations
    for (Int_t mu1 = 0; mu1 < selectedMuons.size(); mu1++){	 
	  for (Int_t mu2 = 0; mu2 < selectedMuons.size(); mu2++){
          
          if (debug) cout <<"making dimuosss   "<< endl;
          
     	    TRootMuon DiMu (TLorentzVector ( selectedMuons[mu1]->Px() + selectedMuons[mu2]->Px() ,selectedMuons[mu1]->Py() + selectedMuons[mu2]->Py()    ,selectedMuons[mu1]->Pz() + selectedMuons[mu2]->Pz() , selectedMuons[mu1]->E() + selectedMuons[mu2]->E() ));
            MSPlot["DiMuon_InvMass"]->Fill( DiMu.M(), datasets[d], true, Luminosity*scaleFactor );  
	  }
	  }
        
        if (debug) cout <<"made dimuons...   "<< endl;
        
	  MSPlot["NbOfIsolatedMuons"]->Fill(selectedMuons.size(), datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["NbOfIsolatedElectrons"]->Fill(selectedElectrons.size(), datasets[d], true, Luminosity*scaleFactor);

        if (debug) cout <<"filling muid   "<< selectedMuons[0]->Pt()<<endl;
        
	  //Fill Muon ID plots, filling these plots is causing a crash.. investoigating...
  
      /* 
      double muzPVz = selectedMuons[0]->vz() - vertex[0]->Z();
      double Mtrans =  (*selectedMuons[0] + *mets[0]).Mt();
	  MSPlot["MuonRelIsolation"]->Fill(selectedMuons[0]->relativePfIso03(), datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["MuonPt"]->Fill(selectedMuons[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["MuonEta"]->Fill(selectedMuons[0]->Eta(), datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["MuonPhi"]->Fill(selectedMuons[0]->Phi(), datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["MuonNValidHits"]->Fill(selectedMuons[0]->nofValidHits(), datasets[d], true, Luminosity*scaleFactor);	 
	  MSPlot["Muond0"]->Fill(selectedMuons[0]->d0(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["MuonDistVzPVz"]->Fill(muzPVz, datasets[d], true, Luminosity*scaleFactor );
      MSPlot["MuonNMatchedStations"]->Fill(selectedMuons[0]->nofMatchedStations(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["MuonTrackerLayersWithMeasurement"]->Fill(selectedMuons[0]->nofTrackerLayersWithMeasurement(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["WMt"]->Fill(Mtrans,datasets[d], true, Luminosity*scaleFactor);
      histo2D[("RelIso_vs_MET_"+datasets[d]->Name()).c_str()]->Fill(mets[0]->Et(),selectedMuons_NoIso[0]->relativePfIso03(), Luminosity*scaleFactor );
*/
         
      if (debug) cout <<"filled all muID plots  .."<<endl;

   
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////Jets
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  MSPlot["NbOfSelectedJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["NbOfSelectedLightJets"]->Fill(selectedLightJets.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["NbOfSelectedBJets"]->Fill(selectedBJets.size(), datasets[d], true, Luminosity*scaleFactor);
	//  MSPlot["RhoCorrection"]->Fill(event->kt6PFJetsPF2PAT_rho(), datasets[d], true, Luminosity*scaleFactor);
	  if (debug) cout <<"per jet plots.."<<endl;

   	  double ljetpt;
	  double bjetpt;
      double jetpt;
      double HT =0;
	  double MHT =0;
	  double sumpx =0, sumpy=0, sumpz=0, sume=0; 

	//plots to to inspire staggered Jet Pt selection
      if (selectedJets.size()>=4) MSPlot["4thJetPt"]->Fill(selectedJets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      if (selectedJets.size()>=5) MSPlot["5thJetPt"]->Fill(selectedJets[4]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      if (selectedJets.size()>=6) MSPlot["6thJetPt"]->Fill(selectedJets[5]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      if (selectedJets.size()>=7) MSPlot["7thJetPt"]->Fill(selectedJets[6]->Pt(), datasets[d], true, Luminosity*scaleFactor);
	
	for (Int_t seljet1 =0; seljet1 < selectedLightJets.size(); seljet1++ ){
		  MSPlot["SelectedJetPt_light"]->Fill(  selectedLightJets[seljet1]->Pt()  , datasets[d], true, Luminosity*scaleFactor);
		  MSPlot["BdiscBJetCand_CSV"]->Fill(selectedJets[seljet1]->btag_combinedSecondaryVertexBJetTags(),datasets[d], true, Luminosity*scaleFactor);
		  MSPlot["SelectedJetPt"]->Fill(jetpt, datasets[d], true, Luminosity*scaleFactor);
		  MSPlot["JetEta"]->Fill(selectedJets[seljet1]->Eta() , datasets[d], true, Luminosity*scaleFactor);
		  MSPlot["JetPhi"]->Fill(selectedJets[seljet1]->Phi() , datasets[d], true, Luminosity*scaleFactor);

		  //Event-level variables
          jetpt = selectedJets[seljet1]->Pt();
		  HT = HT + jetpt;
		  sumpx = sumpx + selectedJets[seljet1]->Px();
		  sumpy = sumpy + selectedJets[seljet1]->Py();
		  sumpz = sumpz + selectedJets[seljet1]->Pz();
		  sume = sume + selectedJets[seljet1]->E();
	}

	    TRootJet sumjet (TLorentzVector (sumpx, sumpy, sumpz,sume ));
        MHT =  sumjet.Pt();
     	double MHTSig = MHT/sqrt(HT);
	    double STJet = HT + selectedMuons[0]->Pt(); //need to add met back in when available
        double EventMass = sumjet.M();
        
        MSPlot["MHT_SelectedJets"]->Fill(MHT, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["HT_SelectedJets"]->Fill(HT, datasets[d], true, Luminosity*scaleFactor);
	    MSPlot["MHTSig_SelectedJets"]->Fill(MHTSig, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["STJet"]->Fill( STJet , datasets[d], true, Luminosity*scaleFactor);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////Jet combinations, Mass Reconstructions...
	//// Logic: 1) Reconstruct best hadronic top candidate (simulataneously minimize Mass_{Topcandidate} - Mass_{Top} && Mass_{Wcandidate} - Mass_{W})
	////        2) Reconstruct best leptonic top candidate without considering b-jet used in hadronic candidate.
	////        3) Reconstruct best hadronic top candidate as in 1) but only use jets not already used.
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if (debug) cout <<"jet comb."<<endl;
	int firstDiMass = 0;
	int secondDiMass = 0;

	int index_j1_fir, index_j2_fir, index_b_fir, index_j1_sec, index_j2_sec, index_b_sec,index_j1_thi, index_j2_thi, index_b_thi ;
	double WMass_had = 87.04;
	double TopMass_had = 180.68;
    double WMass_sigma_had = 11.19;
    double TopMass_sigma_had = 17.3;
    
        
    double chi2_g = 1000000;
	double candchi2 ;
	double firstchi2 =1000;
	double secondchi2 =1000;
	double lepchi2 =1000;
    double candchi2_g ;
    double firstchi2_g =1000;
    double secondchi2_g =1000;
    double mass_diff_chi = 1000;
        
    double fir_top_Pt_chosen_g;
    double fir_top_mass_chosen_g;
    double fir_w_mass_chosen_g;
    double sec_top_Pt_chosen_g;
    double sec_top_mass_chosen_g;
    double sec_w_mass_chosen_g;
    double lep_top_Pt_chosen_g;
    double lep_top_mass_chosen_g;
    double lep_w_mass_chosen_g;
    double  global_chi2 = 999.;
        
        //////////////////////////////////////////////////////
        // Had top candidate
        //////////////////////////////////////////////////////
        
        for (Int_t jet1_g =0; jet1_g < selectedLightJets.size(); jet1_g++ ){
            
            if (debug) cout <<"making global candidate..."<<endl;
        for (Int_t jet2_g =0; jet2_g < selectedLightJets.size(); jet2_g++ ){
            if (jet1_g == jet2_g) continue;
                    for (Int_t bjet1_g =0; bjet1_g < selectedBJets.size(); bjet1_g++ ){
            
                            TRootJet* j1_g = (TRootJet*) selectedLightJets[jet1_g];
                            TRootJet* j2_g = (TRootJet*) selectedLightJets[jet2_g];
                            TRootJet* bj1_g = (TRootJet*) selectedBJets[bjet1_g];
                            
                            double fir_top_mass_g  = (*j1_g + *j2_g + *bj1_g).M(); 
                            double fir_w_mass_g  = (*j1_g + *j2_g).M(); 
                            double fir_top_Pt_g= (*j1_g + *j2_g + *bj1_g).Pt(); 
                            double fir_top_HTRat_g = (HT)/((*j1_g).Pt() + (*j2_g).Pt() + (*bj1_g).Pt()  );
                            
                            double lep_top_mass_g = 0;
                            double lep_w_mass_g = 0;
                            double lep_top_Pt_chosen_g = 0;
                            
                            
                            //delr variable, not used at the moment
                            double fir_top_delr_g =  sqrt(   pow(j1_g->Eta() - j2_g->Eta(),2) + pow(j1_g->Phi() - j2_g->Phi(),2 )   )  +  sqrt(   pow(j2_g->Eta() - bj1_g->Eta(),2) + pow(j2_g->Phi() - bj1_g->Phi(),2 )   )  + sqrt(   pow(bj1_g->Eta() - j1_g->Eta(),2) + pow(bj1_g->Phi() - j1_g->Phi(),2 )) ; 
                            
                                                        
                            //chi2 for hadronic top
                            chi2_g = (pow(fir_w_mass_g  - WMass_had,2)/pow(WMass_sigma_had,2)) + (pow(fir_top_mass_g - TopMass_had,2)/pow(TopMass_sigma_had,2))  ;
                            
                            
                            double chi2_extended_g = chi2_g ;
                            
                            if ( chi2_extended_g < firstchi2_g   ){
                             
                                firstchi2_g = chi2_extended_g;
                                fir_top_mass_chosen_g = fir_top_mass_g;
                                fir_top_Pt_chosen_g =  fir_top_Pt_g;
                                fir_w_mass_chosen_g = fir_w_mass_g;
                                                             
                                //set indices of best comb
                                index_j1_fir = jet1_g;
                                index_j2_fir = jet2_g;
                             
                            }
                        }
            }
        }
                
        if (debug) cout <<" jet combs done"<<endl;
        
    
        //1st hadronic global cand top/W
        MSPlot["FirstDiJetMass_g"]->Fill(fir_w_mass_chosen_g,  datasets[d], true, Luminosity*scaleFactor );
        MSPlot["FirstTriJetMass_1BJet_g"]->Fill(fir_top_mass_chosen_g,  datasets[d], true, Luminosity*scaleFactor );
        MSPlot["FirstTriJetPt_1BJet_g"]->Fill(fir_top_Pt_chosen_g,  datasets[d], true, Luminosity*scaleFactor );
         
        if (debug) cout <<" second had top cand done"<<endl;


	if((dataSetName.find("Data")<=0 || dataSetName.find("data")<=0 || dataSetName.find("DATA")<=0)&& ( selectedBJets.size() > 3 )){
	  cout<<"  "<<endl;
          cout <<"Run: " <<event->runId() << "  Lumi block: " << event->lumiBlockId() <<"  Event: " << event->eventId()   <<endl;
	  cout <<"Interesting event!  1st  DiJet Mass = "<< fir_w_mass_chosen_g << "  1st trijet mass =  "  <<  fir_top_mass_chosen_g  <<endl;
	  cout<<" N Jets: " << selectedJets.size()  <<"  N Tags: " << selectedBJets.size() <<endl;
	  cout<<"  "<<endl;
	}
        
        
     //    cout <<"NJets = " << selectedJets.size()  <<"  NBJets = " << selectedBJets.size() << "    Jet com indices : "<< " j1_fir = "  <<  index_j1_fir << "  j2_fir = "    <<  index_j2_fir   <<" b_fir = "    << index_b_fir   <<" b_lep  = "  << index_b_lep << " j1_sec  = "<< index_j1_sec << "j2_sec = "<< index_j2_sec << "b_sec = "<< index_b_sec << "j1_thi =   "  << index_j1_thi <<" j2_thi = " << index_j2_thi <<"b_thi =  " << index_b_thi <<endl;
        
        
//////////////////////////////////////////////////////////////////////////////////
/////          Calculate derived variables for rejecting ttbar + X          //////
/////                                                                       //////
///// 1. HT, 2. HTX, 3.EventMass, 4.EventMassX                              //////
//////////////////////////////////////////////////////////////////////////////////
        
    double HTX = 0, sumpx_X = 0, sumpy_X= 0, sumpz_X =0, sume_X= 0;
    //Calculate HTX: first caculate HT_light
    for (Int_t seljet1 =0; seljet1 < selectedLightJets.size(); seljet1++ ){
            
        //disregard any previously used light jets
        if ( seljet1 == index_j1_fir ||  seljet1 == index_j2_fir  ) continue;
            //Event-level variables
            double ljetpt = selectedLightJets[seljet1]->Pt();
            HTX = HTX + ljetpt;
            sumpx_X = sumpx_X + selectedLightJets[seljet1]->Px();
            sumpy_X = sumpy_X + selectedLightJets[seljet1]->Py();
            sumpz_X = sumpz_X + selectedLightJets[seljet1]->Pz();
            sume_X = sume_X + selectedLightJets[seljet1]->E();
        }
        
        
    //Calculate HTX: second caculate HT_b
    for (Int_t seljet1 =0; seljet1 < selectedBJets.size(); seljet1++ ){
            
        //disregard any previously used b jets
        if ( seljet1 == index_b_fir  ) continue;
            
            double  bjetpt = selectedBJets[seljet1]->Pt();
            HTX = HTX + bjetpt;
            sumpx_X = sumpx_X + selectedBJets[seljet1]->Px();
            sumpy_X = sumpy_X + selectedBJets[seljet1]->Py();
            sumpz_X = sumpz_X + selectedBJets[seljet1]->Pz();
            sume_X =  sume_X  + selectedBJets[seljet1]->E();
        
        }
        
        
        TRootJet sumjet_X (TLorentzVector (sumpx_X, sumpy_X, sumpz_X,sume_X ));        
        double EventMassX = sumjet_X.M();
        
        
        if (debug) cout <<"arrays  "<<endl;

//////////////////////////////////////////////////////////////////////////////////
/////          Filling NJet Vs NBJet arrays for                             //////
/////                                                                       //////
///// 1. HT, 2. STLep, 3.STJet, 4.HT(w/o ttbar), 5.InvMass (w/o) ttbar      //////
//////////////////////////////////////////////////////////////////////////////////

  for (Int_t b = minNJets; b <= maxNJets; b++){
      for (Int_t c = minNBJets; c<= maxNBJets; c++){
          string NJets_str = static_cast<ostringstream*>( &(ostringstream() << b) )->str();
          string NBJets_str = static_cast<ostringstream*>( &(ostringstream() << c) )->str();
          //string HT_Name = "HT_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
	  string MET_Name = "MET_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
	  string inclMET_Name = "InclMET_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;	  
          //string HTX_Name = "HTX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          //string EventMass_Name = "EventMass_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          //string EventMassX_Name = "EventMassX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          //string FirstTriJetMass_1BJet_g_Name = "FirstTriJetMass_1BJet_g"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          //string FirstDiJetMass_1BJet_g_Name = "FirstDiJetMass_1BJet_g"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          //string SecondDiJetMass_1BJet_g_Name = "SecondDiJetMass_1BJet_g"+NJets_str+"Jets_"+NBJets_str+"Tags" ;

          //string MuMetBMasses_g_Name = "MuMetBMasses_g_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          //string MuMetMasses_g_Name = "MuMetMasses_g_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;

          //string SecondTriJetMass_1BJet_g_Name = "SecondTriJetMass_1BJet_g"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          //string SecondTriJetMass_1BJet_g_chi2_Name = "SecondTriJetMass_1BJet_g_chi2_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          
          
          
        if(b<4 && c<2){
            if(selectedJets.size() == b && selectedBJets.size() == c  ) {
             //   MSPlot[HT_Name.c_str() ]->Fill(HT,datasets[d], true, Luminosity*scaleFactor);
		MSPlot[MET_Name.c_str() ]->Fill(mets[0]->Et(),datasets[d], true, Luminosity*scaleFactor);
               // MSPlot[HTX_Name.c_str() ]->Fill(HTX,datasets[d], true, Luminosity*scaleFactor);
               // MSPlot[EventMass_Name.c_str() ]->Fill(EventMass,datasets[d], true, Luminosity*scaleFactor);
               // MSPlot[EventMassX_Name.c_str() ]->Fill(EventMassX,datasets[d], true, Luminosity*scaleFactor);
               // MSPlot[FirstTriJetMass_1BJet_g_Name.c_str()]->Fill(fir_top_mass_chosen_g,  datasets[d], true, Luminosity*scaleFactor);
                //MSPlot[FirstDiJetMass_1BJet_g_Name.c_str()]->Fill(fir_w_mass_chosen_g,  datasets[d], true, Luminosity*scaleFactor);
                //MSPlot[SecondDiJetMass_1BJet_g_Name.c_str()]->Fill(sec_w_mass_chosen_g,  datasets[d], true, Luminosity*scaleFactor);

              //  MSPlot[SecondTriJetMass_1BJet_g_Name.c_str()]->Fill(sec_top_mass_chosen,  datasets[d], true, Luminosity*scaleFactor );
               // MSPlot[SecondTriJetMass_1BJet_g_chi2_Name.c_str()]->Fill(firstchi2_g ,  datasets[d], true, Luminosity*scaleFactor );
                
            

                      }
                        }
        else{
            if(selectedJets.size() >= b && selectedBJets.size() >= c  ) {
                //MSPlot[HT_Name.c_str() ]->Fill(HT,datasets[d], true, Luminosity*scaleFactor);
		MSPlot[MET_Name.c_str() ]->Fill(mets[0]->Et(),datasets[d], true, Luminosity*scaleFactor);
                //MSPlot[HTX_Name.c_str() ]->Fill(HTX,datasets[d], true, Luminosity*scaleFactor);
                //MSPlot[EventMass_Name.c_str() ]->Fill(EventMass,datasets[d], true, Luminosity*scaleFactor);
                //MSPlot[EventMassX_Name.c_str() ]->Fill(EventMassX,datasets[d], true, Luminosity*scaleFactor);
                //MSPlot[FirstTriJetMass_1BJet_g_Name.c_str()]->Fill(fir_top_mass_chosen_g,  datasets[d], true, Luminosity*scaleFactor);
                //MSPlot[FirstDiJetMass_1BJet_g_Name.c_str()]->Fill(fir_w_mass_chosen_g,  datasets[d], true, Luminosity*scaleFactor);
               
                //MSPlot[SecondTriJetMass_1BJet_g_chi2_Name.c_str()]->Fill(firstchi2_g ,  datasets[d], true, Luminosity*scaleFactor );
           


            }
        }
}
}
 
    
        histo2D[("MassChi2_vs_HT"+datasets[d]->Name()).c_str()]->Fill(firstchi2,HT );
        histo2D[("EventMassX_vs_HT"+datasets[d]->Name()).c_str()]->Fill(EventMassX, HT);
        histo2D[("EventMassX_vs_HTX"+datasets[d]->Name()).c_str()]->Fill(EventMassX, HTX);
        
	}
	else if(Electron){
	  //	MSPlot["1stLeadingElectronRelIsolation"]->Fill(selectedElectrons_NoIso[0]->relativePfIso(), datasets[d], true, Luminosity*scaleFactor);

	}
        
    }//loop on events
    
    if (debug)cout <<"N BBar = = " << nBBBar <<"  N CCBar = = " << nCCBar <<"  N LLBar = =  " << nLLBar << endl;

      
    if(jetTools) delete jetTools;
    
    
      //important: free memory
      treeLoader.UnLoadDataset();
      
      
  } //loop on datasets
    
  //Once everything is filled ...
  cout << " We ran over all the data ;-)" << endl;
  
  ///////////////////
  // Writing
  //////////////////
  cout << " - Writing outputs to the files ..." << endl;

  //Selection tables
  if(Muon){ 
	//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
	//selecTableMu.TableCalculator(  false, true, true, true, true);


  //Options : WithError (false), writeMerged (true), useBookTabs (false), addRawNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false)
  //	selecTableMu.Write(  "FourTop"+postfix+"Table_Mu.tex",    true,true,true,true,false,false,true);
  }
    else if(Electron){
	//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
	selecTableEl.TableCalculator(  false, true, true, true, true);
   //Options : WithError (false), writeMerged (true), useBookTabs (false), addRawNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false)
    selecTableEl.Write(  "FourTop"+postfix+"Table_El.tex",  true,true,true,true,false,false,true);

  }
 cout <<" fouting..."<<endl;
 
  fout->cd();
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
        
      
	MultiSamplePlot *temp = it->second;
	TH1F *tempHisto_data;
	TH1F *tempHisto_TTTT;
	//	temp->addText("CMS preliminary");
	string name = it->first;
	temp->Draw(false, name, true, true, true, true, true,1,false); // merge TT/QCD/W/Z/ST/
	//Draw(bool addRandomPseudoData = false, string label = string("CMSPlot"), bool mergeTT = false, bool mergeQCD = false, bool mergeW = false, bool mergeZ = false, bool mergeST = false, int scaleNPSignal = 1, bool addRatio = false, bool mergeVV = false, bool mergeTTV = false);
      
      cout <<" looping plots..., name ... "<< name<<endl;
        
        temp->Write(fout, name, true, pathPNG, "pdf");
        cout <<" written...."<<endl;

  }

    cout <<" fouting.   2.."<<endl;

  cout <<"1D  "<< histo1D.size()  <<"2D   "  <<  histo2D.size() <<endl;

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

