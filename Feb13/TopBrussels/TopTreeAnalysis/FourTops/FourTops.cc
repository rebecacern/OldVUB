 // rebeca@cern.ch
//
// Fall11 round

#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "../Selection/interface/SelectionTable.h"
#include "../Tools/interface/PlottingTools.h"
#include "../Tools/interface/JetTools.h"
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
#include "../Tools/interface/TopologyWorker.h"
#include "../MCInformation/interface/ResolutionFit.h"
#include "../Reconstruction/interface/JetCorrectorParameters.h"
#include "../Reconstruction/interface/JetCorrectionUncertainty.h"
#include "../MCInformation/interface/Lumi3DReWeighting.h"
#include "../macros/Style.C"

using namespace std;
using namespace TopTree;


int main(int argc, char* argv[]) {
  
  setMyStyle();
  
  clock_t start = clock();
  
  //modes: 0 emu, 1mumu, 2ee 
  int  mode = 0; 
  bool isRAW = false;
  
  //Pile-up reweighting
  bool reweightPU = true;
  
  //b-tag scale factor
  bool scaleFactor = true;
  double SFval = 0.95; 
 
 // Met Type I correction	
  bool metTypeI = false; 
 
  //Systematic calculations 
  bool JESPlus= false;
  bool JESMinus= false;

  bool JERPlus= false;
  bool JERMinus= false;

  bool SFplus = false;
  bool SFminus = false;
  
  bool unclusteredUp = false;
  bool unclusteredDown = false;
  
  bool PUsysUp = false;
  bool PUsysDown = false;
  
  // Special Samples
  bool SystSamples = false;
  bool Spring11 = false;
  bool Special = false;
  
  string xmlfile ="twemu.xml";
  
  // Arguments
  for(int iarg = 0; iarg < argc && argc>1 ; iarg++){
    std::string argval=argv[iarg];
    if(argval=="--help" || argval =="--h"){
      cout << "--ee: Di-Electron" << endl;
      cout << "--emu: Electron-Muon" << endl;
      cout << "--mumu: Di-Muon" << endl;
      cout << "--JESplus: JES sys +1 sigma MET included" << endl;
      cout << "--JESminus: JES sys -1 sigma MET included" << endl;
      cout << "--JERplus: JER +" << endl;
      cout << "--JERminus: JER -" << endl;
      cout << "--SFplus: SF up +10% syst" << endl;
      cout << "--SFminus: SF down -10% syst" << endl;
      cout << "--PUup: PU reweghting scaled up " << endl;
      cout << "--PUdown: PU reweghting scaled down " << endl;
      cout << "--uncMETup: Unclustered MET syst. Up " << endl;
      cout << "--uncMETdown: Unclustered MET syst. Down " << endl;
      cout << "--NoPU: Do not apply pileup re-weighting" << endl;
      cout << "--NoSF: Do not apply b-tag scale factor" << endl;
      cout << "--RAW: Do not apply pileup re-weighting or b-tag scale factor" << endl;
      return 0;
    }
    if (argval=="--ee") mode = 2;
    if (argval=="--emu") mode = 0;
    if (argval=="--mumu") mode = 1;
    if (argval=="--uncMETup") unclusteredUp = true;
    if (argval=="--uncMETdown") unclusteredDown = true;
    if (argval=="--PUup" ){PUsysUp = true;}
    if (argval=="--PUdown" ){PUsysDown = true;}
    if (argval=="--JESplus") {JESPlus = true;}
    if (argval=="--JESminus") {JESMinus = true;}
    if (argval=="--JERplus") {JERPlus = true;}
    if (argval=="--JERminus") {JERMinus = true;}
    if (argval=="--SFplus") SFplus = true;
    if (argval=="--SFminus") SFminus = true;
    if (argval=="--NoPU") reweightPU = false;
    if (argval=="--NoSF") scaleFactor = false;
    if (argval=="--RAW") {reweightPU = false; scaleFactor = false; isRAW = true;}
  }   
  
  // Luminosity and xml files
  double lumi = 0;
  if      (mode == 0){ 	 lumi = 4626.297;	 xmlfile ="twemu.xml";}
  else if (mode == 1){	 lumi = 4534.871;	 xmlfile = "twmumu.xml";}
  else if (mode == 2){	 lumi = 4593.348;	 xmlfile = "twee.xml";}
  
  
  // Analysis environment
  TTree *configTree = new TTree("configTree","configuration Tree");
  TClonesArray* tcdatasets = new TClonesArray("Dataset",1000);
  configTree->Branch("Datasets","TClonesArray",&tcdatasets);
  TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000);
  configTree->Branch("AnaEnv","TClonesArray",&tcAnaEnv);
  
  AnalysisEnvironment anaEnv;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile.c_str());
  new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);
  
  TTreeLoader treeLoader;
  vector < Dataset* > datasets;
  treeLoader.LoadDatasets (datasets, xmlfile.c_str());
  
  // Start of the analysis
  for (unsigned int d = 0; d < datasets.size (); d++)
    {
      treeLoader.LoadDataset (datasets[d], anaEnv);
      string dataSetName = datasets[d]->Name();
      
      bool isData = false;
      char name[100];
      double xlweight;
      
      // cross sections and weights
      if (dataSetName == "data"){		sprintf(name, "data");  	xlweight = 1; isData = true;}
      else if (dataSetName == "tt"){            sprintf(name, "tt");            xlweight = lumi*163/3160707;} //
      else if (dataSetName == "tt2l"){		sprintf(name, "tt2l");  	xlweight = lumi*17.10/8576584;} //
      else if (dataSetName == "twdr"){      	sprintf(name, "tw_dr");		xlweight = lumi*7.87/813743;} //
      else if (dataSetName == "atwdr"){ 	sprintf(name, "atw_dr");        xlweight = lumi*7.87/689462;} //
      else if (dataSetName == "twds"){ 	     	sprintf(name, "tw_ds");	        xlweight = lumi*7.87/794802;} //
      else if (dataSetName == "atwds"){    	sprintf(name, "atw_ds");        xlweight = lumi*7.87/784764;} //
      else if (dataSetName == "t"){    		sprintf(name, "t");        	xlweight = lumi*41.92/3337875;} //
      else if (dataSetName == "at"){    	sprintf(name, "at");        	xlweight = lumi*22.65/1943627;} //
      else if (dataSetName == "ts"){    	sprintf(name, "ts");        	xlweight = lumi*3.19/259777;} //
      else if (dataSetName == "ats"){    	sprintf(name, "ats");        	xlweight = lumi*1.44/137889;} //
      else if (dataSetName == "wjets"){     	sprintf(name, "wjets");	        xlweight = lumi*31314/49708092;} //
      else if (dataSetName == "zjets"){         sprintf(name, "zjets");         xlweight = lumi*3048/26523984;} //
 
      else if (dataSetName == "dymumu"){	sprintf(name, "dymumu");	xlweight = lumi*1666/1;}
      else if (dataSetName == "dyee"){		sprintf(name, "dyee");		xlweight = lumi*1666/1;}
      else if (dataSetName == "dytautau"){	sprintf(name, "dytautau");	xlweight = lumi*1666/1;}
      else if (dataSetName == "ww"){		sprintf(name, "ww");		xlweight = lumi*42.9/4223785;} //
      else if (dataSetName == "wz"){		sprintf(name, "wz");		xlweight = lumi*18.3/3863081;} //
      else if (dataSetName == "zz"){		sprintf(name, "zz");		xlweight = lumi*7.67/4188624;} //
      else if (dataSetName == "qcd_mu"){    	sprintf(name, "qcd_mu");	xlweight = lumi*84679.3/25079892; Special = true;} //OLD 
      
      //special files
      else if (dataSetName == "t_sup"){      		sprintf(name, "t_sup");		        xlweight = lumi*7.87/437736;}
      else if (dataSetName == "tbar_sup"){      	sprintf(name, "tbar_sup");		xlweight = lumi*7.87/437794;}
      else if (dataSetName == "tbar_sdo"){      	sprintf(name, "tbar_sdo");		xlweight = lumi*7.87/437863;}
      
      else if (dataSetName == "tt_largeISR"){   	sprintf(name, "tt_largeISR");   	xlweight = lumi*163/1164194; reweightPU = false;}
      else if (dataSetName == "tt_smallISR"){   	sprintf(name, "tt_smallISR");   	xlweight = lumi*163/1221491; reweightPU = false;}
      else if (dataSetName == "tt_matchingup"){   	sprintf(name, "tt_matchingup");   	xlweight = lumi*163/1036347; reweightPU = false;}
      else if (dataSetName == "tt_matchingdown"){   	sprintf(name, "tt_matchingdown");   	xlweight = lumi*163/937882;  reweightPU = false;}
      else if (dataSetName == "tt_scaleup"){   		sprintf(name, "tt_scaleup");   		xlweight = lumi*163/1010958; reweightPU = false;}
      else if (dataSetName == "tt_scaledown"){   	sprintf(name, "tt_scaledown");   	xlweight = lumi*163/1038835; reweightPU = false;}
      
      //Test file
      else{    	                                sprintf(name, "test");	        xlweight = 1;} 
      
      //Name the output file
      char rootFileName[100];
    
      if  (isRAW) sprintf(rootFileName,"outputs/naked_%d_%s.root", mode, name);
      else if(!isData && JESPlus) sprintf(rootFileName,"outputs/JESsysUp_%d_%s.root", mode, name);
      else if (!isData && JESMinus) sprintf(rootFileName,"outputs/JESsysDown_%d_%s.root", mode, name);
      else if (!isData && JERMinus) sprintf(rootFileName,"outputs/JERsysDown_%d_%s.root", mode, name);
      else if (!isData && JERPlus) sprintf(rootFileName,"outputs/JERsysUp_%d_%s.root", mode, name);
      else if (!isData && SFplus) sprintf(rootFileName,"outputs/SFsysUp_%d_%s.root", mode, name);
      else if (!isData && SFminus) sprintf(rootFileName,"outputs/SFsysDown_%d_%s.root", mode, name);
      else if (!isData && unclusteredUp) sprintf(rootFileName,"outputs/METsysUp_%d_%s.root", mode, name);
      else if (!isData && unclusteredDown) sprintf(rootFileName,"outputs/METsysDown_%d_%s.root", mode, name);
      else if (!isData && PUsysUp) sprintf(rootFileName,"outputs/PUsysUp_%d_%s.root", mode, name);
      else if (!isData && PUsysDown) sprintf(rootFileName,"outputs/PUsysDown_%d_%s.root", mode, name);
      else if (!isData && !reweightPU) sprintf(rootFileName,"outputs/noPU_%d_%s.root", mode, name);
      else sprintf(rootFileName,"outputs/out_%d_%s.root", mode, name);
      
      // Objects
      vector < TRootVertex* > vertex;
      vector < TRootMuon* > init_muons;
      vector < TRootElectron* > init_electrons;
      vector < TRootJet* > init_jets_corrected;
      vector < TRootJet* > init_jets;
      vector < TRootMET* > mets;
      vector < TLorentzVector > allForTopoCalc;
      vector<TRootGenJet*> genjets;
      
      TFile *fout = new TFile (rootFileName, "RECREATE");
      
      TRootEvent* event = 0;
      
      
      // Histograms
      map<string,TH1F*> histo1D;
      map<string,TH2F*> histo2D;
      
      TH1F* cutflow = new TH1F("cutflow", " ", 31,  -0.5, 30.5 );
      TH1F* cutflow_raw = new TH1F("cutflow_raw", " ", 31,  -0.5, 30.5 );
      TH1F* R = new TH1F( "R", " ", 40,  0, 40 );
      
      cutflow->Sumw2();
      cutflow_raw->Sumw2();
      R->Sumw2();
      
      
      // Branches of the output Tree
      double xlWeight; 
      double lum;
      
      int npu;
      int nvertex;
      
      double metPt;
      double metPx;
      double metPy;
      
      std::vector<double> *ptLepton;
      std::vector<double> *pxLepton;
      std::vector<double> *pyLepton;
      std::vector<double> *pzLepton;
      std::vector<double> *eLepton;
      std::vector<double> *qLepton;
      
      std::vector<double> *ptJet;
      std::vector<double> *pxJet;
      std::vector<double> *pyJet;
      std::vector<double> *pzJet;
      std::vector<double> *eJet;
      std::vector<double> *qJet;
      std::vector<double> *btTCHPJet;
      std::vector<double> *btTCHEJet;
      std::vector<double> *btSSVHPJet;
      std::vector<double> *btSSVHEJet;
      
      // Output Tree
      TTree* myTree = new TTree("myTree", "   ");
      
      myTree->Branch("xlWeight", &xlWeight, "xlWeight/D");
      myTree->Branch("lum", &lum, "lum/D");
      
      myTree->Branch("npu", &npu, "npu/I");
      myTree->Branch("nvertex", &nvertex, "nvertex/I");
      
      myTree->Branch("metPt", &metPt, "metPt/D");
      myTree->Branch("metPx", &metPx, "metPx/D");
      myTree->Branch("metPy", &metPy, "metPy/D");
      
      myTree->Branch("ptLepton","std::vector<double>",&ptLepton);
      myTree->Branch("pxLepton","std::vector<double>",&pxLepton);
      myTree->Branch("pyLepton","std::vector<double>",&pyLepton);
      myTree->Branch("pzLepton","std::vector<double>",&pzLepton);
      myTree->Branch("eLepton","std::vector<double>",&eLepton);
      myTree->Branch("qLepton","std::vector<double>",&qLepton);
      
      myTree->Branch("ptJet","std::vector<double>",&ptJet);
      myTree->Branch("pxJet","std::vector<double>",&pxJet);
      myTree->Branch("pyJet","std::vector<double>",&pyJet);
      myTree->Branch("pzJet","std::vector<double>",&pzJet);
      myTree->Branch("eJet","std::vector<double>",&eJet);
      myTree->Branch("qJet","std::vector<double>",&qJet);
      myTree->Branch("btTCHPJet","std::vector<double>",&btTCHPJet);
      myTree->Branch("btTCHEJet","std::vector<double>",&btTCHEJet);
      myTree->Branch("btSSVHPJet","std::vector<double>",&btSSVHPJet);
      myTree->Branch("btSSVHEJet","std::vector<double>",&btSSVHEJet);
      
      
      //Pile-Up reweighting  
      Lumi3DReWeighting Lumi3DWeights = Lumi3DReWeighting("../macros/PileUpReweighting/pileup_MC_Fall11.root","../macros/PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root", "pileup", "pileup");
     
      if(PUsysDown) Lumi3DWeights.weight3D_init(0.92);	
      else if(PUsysUp) Lumi3DWeights.weight3D_init(1.08);
      else Lumi3DWeights.weight3D_init(1.0);

      // Initialize JEC factors
      vector<JetCorrectorParameters> vCorrParam;
      
      // Create the JetCorrectorParameter objects, the order does not matter.
      // YYYY is the first part of the txt files: usually the global tag from which they are retrieved
      JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters("../macros/JECFiles/START42_V17_AK5PFchs_L3Absolute.txt");
      JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters("../macros/JECFiles/START42_V17_AK5PFchs_L2Relative.txt");
      JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters("../macros/JECFiles/START42_V17_AK5PFchs_L1FastJet.txt");
      
      //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
      vCorrParam.push_back(*L1JetPar);
      vCorrParam.push_back(*L2JetPar);
      vCorrParam.push_back(*L3JetPar);
      
      if(isData){
	JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("../macros/JECFiles/START42_V17_AK5PFchs_L2L3Residual.txt");
	vCorrParam.push_back(*ResJetCorPar);
      }
	  
      //OLD
      //JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("../macros/JECFiles/Jec11V2_db_AK5PFchs_Uncertainty.txt");
      //NEW
      JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("../macros/JECFiles/START42_V17_AK5PFchs_Uncertainty.txt");
      
      // true means redo also the L1
      JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);

     
      /// Start message
      
      cout << endl;
      cout << "[Info:] output rootfile named " << rootFileName << endl; 
      cout << "[Info:] mode = " << mode << ", lumi: " <<  lumi << " pb, sample: " << name << ", base weight: " << xlweight << endl;
      
      if (JERPlus ||JERMinus || JESPlus || JESMinus ||  SFplus || SFminus || unclusteredUp || unclusteredDown || SystSamples || Spring11 
	  || !reweightPU || !scaleFactor || PUsysUp || PUsysDown) {
	cout << "[Warning:] Non-standard options, ignore if you did it conciously" << endl;
	if (JERPlus) cout << "[Warning:] JER systematics on, plus" << endl;
	if (JERMinus) cout << "[Warning:] JER systematics on, minus" << endl;
	if (JESPlus) cout << "[Warning:] JES systematics on, plus" << endl;
	if (JESMinus) cout << "[Warning:] JES systematics on, minus" << endl;
	if (SFplus) cout <<"[Warning:] SF up 10% " << endl;
	if (SFminus) cout <<"[Warning:]  SF down 10% " << endl;
	if (unclusteredUp) cout <<"[Warning:] unclustered MET up 10% " << endl;
	if (unclusteredDown) cout <<"[Warning:] unclustered MET down 10% " << endl;
	if (SystSamples) cout <<"[Warning:] Systematic sample " << endl;
	if (Spring11) cout <<"[Warning:] Spring11 sample " << endl;
	if (!reweightPU && !isData) cout << "[Warning:] You are NOT applying PU re-weighting " << endl;
	if (!scaleFactor && !isData) cout << "[Warning:] You are NOT applying the b-tagging SF " << endl;
	if (PUsysUp) cout <<"[Warning:] PU up " << endl;
	if (PUsysDown) cout <<"[Warning:] PU down " << endl;
      } else cout << "[Info:] Standard setup " << endl;
      cout << "[Info:] " << datasets[d]->NofEvtsToRunOver() << " total events" << endl;
      
      
      for (int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
	{
	  
	  if(ievt%500 == 0) std::cout<<"Processing the "<<ievt<<"th event" <<flush<<"\r";
	  event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);
	  if(!isData)
	    {
	      genjets = treeLoader.LoadGenJet(ievt,false);
	      sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
	    }
	  // Weight given by the theoretical cross-section and lumi
	  double weight = xlweight;
          
	  // Pile-Up re-weighting
	  if (reweightPU && !isData && !Special){
	    double lumiWeight3D = 1.0;
	    lumiWeight3D = Lumi3DWeights.weight3D(event->nPu(-1),event->nPu(0),event->nPu(+1));
	    weight *= lumiWeight3D;
	  }
	   
	  //Trigger
	 
          //No trigger after first test
	  bool trigged = true;
	  
	  /*
	  int currentRun = event->runId();
	  bool itrigger = false;
	  bool isecondtrigger = false;
	  if(isData) { 
	    if (mode == 0){
	      if(currentRun >= 150000 && currentRun <= 161176){
		itrigger = treeLoader.iTrigger ("HLT_Mu17_Ele8_CaloIdL_v1", currentRun);
		isecondtrigger = treeLoader.iTrigger ("HLT_Mu8_Ele17_CaloIdL_v1", currentRun);
	      }else if(currentRun >= 161179 && currentRun <= 163261){
		itrigger = treeLoader.iTrigger ("HLT_Mu17_Ele8_CaloIdL_v2", currentRun);
		isecondtrigger = treeLoader.iTrigger ("HLT_Mu8_Ele17_CaloIdL_v2", currentRun);
	      }else if(currentRun >= 163262 && currentRun <= 164237){
		itrigger = treeLoader.iTrigger ("HLT_Mu17_Ele8_CaloIdL_v3", currentRun);
		isecondtrigger = treeLoader.iTrigger ("HLT_Mu8_Ele17_CaloIdL_v3", currentRun);
	      }else if(currentRun >= 165085 && currentRun <= 165888){
		itrigger = treeLoader.iTrigger ("HLT_Mu17_Ele8_CaloIdL_v4", currentRun);
		isecondtrigger = treeLoader.iTrigger ("HLT_Mu8_Ele17_CaloIdL_v4", currentRun);
	      }else if(currentRun >= 165900 && currentRun <= 166967){
		itrigger = treeLoader.iTrigger ("HLT_Mu17_Ele8_CaloIdL_v5", currentRun);
		isecondtrigger = treeLoader.iTrigger ("HLT_Mu8_Ele17_CaloIdL_v5", currentRun);
	      }else if(currentRun >= 166968 && currentRun <= 170053){
		itrigger = treeLoader.iTrigger ("HLT_Mu17_Ele8_CaloIdL_v6", currentRun);
		isecondtrigger = treeLoader.iTrigger ("HLT_Mu8_Ele17_CaloIdL_v6", currentRun);
	      }else if(currentRun >= 170054 && currentRun <= 173198){
		itrigger = treeLoader.iTrigger ("HLT_Mu17_Ele8_CaloIdL_v8", currentRun);
		isecondtrigger = treeLoader.iTrigger ("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v3", currentRun);
	      }else if(currentRun >= 173199 && currentRun <= 178380){
		itrigger = treeLoader.iTrigger ("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v4", currentRun);
		isecondtrigger = treeLoader.iTrigger ("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v4", currentRun);
	      }else if(currentRun >= 178381 && currentRun <= 999999){
		itrigger = treeLoader.iTrigger ("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v7", currentRun);
		isecondtrigger = treeLoader.iTrigger ("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v7", currentRun);
	      }
	    } else if (mode == 1){
	      if(currentRun >= 150000 && currentRun <= 161176){
		itrigger = treeLoader.iTrigger ("HLT_DoubleMu7_v1", currentRun);
	      }else if(currentRun >= 161179 && currentRun <= 163261){
		itrigger = treeLoader.iTrigger ("HLT_DoubleMu7_v1", currentRun);
	      }else if(currentRun >= 163262 && currentRun <= 164237){
		itrigger = treeLoader.iTrigger ("HLT_DoubleMu7_v2", currentRun);
	      }else if(currentRun >= 165085 && currentRun <= 165888){
		itrigger = treeLoader.iTrigger ("HLT_Mu13_Mu8_v2", currentRun);
	      }else if(currentRun >= 165900 && currentRun <= 167043){
		itrigger = treeLoader.iTrigger ("HLT_Mu13_Mu8_v2", currentRun);
	      }else if(currentRun >= 167044 && currentRun <= 170053){
		itrigger = treeLoader.iTrigger ("HLT_Mu13_Mu8_v4", currentRun);
	      }else if(currentRun >= 170054 && currentRun <= 173198){
		itrigger = treeLoader.iTrigger ("HLT_Mu13_Mu8_v6", currentRun);
	      }else if(currentRun >= 173199 && currentRun <= 178380){
		itrigger = treeLoader.iTrigger ("HLT_Mu13_Mu8_v7", currentRun);
	      }else if(currentRun >= 178381 && currentRun <= 999999){
		itrigger = treeLoader.iTrigger ("HLT_Mu17_Mu8_v10", currentRun);
              }
	    } else if (mode == 2){
	      if(currentRun >= 150000 && currentRun <= 161176){
		itrigger = treeLoader.iTrigger ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1", currentRun);
	      }else if(currentRun >= 161179 && currentRun <= 163261){
		itrigger = treeLoader.iTrigger ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2", currentRun);
	      }else if(currentRun >= 163262 && currentRun <= 164237){
		itrigger = treeLoader.iTrigger ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3", currentRun);
	      }else if(currentRun >= 165085 && currentRun <= 165888){
		itrigger = treeLoader.iTrigger ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4", currentRun);
	      }else if(currentRun >= 165900 && currentRun <= 167043){
		itrigger = treeLoader.iTrigger ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5", currentRun);
	      }else if(currentRun >= 167044 && currentRun <= 170053){
		itrigger = treeLoader.iTrigger ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6", currentRun);
	      }else if(currentRun >= 170054 && currentRun <= 170759){
		itrigger = treeLoader.iTrigger ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6", currentRun);
	      }else if(currentRun >= 170760 && currentRun <= 173198){
		itrigger = treeLoader.iTrigger ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7", currentRun);
	      }else if(currentRun >= 173199 && currentRun <= 178380){
		itrigger = treeLoader.iTrigger ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8", currentRun);
	      }else if(currentRun >= 178381 && currentRun <= 999999){
		itrigger = treeLoader.iTrigger ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9", currentRun);
              }
	    }
	    
	   //No trigger for quicker tests
	   // itrigger = true;
	   // isecondtrigger = true;
	    
	  } else {
	    // No trigger in MC
	    itrigger = true;
	    isecondtrigger = true;
	  }
          
	  
	  if (itrigger || isecondtrigger) trigged = true;
	 */
	 
	 
	  // JES CORRECTION
	  // Apply Jet Corrections on-the-fly
	  jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho());

	  // Correct MET Type I
	  if (metTypeI && !Special) jetTools->correctMETTypeOne(init_jets,mets[0]);  //Size of mets is never larger than 1 !!
	  

	  // Systematics
	  //JES and JER
	  //Special = true;
	  if (!Special && !isData){
	   
	    // vector<TRootGenJet*> genjets = treeLoader.LoadGenJet(ievt);

            if(JERPlus)  jetTools->correctJetJER(init_jets, genjets, "plus");
            else if(JERMinus) jetTools->correctJetJER(init_jets, genjets, "minus");
            else jetTools->correctJetJER(init_jets, genjets, "nominal");

            if (JESPlus) jetTools->correctJetJESUnc(init_jets, mets[0], "plus");
	    else if (JESMinus) jetTools->correctJetJESUnc(init_jets, mets[0], "minus");
	  }  
	   
	   
	   
	  //Start selection
	  Selection selection(init_jets, init_muons, init_electrons, mets);
	  
	  // PV cut (useless)
	  //bool isGoodPV = isGoodPV = selection.isPVSelected(vertex, 4,24,2.);  
	  bool isGoodPV = true;
	  
	  // Set up the unclustered MET systematic
	  double uncmet_px = mets[0]->Px();
	  double uncmet_py = mets[0]->Py();
	  for(unsigned int i=0; i<init_jets.size(); i++){
	    uncmet_px += init_jets[i]->Px();
	    uncmet_py += init_jets[i]->Py();
	  }
	  for(unsigned int i=0; i<init_muons.size(); i++){
	    uncmet_px += init_muons[i]->Px();
	    uncmet_py += init_muons[i]->Py();
	  }	
	  for(unsigned int i=0; i<init_electrons.size(); i++){
	    uncmet_px += init_electrons[i]->Px();
	    uncmet_py += init_electrons[i]->Py();
	  }	
	  
	  double met_px = mets[0]->Px();
	  double met_py = mets[0]->Py();
	  
	  if(unclusteredUp){
	    met_px += uncmet_px*0.1;
	    met_py += uncmet_py*0.1;
	  } if(unclusteredDown){
	    met_px -= uncmet_px*0.1;
	    met_py -= uncmet_py*0.1;
	  }
	  
	  double met_pt = sqrt(met_px*met_px + met_py*met_py);

	  // Cut Flow Starts
	  cutflow->Fill(1, weight);
	  cutflow_raw->Fill(1);
	  if(trigged){
	    cutflow->Fill(2, weight);
	    cutflow_raw->Fill(2);
	    if(isGoodPV){
	      cutflow->Fill(3, weight);
	      cutflow_raw->Fill(3);
            
              // Select Objects -> Cuts
	      selection.setJetCuts(20.,2.4,0.01,1.,0.98,0.3,0.1);
	      selection.setDiElectronCuts(20,2.5,0.15,0.02,1);
	      selection.setLooseElectronCuts(15,2.5,0.2);
	      selection.setDiMuonCuts(20,2.4,0.15,10,0.02);
	      selection.setLooseMuonCuts(10,2.5,0.2);
	  
	      //Select Objects 
	      vector<TRootElectron*> selectedElectrons = selection.GetSelectedDiElectrons(vertex[0]);
	      vector<TRootMuon*> selectedMuons = selection.GetSelectedDiMuons();
	      vector<TRootElectron*> looseElectrons = selection.GetSelectedLooseElectrons(true);
	      vector<TRootMuon*> looseMuons = selection.GetSelectedLooseMuons();
	      vector<TRootJet*> selectedJets = selection.GetSelectedJets(20,2.4,true);                    
	  
	      // Tight lepton selection
	      bool leptonSelection = false;
	      if 	(mode == 0 && selectedElectrons.size()== 1 && selectedMuons.size()== 1) leptonSelection = true;
	      else if 	(mode == 1 && selectedElectrons.size()== 0 && selectedMuons.size()== 2) leptonSelection = true;
	      else if 	(mode == 2 && selectedElectrons.size()== 2 && selectedMuons.size()== 0) leptonSelection = true;
              
	      if (leptonSelection) {
	  
		bool charge = false;
		double q0, q1;
		TLorentzVector lepton0, lepton1;
		if (mode == 0){
		  TRootElectron* electron = (TRootElectron*) selectedElectrons[0];
		  TRootMuon* muon = (TRootMuon*) selectedMuons[0];
		  if (electron->charge()*muon->charge() < 0) charge = true;
		  if (electron->Pt() > muon->Pt()){
		    lepton0.SetPxPyPzE(electron->Px(), electron->Py(), electron->Pz(), electron->Energy());
		    lepton1.SetPxPyPzE(muon->Px(), muon->Py(), muon->Pz(), muon->Energy());
		    q0 = electron->charge();
		    q1 = muon->charge();
		  } else {
		    lepton0.SetPxPyPzE(muon->Px(), muon->Py(), muon->Pz(), muon->Energy());
		    lepton1.SetPxPyPzE(electron->Px(), electron->Py(), electron->Pz(), electron->Energy());
		    q0 = muon->charge();
		    q1 = electron->charge();
		  }
		} else if (mode == 1){
		  TRootMuon* muon0 = (TRootMuon*) selectedMuons[0];
		  TRootMuon* muon1 = (TRootMuon*) selectedMuons[1];
		  if (muon0->charge()*muon1->charge() < 0) charge = true;
		  lepton0.SetPxPyPzE(muon0->Px(), muon0->Py(), muon0->Pz(), muon0->Energy());
		  lepton1.SetPxPyPzE(muon1->Px(), muon1->Py(), muon1->Pz(), muon1->Energy());
		  q0 = muon0->charge();
		  q1 = muon1->charge();
		} else {
		  TRootElectron* electron0 = (TRootElectron*) selectedElectrons[0];
		  TRootElectron* electron1 = (TRootElectron*) selectedElectrons[1];
		  if (electron0->charge()*electron1->charge() < 0) charge = true;
		  lepton0.SetPxPyPzE(electron0->Px(), electron0->Py(), electron0->Pz(), electron0->Energy());
		  lepton1.SetPxPyPzE(electron1->Px(), electron1->Py(), electron1->Pz(), electron1->Energy());
		  q0 = electron0->charge();
		  q1 = electron1->charge();
		}
                
		if (charge){
		  cutflow->Fill(4, weight);
		  cutflow_raw->Fill(4);
		  // Loose lepton veto
		  bool leptonVeto = false;
		  if 	  (mode == 0 && looseMuons.size()== 1 && looseElectrons.size() == 1) leptonVeto = true;
		  else if (mode == 1 && looseMuons.size()== 2 && looseElectrons.size() == 0) leptonVeto = true;
		  else if (mode == 2 && looseMuons.size()== 0 && looseElectrons.size() == 2) leptonVeto = true;
                  
		  if (leptonVeto) {
		    cutflow->Fill(5, weight);
		    cutflow_raw->Fill(5);
                   
		    // Low mll cut (all final states)
		    TLorentzVector pair = lepton0 + lepton1;   
		    if (pair.M() > 20){
		      
		      //Filling the Tree (at pre-selection level, leptons and mll)
		      lum = lumi;
                      
		      if (isRAW) weight = 1;
		      xlWeight = weight;
		      
		      npu = event->nPu(0);
		      nvertex = vertex.size();
		      
		      metPt = met_pt;
		      metPx = met_px;
		      metPy = met_py;
		      
		      ptLepton = new std::vector<double>; 
		      pxLepton = new std::vector<double>; 
		      pyLepton = new std::vector<double>; 
		      pzLepton = new std::vector<double>; 
		      eLepton = new std::vector<double>; 
		      qLepton = new std::vector<double>;
                      
		      ptJet = new std::vector<double>; 
		      pxJet = new std::vector<double>; 
		      pyJet = new std::vector<double>; 
		      pzJet = new std::vector<double>; 
		      eJet = new std::vector<double>; 
		      qJet = new std::vector<double>; 
		      btTCHPJet = new std::vector<double>; 
		      btTCHEJet = new std::vector<double>; 
		      btSSVHPJet = new std::vector<double>;
		      btSSVHEJet = new std::vector<double>; 
                      
		      ptLepton->push_back(lepton0.Pt());
		      ptLepton->push_back(lepton1.Pt());
                      
		      pxLepton->push_back(lepton0.Px());
		      pxLepton->push_back(lepton1.Px());
                      
		      pyLepton->push_back(lepton0.Py());
		      pyLepton->push_back(lepton1.Py());
                      
		      pzLepton->push_back(lepton0.Pz());
		      pzLepton->push_back(lepton1.Pz());
                      
		      eLepton->push_back(lepton0.Energy());
		      eLepton->push_back(lepton1.Energy());
                      
		      qLepton->push_back(q0);
		      qLepton->push_back(q1);
                      
		      for (unsigned int i =0; i < selectedJets.size(); i ++){
			TRootJet* tempJet = (TRootJet*) selectedJets[i];
			ptJet->push_back(tempJet->Pt());
			pxJet->push_back(tempJet->Px());
			pyJet->push_back(tempJet->Py());
			pzJet->push_back(tempJet->Pz());
			eJet->push_back(tempJet->Energy());
			qJet->push_back(tempJet->charge());
			btTCHPJet->push_back(tempJet->btag_trackCountingHighPurBJetTags() );
			btTCHEJet->push_back(tempJet->btag_trackCountingHighEffBJetTags() );
			btSSVHPJet->push_back(tempJet->btag_simpleSecondaryVertexHighPurBJetTags() );
			btSSVHEJet->push_back(tempJet->btag_simpleSecondaryVertexHighEffBJetTags() );  
		      }
		      // do calculations as far as topology goes:
                      
		      allForTopoCalc.push_back(lepton0);
		      allForTopoCalc.push_back(lepton1);
		      for(unsigned int ijet=0;ijet<selectedJets.size(); ijet++){
			TRootJet* tempJet = (TRootJet*) selectedJets[ijet];
			TLorentzVector tempLV(tempJet->Px(), tempJet->Py(), tempJet->Pz(), tempJet->E());
			allForTopoCalc.push_back(tempLV);
		      }
		      // std::cout << "creating topology worker " << std::endl;
		      TopologyWorker::TopologyWorker topoWorkerNoBoost(false); // set the bool to true to un-boost events
		      topoWorkerNoBoost.setPartList(allForTopoCalc,allForTopoCalc);
		      topoWorkerNoBoost.setVerbose(true);
		      // see for other methods: header file of TopologyWorker
                      
		      /*
			std::cout << "summary: " << std::endl;
			std::cout << "oblateness: " << topoWorkerNoBoost.oblateness() << std::endl;
			std::cout << "sphericity: " << topoWorkerNoBoost.get_sphericity() << std::endl;
			std::cout << "aplanarity: " << topoWorkerNoBoost.get_aplanarity() << std::endl;
			std::cout << "njetW: " << topoWorkerNoBoost.get_njetW() <<std::endl;
			std::cout << "sqrts: " << topoWorkerNoBoost.get_sqrts() << std::endl; 
		      */
		      
		      // 
		      allForTopoCalc.clear();
		      allForTopoCalc.resize(0);
                      
		      myTree->Fill();
                      
		      delete ptLepton;
		      delete pxLepton;
		      delete pyLepton;
		      delete pzLepton;
		      delete eLepton;
		      delete qLepton;
                      
		      delete ptJet;
		      delete pxJet;
		      delete pyJet;
		      delete pzJet;
		      delete eJet;
		      delete qJet;
		      delete btTCHPJet;
		      delete btTCHEJet;
		      delete btSSVHPJet;
		      delete btSSVHEJet;
                      
		      // Invariant mass in ee and mumu
		      if (pair.M() > 101 || pair.M() < 81 || mode == 0){
			cutflow->Fill(6, weight);
			cutflow_raw->Fill(6);
			// MET in ee and mumu
			if (met_pt > 30 || mode == 0){
			  cutflow->Fill(7, weight);
			  cutflow_raw->Fill(7);
                          
			  //Jet and b-tag selection
			  int nJetsBT = 0;
			  int nJets = 0;
			  bool bTagged = false;
			  int iJet = -5;
			  int iSF;
			  double tempSF = SFval;
			  if (SFminus) 	tempSF = SFval - SFval*10/100;
			  if (SFplus) 	tempSF = SFval + SFval*10/100;
			  int SFvalue = int(tempSF*100 + 1);
			  
			  if (isData || !scaleFactor){
			    for (unsigned int i =0; i < selectedJets.size(); i ++){
			      TRootJet* tempJet = (TRootJet*) selectedJets[i];
			      TLorentzVector tJet(tempJet->Px(), tempJet->Py(), tempJet->Pz(), tempJet->Energy());
			      if (tempJet->Pt() > 30 && TMath::Min(fabs(lepton0.DeltaR(tJet)), fabs(lepton1.DeltaR(tJet))) > 0.3) {
				nJets++;
				//if (iJet == -5) iJet = i;
				iJet = i;
				if (tempJet->btag_simpleSecondaryVertexHighEffBJetTags() > 1.74){
				  bTagged = true;
				  nJetsBT++;
				} 
			      } else if (tempJet->btag_simpleSecondaryVertexHighEffBJetTags() > 1.74) nJetsBT++;
			    }
			  } else {
			    //// Regular SF
			    if (SFvalue < 101){
			      for (unsigned int i =0; i < selectedJets.size(); i ++){
				TRootJet* tempJet = (TRootJet*) selectedJets[i];
				TLorentzVector tJet(tempJet->Px(), tempJet->Py(), tempJet->Pz(), tempJet->Energy());
				if (tempJet->Pt() > 30 && TMath::Min(fabs(lepton0.DeltaR(tJet)), fabs(lepton1.DeltaR(tJet))) > 0.3) {
				  nJets++;
				  //if (iJet == -5) iJet = i;
				  iJet = i;
				  if (tempJet->btag_simpleSecondaryVertexHighEffBJetTags() > 1.74){
				    iSF = rand() % 101;
				    if (iSF < SFvalue ){
				      bTagged = true;
				      nJetsBT++;
				    } 
				  } 
				} else if (tempJet->btag_simpleSecondaryVertexHighEffBJetTags() > 1.74){
				  iSF = rand() % 101;
				  if (iSF < SFvalue ) nJetsBT++;
				}
			      }
			    } else {
			      //// Large SF
			      for (unsigned int i =0; i < selectedJets.size(); i ++){
				TRootJet* tempJet = (TRootJet*) selectedJets[i];
				TLorentzVector tJet(tempJet->Px(), tempJet->Py(), tempJet->Pz(), tempJet->Energy());
				if (tempJet->Pt() > 30 && TMath::Min(fabs(lepton0.DeltaR(tJet)), fabs(lepton1.DeltaR(tJet))) > 0.3) {
				  nJets++;
				  //if (iJet == -5) iJet = i;
				  iJet = i;
				  if (tempJet->btag_simpleSecondaryVertexHighEffBJetTags() > 1.74){
				    bTagged = true;
				    nJetsBT++;
				  } else {
				    iSF = rand() % 101;
				    if (iSF < abs(100 - SFvalue)){
				      nJetsBT++;
				      bTagged = true;
				    }
				  }
				} else if (tempJet->btag_simpleSecondaryVertexHighEffBJetTags() > 1.74){
				  nJetsBT++;
				} else {
				  iSF = rand() % 101;
				  if (iSF < abs(100 - SFvalue)) nJetsBT++;
				}
			      }
			    }
			  }
			  
			  
			  // Filling all the regions
			  if (nJets !=0){
			    TRootJet* jet = (TRootJet*) selectedJets[iJet];
                  	    double ptSysPx = lepton0.Px() + lepton1.Px() + jet->Px() + met_px;
			    double ptSysPy = lepton0.Py() + lepton1.Py() + jet->Py() + met_py;
			    double ptSys = sqrt(ptSysPx*ptSysPx + ptSysPy*ptSysPy);
			    double Ht = lepton0.Pt() + lepton1.Pt() + jet->Pt() + met_pt; 
			    if (ptSys < 60){
			      if (Ht > 160 || mode != 0){
				if (nJets == 1 && nJetsBT == 1 && bTagged)R->Fill(1, weight);
				if (nJets == 1 && nJetsBT == 2)  R->Fill(2, weight);
				if (nJets == 1 && nJetsBT > 0)  R->Fill(3, weight);
				if (nJets == 1 && nJetsBT > 1)  R->Fill(4, weight);
				if (nJets == 2 && nJetsBT == 0)  R->Fill(5, weight);
				if (nJets == 2 && nJetsBT == 1)  R->Fill(6, weight);
				if (nJets == 2 && nJetsBT == 2)  R->Fill(7, weight);
				if (nJets == 2 && nJetsBT > 0)  R->Fill(8, weight);
				if (nJets == 2 && nJetsBT > 1)  R->Fill(9, weight);
				if (nJets > 1 && nJetsBT == 0)  R->Fill(10, weight);
				if (nJets > 1 && nJetsBT == 1)  R->Fill(11, weight);
				if (nJets > 1 && nJetsBT == 2)  R->Fill(12, weight);
				if (nJets > 1 && nJetsBT !=0 )  R->Fill(13, weight);
				if (nJets > 1 && nJetsBT > 1 )  R->Fill(14, weight);
				if (nJets == 3 && nJetsBT ==0 )  R->Fill(15, weight);
				if (nJets == 3 && nJetsBT ==1 )  R->Fill(16, weight);
				if (nJets == 3 && nJetsBT ==2 )  R->Fill(17, weight);
				if (nJets == 3 && nJetsBT ==3 )  R->Fill(18, weight);
			      }
			    }
			  }
			  
			  // Filling the signal region
			  if(nJets == 1){
			    TRootJet* jet = (TRootJet*) selectedJets[iJet];
			    cutflow->Fill(8, weight);
			    cutflow_raw->Fill(8);
			    if (bTagged && nJetsBT == 1){
			      cutflow->Fill(9,weight);
			      cutflow_raw->Fill(9);
                              
			      double ptSysPx = lepton0.Px() + lepton1.Px() + jet->Px() + met_px;
			      double ptSysPy = lepton0.Py() + lepton1.Py() + jet->Py() + met_py;
			      double ptSys = sqrt(ptSysPx*ptSysPx + ptSysPy*ptSysPy);
			      double Ht = lepton0.Pt() + lepton1.Pt() + jet->Pt() + met_pt; 
                              
			      if (ptSys < 60){
				cutflow->Fill(10, weight);
				cutflow_raw->Fill(10);
				if (Ht > 160 || mode != 0){
				  cutflow->Fill(11, weight);
				  cutflow_raw->Fill(11);
				}
			      }
			    }
			  }
			  //

			}	      
		      }      
		    }
		  }
		}
	      }
	    }      
	  }
	
	} // event loop
      
      // cleanup
      allForTopoCalc.clear();
      allForTopoCalc.resize(0);
      
      if(jetTools) delete jetTools;   
      
      cout << "--------------------------------------------------" << endl;
      cout << "[Results Normalized:] " <<  endl;
      cout << "All:       " <<  cutflow->GetBinContent(2) << " +/- "  << cutflow->GetBinError(2) << endl;
      cout << "HLT:       " <<  cutflow->GetBinContent(3) << " +/- "  << cutflow->GetBinError(3) << endl;
      cout << "PV:        " <<  cutflow->GetBinContent(4) << " +/- "  << cutflow->GetBinError(4) << endl;
      cout << "Lep. Sel:  " <<  cutflow->GetBinContent(5) << " +/- "  << cutflow->GetBinError(5) << endl;
      cout << "Lep. Veto: " <<  cutflow->GetBinContent(6) << " +/- "  << cutflow->GetBinError(6) << endl;
      cout << "mll:       " <<  cutflow->GetBinContent(7) << " +/- "  << cutflow->GetBinError(7) << endl;
      cout << "MET:       " <<  cutflow->GetBinContent(8) << " +/- "  << cutflow->GetBinError(8) << endl;
      cout << "1 jet:     " <<  cutflow->GetBinContent(9) << " +/- "  << cutflow->GetBinError(9) << endl;
      cout << "1 jet BT:  " <<  cutflow->GetBinContent(10) << " +/- " << cutflow->GetBinError(10) << endl;
      cout << "Pt tW:     " <<  cutflow->GetBinContent(11) << " +/- " << cutflow->GetBinError(11) << endl;
      cout << "Ht:        " <<  cutflow->GetBinContent(12) << " +/- " << cutflow->GetBinError(12) << endl;
      
      cout << "--------------------------------------------------" << endl;
      cout << "[Results Raw:] " <<  endl;
      cout << "All:       " <<  cutflow_raw->GetBinContent(2) << " +/- "  << cutflow_raw->GetBinError(2) << endl;
      cout << "HLT:       " <<  cutflow_raw->GetBinContent(3) << " +/- "  << cutflow_raw->GetBinError(3) << endl;
      cout << "PV:        " <<  cutflow_raw->GetBinContent(4) << " +/- "  << cutflow_raw->GetBinError(4) << endl;
      cout << "Lep. Sel:  " <<  cutflow_raw->GetBinContent(5) << " +/- "  << cutflow_raw->GetBinError(5) << endl;
      cout << "Lep. Veto: " <<  cutflow_raw->GetBinContent(6) << " +/- "  << cutflow_raw->GetBinError(6) << endl;
      cout << "mll:       " <<  cutflow_raw->GetBinContent(7) << " +/- "  << cutflow_raw->GetBinError(7) << endl;
      cout << "MET:       " <<  cutflow_raw->GetBinContent(8) << " +/- "  << cutflow_raw->GetBinError(8) << endl;
      cout << "1 jet:     " <<  cutflow_raw->GetBinContent(9) << " +/- "  << cutflow_raw->GetBinError(9) << endl;
      cout << "1 jet BT:  " <<  cutflow_raw->GetBinContent(10) << " +/- " << cutflow_raw->GetBinError(10) << endl;
      cout << "Pt tW:     " <<  cutflow_raw->GetBinContent(11) << " +/- " << cutflow_raw->GetBinError(11) << endl;
      cout << "Ht:        " <<  cutflow_raw->GetBinContent(12) << " +/- " << cutflow_raw->GetBinError(12) << endl;
      
      cout << "--------------------------------------------------" << endl;
      cout << "[Jet Multiplicity Check:]" << endl;
      cout << "1 jet 1 tag: " << R->GetBinContent(2) << " +/- " << R->GetBinError(2) << endl;
      cout << "2 jet 1 tag: " << R->GetBinContent(7) << " +/- " << R->GetBinError(2) << endl;
      cout << "2 jet 2 tag: " << R->GetBinContent(8) << " +/- " << R->GetBinError(2) << endl;
      cout << "--------------------------------------------------" << endl;
    
      fout->Write();
      fout->Close();
      
    }// dataset loop
  
  cout << "It took you " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
  
}

