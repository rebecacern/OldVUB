// rebeca@cern.ch
//
// Summer12 round
// 8 Tev

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
#include "../MCInformation/interface/ResolutionFit.h"
#include "../Reconstruction/interface/JetCorrectorParameters.h"
#include "../Reconstruction/interface/JetCorrectionUncertainty.h"
#include "../MCInformation/interface/Lumi3DReWeighting.h"
#include "../MCInformation/interface/LumiReWeighting.h"
#include "../macros/Style.C"

using namespace std;
using namespace TopTree;
using namespace reweight;

int main(int argc, char* argv[]) {
  
  setMyStyle();
  
  clock_t start = clock();
  
  //  Modes: 0 emu, 1 mumu, 2 ee 
  // By defaul, emu channel
  
  int  mode = 0; 
  double lumi = 1000;
  string xmlfile ="twemu.xml";
  
  bool reweightPU = true;
  bool Pu3D = false;
  
  //b-tag scale factor
  bool scaleFactor = false;
  
  // Systematic Calculations
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
  
  // Naked option
  bool isRAW = false;
  
  // Run HLT
  bool runHLT = false;
  
  
  cout << "[Info:] start of the 2012 Data analysis with Summer12 Monte Carlo " << endl;
  
  // Arguments
  for(int iarg = 0; iarg < argc && argc>1 ; iarg++){
    std::string argval=argv[iarg];
    if(argval=="--help" || argval =="--h"){
      cout << "--ee:   Di-Electron" << endl;
      cout << "--emu:  Electron-Muon" << endl;
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
      cout << "--3D: 3D Pileup reweighting" << endl;
      return 0;
    }
    if (argval=="--ee") 	mode = 2;
    if (argval=="--emu") 	mode = 0;
    if (argval=="--mumu") 	mode = 1;
    
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
    if (argval=="--3D") Pu3D = true;
    
  }   
  
  // Luminosity and xml files
  if      (mode == 0){ 	 lumi = 5085.246;  xmlfile ="twemu.xml";}
  else if (mode == 1){	 lumi = 1000;  xmlfile = "twmumu.xml";}
  else if (mode == 2){	 lumi = 5103.58;  xmlfile = "twee.xml";}
 
 
   
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
      bool isTop = false;
      char name[100];
      double xlweight;
      
      // cross sections and weights
      if (dataSetName == "data"){		sprintf(name, "data");  	xlweight = 1; 				isData = true;}//
      else if (dataSetName == "tt"){            sprintf(name, "tt");            xlweight = lumi*225.197/6709118; 	isTop = true;} //
      else if (dataSetName == "twdr"){         	sprintf(name, "tw_dr");         xlweight = lumi*11.1/497657; 		} 
      else if (dataSetName == "atwdr"){         sprintf(name, "atw_dr");        xlweight = lumi*11.1/493460; 		} //
      else if (dataSetName == "t"){         	sprintf(name, "t");         	xlweight = lumi*56.4/23777; 		} //
      else if (dataSetName == "at"){         	sprintf(name, "at");         	xlweight = lumi*30.7/1935071; 		} //
      else if (dataSetName == "ww"){         	sprintf(name, "ww");         	xlweight = lumi*57.07/9969958; 		} //
      else if (dataSetName == "wz"){         	sprintf(name, "wz");         	xlweight = lumi*22.44/8080197; 		} //
      else if (dataSetName == "zz"){         	sprintf(name, "zz");         	xlweight = lumi*9.03/9799902; 		} //
      else if (dataSetName == "zjets"){         sprintf(name, "zjets");         xlweight = lumi*3532.8/16080506; 	} //
      else if (dataSetName == "zjets_lowmll"){  sprintf(name, "zjets_lowmll");  xlweight = lumi*860.5/7132214; 	} //
      else if (dataSetName == "wjets"){  	sprintf(name, "wjets");  	xlweight = lumi*37509/18036994; 	} // 
          
	  
	  
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
      else if (!isData && !reweightPU) sprintf(rootFileName,"outputs/out_noPU_%d_%s.root", mode, name);
      else if (!isData && Pu3D) sprintf(rootFileName,"outputs/out_3D_%d_%s.root", mode, name);
      else sprintf(rootFileName,"outputs/out_%d_%s.root", mode, name);
      
      char myTexFile[300];
      sprintf(myTexFile,"lepsel_info_run_lumi_event_%d_%s.txt", mode, name);
      ofstream salida(myTexFile);
      sprintf(myTexFile,"lepveto_run_lumi_event_%d_%s.txt", mode, name);
      ofstream salida2(myTexFile);
      sprintf(myTexFile,"met_run_lumi_event_%d_%s.txt", mode, name);
      ofstream salida3(myTexFile);
      sprintf(myTexFile,"jet_run_lumi_event_%d_%s.txt", mode, name);
      ofstream salida4(myTexFile);
      sprintf(myTexFile,"bt_run_lumi_event_%d_%s.txt", mode, name);
      ofstream salida5(myTexFile);
      
      
      // Objects
      vector < TRootVertex* > vertex;
      vector < TRootMuon* > init_muons;
      vector < TRootElectron* > init_electrons;
      vector < TRootJet* > init_jets;
      vector < TRootJet* > init_jets_corrected;
      vector < TRootMET* > mets;
      vector<TRootGenJet*> genjets;
       
      
      // Jets TO BE CHECKED!!!!
      
      vector<JetCorrectorParameters> vCorrParam;
      
      // Create the JetCorrectorParameter objects, the order does not matter.
      // YYYY is the first part of the txt files: usually the global tag from which they are retrieved
      JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters("JECFiles/Summer12_V3_MC_L3Absolute_AK5PFchs.txt");
      JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters("JECFiles/Summer12_V3_MC_L2Relative_AK5PFchs.txt");
      JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters("JECFiles/Summer12_V3_MC_L1FastJet_AK5PFchs.txt");

      //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
      vCorrParam.push_back(*L1JetPar);
      vCorrParam.push_back(*L2JetPar);
      vCorrParam.push_back(*L3JetPar);

      if(!isData) { // DATA!
	JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("JECFiles/Summer12_V3_DATA_L2L3Residual_AK5PFchs.txt");
	vCorrParam.push_back(*ResJetCorPar);	
      } 
      
      //I think this is not used!
      JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/Summer12_V3_MC_Uncertainty_AK5PFchs.txt");
      // if (isData) JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/Summer12_V3_DATA_Uncertainty_AK5PFchs.txt");
    
      // true means redo also the L1
      JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true); 
       
      
      
      // Lumi re-weighting  TO CHECK 
      LumiReWeighting LumiWeights;
      LumiWeights = LumiReWeighting("pileupHistos/Summer12.root", "pileupHistos/Run2012AB_new.root", "pileup", "pileup");

      //systematic
      reweight::PoissonMeanShifter PShiftDown_ = reweight::PoissonMeanShifter(-0.6);
      reweight::PoissonMeanShifter PShiftUp_ = reweight::PoissonMeanShifter(0.6);



      TFile *fout = new TFile (rootFileName, "RECREATE");
      
      TRootEvent* event = 0;
      
      
      // Histograms
      map<string,TH1F*> histo1D;
      map<string,TH2F*> histo2D;
      
      TH1F* cutflow = new TH1F("cutflow", " ", 31,  -0.5, 30.5 );
      TH1F* cutflow_raw = new TH1F("cutflow_raw", " ", 31,  -0.5, 30.5 );
      TH1F* R = new TH1F( "R", " ", 40,  0, 40 );
      
      TH1F* pileup_weights = new TH1F( "pileup_weights", " ", 1000,  0, 10 );
      TH1F* pileup_weights_3D = new TH1F( "pileup_weights_3D", " ", 1000,  0, 10 );
      
      cutflow->Sumw2();
      cutflow_raw->Sumw2();
      R->Sumw2();
      
      pileup_weights->Sumw2();
      pileup_weights_3D->Sumw2();
      
      
      // Branches of the output Tree
      double xlWeight; 
      double puweight;
      double rawWeight;
      
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
      std::vector<double> *btJPBJet;
      std::vector<double> *btBJPBJet;
      std::vector<double> *btCSVBJet;
      std::vector<double> *btCSVBmvaJet;
      
      // Output Tree
      TTree* myTree = new TTree("myTree", "   ");
      
      myTree->Branch("xlWeight", &xlWeight, "xlWeight/D");
      myTree->Branch("puweight", &puweight, "puweight/D");
      myTree->Branch("rawWeight", &rawWeight, "rawWeight/D");
      
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
      myTree->Branch("btJPBJet","std::vector<double>",&btJPBJet);
      myTree->Branch("btBJPBJet","std::vector<double>",&btBJPBJet);
      myTree->Branch("btCSVBJet","std::vector<double>",&btCSVBJet);
      myTree->Branch("btCSVBmvaJet","std::vector<double>",&btCSVBmvaJet);
      
      cout << "[Info:] output rootfile named " << rootFileName << endl; 
      cout << "[Info:] mode = " << mode << ", lumi: " <<  lumi << " pb, sample: " << name << ", base weight: " << xlweight << endl;
      
      if (JERPlus ||JERMinus || JESPlus || JESMinus ||  SFplus || SFminus || unclusteredUp || unclusteredDown 
	  || !reweightPU || !scaleFactor || PUsysUp || PUsysDown || Pu3D) {
	cout << "[Warning:] Non-standard options, ignore if you did it conciously" << endl;
	if (JERPlus) cout << "[Warning:] JER systematics on, plus" << endl;
	if (JERMinus) cout << "[Warning:] JER systematics on, minus" << endl;
	if (JESPlus) cout << "[Warning:] JES systematics on, plus" << endl;
	if (JESMinus) cout << "[Warning:] JES systematics on, minus" << endl;
	if (SFplus) cout <<"[Warning:] SF up 10% " << endl;
	if (SFminus) cout <<"[Warning:]  SF down 10% " << endl;
	if (unclusteredUp) cout <<"[Warning:] unclustered MET up 10% " << endl;
	if (unclusteredDown) cout <<"[Warning:] unclustered MET down 10% " << endl;
	if (!reweightPU && !isData) cout << "[Warning:] You are NOT applying PU re-weighting, we keep it this way for Summer12 " << endl;
	if (!scaleFactor && !isData) cout << "[Warning:] You are NOT applying the b-tagging SF " << endl;
	if (PUsysUp) cout <<"[Warning:] PU up " << endl;
	if (PUsysDown) cout <<"[Warning:] PU down " << endl;
	
      } else cout << "[Info:] Standard setup " << endl;
      cout << "[Info:] " << datasets[d]->NofEvtsToRunOver() << " total events" << endl;
      if (runHLT) cout << "[Info:] You have the HLT activated, this might be slower than the usual. " << endl;
      
      for (int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
	{
	  
	  if(ievt%500 == 0) std::cout<<"Processing the "<<ievt<<"th event" <<flush<<"\r";
	  event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets_corrected, mets);

	  double weight = xlweight;
	
	  if(!isData) {
	    genjets = treeLoader.LoadGenJet(ievt,false);
	    sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
	  }
	 
	  // Applying Pile-Up re-weighting
	  double lumiWeight = 1.0;
	  if (!isData ){
	    lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );
	    if(PUsysDown) lumiWeight = lumiWeight*PShiftDown_.ShiftWeight( event->nPu(0) );
            else if(PUsysUp) lumiWeight = lumiWeight*PShiftUp_.ShiftWeight( event->nPu(0) );
	    pileup_weights->Fill(lumiWeight);
	    if (reweightPU) weight *=lumiWeight;
	    
	  }
	   
	  // JER (All MC)
	  if (JERMinus)jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "minus",false); //false means don't use old numbers but newer ones...
	  else if (JERPlus) jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "plus",false);
	  else jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "nominal",false);

          // JES (Systematic only)
	  if (JESPlus) jetTools->correctJetJESUnc(init_jets_corrected, "minus",1);
	  else if (JESMinus) jetTools->correctJetJESUnc(init_jets_corrected, "plus",1);
	   
	  
	  //Trigger
	  bool trigged = false;
	  if (runHLT){
	    int currentRun = event->runId();
	    bool itrigger = false;
	    bool isecondtrigger = false;
	    if(isData) { 
	      if (mode == 0){
		itrigger = treeLoader.iTrigger ("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6", currentRun);
		isecondtrigger = treeLoader.iTrigger ("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6", currentRun);
	      } else if (mode == 1){
		itrigger = treeLoader.iTrigger ("HLT_Mu17_Mu8_v16", currentRun);
		isecondtrigger = treeLoader.iTrigger ("HLT_Mu17_TkMu8_v9", currentRun);
	      } else if (mode == 2){
		itrigger = treeLoader.iTrigger ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17", currentRun);
	      }
	    } else {
	      // No trigger in MC
	      itrigger = true;
	      isecondtrigger = true;
	    }
	    if (itrigger || isecondtrigger) trigged = true;
	  } else trigged = true;
	    
	  //Start selection
	  Selection selection(init_jets_corrected, init_muons, init_electrons, mets);
	  
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
	      selection.setJetCuts(20.,5.,0.01,1.,0.98,0.3,0.1);
              selection.setDiMuonCuts(20,2.4,0.20,999.);
              selection.setDiElectronCuts(20,2.5,0.15,0.04,0.,1,0.3,1);
              selection.setLooseMuonCuts(10,2.5,0.2);
              selection.setLooseDiElectronCuts(15,2.5,0.2); 
		
	      //Select Objects 
	      vector<TRootJet*> selectedJets = selection.GetSelectedJets(true);
	      vector<TRootMuon*> selectedMuons = selection.GetSelectedDiMuons();
	      vector<TRootMuon*> looseMuons = selection.GetSelectedLooseMuons();
	      vector<TRootElectron*> selectedElectrons = selection.GetSelectedDiElectrons();
	      vector<TRootElectron*> looseElectrons = selection.GetSelectedLooseDiElectrons();
	   
	
	      /*
		if (  event->lumiBlockId() == 1104  && event->eventId() == 331025 ){
		//  cout << event->runId() << "\t" << event->lumiBlockId() << "\t" << event->eventId() << endl;
		cout << "HERE I AM! " << selectedMuons.size() << ", " << selectedElectrons.size() << ", " <<looseMuons.size() << ", " << looseElectrons.size() << endl; 
		TRootElectron* muon1 = (TRootElectron*) selectedElectrons[0];
		TRootElectron* muon2 = (TRootElectron*) selectedElectrons[1];
		float reliso = (muon1->chargedHadronIso() + max( 0.0, muon1->neutralHadronIso() + muon1->photonIso() - 0.5*muon1->puChargedHadronIso() ) ) / muon1->Pt(); // dBeta corrected
		cout << "pt \t eta \t reliso " << endl;
		cout << muon1->Pt() << ", " << muon1->Eta() << ", " << reliso << endl;
		reliso = (muon2->chargedHadronIso() + max( 0.0, muon2->neutralHadronIso() + muon2->photonIso() - 0.5*muon2->puChargedHadronIso() ) ) / muon2->Pt(); // dBeta corrected
		     
		cout << muon2->Pt() << ", " << muon2->Eta() << ", " << reliso << endl;
		     
		cout << "Electrons" << endl;
		for(unsigned int i=0; i<init_electrons.size(); i++){
		float reliso = (init_electrons[i]->chargedHadronIso() + max( 0.0, init_electrons[i]->neutralHadronIso() + init_electrons[i]->photonIso() - 0.5*init_electrons[i]->puChargedHadronIso() ) ) / init_electrons[i]->Pt();
		cout << i  << ": " << init_electrons[i]->Pt() << ", " << init_electrons[i]->Eta() << ", " << reliso << ", "
		<<  init_electrons[i]->d0() << ", " <<  init_electrons[i]->passConversion() << ", " << init_electrons[i]->mvaTrigId() << ", " << init_electrons[i]->missingHits() <<  endl;

		}	
		cout << "Muons" << endl;

		for(unsigned int i=0; i<init_muons.size(); i++){
		float reliso = (init_muons[i]->chargedHadronIso() + max( 0.0, init_muons[i]->neutralHadronIso() + init_muons[i]->photonIso() - 0.5*init_muons[i]->puChargedHadronIso() ) ) / init_muons[i]->Pt();
		cout << i  << ": " << init_muons[i]->Pt() << ", " << init_muons[i]->Eta() << ", " << reliso << endl;

		}	
		  
		}
		
	      */
	   
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
		  
		  salida << event->runId() << "\t" << event->lumiBlockId() << "\t" << event->eventId() << endl;
		 
		
		
		  
		  if (leptonVeto) {
		    cutflow->Fill(5, weight);
		    cutflow_raw->Fill(5);

		    salida2 << event->runId() << "\t" << event->lumiBlockId() << "\t" << event->eventId() << endl;
		    // Low mll cut (all final states)
		    TLorentzVector pair = lepton0 + lepton1;   
		    if (pair.M() > 20){
			
		      //Filling the Tree (at pre-selection level, leptons and mll)
		      lum = lumi;
			
			
		      xlWeight = weight;
			
		      puweight =  lumiWeight;
		      rawWeight = xlweight;
					
		      npu = event->nPu(0);
		      nvertex = vertex.size();
			
		      metPt = mets[0]->Pt();
		      metPx = mets[0]->Px();
		      metPy = mets[0]->Py();
		      
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
		      btJPBJet = new std::vector<double>; 
		      btBJPBJet = new std::vector<double>; 
		      btCSVBJet = new std::vector<double>;
		      btCSVBmvaJet = new std::vector<double>; 
			
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
			btJPBJet->push_back(tempJet->btag_jetProbabilityBJetTags() );
			btBJPBJet->push_back(tempJet->btag_jetBProbabilityBJetTags());
			btCSVBJet->push_back(tempJet->btag_combinedSecondaryVertexBJetTags() );
			btCSVBmvaJet->push_back(tempJet->btag_combinedSecondaryVertexMVABJetTags());  
		      }
			
	
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
		      delete btJPBJet;
		      delete btBJPBJet;
		      delete btCSVBJet;
		      delete btCSVBmvaJet;
			
		      
		      double SFval, SFerror;
		      if (isData || !scaleFactor){
			SFval = 1;
			SFerror = 0;
		      } else if (isTop){
			SFval = 0.95;
			SFerror = 0.03;
		      } else {
			SFval = 0.97;
			SFerror = 0.03;
		      } 
			
		      //Jet and b-tag selection
		      int nJetsBT = 0;
		      int nTightJetsBT = 0;
		      int nJets = 0;
		      bool bTagged = false;
		      int iJet = -5;
		      int iSF;
		      double tempSF = SFval;
		      if (SFminus) 	tempSF = SFval - SFerror;
		      if (SFplus) 	tempSF = SFval + SFerror;
		      int SFvalue = int(tempSF*100);
			
		      for (unsigned int i =0; i < selectedJets.size(); i ++){
			TRootJet* tempJet = (TRootJet*) selectedJets[i];
			TLorentzVector tJet(tempJet->Px(), tempJet->Py(), tempJet->Pz(), tempJet->Energy());
			if (tempJet->Pt() > 30 && fabs(tempJet->Eta()) < 2.5 && TMath::Min(fabs(lepton0.DeltaR(tJet)), fabs(lepton1.DeltaR(tJet))) > 0.3) {
			  nJets++;
			  iJet = i;
			  if (tempJet->btag_combinedSecondaryVertexBJetTags()> 0.679){
			    iSF = rand() % 100;
			    if (iSF < SFvalue || SFval == 1){
			      bTagged = true;
			      nJetsBT++;
			      nTightJetsBT++;
			    } 
			  } 
			} else if (tempJet->btag_combinedSecondaryVertexBJetTags()> 0.679 && fabs(tempJet->Eta()) < 2.5){
			  iSF = rand() % 100;
			  if (iSF < SFvalue  || SFval == 1) nJetsBT++;
			}
		      }
			
		
		      // Invariant mass in ee and mumu
		      if (pair.M() > 101 || pair.M() < 81 || mode == 0){
			cutflow->Fill(6, weight);
			cutflow_raw->Fill(6);
			// MET in ee and mumu
			if (met_pt > 30 || mode == 0){
			  cutflow->Fill(7, weight);
			  cutflow_raw->Fill(7);
			  salida3 << event->runId() << "\t" << event->lumiBlockId() << "\t" << event->eventId() << endl;

			  
			  // Filling all the regions
			  if (nJets !=0){
			    TRootJet* jet = (TRootJet*) selectedJets[iJet];
			    double Ht = lepton0.Pt() + lepton1.Pt() + jet->Pt() + met_pt; 
			    if (Ht > 160 || mode != 0){
			      if (nJets == 1 && nTightJetsBT == 1 && nJetsBT == 1 && bTagged)R->Fill(1, weight);
			      if (nJets == 1 && nTightJetsBT == 2)  R->Fill(2, weight);
			      if (nJets == 1 && nTightJetsBT > 0)  R->Fill(3, weight);
			      if (nJets == 1 && nTightJetsBT > 1)  R->Fill(4, weight);
			      if (nJets == 2 && nTightJetsBT == 0)  R->Fill(5, weight);
			      if (nJets == 2 && nTightJetsBT == 1)  R->Fill(6, weight);
			      if (nJets == 2 && nTightJetsBT == 2)  R->Fill(7, weight);
			      if (nJets == 2 && nTightJetsBT > 0)  R->Fill(8, weight);
			      if (nJets == 2 && nTightJetsBT > 1)  R->Fill(9, weight);
			      if (nJets > 1 && nTightJetsBT == 0)  R->Fill(10, weight);
			      if (nJets > 1 && nTightJetsBT == 1)  R->Fill(11, weight);
			      if (nJets > 1 && nTightJetsBT == 2)  R->Fill(12, weight);
			      if (nJets > 1 && nTightJetsBT !=0 )  R->Fill(13, weight);
			      if (nJets > 1 && nTightJetsBT > 1 )  R->Fill(14, weight);
			      if (nJets == 3 && nTightJetsBT ==0 )  R->Fill(15, weight);
			      if (nJets == 3 && nTightJetsBT ==1 )  R->Fill(16, weight);
			      if (nJets == 3 && nTightJetsBT ==2 )  R->Fill(17, weight);
			      if (nJets == 3 && nTightJetsBT ==3 )  R->Fill(18, weight);
			    }
			  }
			    
			  // Filling the signal region
			  if(nJets == 1){
			    salida4 << event->runId() << "\t" << event->lumiBlockId() << "\t" << event->eventId() << endl;

			    TRootJet* jet = (TRootJet*) selectedJets[iJet];
			    cutflow->Fill(8, weight);
			    cutflow_raw->Fill(8);
			    if (nJets == 1 && nTightJetsBT == 1 && nJetsBT == 1 && bTagged){
			      salida5 << event->runId() << "\t" << event->lumiBlockId() << "\t" << event->eventId() << endl;

			      cutflow->Fill(9,weight);
			      cutflow_raw->Fill(9);
			      double Ht = lepton0.Pt() + lepton1.Pt() + jet->Pt() + met_pt; 
			      if (Ht > 160 || mode != 0){
				cutflow->Fill(10, weight);
				cutflow_raw->Fill(10);
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
	
    
      delete jecUnc;
      delete jetTools;

     
      double scaler1 = cutflow->GetBinContent(2) ;
      if(scaler1<=0.0)
	scaler1=1.;
      cout << "--------------------------------------------------" << endl;
      cout << "[Results Normalized:] " <<  endl;
      cout << "All:       " <<  cutflow->GetBinContent(2) << " +/- "  << cutflow->GetBinError(2) << "\t = " << 100.*cutflow->GetBinContent(2)/scaler1 << " +/- "  << 100.*cutflow->GetBinError(2)/scaler1 << "%" << endl;
      cout << "HLT:       " <<  cutflow->GetBinContent(3) << " +/- "  << cutflow->GetBinError(3) <<  "\t = " << 100.*cutflow->GetBinContent(3)/scaler1 << " +/- "  << 100.*cutflow->GetBinError(3)/scaler1 << "%" <<endl;
      //cout << "PV:        " <<  cutflow->GetBinContent(4) << " +/- "  << cutflow->GetBinError(4) <<  "\t = " << 100.*cutflow->GetBinContent(4)/scaler1 << " +/- "  << 100.*cutflow->GetBinError(4)/scaler1 << "%" <<endl;
      cout << "Lep. Sel:  " <<  cutflow->GetBinContent(5) << " +/- "  << cutflow->GetBinError(5) <<  "\t = " << 100.*cutflow->GetBinContent(5)/scaler1 << " +/- "  << 100.*cutflow->GetBinError(5)/scaler1 << "%" <<endl;
      cout << "Lep. Veto: " <<  cutflow->GetBinContent(6) << " +/- "  << cutflow->GetBinError(6) <<  "\t = " << 100.*cutflow->GetBinContent(6)/scaler1 << " +/- "  << 100.*cutflow->GetBinError(6)/scaler1 << "%" <<endl;
      cout << "mll:       " <<  cutflow->GetBinContent(7) << " +/- "  << cutflow->GetBinError(7) <<  "\t = " << 100.*cutflow->GetBinContent(7)/scaler1 << " +/- "  << 100.*cutflow->GetBinError(7)/scaler1 << "%" << endl;
      cout << "MET:       " <<  cutflow->GetBinContent(8) << " +/- "  << cutflow->GetBinError(8) <<  "\t = " << 100.*cutflow->GetBinContent(8)/scaler1 << " +/- "  << 100.*cutflow->GetBinError(8)/scaler1 << "%" << endl;
      cout << "1 jet:     " <<  cutflow->GetBinContent(9) << " +/- "  << cutflow->GetBinError(9) <<  "\t = " << 100.*cutflow->GetBinContent(9)/scaler1 << " +/- "  << 100.*cutflow->GetBinError(9)/scaler1 << "%" <<endl;
      cout << "1 jet BT:  " <<  cutflow->GetBinContent(10) << " +/- " << cutflow->GetBinError(10) <<  "\t = " << 100.*cutflow->GetBinContent(10)/scaler1 << " +/- "  << 100.*cutflow->GetBinError(10)/scaler1 << "%" <<endl;
      cout << "Ht tW:     " <<  cutflow->GetBinContent(11) << " +/- " << cutflow->GetBinError(11) <<  "\t = " << 100.*cutflow->GetBinContent(11)/scaler1 << " +/- "  << 100.*cutflow->GetBinError(11)/scaler1 << "%" <<endl;
    
      
      cout << "--------------------------------------------------" << endl;
      cout << "[Results Raw:] " <<  endl;
      double scaler2 =  cutflow_raw->GetBinContent(2);
      if(scaler2 <=0.0)
	scaler2=1.;
      cout << "All:       " <<  cutflow_raw->GetBinContent(2) << " +/- "  << cutflow_raw->GetBinError(2) <<  "\t = " << 100.*cutflow_raw->GetBinContent(2)/scaler2 << " +/- "  << 100.*cutflow_raw->GetBinError(2)/scaler2 << "%" << endl;
      cout << "HLT:       " <<  cutflow_raw->GetBinContent(3) << " +/- "  << cutflow_raw->GetBinError(3) <<  "\t = " << 100.*cutflow_raw->GetBinContent(3)/scaler2 << " +/- "  << 100.*cutflow_raw->GetBinError(3)/scaler2 << "%" << endl;
      // cout << "PV:        " <<  cutflow_raw->GetBinContent(4) << " +/- "  << cutflow_raw->GetBinError(4) <<  "\t = " << 100.*cutflow_raw->GetBinContent(4)/scaler2 << " +/- "  << 100.*cutflow_raw->GetBinError(4)/scaler2 << "%" <<endl;
      cout << "Lep. Sel:  " <<  cutflow_raw->GetBinContent(5) << " +/- "  << cutflow_raw->GetBinError(5) <<  "\t = " << 100.*cutflow_raw->GetBinContent(5)/scaler2 << " +/- "  << 100.*cutflow_raw->GetBinError(5)/scaler2 << "%" <<endl;
      cout << "Lep. Veto: " <<  cutflow_raw->GetBinContent(6) << " +/- "  << cutflow_raw->GetBinError(6) <<  "\t = " << 100.*cutflow_raw->GetBinContent(6)/scaler2 << " +/- "  << 100.*cutflow_raw->GetBinError(6)/scaler2 << "%" <<endl;
      cout << "mll:       " <<  cutflow_raw->GetBinContent(7) << " +/- "  << cutflow_raw->GetBinError(7) <<  "\t = " << 100.*cutflow_raw->GetBinContent(7)/scaler2 << " +/- "  << 100.*cutflow_raw->GetBinError(7)/scaler2 << "%" <<endl;
      cout << "MET:       " <<  cutflow_raw->GetBinContent(8) << " +/- "  << cutflow_raw->GetBinError(8) <<  "\t = " << 100.*cutflow_raw->GetBinContent(8)/scaler2 << " +/- "  << 100.*cutflow_raw->GetBinError(8)/scaler2 << "%" <<endl;
      cout << "1 jet:     " <<  cutflow_raw->GetBinContent(9) << " +/- "  << cutflow_raw->GetBinError(9) <<  "\t = " << 100.*cutflow_raw->GetBinContent(9)/scaler2 << " +/- "  << 100.*cutflow_raw->GetBinError(9)/scaler2 << "%" <<endl;
      cout << "1 jet BT:  " <<  cutflow_raw->GetBinContent(10) << " +/- " << cutflow_raw->GetBinError(10) <<  "\t = " << 100.*cutflow_raw->GetBinContent(10)/scaler2 << " +/- "  << 100.*cutflow_raw->GetBinError(10)/scaler2 << "%" <<endl;
      cout << "Ht:        " <<  cutflow_raw->GetBinContent(11) << " +/- " << cutflow_raw->GetBinError(11) <<  "\t = " << 100.*cutflow_raw->GetBinContent(11)/scaler2 << " +/- "  << 100.*cutflow_raw->GetBinError(11)/scaler2 << "%" <<endl;
      
      cout << "--------------------------------------------------" << endl;
      cout << "[Jet Multiplicity Check:]" << endl;
      cout << "1 jet 1 tag: " << R->GetBinContent(2) << " +/- " << R->GetBinError(2) << "\t = " << 100.*R->GetBinContent(2)/scaler1 << " +/- "  << 100.*R->GetBinError(2)/scaler1 << "%" << endl;
      cout << "2 jet 1 tag: " << R->GetBinContent(7) << " +/- " << R->GetBinError(2) << "\t = " << 100.*R->GetBinContent(7)/scaler1 << " +/- "  << 100.*R->GetBinError(7)/scaler1 << "%" << endl;
      cout << "2 jet 2 tag: " << R->GetBinContent(8) << " +/- " << R->GetBinError(2) << "\t = " << 100.*R->GetBinContent(8)/scaler1 << " +/- "  << 100.*R->GetBinError(8)/scaler1 << "%" <<endl;
      cout << "--------------------------------------------------" << endl;
    
      fout->Write();
      fout->Close();
      
    }// dataset loop
  
  cout << "It took you " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
  
}

