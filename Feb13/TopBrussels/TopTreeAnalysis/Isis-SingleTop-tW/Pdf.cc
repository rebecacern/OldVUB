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
  

  //Run A/B separation
  bool RunA = false;
  bool RunB = false;

  string xmlfile ="twemu.xml";
  
  bool useTestXML=false;
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
      cout << "--xml <file.xml>: custom xml file. " << endl;
      cout << "--uncMETdown: Unclustered MET syst. Down " << endl;
      cout << "--NoPU: Do not apply pileup re-weighting" << endl;
      cout << "--NoSF: Do not apply b-tag scale factor" << endl;
      cout << "--RAW: Do not apply pileup re-weighting or b-tag scale factor" << endl;
      cout << "--RunA: Run over RunA only" << endl;
      cout << "--RunB: Run over RunB only" << endl;
      return 0;
    }
    if (argval=="--ee") mode = 2;
    if (argval=="--emu") mode = 0;
    if (argval=="--mumu") mode = 1;
    if (argval=="--xml") {useTestXML=true;xmlfile=argv[iarg+1];continue;}
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
    if (argval=="--RunA") RunA = true;
    if (argval=="--RunB") RunB = true;
  }   
  
  // Luminosity and xml files
  double lumi = 0;
  if      (mode == 0){ 	 lumi = 4904.338;	 if(!useTestXML) xmlfile ="twemu.xml";}
  else if (mode == 1){	 lumi = 4919.924;	 if(!useTestXML) xmlfile = "twmumu.xml";}
  else if (mode == 2){	 lumi = 4895.249;	 if(!useTestXML) xmlfile = "twee.xml";}
  if(useTestXML)
    std::cout << "using file: " << xmlfile << std::endl;
  
  double lumia = 2193.338;
  if (mode == 1) lumia = 2203.924 ;
  if (mode == 2) lumia = 2202.249 ;
  
  double lumib = 2711;
  if (mode == 1) lumib = 2716 ;
  if (mode == 2) lumib = 2693 ;
  

  if (RunA) lumi = lumia;
  else if (RunB) lumi = lumib;
  
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
      if (dataSetName == "data"){		sprintf(name, "data");  	xlweight = 1; 				isData = true;}
      else if (dataSetName == "tt"){            sprintf(name, "tt");            xlweight = lumi*163/3628285; 		isTop = true;} 
      else if (dataSetName == "twdr"){      	sprintf(name, "tw_dr");		xlweight = lumi*7.87/811897;} 
      else if (dataSetName == "atwdr"){ 	sprintf(name, "atw_dr");        xlweight = lumi*7.87/807505;} 
      else if (dataSetName == "twds"){ 	     	sprintf(name, "tw_ds");	        xlweight = lumi*7.87/793031;} 
      else if (dataSetName == "atwds"){    	sprintf(name, "atw_ds");        xlweight = lumi*7.87/807505;} 
      else if (dataSetName == "t"){    		sprintf(name, "t");        	xlweight = lumi*41.92/3888894; 		if (mode == 0) xlweight = lumi*41.92/3689490;} 
      else if (dataSetName == "at"){    	sprintf(name, "at");        	xlweight = lumi*22.65/1938263; 		if (mode == 0) xlweight = lumi*22.65/1934547;} 
      else if (dataSetName == "ts"){    	sprintf(name, "ts");        	xlweight = lumi*3.19/259378;}
      else if (dataSetName == "ats"){    	sprintf(name, "ats");        	xlweight = lumi*1.44/137556;} 
      else if (dataSetName == "wjets"){     	sprintf(name, "wjets");	        xlweight = lumi*31314/68090152; 	if (mode == 2) xlweight = lumi*31314/67987281; 	if (mode == 0) xlweight = lumi*31314/67981303;} 
      else if (dataSetName == "zjets"){         sprintf(name, "zjets");         xlweight = lumi*3048/35799272; 		if (mode == 0) xlweight = lumi*3048/35203706; 	if (mode == 2)  xlweight = lumi*3048/35606900;} 
      else if (dataSetName == "zjets_lowmll"){  sprintf(name, "zjets_lowmll");  xlweight = lumi*11908.83/31193930; 	if (mode == 0) xlweight = lumi*11908.83/30978317;}
      else if (dataSetName == "ww"){		sprintf(name, "ww");		xlweight = lumi*42.9/4223584;} 
      else if (dataSetName == "wz"){		sprintf(name, "wz");		xlweight = lumi*18.3/4242845; 		if (mode == 2) xlweight = lumi*18.3/4036353;} 
      else if (dataSetName == "zz"){		sprintf(name, "zz");		xlweight = lumi*7.67/4188621;} 
      else if (dataSetName == "qcd_mu"){    	sprintf(name, "qcd_mu");	xlweight = lumi*84679.3/24787499; 	if (mode == 0) xlweight = lumi*84679.3/24186153;if (mode == 2) xlweight = lumi*84679.3/24785478; } 
      
      //special files
      else if (dataSetName == "t_sup"){      		sprintf(name, "t_sup");		        xlweight = lumi*7.87/437736;}//
      else if (dataSetName == "tbar_sup"){      	sprintf(name, "tbar_sup");		xlweight = lumi*7.87/436710;}
      else if (dataSetName == "t_sdo"){      		sprintf(name, "t_sdo");			xlweight = lumi*7.87/436971;}
      else if (dataSetName == "tbar_sdo"){      	sprintf(name, "tbar_sdo");		xlweight = lumi*7.87/436991;}
      
      else if (dataSetName == "t_sup_ds"){      	sprintf(name, "t_sup_ds");		xlweight = lumi*7.87/436283;}
      else if (dataSetName == "tbar_sup_ds"){      	sprintf(name, "tbar_sup_ds");		xlweight = lumi*7.87/436378;}
      else if (dataSetName == "t_sdo_ds"){      	sprintf(name, "t_sdo_ds");		xlweight = lumi*7.87/436603;}
      else if (dataSetName == "tbar_sdo_ds"){      	sprintf(name, "tbar_sdo_ds");		xlweight = lumi*7.87/436751;} 
     
      
      else if (dataSetName == "tt_matchingup"){   	sprintf(name, "tt_matchingup");   	xlweight = lumi*163/1001708; isTop = true;}
      else if (dataSetName == "tt_matchingdown"){   	sprintf(name, "tt_matchingdown");   	xlweight = lumi*163/1043769; isTop = true;}
      else if (dataSetName == "tt_scaleup"){   		sprintf(name, "tt_scaleup");   		xlweight = lumi*163/911039; isTop = true;}
      else if (dataSetName == "tt_scaledown"){   	sprintf(name, "tt_scaledown");   	xlweight = lumi*163/926575; isTop = true;}
      
      
      else if (dataSetName == "data1"){		sprintf(name, "data1");  	xlweight = 1; isData = true; RunA = true;}
      else if (dataSetName == "data2"){		sprintf(name, "data2");  	xlweight = 1; isData = true; RunB = true;}
      
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
      else sprintf(rootFileName,"outputs/pdf_2j2t_%d_%s.root", mode, name);
      
      char myFile[300];
      sprintf(myFile,"textfiles/pdf_2j2t_%d_%s.txt", mode, name);
      ofstream salida(myFile); 
  
      
      // Objects
      vector < TRootVertex* > vertex;
      vector < TRootMuon* > init_muons;
      vector < TRootElectron* > init_electrons;
      vector < TRootJet* > init_jets_corrected;
      vector < TRootJet* > init_jets;
      vector < TRootMET* > mets;
      vector < TRootTrackMET* > trackmet;
      vector < TLorentzVector > allForTopoCalc;
      vector<TRootGenJet*> genjets;
       
   
    
      TFile *fout = new TFile (rootFileName, "RECREATE");
      
      TRootEvent* event = 0;
      
     
      // Branches of the output Tree
      double id1,id2, x1, x2, q, xlWeight;
     
      // Output Tree
      TTree* myTree = new TTree("myTree", "   ");
      
      myTree->Branch("xlWeight", &xlWeight, "xlWeight/D");
      myTree->Branch("x1", &x1, "x1/D");
      myTree->Branch("x2", &x2, "x2/D");
      myTree->Branch("id1", &id1, "id1/D");
      myTree->Branch("id2", &id2, "id2/D");
      myTree->Branch("q", &q, "q/D");
      
      //Pile-Up reweighting  
      Lumi3DReWeighting Lumi3DWeights;
      
      if (RunA) Lumi3DWeights = Lumi3DReWeighting("pileupHistos/MC_Fall11.root","pileupHistos/RunA.root", "pileup", "pileup");
      else if (RunB) Lumi3DWeights = Lumi3DReWeighting("pileupHistos/MC_Fall11.root","pileupHistos/RunB.root", "pileup", "pileup");
      else Lumi3DWeights =  Lumi3DReWeighting("pileupHistos/MC_Fall11.root","pileupHistos/RunAB.root", "pileup", "pileup");
      //else Lumi3DWeights =  Lumi3DReWeighting("pileupHistos/MC_Fall11.root","pileupHistos/PromptRecov6.root", "pileup", "pileup");
     

      if(PUsysDown) Lumi3DWeights.weight3D_init(0.92);	
      else if(PUsysUp) Lumi3DWeights.weight3D_init(1.08);
      else Lumi3DWeights.weight3D_init(1.0);

      // Initialize JEC factors
      vector<JetCorrectorParameters> vCorrParam;
      JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("../macros/JECFiles/START42_V17_AK5PFchs_Uncertainty.txt");
      
      // true means redo also the L1
      JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);
     
      /// Start message
      
      cout << endl;
      if (RunA)   cout << "[Info:] Running over 2011 RunA " << endl; 
      if (RunB)   cout << "[Info:] Running over 2011 RunB " << endl; 
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


	  if ((RunA && isData && event->runId() < 175860) || (RunB && isData && event->runId() >= 175860) || (isData && !RunA && !RunB) || !isData){
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
	   
	    // Systematics
	    //JES and JER
	    //Special = true;
	    if (!Special && !isData){
	      if(JERPlus)  jetTools->correctJetJER(init_jets, genjets, "plus");
	      else if(JERMinus) jetTools->correctJetJER(init_jets, genjets, "minus");
	      else jetTools->correctJetJER(init_jets, genjets, "nominal");
	      
	      if (JESPlus) jetTools->correctJetJESUnc(init_jets, mets[0], "plus");
	      else if (JESMinus) jetTools->correctJetJESUnc(init_jets, mets[0], "minus");
	    }  
	    
	    //Start selection
	    Selection selection(init_jets, init_muons, init_electrons, mets);
	    trackmet = treeLoader.LoadTrackMET(ievt);
	    
	   
	    
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
	    if(trigged){
	      if(isGoodPV){
		
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
		    
		    // Loose lepton veto
		    bool leptonVeto = false;
		    if 	  (mode == 0 && looseMuons.size()== 1 && looseElectrons.size() == 1) leptonVeto = true;
		    else if (mode == 1 && looseMuons.size()== 2 && looseElectrons.size() == 0) leptonVeto = true;
		    else if (mode == 2 && looseMuons.size()== 0 && looseElectrons.size() == 2) leptonVeto = true;
		    
		    if (leptonVeto) {
		      
		      // Low mll cut (all final states)
		      TLorentzVector pair = lepton0 + lepton1;   
		      if (pair.M() > 20){
			
			
			
			if (isRAW) weight = 1;
			xlWeight = weight;
			
			
			
			// Invariant mass in ee and mumu
			if (pair.M() > 101 || pair.M() < 81 || mode == 0){
			  if (TMath::Min(met_pt, trackmet[0]->Pt()) >= 30 || mode ==0){
			    double SFval, SFerror;
			    if (isData || !scaleFactor){
			      SFval = 1;
			      SFerror = 0;
			    } else if (isTop){
			      SFval = 0.956;
			      SFerror = 0.030;
			    } else {
			      SFval = 0.96;
			      SFerror = 0.04;
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
			      if (tempJet->Pt() > 30 && TMath::Min(fabs(lepton0.DeltaR(tJet)), fabs(lepton1.DeltaR(tJet))) > 0.3) {
				nJets++;
				iJet = i;
				if (tempJet->btag_simpleSecondaryVertexHighEffBJetTags() > 1.74){
				  iSF = rand() % 100;
				  if (iSF < SFvalue || SFval == 1){
				    bTagged = true;
				    nJetsBT++;
				    nTightJetsBT++;
				  } 
				} 
			      } else if (tempJet->btag_simpleSecondaryVertexHighEffBJetTags() > 1.74){
				iSF = rand() % 100;
				if (iSF < SFvalue  || SFval == 1) nJetsBT++;
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
				if (nJets == 2 && nTightJetsBT == 2){
				
				xlWeight = weight;
				id1 = event->idParton1();
				id2 = event->idParton2();
				x1 = event->xParton1();
				x2 = event->xParton2();
				q = event->factorizationScale();
				myTree->Fill();
				salida << weight << ", " << x1 << ", " << x2 << ", " << q << ", " << id1 << ", " << id2 << ", " << name << ", " << mode << endl;
				
				}
			
				
				/*
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
				  */
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
	    }
	  }
	} // event loop
      
      // cleanup
      allForTopoCalc.clear();
      allForTopoCalc.resize(0);
      
      if(jetTools) delete jetTools;   
       fout->Write();
      fout->Close();
   }  
  cout << "It took you " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
  
}

