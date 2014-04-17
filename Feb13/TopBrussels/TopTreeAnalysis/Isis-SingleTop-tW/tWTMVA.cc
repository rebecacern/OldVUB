// rebeca@cern.ch
//
// Fall11 round
// TMVA flat-tuples
// adding run and event number (fblekman@cern.ch)
//

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
  
  // Train and test separation
  bool isTrain = true;
  bool isTest = false;
  bool controlRegion1 = false;
  bool controlRegion2 = false;
  
  //Pile-up reweighting
  bool reweightPU = true;
  
  //b-tag scale factor
  bool scaleFactor = true;
  double SFval = 0.95; 
  
  // Met Type I correction	
  bool metTypeI = false; 
  
  //Systematic calculations 
  bool JESPlus=false;
  bool JESMinus=false;
  bool JERPlus=false;
  bool JERMinus=false;
  
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
      cout << "--uncMETUp: Unclustered MET syst. Up " << endl;
      cout << "--uncMETDown: Unclustered MET syst. Down " << endl;
      cout << "--JESPlus: JES sys +1 sigma MET included" << endl;
      cout << "--JESMinus: JES sys -1 sigma MET included" << endl;
      cout << "--JERPlus: JER +" << endl;
      cout << "--JERMinus: JER -" << endl;
      cout << "--SFplus: SF up +10% syst" << endl;
      cout << "--SFminus: SF down -10% syst" << endl;
      cout << "--NoPU: Do not apply pileup re-weighting" << endl;
      cout << "--NoSF: Do not apply b-tag scale factor" << endl;
      cout << "--RAW: Do not apply pileup re-weighting or b-tag scale factor" << endl;
      return 0;
    }
    if (argval=="--ee") mode = 2;
    if (argval=="--emu") mode = 0;
    if (argval=="--mumu") mode = 1;
    if (argval=="--uncMETUp") unclusteredUp = true;
    if (argval=="--uncMETDown") unclusteredDown = true;
    if (argval=="--JESPlus") {JESPlus = true;}
    if (argval=="--JESMinus") {JESMinus = true;}
    if (argval=="--JERPlus") {JERPlus = true;}
    if (argval=="--JERMinus") {JERMinus = true;}
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
      bool isSignal = false;
      char name[100];
      double xlweight;
      
      // cross sections and weights
      if (dataSetName == "data"){		sprintf(name, "data");  	xlweight = 1; isData = true;}
      else if (dataSetName == "tt"){            sprintf(name, "tt");            xlweight = lumi*163/3160707;} //
      else if (dataSetName == "tt2l"){		sprintf(name, "tt2l");  	xlweight = lumi*17.10/8576584;} //
      else if (dataSetName == "twdr"){      	sprintf(name, "tw_dr");		xlweight = lumi*7.87/813743; isSignal = true;} //
      else if (dataSetName == "atwdr"){ 	sprintf(name, "atw_dr");        xlweight = lumi*7.87/689462; isSignal = true;} //
      else if (dataSetName == "twds"){ 	     	sprintf(name, "tw_ds");	        xlweight = lumi*7.87/794802; isSignal = true;} //
      else if (dataSetName == "atwds"){    	sprintf(name, "atw_ds");        xlweight = lumi*7.87/784764; isSignal = true;} //
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
      char dname[100];
      if (isTest){
	sprintf(dname, "test");
	if (SFplus) sprintf(dname, "btagup_test");
	if (SFminus) sprintf(dname, "btagdown_test");
      } else {
	sprintf(dname, "train");
	if (SFplus) sprintf(dname, "btagup_train");
	if (SFminus) sprintf(dname, "btagdown_train");
      }
      
      char rootFileName[100];
      sprintf(rootFileName,"outputs/tmva_06_%s_%d_%s.root", dname, mode, name);
      if (controlRegion1) sprintf(rootFileName,"outputs/tmva_06_%s_2j1t_%d_%s.root", dname, mode, name);
      if (controlRegion2) sprintf(rootFileName,"outputs/tmva_06_%s_2j2t_%d_%s.root", dname, mode, name);
      
      
      // Objects
      vector < TRootVertex* > vertex;
      vector < TRootMuon* > init_muons;
      vector < TRootElectron* > init_electrons;
      vector < TRootJet* > init_jets;
      vector < TRootMET* > mets;
      vector < TLorentzVector > allForTopoCalc;
      
      
      TFile *fout = new TFile (rootFileName, "RECREATE");
      
      TRootEvent* event = 0;
      
      
      // Histograms
      map<string,TH1F*> histo1D;
      map<string,TH2F*> histo2D;
      
      
      // Branches of the output Tree
      
      double luminosity;
      double weight, weightnopu, savetheweight;
      
      int num_event, num_lumisec, num_run;

      int npu;
      int nvertex;
      
      double metpt, metpx,  metpy, metpro, metphi;
      
      double lep0pt, lep1pt, lep0eta, lep1eta, lep0phi, lep1phi;
      
      double deltaphileps, deltaetaleps, deltarleps;
      
      double philepmetclose, philepmetfar, rlepmetclose, rlepmetfar;
      double philepjetclose, philepjetfar, rlepjetclose, rlepjetfar;
      double phijetmet, rjetmet;
      
      int njets, njetsbt;
      
      double jetpt, jeteta, jetphi, jetbt;
      
      double oblateness, sphericity, aplanarity, njetw, sqrts;
      
      double mll, ptsys, ht, htnomet, ptsysnomet, metminusptsysnomet;
      
      bool data, signal, background;
      
      int nloosejet, nloosejetbtag, nloosejetbtagloose, nloosejetforward;
      double loosejetPtTot, loosejetPtMax;    
      
      
      TTree* myTree = new TTree("myTree", "   ");
      
      myTree->Branch("num_event",&num_event,"num_event/I");
      myTree->Branch("num_lumisec",&num_lumisec,"num_lumisec/I");
      myTree->Branch("num_run",&num_run,"num_run/I");
 
      myTree->Branch("weight", &weight, "weight/D");
      myTree->Branch("weightnopu", &weightnopu, "weightnopu/D");
      myTree->Branch("savetheweight", &savetheweight, "savetheweight/D");
      myTree->Branch("luminosity", &luminosity, "luminosity/D");
      
      myTree->Branch("data", &data, "data/O");
      myTree->Branch("signal", &signal, "signal/O");
      myTree->Branch("background", &background, "background/O");
      
      myTree->Branch("npu", &npu, "npu/I");
      myTree->Branch("nvertex", &nvertex, "nvertex/I");
      
      myTree->Branch("metpt", &metpt, "metpt/D");
      myTree->Branch("metpx", &metpx, "metpx/D");
      myTree->Branch("metpy", &metpy, "metpy/D");
      myTree->Branch("metpro", &metpro, "metpro/D");
      myTree->Branch("metphi", &metphi, "metphi/D");
      
      myTree->Branch("lep0pt", &lep0pt, "lep0pt/D");
      myTree->Branch("lep0eta", &lep0eta, "lep0eta/D");
      myTree->Branch("lep0phi", &lep0phi, "lep0phi/D");

      myTree->Branch("lep1pt", &lep1pt, "lep1pt/D");
      myTree->Branch("lep1eta", &lep1eta, "lep1eta/D");
      myTree->Branch("lep1phi", &lep1phi, "lep1phi/D");
      
      //myTree->Branch("njets", &njets, "njets/I");
      //myTree->Branch("njetsbt", &njetsbt, "njetsbt/I");

      myTree->Branch("jetpt", &jetpt, "jetpt/D");
      myTree->Branch("jeteta", &jeteta, "jeteta/D");
      myTree->Branch("jetphi", &jetphi, "jetPhi/D");
      myTree->Branch("jetbt", &jetbt, "jetbt/D");
      
      myTree->Branch("mll", &mll, "mll/D");
      myTree->Branch("ptsys", &ptsys, "ptsys/D");
      myTree->Branch("ht", &ht, "ht/D");
      myTree->Branch("htnomet", &htnomet, "htnomet/D");
      myTree->Branch("metminusptsysnomet", &metminusptsysnomet, "metminusptsysnomet/D");
      myTree->Branch("ptsysnomet", &ptsysnomet, "ptsysnomet/D");

      myTree->Branch("oblateness", &oblateness, "oblateness/D");
      myTree->Branch("sphericity", &sphericity, "sphericity/D");
      myTree->Branch("aplanarity", &aplanarity, "aplanarity/D");
      myTree->Branch("njetw", &njetw, "njetw/D");
      myTree->Branch("sqrts", &sqrts, "sqrts/D");
      
      myTree->Branch("deltaphileps", &deltaphileps, "deltaphileps/D");
      myTree->Branch("deltaetaleps", &deltaetaleps, "deltaetaleps/D");
      myTree->Branch("deltarleps", &deltarleps, "deltarleps/D");
      
      myTree->Branch("philepmetclose", &philepmetclose, "philepmetclose/D");
      myTree->Branch("philepmetfar", &philepmetfar, "philepmetfar/D");
      myTree->Branch("rlepmetclose", &rlepmetclose, "rlepmetclose/D");
      myTree->Branch("rlepmetfar", &rlepmetfar, "rlepmetfar/D");

      myTree->Branch("philepjetclose", &philepjetclose, "philepjetclose/D");
      myTree->Branch("philepjetfar", &philepjetfar, "philepjetfar/D");
      myTree->Branch("rlepjetclose", &rlepjetclose, "rlepjetclose/D");
      myTree->Branch("rlepjetfar", &rlepjetfar, "rlepjetfar/D");

      myTree->Branch("phijetmet", &phijetmet, "phijetmet/D");
      myTree->Branch("rjetmet", &rjetmet, "rjetmet/D");
      
      myTree->Branch("nloosejet", &nloosejet, "nloosejet/I");
      myTree->Branch("nloosejetbtag", &nloosejetbtag, "nloosejetbtag/I");
      myTree->Branch("nloosejetbtagloose", &nloosejetbtagloose, "nloosejetbtagloose/I");
      myTree->Branch("nloosejetforward", &nloosejetforward, "nloosejetforward/I");
      
      myTree->Branch("loosejetPtTot", &loosejetPtTot, "loosejetPtTot/D");
      myTree->Branch("loosejetPtMax", &loosejetPtMax, "loosejetPtMax/D");
  
  
      //Pile-Up reweighting  
      Lumi3DReWeighting Lumi3DWeights = Lumi3DReWeighting("../macros/PileUpReweighting/pileup_MC_Fall11.root","../macros/PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root", "pileup", "pileup");
      Lumi3DWeights.weight3D_init(1.0);
      
      // Initialize JEC factors
      vector<JetCorrectorParameters> vCorrParam;
      JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("../macros/JECFiles/Jec11V2_db_AK5PFchs_Uncertainty.txt");
      JetTools *jetTools = new JetTools(vCorrParam, jecUnc, false);
      
      /// Start message
      
      cout << endl;
      cout << "[Info:] output rootfile named " << rootFileName << endl; 
      cout << "[Info:] mode = " << mode << ", lumi: " <<  lumi << " pb, sample: " << name << ", base weight: " << xlweight << endl;
      if (isTest) cout << "[Info:] Test sample" << endl;
      else if (isTrain)  cout << "[Info:] Train sample" << endl;

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
      
      
      for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
	{
	  
	  if(ievt%500 == 0) std::cout<<"Processing the "<<ievt<<"th event" <<flush<<"\r";
	  event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);
	  int val = event->eventId();

	  if ( (isTest && TMath::Even(val)) || (isTrain && TMath::Odd(val)) || isData){
	    
	    // Weight given by the theoretical cross-section and lumi
	    double weight = xlweight;
	    if (!isData) weight *= 2;
	    
	    // Pile-Up re-weighting
	    if (reweightPU && !isData && !Special){
	      double lumiWeight3D = 1.0;
	      lumiWeight3D = Lumi3DWeights.weight3D(event->nPu(-1),event->nPu(0),event->nPu(+1));
	      weight *= lumiWeight3D;
	    }
	    
	    //Trigger
	    int currentRun = event->runId();
	    num_run = currentRun;
	    num_event = event->eventId();
	    num_lumisec = event->nb();

	    bool itrigger = false;
	    bool isecondtrigger = false;
	    
	    bool trigged = true;
	    
	    
	    // Correct MET Type I
	    if (metTypeI && !Special) jetTools->correctMETTypeOne(init_jets,mets[0]);  //Size of mets is never larger than 1 !!
	    
	    
	    // Systematics
	    //JES and JER
	    //Special = true;
	    if (!Special && !isData){
	      if (JESPlus) jetTools->correctJetJESUnc(init_jets, mets[0], "plus");
	      else if (JESMinus) jetTools->correctJetJESUnc(init_jets, mets[0], "minus");
	      
	      vector<TRootGenJet*> genjets = treeLoader.LoadGenJet(ievt);
	      
	      if(JERPlus)  jetTools->correctJetJER(init_jets, genjets, "plus");
	      else if(JERMinus) jetTools->correctJetJER(init_jets, genjets, "minus");
	      else jetTools->correctJetJER(init_jets, genjets, "nominal");
	    }  
	    
	    //Start selection
	    Selection selection(init_jets, init_muons, init_electrons, mets);
	  
	    // PV cut (useless)
	    bool isGoodPV = isGoodPV = selection.isPVSelected(vertex, 4,24,2.);  
	    
	    
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
		selection.setJetCuts(20.,5,0.01,1.,0.98,0.3,0.1);
		selection.setDiElectronCuts(20,2.5,0.15,0.02,1);
		selection.setLooseElectronCuts(15,2.5,0.2);
		selection.setDiMuonCuts(20,2.4,0.15,10,0.02);
		selection.setLooseMuonCuts(10,2.5,0.2);
		
		//Select Objects 
		vector<TRootElectron*> selectedElectrons = selection.GetSelectedDiElectrons(vertex[0]);
		vector<TRootMuon*> selectedMuons = selection.GetSelectedDiMuons();
		vector<TRootElectron*> looseElectrons = selection.GetSelectedLooseElectrons(true);
		vector<TRootMuon*> looseMuons = selection.GetSelectedLooseMuons();
		vector<TRootJet*> selectedJets = selection.GetSelectedJets(20,5,true);                    
		
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
			
			//Jet and b-tag selection
			int nJetsBT = 0;
			int nJets = 0;
			bool bTagged = false;
			int iJet = -5;
			int iSF;
			double tempSF = SFval;
			if (SFminus) 	tempSF = SFval - SFval*10/100;
			if (SFplus) 	tempSF = SFval + SFval*10/100;
			int SFvalue = tempSF*100 + 1;
			
			if (isData || !scaleFactor){
			  for (int i =0; i < selectedJets.size(); i ++){
			    TRootJet* tempJet = (TRootJet*) selectedJets[i];
			    TLorentzVector tJet(tempJet->Px(), tempJet->Py(), tempJet->Pz(), tempJet->Energy());
			    if (tempJet->Pt() > 30 && fabs(tempJet->Eta()) < 2.4 && TMath::Min(fabs(lepton0.DeltaR(tJet)), fabs(lepton1.DeltaR(tJet))) > 0.3) {
			      nJets++;
			      if (iJet == -5) iJet = i;
			      if (tempJet->btag_simpleSecondaryVertexHighEffBJetTags() > 1.74){
				bTagged = true;
				nJetsBT++;
			      } 
			    } else if (fabs(tempJet->Eta()) < 2.4 && tempJet->btag_simpleSecondaryVertexHighEffBJetTags() > 1.74) nJetsBT++;
			  }
			} else {
			  //// Regular SF
			  if (SFvalue < 101){
			    for (int i =0; i < selectedJets.size(); i ++){
			      TRootJet* tempJet = (TRootJet*) selectedJets[i];
			      TLorentzVector tJet(tempJet->Px(), tempJet->Py(), tempJet->Pz(), tempJet->Energy());
			      if (tempJet->Pt() > 30  && fabs(tempJet->Eta()) < 2.4 && TMath::Min(fabs(lepton0.DeltaR(tJet)), fabs(lepton1.DeltaR(tJet))) > 0.3) {
				nJets++;
				if (iJet == -5) iJet = i;
				if (tempJet->btag_simpleSecondaryVertexHighEffBJetTags() > 1.74 && fabs(tempJet->Eta()) < 2.4){
				  iSF = rand() % 101;
				  if (iSF < SFvalue ){
				    bTagged = true;
				    nJetsBT++;
				  } 
				} 
			      } else if (fabs(tempJet->Eta()) < 2.4 && tempJet->btag_simpleSecondaryVertexHighEffBJetTags() > 1.74){
				iSF = rand() % 101;
				if (iSF < SFvalue ) nJetsBT++;
			      }
			    }
			  } else {
			    //// Large SF
			    for (int i =0; i < selectedJets.size(); i ++){
			      TRootJet* tempJet = (TRootJet*) selectedJets[i];
			      TLorentzVector tJet(tempJet->Px(), tempJet->Py(), tempJet->Pz(), tempJet->Energy());
			      if (tempJet->Pt() > 30 && fabs(tempJet->Eta()) < 2.4 && TMath::Min(fabs(lepton0.DeltaR(tJet)), fabs(lepton1.DeltaR(tJet))) > 0.3) {
				nJets++;
				if (iJet == -5) iJet = i;
				if (fabs(tempJet->Eta()) < 2.4 && tempJet->btag_simpleSecondaryVertexHighEffBJetTags() > 1.74){
				  bTagged = true;
				  nJetsBT++;
				} else {
				  iSF = rand() % 101;
				  if (iSF < abs(100 - SFvalue)){
				    nJetsBT++;
				    bTagged = true;
				  }
				}
			      } else if (fabs(tempJet->Eta()) < 2.4 && tempJet->btag_simpleSecondaryVertexHighEffBJetTags() > 1.74){
				nJetsBT++;
			      } else {
				iSF = rand() % 101;
				if (iSF < abs(100 - SFvalue)) nJetsBT++;
			      }
			    }
			  }
			}
			
			
			//
			if ((nJets == 1 && bTagged) || (controlRegion1 && nJets == 2 && nJetsBT == 1) ||  (controlRegion2 && nJets == 2 && nJetsBT == 2) ){
			  
			  
			  nloosejetbtag = nJetsBT -1;
			  nloosejetbtagloose = 0;
			  nloosejetforward = 0;
			  loosejetPtTot  = 0;
			  loosejetPtMax = 0;
			  
			  for (int i =0; i < selectedJets.size(); i ++){
			    TRootJet* tempJet = (TRootJet*) selectedJets[i];
			    TLorentzVector tJet(tempJet->Px(), tempJet->Py(), tempJet->Pz(), tempJet->Energy());
			    if (tempJet->Pt() > 30 && fabs(tempJet->Eta()) < 2.4 && TMath::Min(fabs(lepton0.DeltaR(tJet)), fabs(lepton1.DeltaR(tJet))) > 0.3) {}
			    else{
			      loosejetPtTot +=tempJet->Pt();
			      if ( tempJet->Pt() > loosejetPtMax) loosejetPtMax = tempJet->Pt();
			      if (fabs(tempJet->Eta()) >= 2.4 )nloosejetforward ++;
			      if (fabs(tempJet->Eta()) < 2.4 && tempJet->btag_trackCountingHighEffBJetTags() > 1.7)  nloosejetbtagloose++;
			    }
			  }
			  
			  
			 
			  // if (mode != 1 || isSignal || (mode == 1 && !isSignal && nFill < 800)){
			  TRootJet* jet = (TRootJet*) selectedJets[iJet];
			  double ptSysPx = lepton0.Px() + lepton1.Px() + jet->Px() + met_px;
			  double ptSysPy = lepton0.Py() + lepton1.Py() + jet->Py() + met_py;
			  double ptSys = sqrt(ptSysPx*ptSysPx + ptSysPy*ptSysPy);
			  double Ht = lepton0.Pt() + lepton1.Pt() + jet->Pt() + met_pt; 
			  
			  ptsysnomet = sqrt((lepton0.Px()+lepton1.Px()+jet->Px())*(lepton0.Px()+lepton1.Px()+jet->Px()) + 
					    (lepton0.Py()+lepton1.Py()+jet->Py())*(lepton0.Py()+lepton1.Py()+jet->Py()) );
			  
			  // do calculations as far as topology goes:
			  allForTopoCalc.push_back(lepton0);
			  allForTopoCalc.push_back(lepton1);
			  for(int ijet=0;ijet<selectedJets.size(); ijet++){
			    TRootJet* tempJet = (TRootJet*) selectedJets[ijet];
			    TLorentzVector tempLV(tempJet->Px(), tempJet->Py(), tempJet->Pz(), tempJet->E());
			    allForTopoCalc.push_back(tempLV);
			  }
			  // std::cout << "creating topology worker " << std::endl;
			  TopologyWorker::TopologyWorker topoWorkerNoBoost(false); // set the bool to true to un-boost events
			  topoWorkerNoBoost.setPartList(allForTopoCalc,allForTopoCalc);
			  topoWorkerNoBoost.setVerbose(true);
			  // see for other methods: header file of TopologyWorker
			  
			  double phipairmet_t = 0;
			  double pi_m = 3.1416/2;
			  phipairmet_t = pi_m;
			  
			  TVector3 vmet(met_px, met_py, 0);
			  TVector3 vjet(jet->Px(), jet->Py(), jet->Pz());
			  TVector3 l0(lepton0.Px(), lepton0.Py(),lepton0.Pz());
			  TVector3 l1(lepton1.Px(), lepton1.Py(),lepton1.Pz()); 
			  if (TMath::Min(fabs(l0.DeltaPhi(vmet)), fabs(l1.DeltaPhi(vmet))) < phipairmet_t) phipairmet_t = TMath::Min(fabs(l0.DeltaPhi(vmet)), fabs(l1.DeltaPhi(vmet)));
			  
			  double promet = met_pt*sin(phipairmet_t);
			  if (phipairmet_t == pi_m) promet = met_pt;
			  
			  metphi = vmet.Phi();
			  metminusptsysnomet = met_pt - ptsysnomet;
			  
			  //Filling the Tree
			  luminosity = lumi;
			  
			  if (isRAW) weight = 1;
			  weightnopu = xlweight;
			  savetheweight = weight;
			  
			  data = isData;
			  signal = isSignal;
			  background = (!isSignal)*(!isData);
			  
			  npu = event->nPu(0);
			  nvertex = vertex.size();
			  
			  metpt = met_pt;
			  metpx = met_px;
			  metpy = met_py;
			  metpro = promet;
			  
			  lep0pt = lepton0.Pt();
			  lep1pt = lepton1.Pt();
			  
			  lep0eta = lepton0.Eta();
			  lep1eta = lepton1.Eta();
			  
			  lep0phi = lepton0.Phi();
			  lep1phi = lepton1.Phi();
			  
			  njets = nJets;
			  njetsbt = nJetsBT;
			  
			  jetpt = jet->Pt();
			  jeteta = jet->Eta();
			  jetphi = jet->Phi();
			  jetbt = jet->btag_simpleSecondaryVertexHighEffBJetTags();
			  
			  ptsys =  ptSys;
			  ht = Ht; 
			  htnomet = Ht - met_pt;
			  mll = pair.M();
			  
			  deltaphileps = lepton0.DeltaPhi(lepton1);
			  deltaetaleps = lepton0.Eta() - lepton1.Eta();
			  deltarleps = lepton0.DeltaR(lepton1);
			  
			  philepmetclose = TMath::Min(fabs(l0.DeltaPhi(vmet)), fabs(l1.DeltaPhi(vmet)));
			  philepmetfar =  TMath::Max(fabs(l0.DeltaPhi(vmet)), fabs(l1.DeltaPhi(vmet)));
			  
			  rlepmetclose = TMath::Min(fabs(l0.DeltaR(vmet)), fabs(l1.DeltaR(vmet)));
			  rlepmetfar =  TMath::Max(fabs(l0.DeltaR(vmet)), fabs(l1.DeltaR(vmet)));
			  
			  philepjetclose = TMath::Min(fabs(l0.DeltaPhi(vjet)), fabs(l1.DeltaPhi(vjet)));
			  philepjetfar =  TMath::Max(fabs(l0.DeltaPhi(vjet)), fabs(l1.DeltaPhi(vjet)));
			  
			  rlepjetclose = TMath::Min(fabs(l0.DeltaR(vjet)), fabs(l1.DeltaR(vjet)));
			  rlepjetfar =  TMath::Max(fabs(l0.DeltaR(vjet)), fabs(l1.DeltaR(vjet)));
			  
			  phijetmet = vmet.DeltaPhi(vjet);
			  rjetmet = vmet.DeltaR(vjet);
			  
			  oblateness = topoWorkerNoBoost.oblateness();
			  sphericity = topoWorkerNoBoost.get_sphericity();
			  aplanarity = topoWorkerNoBoost.get_aplanarity();
			  njetw = topoWorkerNoBoost.get_njetW();
			  sqrts = topoWorkerNoBoost.get_sqrts(); 
			  
			  // 
			  allForTopoCalc.clear();
			  allForTopoCalc.resize(0);
			  
			 
			  
			  myTree->Fill();
			  
			}
			///
		        
		      }
		    }
		  }
		}
	      }      
	    }
	  } //train/test loop
	} // event loop
      
      // cleanup
      allForTopoCalc.clear();
      allForTopoCalc.resize(0);
      
      if(jetTools) delete jetTools;   
      
      fout->Write();
      fout->Close();
      
    }// dataset loop
  
  cout << "It took you " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
  
}
