// rebeca@cern.ch
//
//Ready for the next

#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

//user code
#include "../TopBrussels/TopTreeProducer/interface/TRootRun.h"
#include "../TopBrussels/TopTreeProducer/interface/TRootEvent.h"
#include "../TopTreeAnalysis/Selection/interface/SelectionTable.h"
#include "../TopTreeAnalysis/Tools/interface/PlottingTools.h"
#include "../TopTreeAnalysis/Tools/interface/JetTools.h"
#include "../TopTreeAnalysis/Tools/interface/MultiSamplePlot.h"
#include "../TopTreeAnalysis/Tools/interface/TTreeLoader.h"
#include "../TopTreeAnalysis/Tools/interface/AnalysisEnvironmentLoader.h"
#include "../TopTreeAnalysis/Content/interface/AnalysisEnvironment.h"
#include "../TopTreeAnalysis/Content/interface/Dataset.h"
#include "../TopTreeAnalysis/MCInformation/interface/MCWeighter.h"
#include "../TopTreeAnalysis/Selection/interface/ElectronPlotter.h"
#include "../TopTreeAnalysis/Selection/interface/MuonPlotter.h"
#include "../TopTreeAnalysis/Selection/interface/JetPlotter.h"
#include "../TopTreeAnalysis/Selection/interface/VertexPlotter.h"
#include "../TopTreeAnalysis/Tools/interface/MVATrainer.h"
#include "../TopTreeAnalysis/Tools/interface/MVAComputer.h"
#include "../TopTreeAnalysis/Tools/interface/TopologyWorker.h"
#include "../TopTreeAnalysis/MCInformation/interface/ResolutionFit.h"
#include "../TopTreeAnalysis/Reconstruction/interface/JetCorrectorParameters.h"
#include "../TopTreeAnalysis/Reconstruction/interface/JetCorrectionUncertainty.h"
#include "../TopTreeAnalysis/MCInformation/interface/LumiReWeighting.h"
#include "../TopTreeAnalysis/macros/Style.C"

using namespace std;
using namespace TopTree;
using namespace reweight;

int main(int argc, char* argv[]) {
  
  setMyStyle();
  
  clock_t start = clock();
  
  //modes: 0 emu, 1mumu, 2ee 
  int  mode = 0; 
  
  bool reweightPU = true;
  
  bool scaleFactor = true;
  double SFval = 0.95; 
  
  bool isRAW = false;
  
  bool doJESsysts=false;
  bool doMETJESToo=false;
  float jetenergyscaleshift=0;
  
  bool doJERsysts=false;
  float jetenergyresoshift=0;
  
  bool SFplus = false;
  bool SFminus = false;
  
  bool unclusteredUp = false;
  bool unclusteredDown = false;
  
  bool SystSamples = false;
  bool Spring11 = false;
  
  bool PUsysUp = false;
  bool PUsysDown = false;
  
  string xmlfile ="tWconfig.xml";
  std::string jetenergyscaleshiftstring="none";
  std::string jetenergyscalelocation="JECFiles/Jec11V2_db_AK5PFchs_Uncertainty.txt";
  
  for(int iarg = 0; iarg < argc && argc>1 ; iarg++){
    std::string argval=argv[iarg];
    if(argval=="--help" || argval =="--h"){
      cout << "--ee: Di-Electron" << endl;
      cout << "--emu: Electron-Muon" << endl;
      cout << "--mumu: Di-Muon" << endl;
      cout << "--uncMETUp: Unclustered MET syst. Up " << endl;
      cout << "--uncMETDown: Unclustered MET syst. Down " << endl;
      cout << "--JESplus: JES sys +1 sigma" << endl;
      cout << "--JESminus: JES sys -1 sigma" << endl;
      cout << "--JESandMETplus: JES sys +1 sigma MET included" << endl;
      cout << "--JESandMETminus: JES sys -1 sigma MET included" << endl;
      cout << "--JER0: JER 0.0" << endl;
      cout << "--JER1: JER 0.1" << endl;
      cout << "--JER2: JER 0.2" << endl;
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
    if (argval=="--JESplus") {doJESsysts=true;	jetenergyscaleshiftstring="plus";}
    if (argval=="--JESminus") {doJESsysts=true;	jetenergyscaleshiftstring="minus";}
    if (argval=="--JESandMETplus") {doJESsysts=true;	doMETJESToo=true;	jetenergyscaleshiftstring="plus";}
    if (argval=="--JESandMETminus") {doJESsysts=true;	doMETJESToo=true;	jetenergyscaleshiftstring="minus";} 
    if (argval=="--JER0") { doJERsysts= true; jetenergyresoshift = 0;}
    if (argval=="--JER1") { doJERsysts= true; jetenergyresoshift = 1;}
    if (argval=="--JER2") { doJERsysts= true; jetenergyresoshift = 2;}
    if (argval=="--SFplus") SFplus = true;
    if (argval=="--SFminus") SFminus = true;
    if (argval=="--NoPU") reweightPU = false;
    if (argval=="--NoSF") scaleFactor = false;
    if (argval=="--RAW") {reweightPU = false; scaleFactor = false; isRAW = true;}
  }   
  
  
  double lumi = 0;
  if      (mode == 0){ 	 lumi = 2121.307;	 xmlfile ="twemu.xml";}
  else if (mode == 1){	 lumi = 2110.25;	 xmlfile = "twmumu.xml";}
  else if (mode == 2){	 lumi = 2096.434;	 xmlfile = "twee.xml";}
  
  
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
  
  /// Start the analyisis
  for (unsigned int d = 0; d < datasets.size (); d++)
    {
      treeLoader.LoadDataset (datasets[d], anaEnv);
      string dataSetName = datasets[d]->Name();
      
      bool isData = false;
      char name[100];
      double xlweight;
      
      // cross sections and weights
      if (dataSetName == "data"){		sprintf(name, "data");  	xlweight = 1; isData = true;}
      else if (dataSetName == "tt"){		sprintf(name, "tt");  		xlweight = lumi*163/2976314;}
      else if (dataSetName == "twdr"){      	sprintf(name, "tw_dr");		xlweight = lumi*7.87/814385;}
      else if (dataSetName == "atwdr"){ 	sprintf(name, "atw_dr");        xlweight = lumi*7.87/809978 ;}
      else if (dataSetName == "twds"){ 	     	sprintf(name, "tw_ds");	        xlweight = lumi*7.87/795373;}
      else if (dataSetName == "atwds"){    	sprintf(name, "atw_ds");        xlweight = lumi*7.87/787627 ;}
      else if (dataSetName == "t"){    		sprintf(name, "t");        	xlweight = lumi*41.92/3900054;}
      else if (dataSetName == "at"){    	sprintf(name, "at");        	xlweight = lumi*22.65/1944782 ;}
      else if (dataSetName == "ts"){    	sprintf(name, "ts");        	xlweight = lumi*3.19/259968;}
      else if (dataSetName == "ats"){    	sprintf(name, "ats");        	xlweight = lumi*1.44/137979 ;}
      else if (dataSetName == "wjets"){     	sprintf(name, "wjets");	        xlweight = lumi*31314/80506476;}
      else if (dataSetName == "zjets"){     	sprintf(name, "zjets");	        xlweight = lumi*3048/35023307;} 
      else if (dataSetName == "dymumu"){	sprintf(name, "dymumu");	xlweight = lumi*1666/2146121;}
      else if (dataSetName == "dyee"){		sprintf(name, "dyee");		xlweight = lumi*1666/2242963;}
      else if (dataSetName == "dytautau"){	sprintf(name, "dytautau");	xlweight = lumi*1666/2030174;}
      else if (dataSetName == "ww"){		sprintf(name, "ww");		xlweight = lumi*42.9/4208152;}
      else if (dataSetName == "wz"){		sprintf(name, "wz");		xlweight = lumi*18.3/4157378;}
      else if (dataSetName == "zz"){		sprintf(name, "zz");		xlweight = lumi*7.67/4164379;}
      else if (dataSetName == "qcd_mu"){    	sprintf(name, "qcd_mu");	xlweight = lumi*84679.3/25079892;} 
      
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
      
      char rootFileName[100];
      sprintf(rootFileName,"outputs/var_%d_%s.root", mode, name);

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
      
      TH1F* R = new TH1F( "R", " ", 40,  0, 40 );;
      R->Sumw2();
      
      
      TH1F* cutflow = new TH1F("cutflow", " ", 31,  -0.5, 30.5 );
      TH1F* cutflow_raw = new TH1F("cutflow_raw", " ", 31,  -0.5, 30.5 );
      
      TH1F* histo_met = new TH1F("histo_met", " ", 100,  0, 200 );
      TH1F* histo_promet = new TH1F( "histo_promet", " ", 100,  0, 200 );
      TH1F* histo_mll = new TH1F( "histo_mll", " ", 200,  0, 400 );
      TH1F* histo_njets = new TH1F( "histo_njets", " ", 10,  -0.5, 9.5 );
      TH1F* histo_njetsbt = new TH1F( "histo_njetsbt", " ", 10,  -0.5, 9.5 );
      TH1F* histo_ptsys = new TH1F( "histo_ptsys", " ", 100,  0, 200 );
      TH1F* histo_ht = new TH1F( "histo_ht", " ", 300,  0, 600 );
      TH1F* histo_ht_nomet = new TH1F( "histo_ht_nomet", " ", 300,  0, 600 );
      TH1F* histo_pt_max = new TH1F( "histo_pt_max", " ", 100,  0, 200 );
      TH1F* histo_pt_min = new TH1F( "histo_pt_min", " ", 100,  0, 200 );
      TH1F* histo_pt_leading = new TH1F( "histo_pt_leading", " ", 100,  0, 200 );
      TH1F* histo_oblateness = new TH1F( "histo_oblateness", " ", 100,  0, 1 );
      TH1F* histo_sphericity = new TH1F( "histo_sphericity", " ", 100,  0, 1 );
      TH1F* histo_aplanarity = new TH1F( "histo_aplanarity", " ", 100,  0, 0.5 );
      TH1F* histo_njetW = new TH1F( "histo_njetW", " ", 100,  0, 10 );
      TH1F* histo_sqrts = new TH1F( "histo_sqrts", " ", 300,  0, 600 );
      TH1F* histo_deltaphi = new TH1F( "histo_deltaphi", " ", 100,  0, 4 );
      TH1F* histo_deltaR = new TH1F( "histo_deltaR", " ", 100,  0, 4 );
      TH1F* histo_deltaeta = new TH1F( "histo_deltaeta", " ", 100,  0, 4 );
      TH1F* histo_deltaphiclosemet = new TH1F( "histo_deltaphiclosemet", " ", 100,  0, 4 );
      TH1F* histo_deltaphiclosejet = new TH1F( "histo_deltaphiclosejet", " ", 100,  0, 4 );
      TH1F* histo_deltaphijetmet = new TH1F( "histo_deltaphijetmet", " ", 100,  0, 4 );
      TH1F* histo_eta_ptmax = new TH1F( "histo_eta_ptmax", " ", 200,  -5, 5 );
      TH1F* histo_eta_ptmin = new TH1F( "histo_eta_ptmin", " ", 200,  -5, 5 );
      TH1F* histo_eta_all = new TH1F( "histo_eta_all", " ", 200,  -5, 5 );
      TH1F* histo_eta_jet = new TH1F( "histo_eta_jet", " ", 200,  -5, 5 );
      TH1F* histo_phi_jet = new TH1F( "histo_phi_jet", " ", 200,  -5, 5 );
      TH1F* histo_phi_met = new TH1F( "histo_phi_met", " ", 200,  -5, 5 );
      TH1F* histo_total_pt = new TH1F( "histo_total_pt", " ", 100,  0, 200 );
      TH1F* histo_total_eta = new TH1F( "histo_total_eta", " ", 200,  -5, 5 );
      TH1F* histo_total_phi = new TH1F( "histo_total_phi", " ", 200,  -5, 5 );
      
      TH1F* histo_njetsextra = new TH1F( "histo_njetsextra", " ", 10,  -0.5, 9.5 );
      TH1F* histo_ptjetsextra = new TH1F( "histo_ptjetsextra", " ", 100,  0, 200 );
      TH1F* histo_etajetsextra = new TH1F( "histo_etajetsextra", " ", 200,  -5, 5 );
      TH2F* histo_extrajets = new TH2F("histo_extrajets", " ", 100,  0, 200, 200,  -5, 5 );
      TH1F* histo_mll_lepclosejet = new TH1F( "histo_mll_lepclosejet", " ", 200,  0, 400 );
      TH1F* histo_mll_lepfarjet = new TH1F( "histo_mll_lepfarjet", " ", 200,  0, 400 );
      
      TH1F* histo_metminusht = new TH1F( "histo_metminusht", " ", 200,  -100, 100 );
      TH1F* histo_metminuspt = new TH1F( "histo_metminuspt", " ", 200,  -100, 100 );
      
      cutflow->Sumw2();
      cutflow_raw->Sumw2();
      
      histo_met->Sumw2();
      histo_promet->Sumw2();
      histo_mll->Sumw2();
      histo_njets->Sumw2();
      histo_njetsbt->Sumw2();
      histo_ptsys->Sumw2();
      histo_ht->Sumw2();
      histo_ht_nomet->Sumw2();
      histo_pt_max->Sumw2();
      histo_pt_min->Sumw2();
      histo_pt_leading->Sumw2();
      histo_oblateness->Sumw2();
      histo_sphericity->Sumw2();
      histo_aplanarity->Sumw2();
      histo_njetW->Sumw2();
      histo_sqrts->Sumw2();
      histo_deltaphi->Sumw2();
      histo_deltaeta->Sumw2();
      histo_deltaR->Sumw2();
      histo_deltaphiclosemet->Sumw2();
      histo_deltaphiclosejet->Sumw2();
      histo_deltaphijetmet->Sumw2();
      
      histo_eta_ptmax->Sumw2();
      histo_eta_ptmin->Sumw2(); 
      histo_eta_all->Sumw2();
      
      histo_eta_jet->Sumw2();
      
      histo_phi_jet->Sumw2();
      histo_phi_met->Sumw2();
      
      histo_total_pt->Sumw2();
      histo_total_eta->Sumw2();
      histo_total_phi->Sumw2();
      
      histo_njetsextra->Sumw2();
      histo_ptjetsextra->Sumw2();
      histo_etajetsextra->Sumw2();
      histo_extrajets->Sumw2();
      
      histo_mll_lepclosejet->Sumw2();
      histo_mll_lepfarjet->Sumw2();
      
      histo_metminusht->Sumw2();
      histo_metminuspt->Sumw2();
        
     
      //Pile-Up reweighting  
      LumiReWeighting LumiWeights = LumiReWeighting("PileUpReweighting/pileup_WJets_36bins.root", "PileUpReweighting/pileup_2011Data_UpToRun177515.root", "pileup2", "pileup"); 
    
      PoissonMeanShifter PShiftUp_ = PoissonMeanShifter(0.6); // PU-systematic
      PoissonMeanShifter PShiftDown_ = PoissonMeanShifter(-0.6); // PU-systematic  
      
      // Initialize JEC factor uncertainties
      
      vector<JetCorrectorParameters> vCorrParam;  
      if(name == "data"){// data
	JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("JECFiles/Jec11V2_db_AK5PFchs_L2L3Residual.txt"); // hardcoded here!!
	vCorrParam.push_back(*ResJetCorPar);
      }
      JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(jetenergyscalelocation); // hardcoded at start of main!!!
      JetTools *jetTools = new JetTools(vCorrParam, jecUnc, false);   
      
      /// Start message
      
      cout << endl;
      cout << "[Info:] output rootfile named " << rootFileName << endl; 
      cout << "[Info:] mode = " << mode << ", lumi: " <<  lumi << " pb, sample: " << name << ", base weight: " << xlweight << endl;
      
      if (doJERsysts || doJESsysts || doMETJESToo ||  SFplus || SFminus || unclusteredUp || unclusteredDown || SystSamples || Spring11 || !reweightPU || !scaleFactor || PUsysUp || PUsysDown) 
	{
	  cout << "[Warning:] Non-standard options, ignore if you did it conciously" << endl;
	  if (doJERsysts) cout << "[Warning:] JER systematics on, shift" << jetenergyresoshift << endl;
	  if (doJESsysts){ cout << "[Warning:] JES systematics on"; if (doMETJESToo) cout << " with MET " << endl; else cout << endl;}
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
          
	  double weight = xlweight;
          
	  // Pile-Up re-weighting
	  if (reweightPU && !isData){
	    if(!Spring11){
	      float avPU = ( (float)event->nPu(-1) + (float)event->nPu(0) + (float)event->nPu(+1) ) / 3.; // average in 3 BX!!!, as recommended
	      weight *= LumiWeights.ITweight(avPU);
	      if (PUsysUp) weight *= PShiftUp_.ShiftWeight(avPU);
	      else if (PUsysDown) weight *= PShiftDown_.ShiftWeight(avPU);
	    } else {
	      weight *= LumiWeights.ITweight((float)event->nPu(0));
	      if (PUsysUp) weight *= PShiftUp_.ShiftWeight(event->nPu(0));
	      else if (PUsysDown) weight *= PShiftDown_.ShiftWeight(event->nPu(0));
	    }
	  }
	  
	  //Trigger
          bool trigged = true;
          
          //Start selection
	  Selection selection(init_jets, init_muons, init_electrons, mets);
	  
	  // Systematics
	  if(doJESsysts){
	    if(doMETJESToo)
	      jetTools->correctJetJESUnc(init_jets, mets[0] ,jetenergyscaleshiftstring);
	    else
	      jetTools->correctJetJESUnc(init_jets, jetenergyscaleshiftstring);
	  }
	  if(doJERsysts){
	    // match the jets to generator jets and shift <jetenergyresoshift> sigma according to definitions on https://twiki.cern.ch/twiki/bin/view/CMS/TopLeptonPlusJets2010Systematics    
	    if( name != "data")// NOT data so MC ;-)
	      {
		vector<TRootGenJet*> genjets = treeLoader.LoadGenJet(ievt);
		if(genjets.size()==0)
		  cerr << "!!!!trying JER correction but no genjets in event!!!!!"  << endl;
		// loop over all genjets and match to reco jets (in deltar<0.4)
		for(size_t irecojet=0; irecojet<init_jets.size(); irecojet++){ 
		  if(init_jets[irecojet]->Pt()<10)// don't bother with low-pT jets.
		    continue;
		  float minDR=1000;
		  size_t genmatchedjet=0;
		  for(size_t igenjet=0; igenjet<genjets.size(); igenjet++){
		    if(genjets[igenjet]->Pt()<15)
		      continue;
		    if(init_jets[irecojet]->DeltaR(*genjets[igenjet])<0.4 && init_jets[irecojet]->DeltaR(*genjets[igenjet])<minDR){
		      minDR=init_jets[irecojet]->DeltaR(*genjets[igenjet]);
		      genmatchedjet=igenjet;
		    }
		  }
		  if(minDR<=0.4){// found match
		    // officially should also correct MET, but not implemented yet:
		    // now calculating how much the jet should be shifted.
		    float jercorrfac=0.1;// basic correction
		    float extracorr = 0.1 /*value*/ * jetenergyresoshift;
		    float ptdiff = (init_jets[irecojet]->Pt() - genjets[genmatchedjet]->Pt())*(jercorrfac+extracorr);
		    float ptscale=max(0.0,(init_jets[irecojet]->Pt()+ptdiff)/init_jets[irecojet]->Pt());
		    // and scale the jet:
		    //                            cout << "scaling the jet by " << ptscale << " !!!" ;
		    //                            cout << "old value : " << init_jets[irecojet]->Pt() << " ";
		    jetTools->scaleJet(init_jets[irecojet], ptscale);
		    //                            cout << "new value: " << init_jets[irecojet]->Pt() << endl;
		  }
		}
	      }
	  }
	  
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
	  
	  
	  // Cut Flow
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
	      vector<TRootJet*> selectedJets = selection.GetSelectedJets(20,2.4,selectedElectrons,true);                    
	      
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
                  
		  bool leptonVeto = false;
		  if 	  (mode == 0 && looseMuons.size()== 1 && looseElectrons.size() == 1) leptonVeto = true;
		  else if (mode == 1 && looseMuons.size()== 2 && looseElectrons.size() == 0) leptonVeto = true;
		  else if (mode == 2 && looseMuons.size()== 0 && looseElectrons.size() == 2) leptonVeto = true;
                  
		  if (leptonVeto) {
		    cutflow->Fill(5, weight);
		    cutflow_raw->Fill(5);
                    
		    TLorentzVector pair = lepton0 + lepton1;   
		    if (pair.M() > 20){
		      
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
			   if (tempJet->Pt() > 30 && TMath::Min(fabs(lepton0.DeltaR(tJet)), fabs(lepton1.DeltaR(tJet))) > 0.3) {
			     nJets++;
			     if (iJet == -5) iJet = i;
			     if (tempJet->btag_simpleSecondaryVertexHighEffBJetTags() > 1.74){
			       bTagged = true;
			       nJetsBT++;
			     } 
			   } else if (tempJet->btag_simpleSecondaryVertexHighEffBJetTags() > 1.74) nJetsBT++;
			 }
		       } else {
			 //// Regular SF
			 if (SFvalue < 101){
			   for (int i =0; i < selectedJets.size(); i ++){
			     TRootJet* tempJet = (TRootJet*) selectedJets[i];
			     TLorentzVector tJet(tempJet->Px(), tempJet->Py(), tempJet->Pz(), tempJet->Energy());
			     if (tempJet->Pt() > 30 && TMath::Min(fabs(lepton0.DeltaR(tJet)), fabs(lepton1.DeltaR(tJet))) > 0.3) {
			       nJets++;
			       if (iJet == -5) iJet = i;
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
			   for (int i =0; i < selectedJets.size(); i ++){
			     TRootJet* tempJet = (TRootJet*) selectedJets[i];
			     TLorentzVector tJet(tempJet->Px(), tempJet->Py(), tempJet->Pz(), tempJet->Energy());
			     if (tempJet->Pt() > 30 && TMath::Min(fabs(lepton0.DeltaR(tJet)), fabs(lepton1.DeltaR(tJet))) > 0.3) {
			       nJets++;
			       if (iJet == -5) iJet = i;
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
                      
		     
		       if (nJets!=0){
			 TRootJet* jet = (TRootJet*) selectedJets[iJet];
			 TLorentzVector lvjet(jet->Px(), jet->Py(), jet->Pz(), jet->Energy());
			 double ptSysPx = lepton0.Px() + lepton1.Px() + jet->Px() + met_px;
			 double ptSysPy = lepton0.Py() + lepton1.Py() + jet->Py() + met_py;
			 double ptSys = sqrt(ptSysPx*ptSysPx + ptSysPy*ptSysPy);
			 double Ht = lepton0.Pt() + lepton1.Pt() + jet->Pt() + met_pt; 
			 double Ht_nomet = lepton0.Pt() + lepton1.Pt() + jet->Pt(); 
			 
			 double phipairmet_t = 0;
			 double pi_m = 3.1416/2;
			 phipairmet_t = pi_m;
			 
			 TVector3 vmet(met_px, met_py, 0);
			 TVector3 vjet(jet->Px(), jet->Py(), jet->Pz());
			 
			 TVector3 m0(lepton0.Px(),lepton0.Py(), lepton0.Pz());
			 TVector3 m1(lepton1.Px(),lepton1.Py(), lepton1.Pz()); 
			 if (fabs(m0.DeltaPhi(vmet)) < phipairmet_t) phipairmet_t = fabs(m0.DeltaPhi(vmet));
			 if (fabs(m1.DeltaPhi(vmet)) < phipairmet_t) phipairmet_t = fabs(m1.DeltaPhi(vmet));
			 
			 double promet = met_pt*sin(phipairmet_t);
			 if (phipairmet_t == pi_m) promet = met_pt;
			 
			 histo_met->Fill(met_pt, weight);
			 histo_promet->Fill(promet , weight);
			 histo_mll->Fill(pair.M() , weight);
			 histo_njets->Fill(nJets , weight);
			 histo_njetsbt->Fill(nJetsBT , weight);
			 histo_ptsys->Fill(ptSys , weight);
			 histo_ht->Fill(Ht , weight);
			 histo_ht_nomet->Fill(Ht_nomet , weight);
			 histo_pt_max->Fill(TMath::Max(lepton0.Pt(), lepton1.Pt()) , weight);
			 histo_pt_min->Fill(TMath::Min(lepton0.Pt(), lepton1.Pt()) , weight);
			 histo_pt_leading->Fill(jet->Pt() , weight);  
			 histo_deltaphi->Fill(fabs(lepton0.DeltaPhi(lepton1)), weight);
			 histo_deltaeta->Fill(fabs(lepton0.Eta() - lepton1.Eta()), weight);
			 
			 histo_deltaR->Fill(fabs(lepton0.DeltaR(lepton1)), weight);
			 histo_deltaphiclosemet->Fill(TMath::Min(fabs(m0.DeltaPhi(vmet)), fabs(m1.DeltaPhi(vmet))) , weight);
			 histo_deltaphiclosejet->Fill(TMath::Min(fabs(m0.DeltaPhi(vjet)), fabs(m1.DeltaPhi(vjet))) , weight);
			 histo_deltaphijetmet->Fill(fabs(vjet.DeltaPhi(vmet)), weight);
			 histo_eta_ptmax->Fill(lepton0.Eta(), weight);
			 histo_eta_ptmin->Fill(lepton1.Eta(), weight);
			 histo_eta_all->Fill(lepton0.Eta(), weight);
			 histo_eta_all->Fill(lepton1.Eta(), weight);
			 
			 histo_eta_jet->Fill(vjet.Eta(), weight);
			 
			 histo_phi_met->Fill(vmet.Phi(), weight);
			 histo_phi_jet->Fill(vjet.Phi(), weight);
			 
			 TLorentzVector total = lepton0 + lepton1 + lvjet;
			 
			 histo_total_pt->Fill(total.Pt(), weight);
			 histo_total_eta->Fill(total.Eta(), weight);
			 histo_total_phi->Fill(total.Phi(), weight);
			 
			 int njetsextra = 0;
			 for(unsigned int i=0; i<init_jets.size(); i++){
			   TRootJet* tempJet = (TRootJet*) init_jets[i];
			   if (tempJet != jet) {
			     njetsextra++;
			     histo_ptjetsextra->Fill(tempJet->Pt(), weight);
			     histo_etajetsextra->Fill(tempJet->Eta(), weight);
			     histo_extrajets->Fill(tempJet->Pt(), tempJet->Eta(), weight);
			   }
			   
			 }
			 histo_njetsextra->Fill(njetsextra, weight);
			 
			 TLorentzVector pair0 = lepton0 + lvjet;
			 TLorentzVector pair1 = lepton1 + lvjet;
			 
			 if (fabs(m0.DeltaPhi(vjet)) <  fabs(m1.DeltaPhi(vjet))){
			   histo_mll_lepclosejet->Fill(pair0.M(), weight);
			   histo_mll_lepfarjet->Fill(pair1.M(), weight);
			 } else {
			   histo_mll_lepclosejet->Fill(pair1.M(), weight);
			   histo_mll_lepfarjet->Fill(pair0.M(), weight);
			 }
			 
			 histo_metminusht->Fill(met_pt - Ht_nomet, weight);
			 histo_metminuspt->Fill(met_pt - total.Pt(),weight);
			 
			 allForTopoCalc.push_back(lepton0);
			 allForTopoCalc.push_back(lepton1);
			 for(int ijet=0;ijet<selectedJets.size(); ijet++){
			   TRootJet* tempJet = (TRootJet*) selectedJets[ijet];
			   TLorentzVector tempLV(tempJet->Px(), tempJet->Py(), tempJet->Pz(), tempJet->E());
			   allForTopoCalc.push_back(tempLV);
			 }
			 
			 TopologyWorker::TopologyWorker topoWorkerNoBoost(false); // set the bool to true to un-boost events
			 topoWorkerNoBoost.setPartList(allForTopoCalc,allForTopoCalc);
			 topoWorkerNoBoost.setVerbose(true);
			 
			 histo_oblateness->Fill(topoWorkerNoBoost.oblateness(), weight);
			 histo_sphericity->Fill(topoWorkerNoBoost.get_sphericity(), weight);
			 histo_aplanarity->Fill(topoWorkerNoBoost.get_aplanarity(), weight);
			 histo_njetW->Fill(topoWorkerNoBoost.get_njetW(), weight);
			 histo_sqrts->Fill(topoWorkerNoBoost.get_sqrts(), weight);
			 allForTopoCalc.clear();
			 allForTopoCalc.resize(0);
		       }
		       
		       
		    }
		  }
		}
	      }
	    }      
	  }
	} // event loop
      
      // cleanup
      
      if(jetTools) delete jetTools;   
      fout->Write();
      fout->Close();
      
    }// dataset loop
  
  cout << "It took you " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
  
}

