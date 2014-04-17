///////////////////////////
///// TODO & COMMENTS /////
/////////////////////////// 

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
#include "../JESMeasurement/interface/FullKinFit.h"
#include "../JESMeasurement/interface/JetCombiner.h"
#include "Style.C"

using namespace std;
using namespace TopTree;

int main() {

  clock_t start = clock();

  cout << "*************************************************" << endl;
  cout << " Running The Lepton+Jets RefSel for the SyncEx ! " << endl;
  cout << "*************************************************" << endl;
  
  //SetStyle if needed
  //setTDRStyle(); 
  setMyStyle();

  /////////////////////
  // Configuration
  /////////////////////

    unsigned int lepton = 1; // 0 = muon+jets 1= electron +jets

  //xml file
  string xmlfile ="../config/myRefSelconfig.xml";

  //Output ROOT file
  string rootFileName ("Esel.root");

  //Configuration output format
  TTree *configTree = new TTree("configTree","configuration Tree");
  TClonesArray* tcdatasets = new TClonesArray("Dataset",1000);
  configTree->Branch("Datasets","TClonesArray",&tcdatasets);
  TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000);
  configTree->Branch("AnaEnv","TClonesArray",&tcAnaEnv);

  ////////////////////////////////////
  /// AnalysisEnvironment  
  ////////////////////////////////////

  AnalysisEnvironment anaEnv;
  cout<<"Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile.c_str());
  new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);
  int verbose = anaEnv.Verbose;
  
  /////////////////////
  // Load Datasets
  /////////////////////

  TTreeLoader treeLoader;
  cout << " - Load datasets ..." << endl;
  vector < Dataset* > datasets;
  treeLoader.LoadDatasets (datasets, xmlfile.c_str());
  for(unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);

  float oldLuminosity = anaEnv.Luminosity;	// in 1/pb
  float Luminosity = oldLuminosity;
  for (unsigned int d = 0; d < datasets.size (); d++)
    {
      //cout << "luminosity of dataset "<< d << " is " << datasets[d]->EquivalentLumi() << endl;
      if(Luminosity > datasets[d]->EquivalentLumi() ) Luminosity = datasets[d]->EquivalentLumi();
      string dataSetName = datasets[d]->Name();
      if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
	{
	  Luminosity = datasets[d]->EquivalentLumi();
	}
    }
  if(Luminosity != oldLuminosity) cout << "changed analysis environment luminosity to "<< Luminosity << endl;
  

  //vector of objects
  cout << " - Variable declaration ..." << endl;
  vector < TRootVertex* > vertex;
  vector < TRootMuon* > init_muons;
  vector < TRootElectron* > init_electrons;
  vector < TRootJet* > init_jets;
  vector < TRootMET* > mets;
      vector < TRootGenJet* > genjets;

  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");
  //Global variable
  TRootEvent* event = 0;
  
  ////////////////////////////////////
  /// Normal Plots (TH1F* and TH2F*)
  ////////////////////////////////////

  map<string,TH1F*> histo1D;
  map<string,TH2F*> histo2D;

   ////////////////////////////////////
  /// Selection table
  ////////////////////////////////////

  vector<string> CutsSelecTable;

  if (lepton == 0) {
    cout << " - Preparing the cutflow table for mu+jets refSelV4" << endl;
    CutsSelecTable.push_back(string("initial"));
    CutsSelecTable.push_back(string("trigged"));
    CutsSelecTable.push_back(string("Good PV"));
    CutsSelecTable.push_back(string("1 selected muon"));
    CutsSelecTable.push_back(string("Veto 2nd muon"));
    CutsSelecTable.push_back(string("Veto electron"));
    char LabelNJets[100];
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-3);
    CutsSelecTable.push_back(string(LabelNJets));
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-2);
    CutsSelecTable.push_back(string(LabelNJets));
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-1);
    CutsSelecTable.push_back(string(LabelNJets));
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets);
    CutsSelecTable.push_back(string(LabelNJets));

    CutsSelecTable.push_back("$\\geq$ 1 CSVM tag");
  } else if (lepton == 1){
    cout << " - Preparing the cutflow table for e+jets refSelV4" << endl;
    CutsSelecTable.push_back(string("initial"));
    CutsSelecTable.push_back(string("trigged"));
    CutsSelecTable.push_back(string("Good PV"));
    CutsSelecTable.push_back(string("1 isolated electron"));
    CutsSelecTable.push_back(string("Loose muon veto"));
    CutsSelecTable.push_back(string("Dilepton veto"));
    CutsSelecTable.push_back(string("Conversion rejection"));
    CutsSelecTable.push_back(string("Partnertrack veto"));
    char LabelNJets[100];
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-3);
    CutsSelecTable.push_back(string(LabelNJets));
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-2);
    CutsSelecTable.push_back(string(LabelNJets));
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-1);
    CutsSelecTable.push_back(string(LabelNJets));
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets);
    CutsSelecTable.push_back(string(LabelNJets));

    CutsSelecTable.push_back("$\\geq$ 1 CSVM tag");

  } else {

    cout << "!!!!Lepton should equal 0 or 1, exiting" << endl;

    exit(1);

  }

  if (verbose > 0)
    cout << " - CutsSelectionTable instantiated ... size: " << CutsSelecTable.size() <<endl;
  SelectionTable selecTable(CutsSelecTable, datasets);
  selecTable.SetLuminosity(Luminosity);
  if (verbose > 0)
    cout << " - SelectionTable instantiated ..." << endl;


  vector<JetCorrectorParameters> vCorrParam;

    JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters("JECFiles/START42_V17_AK5PFchs_L3Absolute.txt");
    JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters("JECFiles/START42_V17_AK5PFchs_L2Relative.txt");
    JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters("JECFiles/START42_V17_AK5PFchs_L1FastJet.txt");

    //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
    vCorrParam.push_back(*L1JetPar);
    vCorrParam.push_back(*L2JetPar);
    vCorrParam.push_back(*L3JetPar);

 
  //vCorrParam.push_back(new JetCorrectorParameters());

    JetCorrectionUncertainty *jecUnc = NULL;

    // true means redo also the L1
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);
    
    ////////////////////////////////////////////////////////////
    // CREATE OUTPUT FILE AND TTREE FOR STORAGE OF THE NTUPLE //
    ////////////////////////////////////////////////////////////

  ////////////////////////////////////
  //	Loop on datasets
  ////////////////////////////////////
  
  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
  for (unsigned int d = 0; d < datasets.size (); d++)
  {
    if (verbose > 1)
      cout << "   Dataset " << d << ": " << datasets[d]->Name () << "/ title : " << datasets[d]->Title () << endl;
    if (verbose > 1)
      std::cout<<"      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
		
		//open files and load
    cout<<"LoadEvent"<<endl;
    treeLoader.LoadDataset (datasets[d], anaEnv);

    //seleTable.Fill(d,0, datasets[d]->Xsection() * datasets[d]->EquivalentLumi() );

    string dataSetName = datasets[d]->Name();

    ////////////////////////////////////
    //	Loop on events
    ////////////////////////////////////
    
    //nEvents[d] = 0;
    int itrigger = -1, previousRun = -1;
    if (verbose > 1)
      cout << "	Loop over events " << endl;
    for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
//    for (unsigned int ievt = 0; ievt < 50000; ievt++)
    {
      
      //nEvents[d]++;

      //cout << "anaEnv.JetType -> " << anaEnv.JetType << " --- " <<  endl;
      //cout << "anaEnv.METType -> " << anaEnv.METType << " --- " <<  endl;


      //if(ievt%500 == 0)
      //std::cout<<"Processing the "<<ievt<<"th event" <<flush<<"\r";
        
      //load event

      event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);

      //jetTools->unCorrectJets(init_jets,false);

      //cout << mets.size() << endl;

      // check with genEvent which ttbar channel it is
      
      /*if(dataSetName.find("TTbarJets") == 0)  {
	//cout << "Loading GenEvent" << endl;
        TRootGenEvent* genEvt = treeLoader.LoadGenEvent(ievt,false);
        if( ! genEvt->isSemiLeptonic(TRootGenEvent::kMuon) ) {
	  continue;
	}
	}*/

      selecTable.Fill(d,0,1.);

      /////////////
      // TRIGGER //
      /////////////
      
      bool trigged=false;
      
      if (lepton == 0) {
	//if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") {
	
	int currentRun = event->runId();
	//if(previousRun != currentRun) {
	//  previousRun = currentRun;
	//cout << currentRun << endl;
	itrigger = treeLoader.iTrigger ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2", currentRun,0);
	//itrigger = treeLoader.iTrigger ("HLT_IsoMu17_eta2p1_TriCentralPFJet30_v2", currentRun,0);
	//itrigger = treeLoader.iTrigger ("HLT_IsoMu17_eta2p1_TriCentralPFJet30_v4", currentRun);
	//}

	trigged = treeLoader.EventTrigged (itrigger);
	
	//} else 
	//trigged=true;
      }
      
      else if (lepton == 1) {

	int currentRun = event->runId();
	
	itrigger = treeLoader.iTrigger ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v8", currentRun,0);
	  
	trigged = treeLoader.EventTrigged (itrigger);

	//}
	//else*/
	//  trigged=true;
      }
      
      //jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho(),false); //last boolean: isData (needed for L2L3Residual...)

      if(! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) ) {
        genjets = treeLoader.LoadGenJet(ievt,false);
        sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
      }

      //jetTools->correctJetJER(init_jets, genjets, mets[0], "nominal",false);

      /////////////////////////////
      //   Selection
      /////////////////////////////
      
      //Declare selection instance    
      Selection selection(init_jets, init_muons, init_electrons, mets);
      
      bool isGoodPV = isGoodPV = selection.isPVSelected(vertex, 4,24,2.);  

      //*******************************//
      //***** Muon+Jets RefSel V4 *****//
      //*******************************//

      if (lepton == 0) {
    
	// FILL THE SELECTION TABLE //
	
	if(trigged){
	  selecTable.Fill(d,1,1.);
 	  
	  if(isGoodPV){
	    selecTable.Fill(d,2,1.);
	    
	    vector<TRootJet*> selectedJets;
	    vector<TRootMuon*> selectedMuons;
	    
	    selection.setJetCuts(35.,2.5,0.01,1.,0.98,0.3,0.1);
	    selection.setMuonCuts(20,2.1,0.125,0,0.02,0.3,1,1,5);
	    selection.setLooseMuonCuts(10,2.5,0.2);
	    selection.setLooseElectronCuts(20,2.5,0.2,0.); // semiMu looseElectron cuts

	    selectedJets = selection.GetSelectedJets(true);
	    selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);

	    vector<TRootMuon*> vetoMuons = selection.GetSelectedLooseMuons();
	    vector<TRootElectron*> vetoElectrons = selection.GetSelectedLooseElectrons(false);

	    //cout << event->runId() << " : " << event->eventId() << " : " << event->lumiBlockId() << endl;

	    //cout << "" << endl;
	    //cout << "Event " << ievt << " nMuons: " << init_muons.size() << endl;
	    //cout << "" << endl;
	    //cout << "muon reliso " << 

	    //cout << event->eventId() << " : " << event->runId() << " : " << event->lumiBlockId() << endl;
	    /*if (event->eventId() == 95595184 && event->lumiBlockId() == 318715) {
	      cout << "Special event " << ievt << " -> vetoMuons.size() = " << vetoMuons.size() << endl;
	      cout << "Special event " << ievt << " -> initMuons.size() = " << init_muons.size() << endl;
	      for (unsigned int i=0; i<init_muons.size(); i++) {
		if (init_muons[i]->puChargedHadronIso() == -9999.) init_muons[i]->setPuChargedHadronIso(0);

		float reliso = (init_muons[i]->chargedHadronIso()+init_muons[i]->neutralHadronIso()+init_muons[i]->photonIso())/init_muons[i]->Pt();
		float dZ = fabs(init_muons[i]->vz() - vertex[0]->Z());

		cout << "Muon " << i << " pT " << init_muons[i]->Pt() << " " << " relIso " << reliso << " dZ " << dZ <<  endl;
	      }
	      }*/

	    //if (event->eventId() == 64704220) {
	    if (event->eventId() == 38694703) {
	      cout << "Special event # " << ievt << " ID " << event->eventId() << " -> jets.size() = " << init_jets.size() << endl;
	      for (unsigned int i=0; i<init_jets.size(); i++) {
		TRootJet* raw = (TRootJet*) init_jets[i]->Clone();
		jetTools->unCorrectJet(raw,false);
		cout << "jet " << i << " -- raw pt " << raw->Pt() << " -- pt " << init_jets[i]->Pt() << " -- eta " << init_jets[i]->Eta() << " -- CSV value " << init_jets[i]->btag_combinedSecondaryVertexBJetTags() << endl;
	      }
	    }

	    if (event->eventId() == 96660436) {
	      cout << "Special event # " << ievt << " ID " << event->eventId() << " -> electrons.size() = " << init_electrons.size() << endl;
	      cout << endl << "Initial electrons in this event: " << endl;
	      for (unsigned int i=0; i<init_electrons.size(); i++) {
		float RelIso = (init_electrons[i]->chargedHadronIso()+init_electrons[i]->neutralHadronIso()+init_electrons[i]->photonIso())/init_electrons[i]->Pt();
		cout << "electron " << i << " -- pt " << init_electrons[i]->Pt() << " -- eta " << init_electrons[i]->Eta() << " -- RelIso " << RelIso << endl;
	      }

	      cout << endl << "Loose electrons in this event: " << endl;
	      for (unsigned int i=0; i<vetoElectrons.size(); i++) {
		float RelIso = (vetoElectrons[i]->chargedHadronIso()+vetoElectrons[i]->neutralHadronIso()+vetoElectrons[i]->photonIso())/vetoElectrons[i]->Pt();
		cout << "electron " << i << " -- pt " << vetoElectrons[i]->Pt() << " -- eta " << vetoElectrons[i]->Eta() << " -- RelIso " << RelIso <<  endl;
	      }

	    }
	    
	    if (event->lumiBlockId()==129009 && init_muons.size() > 0) {

	      cout << "Run " <<event->runId() << " Event " << event->eventId() << " Lumi " << event->lumiBlockId();

	      for (unsigned int i=0; i<init_muons.size(); i++) {
		float reliso = (init_muons[i]->chargedHadronIso()+init_muons[i]->neutralHadronIso()+init_muons[i]->photonIso())/init_muons[i]->Pt();
		float reliso2 = (init_muons[i]->chargedHadronIso() + max( 0.0, init_muons[i]->neutralHadronIso() + init_muons[i]->photonIso() - 0.5*init_muons[i]->puChargedHadronIso() ) ) / init_muons[i]->Pt();

		cout << " ** Muon " << i << " pT " << init_muons[i]->Pt() << " " << " relIso " << reliso << " relIso(dBCorr) " << reliso2;
	      }
	      
	      cout << endl;

	    }
	    
	    if(selectedMuons.size()==1){
	      selecTable.Fill(d,3,1.);

	      if(vetoMuons.size()==1){
		selecTable.Fill(d,4,1.);
		
		//cout << event->runId() << " : " << event->eventId() << " : " << event->lumiBlockId() << endl;

		if(vetoElectrons.size()==0) {
		  selecTable.Fill(d,5,1.);

		  if (event->eventId() == 38694853) {
		    cout << "After electron veto event # " << ievt << " ID " << event->eventId() << " -> jets.size() = " << init_jets.size() << endl;
		    for (unsigned int i=0; i<init_jets.size(); i++) {
		      TRootJet* raw = (TRootJet*) init_jets[i]->Clone();
		      jetTools->unCorrectJet(raw,false);
		      cout << "jet " << i << " -- raw pt " << raw->Pt() << " -- pt " << init_jets[i]->Pt() << " -- eta " << init_jets[i]->Eta() << endl;
		    }
		  }
		  
		  
		  if(selectedJets.size()>=1) {
		    selecTable.Fill(d,6,1.);
		    
		    if(selectedJets.size()>=2) {
		      selecTable.Fill(d,7,1.);
		      
		      if(selectedJets.size()>=3)
			{
			  selecTable.Fill(d,8,1.);
			  
			  if(selectedJets.size()>=4)
			    {
			      selecTable.Fill(d,9,1.);

			      int nTags = 0;

			      for (int j=0; j<selectedJets.size(); j++) {
				if (selectedJets[j]->btag_combinedSecondaryVertexBJetTags() > 0.679) {
				  nTags++;
				}
			      }
			      
			      if (nTags >= 1)
				selecTable.Fill(d,10,1.);
			      // EVENT IS SELECTED!!!!
 
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

      //***********************************//
      //***** ELECTRON+Jets RefSel V4 *****//
      //***********************************//
      
      else if (lepton == 1) {
	
	// FILL THE SELECTION TABLE //
	
	if(trigged){
	  selecTable.Fill(d,1,1.);
	  
	  if(isGoodPV){
	    selecTable.Fill(d,2,1.);

	    selection.setJetCuts(35.,2.5,0.01,1.,0.98,0.3,0.1);
	    selection.setLooseElectronCuts(35,2.5,0.2,0.);
		    
	    vector<TRootElectron*> selectedElectrons = selection.GetSelectedElectrons();
	    vector<TRootJet*> selectedJets = selection.GetSelectedJets(true);


	    vector<TRootElectron*> vetoElectrons = selection.GetSelectedLooseElectrons(true);
	    vector<TRootMuon*> vetoMuons = selection.GetSelectedLooseMuons();

	    if (event->eventId() == 65253633 || event->eventId() == 38694968) {
	      cout << "Special event # " << ievt << " ID " << event->eventId() << " -> electrons.size() = " << init_electrons.size() << endl;
	      cout << endl << "Initial electrons in this event: " << endl;
	      for (unsigned int i=0; i<init_electrons.size(); i++) {
		//float RelIso = (init_electrons[i]->chargedHadronIso()+init_electrons[i]->neutralHadronIso()+init_electrons[i]->photonIso())/init_electrons[i]->Pt();
		float RelIso = (init_electrons[i]->chargedHadronIso() + max( 0.0, init_electrons[i]->neutralHadronIso() + init_electrons[i]->photonIso() - 0.5*init_electrons[i]->puChargedHadronIso() ) ) / init_electrons[i]->Pt();
		cout << "electron " << i;
		cout << " -- pt " << init_electrons[i]->Pt();
		cout << " -- eta " << init_electrons[i]->Eta();
		cout << " -- RelIso " << RelIso;
		cout << " -- dXY " << fabs(init_electrons[i]->d0());
		cout << " -- mvaTrigV0 " << init_electrons[i]->mvaTrigId();
		cout << " -- mvaNonTrigV0 " << init_electrons[i]->mvaNonTrigId();
		cout << " -- conv " << init_electrons[i]->passConversion();
		cout << endl;


		cout << "electron " << i;
		cout << " -- chargedHadronIso " << init_electrons[i]->chargedHadronIso();
		cout << " -- puChargedHadronIso " << init_electrons[i]->puChargedHadronIso();
		cout << " -- neutralHadronIso " << init_electrons[i]->neutralHadronIso();
		cout << " -- photonIso " << init_electrons[i]->photonIso();

		cout << endl;

	      }

	    }
	    
	    if (selectedElectrons.size()==1) {
	      selecTable.Fill(d,3,1.);

	      //cout << event->runId() << " : " << event->eventId() << " : " << event->lumiBlockId() << endl;
	      
	      TRootElectron* electron = (TRootElectron*) selectedElectrons[0];
	      
	      if (histo1D.find("Et") == histo1D.end())
		histo1D["Et"] = new TH1F("Et","Et",300,0,300);
	      
	      histo1D["Et"]->Fill(electron->Et());
	      
	      if (histo1D.find("Eta") == histo1D.end())
		histo1D["Eta"] = new TH1F("Eta","Eta",60,0,6);
	      
	      histo1D["Eta"]->Fill(fabs(electron->Eta()));

	      if (histo1D.find("d0") == histo1D.end())
		histo1D["d0"] = new TH1F("d0","d0",50,0,0.5);
		
	      histo1D["d0"]->Fill(fabs(electron->d0()));
	      
	      if (histo1D.find("RelIso") == histo1D.end())
		histo1D["RelIso"] = new TH1F("RelIso","RelIso",50,0,0.5);
	      
	      float relISO = (electron->caloIso(3)+electron->trackerIso(3)) / electron->Et();
	      histo1D["RelIso"]->Fill(relISO);
	      
	      if (vetoMuons.size()==0) {
		selecTable.Fill(d,4,1.);
		
		if (vetoElectrons.size() == 1) {
		  selecTable.Fill(d,5,1.);
		  
		  if (selectedElectrons[0]->missingHits() == 0) {
		    selecTable.Fill(d,6,1.);
		    		    
		    if (fabs(selectedElectrons[0]->Dist()) >= 0.02 || fabs(selectedElectrons[0]->DCot()) >= 0.02) {
		      selecTable.Fill(d,7,1.);
			
		      if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-3) {
			selecTable.Fill(d,8,1.);
			
			if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-2) {
			  selecTable.Fill(d,9,1.);
			  
			  if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-1) {
			    selecTable.Fill(d,10,1.);
			    
			    if(selectedJets.size()>=(unsigned int)anaEnv.NofJets) {
			      selecTable.Fill(d,11,1.);

			      int nTags = 0;

			      for (int j=0; j<selectedJets.size(); j++) {
				if (selectedJets[j]->btag_combinedSecondaryVertexBJetTags() > 0.679) {
				  nTags++;
				}
			      }
			      
			      if (nTags >= 1)
				selecTable.Fill(d,12,1.);
			      // EVENT IS SELECTED!!!!

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
    
  } // dataset loop

  //Selection tables
  selecTable.TableCalculator(false, true, true, true,true);
  string selectiontable = "SelectionTable_SyncEx";
  selectiontable = selectiontable +".tex"; 	
  selecTable.Write(selectiontable.c_str());


  //TFile* fout = new TFile("out.root","RECREATE");

  fout->cd();

  TDirectory* th1dir = fout->mkdir("1D_histograms");

  th1dir->cd();

  fout->cd();

  for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++) {
    TH1F *temp = it->second;
    int N = temp->GetNbinsX();
    //temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
    //temp->SetBinContent(N+1,0);
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( ("plots/"+it->first+".png").c_str() );
  }
  
  fout->Close();

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;

}
 
