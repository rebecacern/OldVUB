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
#include "../Tools/interface/MultiSamplePlot.h"
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

int main(int argc, char *argv[]) {

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

    unsigned int lepton = 0; // 0 = muon+jets 1= electron +jets

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
  cout<<"Loading environment 2..."<<endl;
  new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);
  cout<<"Loading environment 3..."<<endl;
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

  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");
  //Global variable
  TRootEvent* event = 0;
  
  string pathPNG = "ValidationPlots";
  if (argc >= 3){
		string sample=string(argv[2]);
 		pathPNG = pathPNG+"_"+sample;
 	}
  pathPNG = pathPNG +"/"; 

  mkdir(pathPNG.c_str(),0777);
  
  ////////////////////////////////////
  /// Normal Plots (TH1F* and TH2F*) //note: not used for the moment
  ////////////////////////////////////

  map<string,TH1F*> histo1D;
  map<string,TH2F*> histo2D;


  ////////////////////////////////////
  /// Multisample Plots
  ////////////////////////////////////
  map<string,MultiSamplePlot*> MSPlot;
  
  MSPlot["nEventsAfterCuts"] = new MultiSamplePlot(datasets, "nEventsAfterCuts", 23, -0.5, 22.5, "Nr. of events after each cut");
  
  MSPlot["nPrimaryVertices"] = new MultiSamplePlot(datasets, "nPrimaryVertices", 12, -0.5, 11.5, "Nr. of primary vertices");
  MSPlot["nGoodPrimaryVertices"] = new MultiSamplePlot(datasets, "nGoodPrimaryVertices", 12, -0.5, 11.5, "Nr. of good primary vertices");
  // Variables used during event selection
  MSPlot["globalMuPt"] = new MultiSamplePlot(datasets, "globalMuPt", 200, 0, 200, "Global Muon P_{T}");
  MSPlot["globalMuEta"] = new MultiSamplePlot(datasets, "globalMuEta", 50, -2.5, 2.5, "Global Muon #eta");
  MSPlot["globalMuRelIso"] = new MultiSamplePlot(datasets, "globalMuRelIso", 100, -0.0001, 5, "Global Muon RelIso");
  MSPlot["globalMuChi2"] = new MultiSamplePlot(datasets, "globalMuChi2", 100, 0, 25, "Global Muon #Chi^{2}/ndf");
  MSPlot["globalMuNValidHits"] = new MultiSamplePlot(datasets, "globalMuNValidHits", 35, -0.5, 34.5, "Global Muon Number of Valid Hits");
  MSPlot["globalMuMinDRMuJet"] = new MultiSamplePlot(datasets, "globalMuMinDRMuJet", 110, 0, 5.5, "Global Muon Minimal #DeltaR(#mu, jet)");
  MSPlot["globalMud0"] = new MultiSamplePlot(datasets, "globalMud0", 100, -0.1, 0.1, "Global Muon d0");
  
  MSPlot["nSelectedMuons"] = new MultiSamplePlot(datasets, "nSelectedMuons", 5, -0.5, 4.5, "Nr. of selected Muons");

  MSPlot["otherGlobalMuPt"] = new MultiSamplePlot(datasets, "otherGlobalMuPt", 10, 0, 100, "Other Global Muon P_{T}");
  MSPlot["otherGlobalMuEta"] = new MultiSamplePlot(datasets, "otherGlobalMuEta", 7, -2.8, 2.8, "Other Global Muon #eta");
  MSPlot["otherGlobalMuRelIso"] = new MultiSamplePlot(datasets, "otherGlobalMuRelIso", 25, 0, 5, "Other Global Muon RelIso");
  
  MSPlot["nLooseOtherMuons"] = new MultiSamplePlot(datasets, "nLooseOtherMuons", 5, -0.5, 4.5, "Nr. of Loose Other Muons");

  MSPlot["BeforeSelection_ElectronEt"] = new MultiSamplePlot(datasets, "ElectronEt", 15, 0, 75, "Electron E_{T}");
  MSPlot["BeforeSelection_ElectronEta"] = new MultiSamplePlot(datasets, "ElectronEta", 13, -2.6, 2.6, "Electron #eta");
  MSPlot["BeforeSelection_ElectronRelIso"] = new MultiSamplePlot(datasets, "ElectronRelIso", 20, -0.0001, 10, "Electron RelIso");
  MSPlot["ElectronEt"] = new MultiSamplePlot(datasets, "ElectronEt", 15, 0, 75, "Electron E_{T}");
  MSPlot["ElectronEta"] = new MultiSamplePlot(datasets, "ElectronEta", 13, -2.6, 2.6, "Electron #eta");
  MSPlot["ElectronRelIso"] = new MultiSamplePlot(datasets, "ElectronRelIso", 20, -0.0001, 10, "Electron RelIso");
  
  MSPlot["nLooseElectrons"] = new MultiSamplePlot(datasets, "nLooseElectrons", 5, -0.5, 4.5, "Nr. of Loose Electrons");
  
  MSPlot["allJetsPt1"] = new MultiSamplePlot(datasets, "allJetsPt1", 20, 0, 250, "Hardest Jet P_{T}");
  MSPlot["allJetsPt2"] = new MultiSamplePlot(datasets, "allJetsPt2", 20, 0, 200, "2nd Hardest Jet P_{T}");
  MSPlot["allJetsPt3"] = new MultiSamplePlot(datasets, "allJetsPt3", 25, 0, 150, "3rd Hardest Jet P_{T}");
  MSPlot["allJetsPt4"] = new MultiSamplePlot(datasets, "allJetsPt4", 25, 0, 100, "4th Hardest Jet P_{T}");
  MSPlot["allJetsEta"] = new MultiSamplePlot(datasets, "allJetsEta", 25, -2.5, 2.5, "Jet #eta");
  MSPlot["allJetsEMF"] = new MultiSamplePlot(datasets, "allJetsEMF", 20, 0, 1, "Jet EMF");
  MSPlot["allJetsn90Hits"] = new MultiSamplePlot(datasets, "allJetsn90Hits", 60, -0.5, 59.5, "Jet n90Hits");
  MSPlot["allJetsfHPD"] = new MultiSamplePlot(datasets, "allJetsfHPD", 25, 0, 1, "Jet fHPD");
  
  MSPlot["nSelectedJets"] = new MultiSamplePlot(datasets, "nSelectedJets", 10, -0.5, 9.5, "Nr. of Selected Jets");
  
  // Plots of selected events
  MSPlot["selectedEventsMuPt"] = new MultiSamplePlot(datasets, "selectedEventsMuPt", 8, 0, 160, "Selected Muon P_{T}");
  MSPlot["selectedEventsMuEta"] = new MultiSamplePlot(datasets, "selectedEventsMuEta", 11, -2.2, 2.2, "Selected Muon #eta");
  MSPlot["selectedEventsMuRelIso"] = new MultiSamplePlot(datasets, "selectedEventsMuRelIso", 10, -0.0001, 0.05, "Selected Muon RelIso");
  MSPlot["selectedEventsMuChi2"] = new MultiSamplePlot(datasets, "selectedEventsMuChi2", 20, 0, 10, "Selected Muon #Chi^{2}/ndf");
  MSPlot["selectedEventsMuNValidHits"] = new MultiSamplePlot(datasets, "selectedEventsMuNValidHits", 25, 9.5, 34.5, "Selected Muon Number of Valid Hits");
  MSPlot["selectedEventsMuMinDRMuJet"] = new MultiSamplePlot(datasets, "selectedEventsMuMinDRMuJet", 10, 0, 4, "Selected Muon Minimal #DeltaR(#mu, jet)");    
  MSPlot["selectedEventsMud0"] = new MultiSamplePlot(datasets, "selectedEventsMud0", 10, -0.02, 0.02, "Selected Muon d0");
  
  MSPlot["selectedEventsMET"] = new MultiSamplePlot(datasets, "selectedEventsMET", 10, 0, 200, "Selected Events MET");
  MSPlot["selectedEventsM3"] = new MultiSamplePlot(datasets, "selectedEventsM3", 20, 0, 1000, "Selected Events M3");
  MSPlot["selectedEventsWMT"] = new MultiSamplePlot(datasets, "selectedEventsWMT", 10, 0, 400, "Selected Events W M_{T}");

  MSPlot["selectedEventsJetsPt"] = new MultiSamplePlot(datasets, "selectedEventsJetsPt", 20, 0, 200, "Selected Jets P_{T}");
  MSPlot["selectedEventsJetsEta"] = new MultiSamplePlot(datasets, "selectedEventsJetsEta", 13, -2.6, 2.6, "Jet #eta");

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
  } else if (lepton == 1){
    cout << " - Preparing the cutflow table for e+jets refSelV4" << endl;
    CutsSelecTable.push_back(string("initial"));
    CutsSelecTable.push_back(string("trigged"));
    CutsSelecTable.push_back(string("Good PV"));
    CutsSelecTable.push_back(string("1 isolated electron"));
    CutsSelecTable.push_back(string("Loose muon veto"));
    CutsSelecTable.push_back(string("Z veto"));
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

  } else {

    cout << "!!!!Lepton should equal 0 or 1, exiting" << endl;

    exit(1);

  }

  if (verbose > 0)
    cout << " - CutsSelectionTable instantiated ... size: " << CutsSelecTable.size() <<endl;
  SelectionTable selecTable(CutsSelecTable, datasets);
  //selecTable.SetLuminosity(Luminosity);
  if (verbose > 0)
    cout << " - SelectionTable instantiated ..." << endl;


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

    //selecTable.Fill(d,0, datasets[d]->Xsection() * datasets[d]->EquivalentLumi() );

    string dataSetName = datasets[d]->Name();
    string previousFilename = "";
    int iFile = -1;

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

      selecTable.Fill(d,0,1.);
      
      // scale factor for the event
      float scaleFactor = 1.;
      
      //nEvents[d]++;

      //cout << "anaEnv.JetType -> " << anaEnv.JetType << " --- " <<  endl;
      //cout << "anaEnv.METType -> " << anaEnv.METType << " --- " <<  endl;


      if(ievt%500 == 0)
        std::cout<<"Processing the "<<ievt<<"th event" <<flush<<"\r";
        
      //cout << "load event" << endl;

      event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);
      //cout << "load event - done" << endl;

      /////////////
      // TRIGGER //
      /////////////
      
      bool trigged=false;
      
      if (lepton == 0) {
	//if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") {
/*
        // check which file in the dataset it is to have the HLTInfo right
        string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
        if(previousFilename != currentFilename)
        {
      	  previousFilename = currentFilename;
	      iFile++;
	      cout<<"File changed!!! => iFile = "<<iFile<<endl;
        }
      
        int currentRun = event->runId();
        if(previousRun != currentRun)
        {
          previousRun = currentRun;
          if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
          {
          //cout << "data" << endl;
	        if (event->runId() < 147196)
	          itrigger = treeLoader.iTrigger (string ("HLT_Mu9"), currentRun, iFile);
	        else if (event->runId() >= 147196)
	          itrigger = treeLoader.iTrigger (string ("HLT_Mu15_v1"), currentRun, iFile);
	        }
          else
          {
            itrigger = treeLoader.iTrigger (string ("HLT_Mu9"), currentRun);
          }
        }
	trigged = treeLoader.EventTrigged (itrigger);	
*/
/*	int currentRun = event->runId();
	//if(previousRun != currentRun) {
	//  previousRun = currentRun;
	//cout << currentRun << endl;
	itrigger = treeLoader.iTrigger ("HLT_Mu9", currentRun);
	//}
	
	trigged = treeLoader.EventTrigged (itrigger);
	
	//} else 
*/
	trigged=true;//putting to true until I understand how the trigger info should be used
      }
      
      else if (lepton == 1) {

	trigged=false;
	
	if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") {
	  
	  // SHOULD BE CHECKED ?!
	  std::map<std::string,std::vector<TopTree::triggeredObject> > filters= event->getTriggerFilters(); //previously: double instead of TopTree::triggeredObject
	  
	  if (event->runId() < 140041 && filters.find("hltL1NonIsoHLTNonIsoSingleElectronLWEt10PixelMatchFilter") != filters.end())
	    trigged=true;
	  else if (event->runId() >= 140041 && event->runId() <= 143962 && filters.find("hltL1NonIsoHLTNonIsoSingleElectronEt15PixelMatchFilter") != filters.end())
	    trigged=true;
	  else if (event->runId() >= 143963 && event->runId() <= 146427 && filters.find("hltL1NonIsoHLTNonIsoSingleElectronEt15CaloEleIdPixelMatchFilter") != filters.end())
	    trigged=true;
	  else if (event->runId() >= 146428 && event->runId() <= 147116 && filters.find("hltL1NonIsoHLTNonIsoSingleElectronEt17CaloEleIdPixelMatchFilter") != filters.end())
	    trigged=true;
	  else if (event->runId() >= 147117 && filters.find("hltL1NonIsoHLTNonIsoSingleElectronEt17CaloEleIdPixelMatchFilter") != filters.end())
	    trigged=true;
	  else
	    trigged=false;
	  
	  //trigged=true; //putting to true until I understand how the trigger info should be used
	}
	else
	  trigged=true;
      }
 
      //cout << "passed trigger section (switched off for muons)" << endl;

      /////////////////////////////
      //   Selection
      /////////////////////////////
      
      //Declare selection instance    
      Selection selection(init_jets, init_muons, init_electrons, mets);
      //cout << "declared selection instance" << endl;
      
      bool isGoodPV = isGoodPV = selection.isPVSelected(vertex, 4,24,2.);  
    
	MSPlot["nPrimaryVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
	MSPlot["nEventsAfterCuts"]->Fill(0, datasets[d], true, Luminosity*scaleFactor);
	for(unsigned int i=0; i<init_electrons.size(); i++)
        {
            MSPlot["BeforeSelection_ElectronEt"]->Fill(init_electrons[i]->Et(), datasets[d], true, Luminosity*scaleFactor);
            MSPlot["BeforeSelection_ElectronEta"]->Fill(init_electrons[i]->Eta(), datasets[d], true, Luminosity*scaleFactor);
            MSPlot["BeforeSelection_ElectronRelIso"]->Fill(init_electrons[i]->combinedIso(3,3,3)/init_electrons[i]->Et(), datasets[d], true, Luminosity*scaleFactor);
        }
	// FILL THE SELECTION TABLE //
    //  cout << "fill the selection table" << endl;
	selecTable.Fill(d,0,1.);

      //*******************************//
      //***** Muon+Jets RefSel V4 *****//
      //*******************************//

      if (lepton == 0) {
	if(trigged){
	  MSPlot["nEventsAfterCuts"]->Fill(1, datasets[d], true, Luminosity*scaleFactor);
	  selecTable.Fill(d,1,1.);
	  	  
	  if(isGoodPV){
	    MSPlot["nEventsAfterCuts"]->Fill(2, datasets[d], true, Luminosity*scaleFactor);
	    MSPlot["nGoodPrimaryVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
	    selecTable.Fill(d,2,1.);
	    
	    vector<TRootJet*> selectedJets;
	    vector<TRootMuon*> selectedMuons;
	   
	    bool doPF2PAT = false; // not supported/synced atm...

	    if (init_jets.size() > 0) {
	      if (init_jets[0]->jetType() == 1 || doPF2PAT) { // calojets
		//cout << "Selecting for caloJets" << endl;
		selectedJets = selection.GetSelectedJets(true);
		selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
	      }	else {
		//cout << "Selecting for PF/JPT jets" << endl;
		vector<TRootMuon*> overlapMuons = selection.GetSelectedMuons(vertex[0]);
		//selection.setJetCuts(30.,2.4,0.01,1.,0.98,0.3,0.1); // refSelV4 values
		selectedJets = selection.GetSelectedJets(overlapMuons,true);
		selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);	
	      }
	    }
	    
     // cout << "will start the real selection" << endl;
            for(unsigned int i=0; i<init_muons.size(); i++)
            {
              if(init_muons[i]->isGlobalMuon() && init_muons[i]->Pt() > 18)
              {
                MSPlot["globalMuPt"]->Fill(init_muons[i]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["globalMuEta"]->Fill(init_muons[i]->Eta(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["globalMuChi2"]->Fill(init_muons[i]->chi2(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["globalMuNValidHits"]->Fill(init_muons[i]->nofValidHits(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["globalMud0"]->Fill(init_muons[i]->d0(), datasets[d], true, Luminosity*scaleFactor);
              
//              if()
                {
                  MSPlot["globalMuRelIso"]->Fill(init_muons[i]->relativeIso03(), datasets[d], true, Luminosity*scaleFactor);

                  float mindRMuJet = 999.;
                  for(unsigned int j=0;j<selectedJets.size();j++)
                  {
                    float dRMuJet = init_muons[i]->DeltaR(*selectedJets[j]);
                    if(dRMuJet < mindRMuJet) mindRMuJet = dRMuJet;
                  }
                  if(mindRMuJet != 999.)
                    MSPlot["globalMuMinDRMuJet"]->Fill(mindRMuJet, datasets[d], true, Luminosity*scaleFactor);
                }
              }
            }
	    MSPlot["nSelectedMuons"]->Fill(selectedMuons.size(), datasets[d], true, Luminosity*scaleFactor);


	    vector<TRootMuon*> vetoMuons = selection.GetSelectedLooseMuons();
	    vector<TRootElectron*> vetoElectrons = selection.GetSelectedLooseElectrons(false);

	    if(selectedMuons.size()==1){
	      MSPlot["nEventsAfterCuts"]->Fill(3, datasets[d], true, Luminosity*scaleFactor);
	      selecTable.Fill(d,3,1.);
	      
	      if(vetoMuons.size()==1){
	        MSPlot["nEventsAfterCuts"]->Fill(4, datasets[d], true, Luminosity*scaleFactor);
		selecTable.Fill(d,4,1.);
		for(unsigned int i=0; i<init_electrons.size(); i++)
                {
                  MSPlot["ElectronEt"]->Fill(init_electrons[i]->Et(), datasets[d], true, Luminosity*scaleFactor);
                  MSPlot["ElectronEta"]->Fill(init_electrons[i]->Eta(), datasets[d], true, Luminosity*scaleFactor);
                  MSPlot["ElectronRelIso"]->Fill(init_electrons[i]->combinedIso(3,3,3)/init_electrons[i]->Et(), datasets[d], true, Luminosity*scaleFactor);
                }
              
                MSPlot["nLooseElectrons"]->Fill(vetoElectrons.size(), datasets[d], true, Luminosity*scaleFactor);
		
		if(vetoElectrons.size()==0) {
		  MSPlot["nEventsAfterCuts"]->Fill(5, datasets[d], true, Luminosity*scaleFactor);
		  MSPlot["nSelectedJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
		  selecTable.Fill(d,5,1.);			    
		  
		  if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-3) {
		    MSPlot["nEventsAfterCuts"]->Fill(6, datasets[d], true, Luminosity*scaleFactor);
		    selecTable.Fill(d,6,1.);
		    
		    if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-2) {
		      MSPlot["nEventsAfterCuts"]->Fill(7, datasets[d], true, Luminosity*scaleFactor);
		      selecTable.Fill(d,7,1.);
		      
		      if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-1)
			{
			  MSPlot["nEventsAfterCuts"]->Fill(8, datasets[d], true, Luminosity*scaleFactor);
			  selecTable.Fill(d,8,1.);
			  
			  if(selectedJets.size()>=(unsigned int)anaEnv.NofJets)
			    {
			      // EVENT IS SELECTED!!!!
      //cout << "event is selected" << endl;
			      MSPlot["nEventsAfterCuts"]->Fill(9, datasets[d], true, Luminosity*scaleFactor);
      //cout << "fill the selection table with the number of events after cuts" << endl;
			      selecTable.Fill(d,9,1.);
			      
			      // plot some properties of the selected events
      //cout << "plot some properties of the selected events" << endl;
                              MSPlot["selectedEventsMuPt"]->Fill(selectedMuons[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                              MSPlot["selectedEventsMuEta"]->Fill(selectedMuons[0]->Eta(), datasets[d], true, Luminosity*scaleFactor);
                              MSPlot["selectedEventsMuRelIso"]->Fill(selectedMuons[0]->relativeIso03(), datasets[d], true, Luminosity*scaleFactor);
                              MSPlot["selectedEventsMuChi2"]->Fill(selectedMuons[0]->chi2(), datasets[d], true, Luminosity*scaleFactor);
                              MSPlot["selectedEventsMuNValidHits"]->Fill(selectedMuons[0]->nofValidHits(), datasets[d], true, Luminosity*scaleFactor);
                              MSPlot["selectedEventsMud0"]->Fill(selectedMuons[0]->d0(), datasets[d], true, Luminosity*scaleFactor);
  			      
			      float mindRMuJet = 999.;
                              for(unsigned int j=0;j<selectedJets.size();j++)
                              {
                                float dRMuJet = selectedMuons[0]->DeltaR(*selectedJets[j]);
                                if(dRMuJet < mindRMuJet) mindRMuJet = dRMuJet;
                              }
                              MSPlot["selectedEventsMuMinDRMuJet"] ->Fill(mindRMuJet, datasets[d], true, Luminosity*scaleFactor);   
                        
                              if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
                                 cout << event->runId() << ":" << event->eventId() << ":" << event->lumiBlockId()  << endl;
                        
                              MSPlot["selectedEventsMET"]->Fill( mets[0]->Et(), datasets[d], true, Luminosity*scaleFactor);
			      float newMT = sqrt( 2 * selectedMuons[0]->Pt() * mets[0]->Et() * ( 1 - cos( selectedMuons[0]->DeltaPhi(*mets[0]) ) ) );
 	  	              MSPlot["selectedEventsWMT"]->Fill( newMT, datasets[d], true, Luminosity*scaleFactor);

                              float M3 = -1, maxPt = -1;
			      int previousJet = 0;
			      for(int i=0;i<selectedJets.size();i++)
			      {
				                  for(int j=0;j<i;j++)
				                  {
					                  for(int k=0;k<j;k++)
					                  {
						                  float combinedPt = (*selectedJets[i] + *selectedJets[j] + *selectedJets[k]).Pt();
						                  if(combinedPt > maxPt)
						                  {
							                  maxPt = combinedPt;
							                  M3 = (*selectedJets[i] + *selectedJets[j] + *selectedJets[k]).M();
						                  }
					                  }
				                  }
				                  
				                  //Fill jet plots inside the first jet-loop
      //cout << "fill jet plots inside the first jet-loop" << endl;
				                  MSPlot["selectedEventsJetsPt"]->Fill(selectedJets[i]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                                                  MSPlot["selectedEventsJetsEta"]->Fill(selectedJets[i]->Eta(), datasets[d], true, Luminosity*scaleFactor);
			                         /*added the following myself...*/
					         if(previousJet == 0){
					           MSPlot["allJetsPt1"]->Fill(selectedJets[i]->Pt(), datasets[d], true, Luminosity*scaleFactor); 
					           previousJet = 1;
						 }
						 else if(previousJet == 1){
					           MSPlot["allJetsPt2"]->Fill(selectedJets[i]->Pt(), datasets[d], true, Luminosity*scaleFactor); 
					           previousJet = 2;
						 }
						 else if(previousJet == 2){
                      				   MSPlot["allJetsPt3"]->Fill(selectedJets[i]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                      				   previousJet = 3;
                    				 }
                    				 else if(previousJet == 3){
                      				   MSPlot["allJetsPt4"]->Fill(selectedJets[i]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                      				   previousJet = 4;
                    				 }	
                    				 MSPlot["allJetsEta"]->Fill(selectedJets[i]->Eta(), datasets[d], true, Luminosity*scaleFactor);
						 
			      }
                              MSPlot["selectedEventsM3"]->Fill(M3, datasets[d], true, Luminosity*scaleFactor);
			    }
			}
		    }
		  }
		}
	      }
	    }
	  }
	}
//cout << "done with refsel muon" << endl;
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
	    
	    vector<TRootElectron*> selectedElectrons = selection.GetSelectedElectrons(vertex[0]);
	    vector<TRootJet*> selectedJets = selection.GetSelectedJets(selectedElectrons,true);
	    vector<TRootElectron*> looseElectrons = selection.GetSelectedLooseElectrons(20,2.5,1.,true);
	    vector<TRootMuon*> vetoMuons = selection.GetSelectedLooseMuons();
	    
	    if (selectedElectrons.size()==1) {
	      selecTable.Fill(d,3,1.);
	      
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

		bool passZVeto = true;
		for (unsigned int e=0;e<looseElectrons.size();e++) {
		  
		  TRootElectron* el = (TRootElectron*) looseElectrons[e];
		  
		  if (fabs(el->superClusterEta()) > 1.5660 || fabs(el->superClusterEta()) < 1.4442) {
		    TLorentzVector Zcand = *looseElectrons[e]+*selectedElectrons[0];
		    
		    //cout << Zcand.M() << endl;
		    if (Zcand.M() > 76 && Zcand.M() < 106)
		      passZVeto = false;
		  }
		}
		
		if (passZVeto) {
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

			      // EVENT IS SELECTED
			      if(verbose>1) cout<<"EVENT SELECTED"<<endl;

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
  selecTable.TableCalculator(false, true, true, true);
  string selectiontable = "SelectionTable_SyncEx";
  selectiontable = selectiontable +".tex"; 	
  selecTable.Write(selectiontable.c_str());


  //TFile* fout = new TFile("out.root","RECREATE");

  //fout->cd();
  mkdir((pathPNG+"MSPlot/").c_str(),0777);
  TDirectory* th1dir = fout->mkdir("1D_histograms");

  /*th1dir->cd();

  fout->cd();

  for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++) {
    cout<<"??????"<<endl;
    TH1F *temp = it->second;
    int N = temp->GetNbinsX();
    //temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
    //temp->SetBinContent(N+1,0);
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( ("plots/"+it->first+".png").c_str() );
  }*/
	
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
    MultiSamplePlot *temp = it->second;
    string name = it->first;
    temp->Draw(false, name, true, true, true, true, true);
    temp->Write(fout, name, true, pathPNG+"MSPlot/");
  }
  
  fout->Close();

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;

}
 
