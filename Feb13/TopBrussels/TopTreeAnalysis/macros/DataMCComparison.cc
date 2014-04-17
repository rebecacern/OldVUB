#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "../Selection/interface/SelectionTable.h"
#include "../Tools/interface/PlottingTools.h"
#include "../Tools/interface/MultiSamplePlot.h"
#include "../Tools/interface/MultiCutPlot.h"
#include "../Tools/interface/CutImpactEvaluation.h"
#include "../Tools/interface/TemplateComparator.h"
#include "../Tools/interface/TTreeLoader.h"
#include "../Tools/interface/AnalysisEnvironmentLoader.h"
#include "../Content/interface/AnalysisEnvironment.h"
#include "../Content/interface/MCExpectation.h"
#include "../Content/interface/MCObsExpectation.h"
#include "../Content/interface/Dataset.h"
#include "../MCInformation/interface/MCWeighter.h"
#include "../Selection/interface/ElectronPlotter.h"
#include "../Selection/interface/MuonPlotter.h"
#include "../Selection/interface/JetPlotter.h"
#include "../Selection/interface/VertexPlotter.h"

#include "Style.C"

using namespace std;
using namespace TopTree;

int main (int argc, char *argv[])
{
  cout << "**************************************************" << endl;
  cout << " Begining of the program for Data-MC Comparison ! " << endl;
  cout << "**************************************************" << endl;

  //SetStyle if needed
  setTDRStyle(); 
  //setMyStyle();

  /////////////////////
  // Configuration
  /////////////////////
  ////////////////////////////////////////

  //xml file
  string xmlFileName ="../config/myconfig.xml";
  if (argc >= 2)
    xmlFileName=string(argv[1]);
  const char *xmlfile = xmlFileName.c_str();
  //Output files
  string rootFileName ("DataMCComparison.root");

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
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);
  int verbose = anaEnv.Verbose;

  /////////////////////
  // Load Datasets
  /////////////////////

  TTreeLoader treeLoader;
  cout << " - Load datasets ..." << endl;
  vector < Dataset > datasets;
  treeLoader.LoadDatasets (datasets, xmlfile);
  for(unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(datasets[i]);
  /////////////////////
  
  //vector of objects
  cout << " - Variable declaration ..." << endl;
  vector < TRootVertex > vertex;
  vector < TRootMuon > init_muons;
  vector < TRootElectron > init_electrons;
  vector < TRootJet > init_jets;
  vector < TRootCaloJet > init_Calojets;
  vector < TRootPFJet > init_PFjets;
  vector < TRootMET > mets;
  ////////////////////////////////////////

  TFile *fout = new TFile (rootFileName.c_str (), "RECREATE");
  //Global variable
  TRootEvent* event = 0;

  //nof selected events
  float NEvtsData = 0;
  Double_t *nEvents = new Double_t[datasets.size()];

  ////////////////////////////////////
  /// MultiSamplePlot
  ////////////////////////////////////
//  MultiSamplePlot* SelMuonPtMSPlot = new MultiSamplePlot(datasets, "SelectedMuonPt", 100, 0, 200, "Muon P_{T}");
//  MultiSamplePlot* AllMuonPtMSPlot = new MultiSamplePlot(datasets, "SelectedMuonPt", 100, 0, 200, "Muon P_{T}");

  map<string,MultiSamplePlot*> MSPlot;

  MSPlot["nEventsAfterCuts"] = new MultiSamplePlot(datasets, "nEventsAfterCuts", 20, -0.5, 19.5, "Nr. of events after each cut");
  MSPlot["nMuonsAfterCuts"] = new MultiSamplePlot(datasets, "nMuonsAfterCuts", 20, -0.5, 19.5, "Nr. of muons after each cut");

  // Muon stuff
  MSPlot["AllMuonsPt"] = new MultiSamplePlot(datasets, "AllMuonsPt", 86, 14, 100, "Muon P_{T}");
  MSPlot["AllMuonsEta"] = new MultiSamplePlot(datasets, "AllMuonsEta", 25, -2.5, 2.5, "Muon #eta");
  MSPlot["AllMuonsPhi"] = new MultiSamplePlot(datasets, "AllMuonsPhi", 34, -3.4, 3.4, "Muon #phi");
  MSPlot["AllMuonsRelIso"] = new MultiSamplePlot(datasets, "AllMuonsRelIso", 50, 0, 1, "Muon RelIso");
  MSPlot["AllMuonsECALIso"] = new MultiSamplePlot(datasets, "AllMuonsECALIso", 60, 0, 30, "Muon ECAL Iso");
  MSPlot["AllMuonsHCALIso"] = new MultiSamplePlot(datasets, "AllMuonsHCALIso", 60, 0, 30, "Muon HCAL Iso");
  MSPlot["AllMuonsTrackerIso"] = new MultiSamplePlot(datasets, "AllMuonsTrackerIso", 80, 0, 40, "Muon Tracker Iso");
  MSPlot["AllMuonsChi2"] = new MultiSamplePlot(datasets, "AllMuonsChi2", 50, 0, 25, "Muon normalized #chi^{2}");
  MSPlot["AllMuonsNHits"] = new MultiSamplePlot(datasets, "AllMuonsNHits", 40, -0.5, 39.5, "Muon Number of Valid Hits");
  MSPlot["AllMuonsd0"] = new MultiSamplePlot(datasets, "AllMuonsd0", 100, -0.2, 0.2, "Muon d0");
  MSPlot["AllMuonsMinDR"] = new MultiSamplePlot(datasets, "AllMuonsMinDR", 100, 0, 5, "Muon Minimal #DeltaR(#mu, loose jet)");
  MSPlot["AllMuonsWMt"] = new MultiSamplePlot(datasets, "AllMuonsWMt", 50, 0, 250, "W Transverse Mass");
  MSPlot["AllMuonsNumber"] = new MultiSamplePlot(datasets, "AllMuonsNumber", 10, -0.5, 9.5, "Nr. of muons");
  MSPlot["AllMuonsDirection"] = new MultiSamplePlot(datasets, "AllMuonsDirection", 9, -4.5, 4.5, "Direction of muons");
  
  MSPlot["LooseMuonsPt"] = new MultiSamplePlot(datasets, "LooseMuonsPt", 86, 14, 100, "Muon P_{T}");
  MSPlot["LooseMuonsEta"] = new MultiSamplePlot(datasets, "LooseMuonsEta", 25, -2.5, 2.5, "Muon #eta");
  MSPlot["LooseMuonsPhi"] = new MultiSamplePlot(datasets, "LooseMuonsPhi", 34, -3.4, 3.4, "Muon #phi");
  MSPlot["LooseMuonsRelIso"] = new MultiSamplePlot(datasets, "LooseMuonsRelIso", 50, 0, 1, "Muon RelIso");
  MSPlot["LooseMuonsECALIso"] = new MultiSamplePlot(datasets, "LooseMuonsECALIso", 60, 0, 30, "Muon ECAL Iso");
  MSPlot["LooseMuonsHCALIso"] = new MultiSamplePlot(datasets, "LooseMuonsHCALIso", 60, 0, 30, "Muon HCAL Iso");
  MSPlot["LooseMuonsTrackerIso"] = new MultiSamplePlot(datasets, "LooseMuonsTrackerIso", 80, 0, 40, "Muon Tracker Iso");
  MSPlot["LooseMuonsChi2"] = new MultiSamplePlot(datasets, "LooseMuonsChi2", 50, 0, 25, "Muon normalized #chi^{2}");
  MSPlot["LooseMuonsNHits"] = new MultiSamplePlot(datasets, "LooseMuonsNHits", 40, -0.5, 39.5, "Muon Number of Valid Hits");
  MSPlot["LooseMuonsd0"] = new MultiSamplePlot(datasets, "LooseMuonsd0", 100, -0.2, 0.2, "Muon d0");
  MSPlot["LooseMuonsMinDR"] = new MultiSamplePlot(datasets, "LooseMuonsMinDR", 100, 0, 5, "Muon Minimal #DeltaR(#mu, loose jet)");
  MSPlot["LooseMuonsWMt"] = new MultiSamplePlot(datasets, "LooseMuonsWMt", 50, 0, 250, "W Transverse Mass");
  MSPlot["LooseMuonsNumber"] = new MultiSamplePlot(datasets, "LooseMuonsNumber", 10, -0.5, 9.5, "Nr. of muons");
  MSPlot["LooseMuonsDirection"] = new MultiSamplePlot(datasets, "LooseMuonsDirection", 9, -4.5, 4.5, "Direction of muons");

  MSPlot["TightMuonsPt"] = new MultiSamplePlot(datasets, "TightMuonsPt", 43, 14, 100, "Muon P_{T}");
  MSPlot["TightMuonsEta"] = new MultiSamplePlot(datasets, "TightMuonsEta", 25, -2.5, 2.5, "Muon #eta");
  MSPlot["TightMuonsPhi"] = new MultiSamplePlot(datasets, "TightMuonsPhi", 34, -3.4, 3.4, "Muon #phi");
  MSPlot["TightMuonsRelIso"] = new MultiSamplePlot(datasets, "TightMuonsRelIso", 50, 0, 0.1, "Muon RelIso");
  MSPlot["TightMuonsECALIso"] = new MultiSamplePlot(datasets, "TightMuonsECALIso", 20, 0, 5, "Muon ECAL Iso");
  MSPlot["TightMuonsHCALIso"] = new MultiSamplePlot(datasets, "TightMuonsHCALIso", 20, 0, 5, "Muon HCAL Iso");
  MSPlot["TightMuonsTrackerIso"] = new MultiSamplePlot(datasets, "TightMuonsTrackerIso", 20, 0, 5, "Muon Tracker Iso");
  MSPlot["TightMuonsChi2"] = new MultiSamplePlot(datasets, "TightMuonsChi2", 40, 0, 10, "Muon normalized #chi^{2}");
  MSPlot["TightMuonsNHits"] = new MultiSamplePlot(datasets, "TightMuonsNHits", 30, 9.5, 39.5, "Muon Number of Valid Hits");
  MSPlot["TightMuonsd0"] = new MultiSamplePlot(datasets, "TightMuonsd0", 42, -0.021, 0.021, "Muon d0");
  MSPlot["TightMuonsMinDR"] = new MultiSamplePlot(datasets, "TightMuonsMinDR", 50, 0, 5, "Muon Minimal #DeltaR(#mu, loose jet)");
  MSPlot["TightMuonsWMt"] = new MultiSamplePlot(datasets, "TightMuonsWMt", 50, 0, 250, "W Transverse Mass");
  MSPlot["TightMuonsNumber"] = new MultiSamplePlot(datasets, "TightMuonsNumber", 10, -0.5, 9.5, "Nr. of muons");
  MSPlot["TightMuonsDirection"] = new MultiSamplePlot(datasets, "TightMuonsDirection", 9, -4.5, 4.5, "Direction of muons");
 
  MSPlot["2TightMuonsMass"] = new MultiSamplePlot(datasets, "2TightMuonsMass", 75, 0, 150, "M^{inv}_{#mu, #mu}");
  MSPlot["2+AllMuonsDR"] = new MultiSamplePlot(datasets, "2+AllMuonsDR", 50, 0, 5, "Muon #DeltaR(#mu, #mu)");
  MSPlot["2+TightMuonsDR"] = new MultiSamplePlot(datasets, "2+TightMuonsDR", 50, 0, 5, "Muon #DeltaR(#mu, #mu)");
  
  MSPlot["1TightMuPt"] = new MultiSamplePlot(datasets, "1TightMuPt", 50, 0, 100, "Muon P_{T}");
  MSPlot["1TightMuEta"] = new MultiSamplePlot(datasets, "1TightMuEta", 25, -2.5, 2.5, "Muon #eta");
  MSPlot["1TightMuPhi"] = new MultiSamplePlot(datasets, "1TightMuPhi", 34, -3.4, 3.4, "Muon #phi");
  MSPlot["1TightMuRelIso"] = new MultiSamplePlot(datasets, "1TightMuRelIso", 25, 0, 0.1, "Muon RelIso");
  MSPlot["1TightMuECALIso"] = new MultiSamplePlot(datasets, "1TightMuECALIso", 20, 0, 5, "Muon ECAL Iso");
  MSPlot["1TightMuHCALIso"] = new MultiSamplePlot(datasets, "1TightMuHCALIso", 20, 0, 5, "Muon HCAL Iso");
  MSPlot["1TightMuTrackerIso"] = new MultiSamplePlot(datasets, "1TightMuTrackerIso", 20, 0, 5, "Muon Tracker Iso");
  MSPlot["1TightMuChi2"] = new MultiSamplePlot(datasets, "1TightMuChi2", 40, 0, 10, "Muon normalized #chi^{2}");
  MSPlot["1TightMuNHits"] = new MultiSamplePlot(datasets, "1TightMuNHits", 30, 9.5, 39.5, "Muon Number of Valid Hits");
  MSPlot["1TightMud0"] = new MultiSamplePlot(datasets, "1TightMud0", 42, -0.021, 0.021, "Muon d0");
  MSPlot["1TightMuMinDR"] = new MultiSamplePlot(datasets, "1TightMuMinDR", 60, 0, 6, "Muon Minimal #DeltaR(#mu, loose jet)");
  MSPlot["1TightMuWMt"] = new MultiSamplePlot(datasets, "1TightMuWMt", 50, 0, 250, "W Transverse Mass");
  MSPlot["1TightMuNumber"] = new MultiSamplePlot(datasets, "1TightMuNumber", 10, -0.5, 9.5, "Nr. of muons");
  MSPlot["1TightMuDirection"] = new MultiSamplePlot(datasets, "1TightMuDirection", 9, -4.5, 4.5, "Direction of muons");
  
  MSPlot["1TightMu3+LooseJetsMuPt"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsMuPt", 50, 0, 100, "Muon P_{T}");
  MSPlot["1TightMu3+LooseJetsMuEta"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsMuEta", 25, -2.5, 2.5, "Muon #eta");
  MSPlot["1TightMu3+LooseJetsMuPhi"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsMuPhi", 17, -3.4, 3.4, "Muon #phi");
  MSPlot["1TightMu3+LooseJetsMuRelIso"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsMuRelIso", 25, 0, 0.1, "Muon RelIso");
  MSPlot["1TightMu3+LooseJetsMuECALIso"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsMuECALIso", 20, 0, 5, "Muon ECAL Iso");
  MSPlot["1TightMu3+LooseJetsMuHCALIso"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsMuHCALIso", 20, 0, 5, "Muon HCAL Iso");
  MSPlot["1TightMu3+LooseJetsMuTrackerIso"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsMuTrackerIso", 20, 0, 5, "Muon Tracker Iso");
  MSPlot["1TightMu3+LooseJetsMuChi2"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsMuChi2", 40, 0, 10, "Muon normalized #chi^{2}");
  MSPlot["1TightMu3+LooseJetsMuNHits"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsMuNHits", 30, 9.5, 39.5, "Muon Number of Valid Hits");
  MSPlot["1TightMu3+LooseJetsMud0"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsMud0", 42, -0.021, 0.021, "Muon d0");
  MSPlot["1TightMu3+LooseJetsMuMinDR"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsMuMinDR", 60, 0, 6, "Muon Minimal #DeltaR(#mu, loose jet)");
  MSPlot["1TightMu3+LooseJetsMuWMt"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsMuWMt", 50, 0, 250, "W Transverse Mass");
  MSPlot["1TightMu3+LooseJetsMuNumber"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsMuNumber", 10, -0.5, 9.5, "Nr. of muons");
  MSPlot["1TightMu3+LooseJetsMuDirection"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsMuDirection", 9, -4.5, 4.5, "Direction of muons");

  MSPlot["1TightMu3+TightJetsMuPt"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsMuPt", 50, 0, 100, "Muon P_{T}");
  MSPlot["1TightMu3+TightJetsMuEta"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsMuEta", 25, -2.5, 2.5, "Muon #eta");
  MSPlot["1TightMu3+TightJetsMuPhi"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsMuPhi", 17, -3.4, 3.4, "Muon #phi");
  MSPlot["1TightMu3+TightJetsMuRelIso"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsMuRelIso", 25, 0, 0.1, "Muon RelIso");
  MSPlot["1TightMu3+TightJetsMuECALIso"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsMuECALIso", 20, 0, 5, "Muon ECAL Iso");
  MSPlot["1TightMu3+TightJetsMuHCALIso"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsMuHCALIso", 20, 0, 5, "Muon HCAL Iso");
  MSPlot["1TightMu3+TightJetsMuTrackerIso"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsMuTrackerIso", 20, 0, 5, "Muon Tracker Iso");
  MSPlot["1TightMu3+TightJetsMuChi2"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsMuChi2", 40, 0, 10, "Muon normalized #chi^{2}");
  MSPlot["1TightMu3+TightJetsMuNHits"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsMuNHits", 30, 9.5, 39.5, "Muon Number of Valid Hits");
  MSPlot["1TightMu3+TightJetsMud0"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsMud0", 42, -0.021, 0.021, "Muon d0");
  MSPlot["1TightMu3+TightJetsMuMinDR"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsMuMinDR", 30, 0, 6, "Muon Minimal #DeltaR(#mu, Tight jet)");
  MSPlot["1TightMu3+TightJetsMuWMt"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsMuWMt", 50, 0, 250, "W Transverse Mass");
  MSPlot["1TightMu3+TightJetsMuNumber"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsMuNumber", 10, -0.5, 9.5, "Nr. of muons");
  MSPlot["1TightMu3+TightJetsMuDirection"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsMuDirection", 9, -4.5, 4.5, "Direction of muons");

  // Jet stuff
  MSPlot["1+LooseMuAllJetsPt"] = new MultiSamplePlot(datasets, "1+LooseMuAllJetsPt", 100, 0, 100, "Jet P_{T}");
  MSPlot["1+LooseMuAllJetsEta"] = new MultiSamplePlot(datasets, "1+LooseMuAllJetsEta", 50, -2.5, 2.5, "Jet #eta");
  MSPlot["1+LooseMuAllJetsPhi"] = new MultiSamplePlot(datasets, "1+LooseMuAllJetsPhi", 34, -3.4, 3.4, "Jet #phi");
  MSPlot["1+LooseMuAllJetsEMF"] = new MultiSamplePlot(datasets, "1+LooseMuAllJetsEMF", 100, 0, 1, "Jet EMF");
  MSPlot["1+LooseMuAllJetsN90Hits"] = new MultiSamplePlot(datasets, "1+LooseMuAllJetsN90Hits", 50, -0.5, 49.5, "Jet n90Hits");
  MSPlot["1+LooseMuAllJetsFHPD"] = new MultiSamplePlot(datasets, "1+LooseMuAllJetsFHPD", 50, 0, 1, "Jet fHPD");
  MSPlot["1+LooseMuAllJetsNumber"] = new MultiSamplePlot(datasets, "1+LooseMuAllJetsNumber", 60, -0.5, 59.5, "Nr. of Jets");
  MSPlot["1+LooseMuAllJetsTrackCountingHighEff"] = new MultiSamplePlot(datasets, "1+LooseMuAllJetsTrackCountingHighEff", 60, -15, 15, "TrackCountingHighEff_btag");
  MSPlot["1+LooseMuAllJetsTrackCountingHighPur"] = new MultiSamplePlot(datasets, "1+LooseMuAllJetsTrackCountingHighPur", 60, -15, 15, "TrackCountingHighPur_btag");
  MSPlot["1+LooseMuAllJetsCombinedSV"] = new MultiSamplePlot(datasets, "1+LooseMuAllJetsCombinedSV", 50, 0, 1, "CombinedSecondaryVertexBJet_btag");
  MSPlot["1+LooseMuAllJetsCombinedSVMVA"] = new MultiSamplePlot(datasets, "1+LooseMuAllJetsCombinedSVMVA", 50, 0, 1, "CombinedSecondaryVertexMVABJet_btag");
  
  MSPlot["1+LooseMuLooseJetsPt"] = new MultiSamplePlot(datasets, "1+LooseMuLooseJetsPt", 100, 0, 200, "Jet P_{T}");
  MSPlot["1+LooseMuLooseJetsEta"] = new MultiSamplePlot(datasets, "1+LooseMuLooseJetsEta", 50, -2.5, 2.5, "Jet #eta");
  MSPlot["1+LooseMuLooseJetsPhi"] = new MultiSamplePlot(datasets, "1+LooseMuLooseJetsPhi", 34, -3.4, 3.4, "Jet #phi");
  MSPlot["1+LooseMuLooseJetsEMF"] = new MultiSamplePlot(datasets, "1+LooseMuLooseJetsEMF", 25, 0, 1, "Jet EMF");
  MSPlot["1+LooseMuLooseJetsN90Hits"] = new MultiSamplePlot(datasets, "1+LooseMuLooseJetsN90Hits", 50, -0.5, 49.5, "Jet n90Hits");
  MSPlot["1+LooseMuLooseJetsFHPD"] = new MultiSamplePlot(datasets, "1+LooseMuLooseJetsFHPD", 50, 0, 1, "Jet fHPD");
  MSPlot["1+LooseMuLooseJetsNumber"] = new MultiSamplePlot(datasets, "1+LooseMuLooseJetsNumber", 10, -0.5, 9.5, "Nr. of Jets");
  MSPlot["1+LooseMuLooseJetsTrackCountingHighEff"] = new MultiSamplePlot(datasets, "1+LooseMuLooseJetsTrackCountingHighEff", 60, -15, 15, "TrackCountingHighEff_btag");
  MSPlot["1+LooseMuLooseJetsTrackCountingHighPur"] = new MultiSamplePlot(datasets, "1+LooseMuLooseJetsTrackCountingHighPur", 60, -15, 15, "TrackCountingHighEff_btag");
  MSPlot["1+LooseMuLooseJetsCombinedSV"] = new MultiSamplePlot(datasets, "1+LooseMuLooseJetsCombinedSV", 50, 0, 1, "CombinedSecondaryVertexBJet_btag");
  MSPlot["1+LooseMuLooseJetsCombinedSVMVA"] = new MultiSamplePlot(datasets, "1+LooseMuLooseJetsCombinedSVMVA", 50, 0, 1, "CombinedSecondaryVertexMVABJet_btag");

  MSPlot["1+LooseMuTightJetsPt"] = new MultiSamplePlot(datasets, "1+LooseMuTightJetsPt", 100, 0, 200, "Jet P_{T}");
  MSPlot["1+LooseMuTightJetsEta"] = new MultiSamplePlot(datasets, "1+LooseMuTightJetsEta", 50, -2.5, 2.5, "Jet #eta");
  MSPlot["1+LooseMuTightJetsPhi"] = new MultiSamplePlot(datasets, "1+LooseMuTightJetsPhi", 34, -3.4, 3.4, "Jet #phi");
  MSPlot["1+LooseMuTightJetsEMF"] = new MultiSamplePlot(datasets, "1+LooseMuTightJetsEMF", 25, 0, 1, "Jet EMF");
  MSPlot["1+LooseMuTightJetsN90Hits"] = new MultiSamplePlot(datasets, "1+LooseMuTightJetsN90Hits", 50, -0.5, 49.5, "Jet n90Hits");
  MSPlot["1+LooseMuTightJetsFHPD"] = new MultiSamplePlot(datasets, "1+LooseMuTightJetsFHPD", 50, 0, 1, "Jet fHPD");
  MSPlot["1+LooseMuTightJetsNumber"] = new MultiSamplePlot(datasets, "1+LooseMuTightJetsNumber", 10, -0.5, 9.5, "Nr. of Jets");
  MSPlot["1+LooseMuTightJetsTrackCountingHighEff"] = new MultiSamplePlot(datasets, "1+LooseMuTightJetsTrackCountingHighEff", 60, -15, 15, "TrackCountingHighEff_btag");
  MSPlot["1+LooseMuTightJetsTrackCountingHighPur"] = new MultiSamplePlot(datasets, "1+LooseMuTightJetsTrackCountingHighPur", 60, -15, 15, "TrackCountingHighEff_btag");
  MSPlot["1+LooseMuTightJetsCombinedSV"] = new MultiSamplePlot(datasets, "1+LooseMuTightJetsCombinedSV", 50, 0, 1, "CombinedSecondaryVertexBJet_btag");
  MSPlot["1+LooseMuTightJetsCombinedSVMVA"] = new MultiSamplePlot(datasets, "1+LooseMuTightJetsCombinedSVMVA", 50, 0, 1, "CombinedSecondaryVertexMVABJet_btag");

  MSPlot["1TightMuLooseJetsPt"] = new MultiSamplePlot(datasets, "1TightMuLooseJetsPt", 100, 0, 200, "Jet P_{T}");
  MSPlot["1TightMuLooseJetsEta"] = new MultiSamplePlot(datasets, "1TightMuLooseJetsEta", 50, -2.5, 2.5, "Jet #eta");
  MSPlot["1TightMuLooseJetsPhi"] = new MultiSamplePlot(datasets, "1TightMuLooseJetsPhi", 34, -3.4, 3.4, "Jet #phi");
  MSPlot["1TightMuLooseJetsEMF"] = new MultiSamplePlot(datasets, "1TightMuLooseJetsEMF", 25, 0, 1, "Jet EMF");
  MSPlot["1TightMuLooseJetsN90Hits"] = new MultiSamplePlot(datasets, "1TightMuLooseJetsN90Hits", 50, -0.5, 49.5, "Jet n90Hits");
  MSPlot["1TightMuLooseJetsFHPD"] = new MultiSamplePlot(datasets, "1TightMuLooseJetsFHPD", 50, 0, 1, "Jet fHPD");
  MSPlot["1TightMuLooseJetsNumber"] = new MultiSamplePlot(datasets, "1TightMuLooseJetsNumber", 10, -0.5, 9.5, "Nr. of Jets");
  MSPlot["1TightMuLooseJetsTrackCountingHighEff"] = new MultiSamplePlot(datasets, "1TightMuLooseJetsTrackCountingHighEff", 60, -15, 15, "TrackCountingHighEff_btag");
  MSPlot["1TightMuLooseJetsTrackCountingHighPur"] = new MultiSamplePlot(datasets, "1TightMuLooseJetsTrackCountingHighPur", 60, -15, 15, "TrackCountingHighEff_btag");
  MSPlot["1TightMuLooseJetsCombinedSV"] = new MultiSamplePlot(datasets, "1TightMuLooseJetsCombinedSV", 50, 0, 1, "CombinedSecondaryVertexBJet_btag");
  MSPlot["1TightMuLooseJetsCombinedSVMVA"] = new MultiSamplePlot(datasets, "1TightMuLooseJetsCombinedSVMVA", 50, 0, 1, "CombinedSecondaryVertexMVABJet_btag");

  MSPlot["1TightMuTightJetsPt"] = new MultiSamplePlot(datasets, "1TightMuTightJetsPt", 100, 0, 200, "Jet P_{T}");
  MSPlot["1TightMuTightJetsEta"] = new MultiSamplePlot(datasets, "1TightMuTightJetsEta", 50, -2.5, 2.5, "Jet #eta");
  MSPlot["1TightMuTightJetsPhi"] = new MultiSamplePlot(datasets, "1TightMuTightJetsPhi", 34, -3.4, 3.4, "Jet #phi");
  MSPlot["1TightMuTightJetsEMF"] = new MultiSamplePlot(datasets, "1TightMuTightJetsEMF", 25, 0, 1, "Jet EMF");
  MSPlot["1TightMuTightJetsN90Hits"] = new MultiSamplePlot(datasets, "1TightMuTightJetsN90Hits", 50, -0.5, 49.5, "Jet n90Hits");
  MSPlot["1TightMuTightJetsFHPD"] = new MultiSamplePlot(datasets, "1TightMuTightJetsFHPD", 50, 0, 1, "Jet fHPD");
  MSPlot["1TightMuTightJetsNumber"] = new MultiSamplePlot(datasets, "1TightMuTightJetsNumber", 10, -0.5, 9.5, "Nr. of Jets");
  MSPlot["1TightMuTightJetsTrackCountingHighEff"] = new MultiSamplePlot(datasets, "1TightMuTightJetsTrackCountingHighEff", 60, -15, 15, "TrackCountingHighEff_btag");
  MSPlot["1TightMuTightJetsTrackCountingHighPur"] = new MultiSamplePlot(datasets, "1TightMuTightJetsTrackCountingHighPur", 60, -15, 15, "TrackCountingHighEff_btag");
  MSPlot["1TightMuTightJetsCombinedSV"] = new MultiSamplePlot(datasets, "1TightMuTightJetsCombinedSV", 50, 0, 1, "CombinedSecondaryVertexBJet_btag");
  MSPlot["1TightMuTightJetsCombinedSVMVA"] = new MultiSamplePlot(datasets, "1TightMuTightJetsCombinedSVMVA", 50, 0, 1, "CombinedSecondaryVertexMVABJet_btag");

  MSPlot["1TightMu3+LooseJetsPt"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsPt", 100, 0, 200, "Jet P_{T}");
  MSPlot["1TightMu3+LooseJetsEta"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsEta", 50, -2.5, 2.5, "Jet #eta");
  MSPlot["1TightMu3+LooseJetsPhi"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsPhi", 34, -3.4, 3.4, "Jet #phi");
  MSPlot["1TightMu3+LooseJetsEMF"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsEMF", 25, 0, 1, "Jet EMF");
  MSPlot["1TightMu3+LooseJetsN90Hits"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsN90Hits", 50, -0.5, 49.5, "Jet n90Hits");
  MSPlot["1TightMu3+LooseJetsFHPD"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsFHPD", 50, 0, 1, "Jet fHPD");
  MSPlot["1TightMu3+LooseJetsNumber"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsNumber", 10, -0.5, 9.5, "Nr. of Jets");
  MSPlot["1TightMu3+LooseJetsM3"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsM3", 60, 0, 300, "Jet P_{T}");
  MSPlot["1TightMu3+LooseJetsTrackCountingHighEff"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsTrackCountingHighEff", 60, -15, 15, "TrackCountingHighEff_btag");
  MSPlot["1TightMu3+LooseJetsTrackCountingHighPur"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsTrackCountingHighPur", 60, -15, 15, "TrackCountingHighEff_btag");
  MSPlot["1TightMu3+LooseJetsCombinedSV"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsCombinedSV", 50, 0, 1, "CombinedSecondaryVertexBJet_btag");
  MSPlot["1TightMu3+LooseJetsCombinedSVMVA"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsCombinedSVMVA", 50, 0, 1, "CombinedSecondaryVertexMVABJet_btag");
  
  MSPlot["1TightMu3+TightJetsPt"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsPt", 100, 0, 200, "Jet P_{T}");
  MSPlot["1TightMu3+TightJetsEta"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsEta", 50, -2.5, 2.5, "Jet #eta");
  MSPlot["1TightMu3+TightJetsPhi"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsPhi", 34, -3.4, 3.4, "Jet #phi");
  MSPlot["1TightMu3+TightJetsEMF"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsEMF", 25, 0, 1, "Jet EMF");
  MSPlot["1TightMu3+TightJetsN90Hits"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsN90Hits", 50, -0.5, 49.5, "Jet n90Hits");
  MSPlot["1TightMu3+TightJetsFHPD"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsFHPD", 50, 0, 1, "Jet fHPD");
  MSPlot["1TightMu3+TightJetsNumber"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsNumber", 10, -0.5, 9.5, "Nr. of Jets");
  MSPlot["1TightMu3+TightJetsM3"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsM3", 60, 0, 300, "Jet P_{T}");
  MSPlot["1TightMu3+TightJetsTrackCountingHighEff"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsTrackCountingHighEff", 60, -15, 15, "TrackCountingHighEff_btag");
  MSPlot["1TightMu3+TightJetsTrackCountingHighPur"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsTrackCountingHighPur", 60, -15, 15, "TrackCountingHighEff_btag");
  MSPlot["1TightMu3+TightJetsCombinedSV"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsCombinedSV", 50, 0, 1, "CombinedSecondaryVertexBJet_btag");
  MSPlot["1TightMu3+TightJetsCombinedSVMVA"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsCombinedSVMVA", 50, 0, 1, "CombinedSecondaryVertexMVABJet_btag");

  //MET stuff
  MSPlot["1+LooseMuMET"] = new MultiSamplePlot(datasets, "1+LooseMuMET", 100, 0, 200, "MET");
  MSPlot["1+LooseMuMETx"] = new MultiSamplePlot(datasets, "1+LooseMuMETx", 200, -200, 200, "MET_{x}");
  MSPlot["1+LooseMuMETy"] = new MultiSamplePlot(datasets, "1+LooseMuMETy", 200, -200, 200, "MET_{y}");  
  
  MSPlot["1TightMuMET"] = new MultiSamplePlot(datasets, "1TightMuMET", 100, 0, 200, "MET");
  MSPlot["1TightMuMETx"] = new MultiSamplePlot(datasets, "1TightMuMETx", 200, -200, 200, "MET_{x}");
  MSPlot["1TightMuMETy"] = new MultiSamplePlot(datasets, "1TightMuMETy", 200, -200, 200, "MET_{y}");  
  
  MSPlot["1TightMu3+LooseJetsMET"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsMET", 100, 0, 200, "MET");
  MSPlot["1TightMu3+LooseJetsMETx"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsMETx", 200, -200, 200, "MET_{x}");
  MSPlot["1TightMu3+LooseJetsMETy"] = new MultiSamplePlot(datasets, "1TightMu3+LooseJetsMETy", 200, -200, 200, "MET_{y}");  
  
  MSPlot["1TightMu3+TightJetsMET"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsMET", 100, 0, 200, "MET");
  MSPlot["1TightMu3+TightJetsMETx"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsMETx", 200, -200, 200, "MET_{x}");
  MSPlot["1TightMu3+TightJetsMETy"] = new MultiSamplePlot(datasets, "1TightMu3+TightJetsMETy", 200, -200, 200, "MET_{y}");  
  
  ////////////////////////////////////
  /// Normal Plots
  ////////////////////////////////////
//  TH1F* h_muonPt = new TH1F("h_MuonPt","Muon P_{T}", 100, 0, 200);
//  TH1F* h_muonPtWeighted = new TH1F("h_MuonPtWeighted","Muon P_{T}", 100, 0, 200);
  
  ////////////////////////////////////
  /// Fill dummy histos for reweighting
  ////////////////////////////////////
//  TH1F* h_DummyMC = new TH1F("h_DummyMC","Dummy MC histo for Reweighting",100,0,200);
//  TH1F* h_DummyData = new TH1F("h_DummyData","Dummy Data histo for Reweighting",100,0,200);

//  for(unsigned int i=0; i<100000; i++)
//  {
//    h_DummyMC->Fill(i*200/100000); // uniform distribution
//    h_DummyData->Fill(0.5*i*200/100000); // unform but half as big
//  }

  ////////////////////////////////////
  /// Make and initiate MCWeighter objects
  ////////////////////////////////////
//  cout << "Initiating MCWeighter" << endl;
//  MCWeighter* weighter = new MCWeighter();

//  weighter->LoadVariable("muonPt", true, (TH1F*) h_DummyMC->Clone(), (TH1F*) h_DummyData->Clone(), 0, 0);
//  weighter->Write(fout);

//  cout << "MCWeighter Initiated" << endl;

  ////////////////////////////////////
  //	Loop on datasets
  ////////////////////////////////////
  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
  for (unsigned int d = 0; d < datasets.size (); d++)
  {
    if (verbose > 1)
      cout << "   Dataset " << d << ": " << datasets[d].Name () << "/ title : " << datasets[d].Title () << endl;

    //open files and load
    cout<<"LoadEvent"<<endl;
    treeLoader.LoadDataset (datasets[d], anaEnv);
    cout<<"LoadEvent"<<endl;

    ///////////////////////////////
    //  Plotters
    ///////////////////////////////
//    ElectronPlotter* allElecPlotter = new ElectronPlotter( ("AllElectrons"+datasets[d].Name()).c_str() );
//    ElectronPlotter* vetoElecPlotter = new ElectronPlotter( ("VetoElectrons"+datasets[d].Name()).c_str() );
//    MuonPlotter* allMuonPlotter = new MuonPlotter( ("AllMuons"+datasets[d].Name()).c_str() );
//    MuonPlotter* selMuonPlotter = new MuonPlotter( ("SelectedMuons"+datasets[d].Name()).c_str() );
//    MuonPlotter* vetoMuonPlotter = new MuonPlotter( ("VetoMuons"+datasets[d].Name()).c_str() );
//    JetPlotter* allJetPlotter = new JetPlotter( ("AllJets"+datasets[d].Name()).c_str() );
//    JetPlotter* selJetPlotter = new JetPlotter( ("SelectedJets"+datasets[d].Name()).c_str() );
//    VertexPlotter* allVertexPlotter = new VertexPlotter( ("AllVertex"+datasets[d].Name()).c_str() );
//    VertexPlotter* selVertexPlotter = new VertexPlotter( ("SelectedVertex"+datasets[d].Name()).c_str() );

    ////////////////////////////////////
    //	Loop on events
    ////////////////////////////////////
    nEvents[d] = 0;
    int itrigger = -1, previousRun = -1;
    if (verbose > 1)
      cout << "	Loop over events " << endl;
    for (unsigned int ievt = 0; ievt < datasets[d].NofEvtsToRunOver(); ievt++)
//    for (unsigned int ievt = 0; ievt < 100; ievt++)
    {
      nEvents[d]++;
      if(ievt%5000 == 0)
        std::cout<<"Processing the "<<ievt<<"th event"<<flush<<"\r";

//      TLorentzVector JESbalance (0, 0, 0, 0);
     
      //load event
      if(anaEnv.JetType == 0 ) event = treeLoader.LoadEvent (ievt, &datasets[d], vertex, init_muons, init_electrons, init_jets, mets);
      if(anaEnv.JetType == 1 ) event = treeLoader.LoadEvent (ievt, &datasets[d], vertex, init_muons, init_electrons, init_Calojets, mets);
      if(anaEnv.JetType == 2 ) event = treeLoader.LoadEvent (ievt, &datasets[d], vertex, init_muons, init_electrons, init_PFjets, mets);
      ////////////////////////////////////////
    
      int currentRun = event->runId();
      if(previousRun != currentRun)
      {
        previousRun = currentRun;
        itrigger = treeLoader.iTrigger (&datasets[d], string ("HLT_L2Mu9"), currentRun);
      }
    
//      if(ievt == 0)
//        itrigger = treeLoader.iTrigger (&datasets[d], string ("HLT_L2Mu9"), event->runId());
//      else if(datasets[d].Name() == "Data" || datasets[d].Name() == "data" || datasets[d].Name() == "DATA")
//        itrigger = treeLoader.iTrigger (&datasets[d], string ("HLT_L2Mu9"), event->runId());
      
      /////////////////////////////
      //   Selection
      /////////////////////////////
      //Declare selection instance    
      Selection selection(anaEnv.JetType, init_jets, init_Calojets, init_PFjets, init_muons, init_electrons, mets);
      selection.SetConfiguration(anaEnv);

      MSPlot["nEventsAfterCuts"]->Fill(0, &datasets[d]);
//      if(datasets[d].Name() == "Data" || datasets[d].Name() == "data" || datasets[d].Name() == "DATA")
//      {
//        bool passedJSON = treeLoader.EventPassedJSON(datasets[d], event->runId(), event->lumiBlockId());
//        if(!passedJSON) continue;
//      }

      MSPlot["nEventsAfterCuts"]->Fill(1, &datasets[d]);

      bool trigged = treeLoader.EventTrigged (itrigger);
      if(!trigged) continue;
      
      MSPlot["nEventsAfterCuts"]->Fill(2, &datasets[d]);
      bool isGoodPV = selection.isPVSelected(vertex, anaEnv.PVertexNdofCut, anaEnv.PVertexZCut, anaEnv.PVertexRhoCut);
      if(!isGoodPV) continue;

      MSPlot["nEventsAfterCuts"]->Fill(3, &datasets[d]);
      bool isBeamBG = true;
      if(event->nTracks() > 10)
      {
        if( ( (float) event->nHighPurityTracks() ) / ( (float) event->nTracks() ) > 0.25 )
          isBeamBG = false;
      }
      if(isBeamBG) continue;

      MSPlot["nEventsAfterCuts"]->Fill(4, &datasets[d]);

      vector<TRootJet> selectedJets = selection.GetSelectedJets(anaEnv.JetsPtCutSR,anaEnv.JetsEtaCutSR, anaEnv.applyJetID);
//      vector<TRootMuon> selectedMuons = selection.GetSelectedMuons(anaEnv.MuonPtCutSR,anaEnv.MuonEtaCutSR,anaEnv.MuonRelIsoCutSR, selectedJets);
//      vector<TRootMuon> vetoMuons = selection.GetSelectedLooseMuons(anaEnv.MuonPtCutVetoSR,anaEnv.MuonEtaCutVetoSR,anaEnv.MuonRelIsoCutVetoSR);
//      vector<TRootElectron> vetoElectrons = selection.GetSelectedLooseElectrons(anaEnv.ElectronPtCut,anaEnv.ElectronEtaCut,anaEnv.ElectronRelIsoCut);
      
      vector<TRootCaloJet> selectedLooseJets, selectedTightJets;
      vector<TRootMuon> selectedTightMuons;

      unsigned int nAllMuons = 0, nLooseMuons = 0, nTightMuons = 0;
      float minDr2Muons = 9999.;
      for(unsigned int i=0; i<init_muons.size(); i++)
      {
        if(i > 0)
        {
          for(int j=i-1; j>=0; j--)
          {
//            cout << "i = " << i << "  j = " << j << endl;
            float tempDr2Mu = init_muons[i].DeltaR(init_muons[j]);
            if( fabs(tempDr2Mu - TMath::Pi()) < fabs(minDr2Muons - TMath::Pi()) ) minDr2Muons = tempDr2Mu;
          }
        }

        if(init_muons[i].Pt() <= 17.5 || fabs(init_muons[i].Eta()) >= 2.1) continue;
        
        nAllMuons++;
        MSPlot["nMuonsAfterCuts"]->Fill(1, &datasets[d], true, 0.253803345);
        
        MSPlot["AllMuonsPt"]->Fill(init_muons[i].Pt(), &datasets[d], true, 0.253803345);
        MSPlot["AllMuonsEta"]->Fill(init_muons[i].Eta(), &datasets[d], true, 0.253803345);
        MSPlot["AllMuonsPhi"]->Fill(init_muons[i].Phi(), &datasets[d], true, 0.253803345);
        MSPlot["AllMuonsRelIso"]->Fill(init_muons[i].relativeIso03(), &datasets[d], true, 0.253803345);
        MSPlot["AllMuonsECALIso"]->Fill(init_muons[i].isoR03_emEt(), &datasets[d], true, 0.253803345);
        MSPlot["AllMuonsHCALIso"]->Fill(init_muons[i].isoR03_hadEt(), &datasets[d], true, 0.253803345);
        MSPlot["AllMuonsTrackerIso"]->Fill(init_muons[i].isoR03_sumPt(), &datasets[d], true, 0.253803345);
        MSPlot["AllMuonsChi2"]->Fill(init_muons[i].chi2(), &datasets[d], true, 0.253803345);
        MSPlot["AllMuonsNHits"]->Fill(init_muons[i].nofValidHits(), &datasets[d], true, 0.253803345);
        MSPlot["AllMuonsd0"]->Fill(init_muons[i].d0(), &datasets[d], true, 0.253803345);
        MSPlot["AllMuonsWMt"]->Fill((init_muons[i]+mets[0]).Mt(), &datasets[d], true, 0.253803345);
        MSPlot["AllMuonsDirection"]->Fill(init_muons[i].direction(), &datasets[d], true, 0.253803345);
        
        float mindRMuJet = 999.;
        for(unsigned int j=0;j<selectedJets.size();j++)
        {
          float dRMuJet = init_muons[i].DeltaR(selectedJets[j]);
          if(dRMuJet < mindRMuJet) mindRMuJet = dRMuJet;
        }
        MSPlot["AllMuonsMinDR"]->Fill(mindRMuJet, &datasets[d], true, 0.253803345);
        
        if( !( init_muons[i].isGlobalMuon() ) ) continue;
        MSPlot["nMuonsAfterCuts"]->Fill(2, &datasets[d], true, 0.253803345);
        if( !( init_muons[i].isTrackerMuon() ) ) continue;
        MSPlot["nMuonsAfterCuts"]->Fill(3, &datasets[d], true, 0.253803345);
        nLooseMuons++;

        MSPlot["LooseMuonsPt"]->Fill(init_muons[i].Pt(), &datasets[d], true, 0.253803345);
        MSPlot["LooseMuonsEta"]->Fill(init_muons[i].Eta(), &datasets[d], true, 0.253803345);
        MSPlot["LooseMuonsPhi"]->Fill(init_muons[i].Phi(), &datasets[d], true, 0.253803345);
        MSPlot["LooseMuonsRelIso"]->Fill(init_muons[i].relativeIso03(), &datasets[d], true, 0.253803345);
        MSPlot["LooseMuonsECALIso"]->Fill(init_muons[i].isoR03_emEt(), &datasets[d], true, 0.253803345);
        MSPlot["LooseMuonsHCALIso"]->Fill(init_muons[i].isoR03_hadEt(), &datasets[d], true, 0.253803345);
        MSPlot["LooseMuonsTrackerIso"]->Fill(init_muons[i].isoR03_sumPt(), &datasets[d], true, 0.253803345);
        MSPlot["LooseMuonsChi2"]->Fill(init_muons[i].chi2(), &datasets[d], true, 0.253803345);
        MSPlot["LooseMuonsNHits"]->Fill(init_muons[i].nofValidHits(), &datasets[d], true, 0.253803345);
        MSPlot["LooseMuonsd0"]->Fill(init_muons[i].d0(), &datasets[d], true, 0.253803345);
        MSPlot["LooseMuonsMinDR"]->Fill(mindRMuJet, &datasets[d], true, 0.253803345);
        MSPlot["LooseMuonsWMt"]->Fill((init_muons[i]+mets[0]).Mt(), &datasets[d], true, 0.253803345);
        MSPlot["LooseMuonsDirection"]->Fill(init_muons[i].direction(), &datasets[d], true, 0.253803345);

        if( !( init_muons[i].idGlobalMuonPromptTight() ) ) continue;
        MSPlot["nMuonsAfterCuts"]->Fill(4, &datasets[d], true, 0.253803345);
        if( init_muons[i].nofValidHits() <= 10 ) continue;
        MSPlot["nMuonsAfterCuts"]->Fill(5, &datasets[d], true, 0.253803345);
        if( fabs(init_muons[i].d0()) >= 0.02 ) continue;
        MSPlot["nMuonsAfterCuts"]->Fill(6, &datasets[d], true, 0.253803345);
        if( mindRMuJet <= 0.3 ) continue;
        MSPlot["nMuonsAfterCuts"]->Fill(7, &datasets[d], true, 0.253803345);
        if( init_muons[i].relativeIso03() >= 0.1 ) continue;
        MSPlot["nMuonsAfterCuts"]->Fill(8, &datasets[d], true, 0.253803345);
        nTightMuons++;

        MSPlot["TightMuonsPt"]->Fill(init_muons[i].Pt(), &datasets[d], true, 0.253803345);
        MSPlot["TightMuonsEta"]->Fill(init_muons[i].Eta(), &datasets[d], true, 0.253803345);
        MSPlot["TightMuonsPhi"]->Fill(init_muons[i].Phi(), &datasets[d], true, 0.253803345);
        MSPlot["TightMuonsRelIso"]->Fill(init_muons[i].relativeIso03(), &datasets[d], true, 0.253803345);
        MSPlot["TightMuonsECALIso"]->Fill(init_muons[i].isoR03_emEt(), &datasets[d], true, 0.253803345);
        MSPlot["TightMuonsHCALIso"]->Fill(init_muons[i].isoR03_hadEt(), &datasets[d], true, 0.253803345);
        MSPlot["TightMuonsTrackerIso"]->Fill(init_muons[i].isoR03_sumPt(), &datasets[d], true, 0.253803345);
        MSPlot["TightMuonsChi2"]->Fill(init_muons[i].chi2(), &datasets[d], true, 0.253803345);
        MSPlot["TightMuonsNHits"]->Fill(init_muons[i].nofValidHits(), &datasets[d], true, 0.253803345);
        MSPlot["TightMuonsd0"]->Fill(init_muons[i].d0(), &datasets[d], true, 0.253803345);
        MSPlot["TightMuonsMinDR"]->Fill(mindRMuJet, &datasets[d], true, 0.253803345);
        MSPlot["TightMuonsWMt"]->Fill((init_muons[i]+mets[0]).Mt(), &datasets[d], true, 0.253803345);
        MSPlot["TightMuonsDirection"]->Fill(init_muons[i].direction(), &datasets[d], true, 0.253803345);
        
        selectedTightMuons.push_back(init_muons[i]);
      }
      
      if(nLooseMuons < 1) continue;

      MSPlot["AllMuonsNumber"]->Fill(nAllMuons, &datasets[d], true, 0.253803345);
      MSPlot["LooseMuonsNumber"]->Fill(nLooseMuons, &datasets[d], true, 0.253803345);
      MSPlot["TightMuonsNumber"]->Fill(nTightMuons, &datasets[d], true, 0.253803345);
      
      if(nTightMuons == 2)
        MSPlot["2TightMuonsMass"]->Fill((selectedTightMuons[0]+selectedTightMuons[1]).M(), &datasets[d], true, 0.253803345);
      
      MSPlot["1+LooseMuMET"]->Fill(mets[0].Pt(), &datasets[d], true, 0.253803345);
      MSPlot["1+LooseMuMETx"]->Fill(mets[0].Px(), &datasets[d], true, 0.253803345);
      MSPlot["1+LooseMuMETy"]->Fill(mets[0].Py(), &datasets[d], true, 0.253803345);
      
      MSPlot["2+AllMuonsDR"]->Fill(minDr2Muons, &datasets[d], true, 0.253803345);
      
      if(nTightMuons > 0)
        MSPlot["2+TightMuonsDR"]->Fill(minDr2Muons, &datasets[d], true, 0.253803345);

      unsigned int nLooseJets = 0, nTightJets = 0;

      for(unsigned int i=0; i<init_Calojets.size(); i++)
      {
        if( fabs(init_Calojets[i].Eta()) >= 2.4 ) continue;
        
        MSPlot["1+LooseMuAllJetsPt"]->Fill(init_Calojets[i].Pt(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuAllJetsEta"]->Fill(init_Calojets[i].Eta(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuAllJetsPhi"]->Fill(init_Calojets[i].Phi(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuAllJetsEMF"]->Fill(init_Calojets[i].ecalEnergyFraction(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuAllJetsN90Hits"]->Fill(init_Calojets[i].n90Hits(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuAllJetsFHPD"]->Fill(init_Calojets[i].fHPD(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuAllJetsTrackCountingHighEff"]->Fill(init_Calojets[i].btag_trackCountingHighEffBJetTags(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuAllJetsTrackCountingHighPur"]->Fill(init_Calojets[i].btag_trackCountingHighPurBJetTags(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuAllJetsCombinedSV"]->Fill(init_Calojets[i].btag_combinedSecondaryVertexBJetTags(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuAllJetsCombinedSVMVA"]->Fill(init_Calojets[i].btag_combinedSecondaryVertexMVABJetTags(), &datasets[d], true, 0.253803345);
        
        if(init_Calojets[i].ecalEnergyFraction() <= 0.01 || init_Calojets[i].n90Hits() <= 1 || init_Calojets[i].fHPD() >= 0.98 ) continue;
        if(init_Calojets[i].Pt() <= 15) continue;
        
        nLooseJets++;
        MSPlot["1+LooseMuLooseJetsPt"]->Fill(init_Calojets[i].Pt(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuLooseJetsEta"]->Fill(init_Calojets[i].Eta(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuLooseJetsPhi"]->Fill(init_Calojets[i].Phi(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuLooseJetsEMF"]->Fill(init_Calojets[i].ecalEnergyFraction(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuLooseJetsN90Hits"]->Fill(init_Calojets[i].n90Hits(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuLooseJetsFHPD"]->Fill(init_Calojets[i].fHPD(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuLooseJetsTrackCountingHighEff"]->Fill(init_Calojets[i].btag_trackCountingHighEffBJetTags(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuLooseJetsTrackCountingHighPur"]->Fill(init_Calojets[i].btag_trackCountingHighPurBJetTags(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuLooseJetsCombinedSV"]->Fill(init_Calojets[i].btag_combinedSecondaryVertexBJetTags(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuLooseJetsCombinedSVMVA"]->Fill(init_Calojets[i].btag_combinedSecondaryVertexMVABJetTags(), &datasets[d], true, 0.253803345);
        
        selectedLooseJets.push_back(init_Calojets[i]);
        
        if(init_Calojets[i].Pt() <= 25) continue;
        
        nTightJets++;
        MSPlot["1+LooseMuTightJetsPt"]->Fill(init_Calojets[i].Pt(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuTightJetsEta"]->Fill(init_Calojets[i].Eta(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuTightJetsPhi"]->Fill(init_Calojets[i].Phi(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuTightJetsEMF"]->Fill(init_Calojets[i].ecalEnergyFraction(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuTightJetsN90Hits"]->Fill(init_Calojets[i].n90Hits(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuTightJetsFHPD"]->Fill(init_Calojets[i].fHPD(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuTightJetsTrackCountingHighEff"]->Fill(init_Calojets[i].btag_trackCountingHighEffBJetTags(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuTightJetsTrackCountingHighPur"]->Fill(init_Calojets[i].btag_trackCountingHighPurBJetTags(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuTightJetsCombinedSV"]->Fill(init_Calojets[i].btag_combinedSecondaryVertexBJetTags(), &datasets[d], true, 0.253803345);
        MSPlot["1+LooseMuTightJetsCombinedSVMVA"]->Fill(init_Calojets[i].btag_combinedSecondaryVertexMVABJetTags(), &datasets[d], true, 0.253803345);
        
        selectedTightJets.push_back(init_Calojets[i]);
      }
      
      MSPlot["1+LooseMuAllJetsNumber"]->Fill(init_Calojets.size(), &datasets[d], true, 0.253803345);
      MSPlot["1+LooseMuLooseJetsNumber"]->Fill(nLooseJets, &datasets[d], true, 0.253803345);
      MSPlot["1+LooseMuTightJetsNumber"]->Fill(nTightJets, &datasets[d], true, 0.253803345);
      
      if(nTightMuons != 1) continue;
      
      MSPlot["1TightMuMET"]->Fill(mets[0].Pt(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMuMETx"]->Fill(mets[0].Px(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMuMETy"]->Fill(mets[0].Py(), &datasets[d], true, 0.253803345);
      
      float mindRMuLooseJet = 9999;
      
      for(unsigned int i=0; i<selectedLooseJets.size(); i++)
      {
        float dRMuJet = selectedTightMuons[0].DeltaR(selectedLooseJets[i]);
        if(dRMuJet < mindRMuLooseJet) mindRMuLooseJet = dRMuJet;

        MSPlot["1TightMuLooseJetsPt"]->Fill(selectedLooseJets[i].Pt(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMuLooseJetsEta"]->Fill(selectedLooseJets[i].Eta(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMuLooseJetsPhi"]->Fill(selectedLooseJets[i].Phi(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMuLooseJetsEMF"]->Fill(selectedLooseJets[i].ecalEnergyFraction(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMuLooseJetsN90Hits"]->Fill(selectedLooseJets[i].n90Hits(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMuLooseJetsFHPD"]->Fill(selectedLooseJets[i].fHPD(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMuLooseJetsTrackCountingHighEff"]->Fill(selectedLooseJets[i].btag_trackCountingHighEffBJetTags(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMuLooseJetsTrackCountingHighPur"]->Fill(selectedLooseJets[i].btag_trackCountingHighPurBJetTags(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMuLooseJetsCombinedSV"]->Fill(selectedLooseJets[i].btag_combinedSecondaryVertexBJetTags(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMuLooseJetsCombinedSVMVA"]->Fill(selectedLooseJets[i].btag_combinedSecondaryVertexMVABJetTags(), &datasets[d], true, 0.253803345);
      }

      MSPlot["1TightMuLooseJetsNumber"]->Fill(selectedLooseJets.size(), &datasets[d], true, 0.253803345);
            
      for(unsigned int i=0; i<selectedTightJets.size(); i++)
      {
        MSPlot["1TightMuTightJetsPt"]->Fill(selectedTightJets[i].Pt(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMuTightJetsEta"]->Fill(selectedTightJets[i].Eta(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMuTightJetsPhi"]->Fill(selectedTightJets[i].Phi(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMuTightJetsEMF"]->Fill(selectedTightJets[i].ecalEnergyFraction(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMuTightJetsN90Hits"]->Fill(selectedTightJets[i].n90Hits(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMuTightJetsFHPD"]->Fill(selectedTightJets[i].fHPD(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMuTightJetsTrackCountingHighEff"]->Fill(selectedTightJets[i].btag_trackCountingHighEffBJetTags(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMuTightJetsTrackCountingHighPur"]->Fill(selectedTightJets[i].btag_trackCountingHighPurBJetTags(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMuTightJetsCombinedSV"]->Fill(selectedTightJets[i].btag_combinedSecondaryVertexBJetTags(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMuTightJetsCombinedSVMVA"]->Fill(selectedTightJets[i].btag_combinedSecondaryVertexMVABJetTags(), &datasets[d], true, 0.253803345);
      }
      
      MSPlot["1TightMuTightJetsNumber"]->Fill(selectedTightJets.size(), &datasets[d], true, 0.253803345);

      MSPlot["1TightMuPt"]->Fill(selectedTightMuons[0].Pt(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMuEta"]->Fill(selectedTightMuons[0].Eta(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMuPhi"]->Fill(selectedTightMuons[0].Phi(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMuRelIso"]->Fill(selectedTightMuons[0].relativeIso03(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMuECALIso"]->Fill(selectedTightMuons[0].isoR03_emEt(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMuHCALIso"]->Fill(selectedTightMuons[0].isoR03_hadEt(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMuTrackerIso"]->Fill(selectedTightMuons[0].isoR03_sumPt(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMuChi2"]->Fill(selectedTightMuons[0].chi2(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMuNHits"]->Fill(selectedTightMuons[0].nofValidHits(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMud0"]->Fill(selectedTightMuons[0].d0(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMuMinDR"]->Fill(mindRMuLooseJet, &datasets[d], true, 0.253803345);
      MSPlot["1TightMuWMt"]->Fill((selectedTightMuons[0]+mets[0]).Mt(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMuDirection"]->Fill(selectedTightMuons[0].direction(), &datasets[d], true, 0.253803345);
      
      MSPlot["1TightMuNumber"]->Fill(1, &datasets[d], true, 0.253803345);

      if(selectedLooseJets.size() < 3) continue;
      
      float mindRMu3pLooseJet = 9999;
      
      MSPlot["1TightMu3+LooseJetsMET"]->Fill(mets[0].Pt(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+LooseJetsMETx"]->Fill(mets[0].Px(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+LooseJetsMETy"]->Fill(mets[0].Py(), &datasets[d], true, 0.253803345);
      
      float M3 = -1, maxPt = -1;
			for(int i=0;i<selectedLooseJets.size();i++)
			{
				for(int j=0;j<i;j++)
				{
					for(int k=0;k<j;k++)
					{
						float combinedPt = (selectedLooseJets[i] + selectedLooseJets[j] + selectedLooseJets[k]).Pt();
						if(combinedPt > maxPt)
						{
							maxPt = combinedPt;
							M3 = (selectedLooseJets[i] + selectedLooseJets[j] + selectedLooseJets[k]).M();
						}
					}
				}
			}
			
			MSPlot["1TightMu3+LooseJetsM3"]->Fill(M3, &datasets[d], true, 0.253803345);

      for(unsigned int i=0; i<selectedLooseJets.size(); i++)
      {
        float dRMuJet = selectedTightMuons[0].DeltaR(selectedLooseJets[i]);
        if(dRMuJet < mindRMu3pLooseJet) mindRMu3pLooseJet = dRMuJet;

        MSPlot["1TightMu3+LooseJetsPt"]->Fill(selectedLooseJets[i].Pt(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMu3+LooseJetsEta"]->Fill(selectedLooseJets[i].Eta(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMu3+LooseJetsPhi"]->Fill(selectedLooseJets[i].Phi(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMu3+LooseJetsEMF"]->Fill(selectedLooseJets[i].ecalEnergyFraction(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMu3+LooseJetsN90Hits"]->Fill(selectedLooseJets[i].n90Hits(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMu3+LooseJetsFHPD"]->Fill(selectedLooseJets[i].fHPD(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMu3+LooseJetsTrackCountingHighEff"]->Fill(selectedLooseJets[i].btag_trackCountingHighEffBJetTags(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMu3+LooseJetsTrackCountingHighPur"]->Fill(selectedLooseJets[i].btag_trackCountingHighPurBJetTags(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMu3+LooseJetsCombinedSV"]->Fill(selectedLooseJets[i].btag_combinedSecondaryVertexBJetTags(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMu3+LooseJetsCombinedSVMVA"]->Fill(selectedLooseJets[i].btag_combinedSecondaryVertexMVABJetTags(), &datasets[d], true, 0.253803345);
        
      }
      
      MSPlot["1TightMu3+LooseJetsNumber"]->Fill(selectedLooseJets.size(), &datasets[d], true, 0.253803345);
      
      MSPlot["1TightMu3+LooseJetsMuPt"]->Fill(selectedTightMuons[0].Pt(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+LooseJetsMuEta"]->Fill(selectedTightMuons[0].Eta(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+LooseJetsMuPhi"]->Fill(selectedTightMuons[0].Phi(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+LooseJetsMuRelIso"]->Fill(selectedTightMuons[0].relativeIso03(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+LooseJetsMuECALIso"]->Fill(selectedTightMuons[0].isoR03_emEt(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+LooseJetsMuHCALIso"]->Fill(selectedTightMuons[0].isoR03_hadEt(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+LooseJetsMuTrackerIso"]->Fill(selectedTightMuons[0].isoR03_sumPt(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+LooseJetsMuChi2"]->Fill(selectedTightMuons[0].chi2(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+LooseJetsMuNHits"]->Fill(selectedTightMuons[0].nofValidHits(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+LooseJetsMud0"]->Fill(selectedTightMuons[0].d0(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+LooseJetsMuMinDR"]->Fill(mindRMu3pLooseJet, &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+LooseJetsMuWMt"]->Fill((selectedTightMuons[0]+mets[0]).Mt(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+LooseJetsMuDirection"]->Fill(selectedTightMuons[0].direction(), &datasets[d], true, 0.253803345);
      
      MSPlot["1TightMu3+LooseJetsMuNumber"]->Fill(1, &datasets[d], true, 0.253803345);

      if(selectedTightJets.size() < 3) continue;
      
      float mindRMuTightJet = 9999;
      
      MSPlot["1TightMu3+TightJetsMET"]->Fill(mets[0].Pt(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+TightJetsMETx"]->Fill(mets[0].Px(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+TightJetsMETy"]->Fill(mets[0].Py(), &datasets[d], true, 0.253803345);
      
      float M3Tight = -1, maxPtTight = -1;
			for(int i=0;i<selectedTightJets.size();i++)
			{
				for(int j=0;j<i;j++)
				{
					for(int k=0;k<j;k++)
					{
						float combinedPt = (selectedTightJets[i] + selectedTightJets[j] + selectedTightJets[k]).Pt();
						if(combinedPt > maxPtTight)
						{
							maxPtTight = combinedPt;
							M3Tight = (selectedTightJets[i] + selectedTightJets[j] + selectedTightJets[k]).M();
						}
					}
				}
			}
			
			MSPlot["1TightMu3+TightJetsM3"]->Fill(M3Tight, &datasets[d], true, 0.253803345);

      for(unsigned int i=0; i<selectedTightJets.size(); i++)
      {
        float dRMuJet = selectedTightMuons[0].DeltaR(selectedTightJets[i]);
        if(dRMuJet < mindRMuTightJet) mindRMuTightJet = dRMuJet;

        MSPlot["1TightMu3+TightJetsPt"]->Fill(selectedTightJets[i].Pt(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMu3+TightJetsEta"]->Fill(selectedTightJets[i].Eta(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMu3+TightJetsPhi"]->Fill(selectedTightJets[i].Phi(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMu3+TightJetsEMF"]->Fill(selectedTightJets[i].ecalEnergyFraction(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMu3+TightJetsN90Hits"]->Fill(selectedTightJets[i].n90Hits(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMu3+TightJetsFHPD"]->Fill(selectedTightJets[i].fHPD(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMu3+TightJetsTrackCountingHighEff"]->Fill(selectedTightJets[i].btag_trackCountingHighEffBJetTags(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMu3+TightJetsTrackCountingHighPur"]->Fill(selectedTightJets[i].btag_trackCountingHighPurBJetTags(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMu3+TightJetsCombinedSV"]->Fill(selectedTightJets[i].btag_combinedSecondaryVertexBJetTags(), &datasets[d], true, 0.253803345);
        MSPlot["1TightMu3+TightJetsCombinedSVMVA"]->Fill(selectedTightJets[i].btag_combinedSecondaryVertexMVABJetTags(), &datasets[d], true, 0.253803345);
        
      }
      
      MSPlot["1TightMu3+TightJetsNumber"]->Fill(selectedTightJets.size(), &datasets[d], true, 0.253803345);
      
      MSPlot["1TightMu3+TightJetsMuPt"]->Fill(selectedTightMuons[0].Pt(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+TightJetsMuEta"]->Fill(selectedTightMuons[0].Eta(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+TightJetsMuPhi"]->Fill(selectedTightMuons[0].Phi(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+TightJetsMuRelIso"]->Fill(selectedTightMuons[0].relativeIso03(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+TightJetsMuECALIso"]->Fill(selectedTightMuons[0].isoR03_emEt(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+TightJetsMuHCALIso"]->Fill(selectedTightMuons[0].isoR03_hadEt(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+TightJetsMuTrackerIso"]->Fill(selectedTightMuons[0].isoR03_sumPt(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+TightJetsMuChi2"]->Fill(selectedTightMuons[0].chi2(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+TightJetsMuNHits"]->Fill(selectedTightMuons[0].nofValidHits(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+TightJetsMud0"]->Fill(selectedTightMuons[0].d0(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+TightJetsMuMinDR"]->Fill(mindRMuTightJet, &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+TightJetsMuWMt"]->Fill((selectedTightMuons[0]+mets[0]).Mt(), &datasets[d], true, 0.253803345);
      MSPlot["1TightMu3+TightJetsMuDirection"]->Fill(selectedTightMuons[0].direction(), &datasets[d], true, 0.253803345);
      
      MSPlot["1TightMu3+TightJetsMuNumber"]->Fill(1, &datasets[d], true, 0.253803345);
     
      
      
      
      
      
//      allElecPlotter->Fill(init_electrons);
//      vetoElecPlotter->Fill(vetoElectrons);
//      allMuonPlotter->Fill(init_muons, selectedJets);
//      selMuonPlotter->Fill(selectedMuons, selectedJets);
//      vetoMuonPlotter->Fill(vetoMuons, selectedJets);
//      if(anaEnv.JetType == 0 ) allJetPlotter->Fill(init_jets);
//      if(anaEnv.JetType == 1 ) allJetPlotter->Fill(init_Calojets);
//      if(anaEnv.JetType == 2 ) allJetPlotter->Fill(init_PFjets);
//      selJetPlotter->Fill(selectedJets);
//      allVertexPlotter->Fill(vertex);
//      if(isGoodPV) selVertexPlotter->Fill(vertex);

//      for(unsigned int i=0; i<selectedMuons.size(); i++)
//      {
//        SelMuonPtMSPlot->Fill(selectedMuons[i].Pt(), datasets[d]);
//        h_muonPt->Fill(selectedMuons[i].Pt());
//        h_muonPtWeighted->Fill(selectedMuons[i].Pt(),weighter->GetWeight("muonPt",selectedMuons[i].Pt()));
//      }

      //////////////////////////////
      
      //delete selection;
    }				//loop on events
    cout<<endl;

//    allElecPlotter->Write(fout);
//    delete allElecPlotter;
//    vetoElecPlotter->Write(fout);
//    delete vetoElecPlotter;
//    allMuonPlotter->Write(fout);
//    delete allMuonPlotter;
//    selMuonPlotter->Write(fout);
//    delete selMuonPlotter;
//    vetoMuonPlotter->Write(fout);
//    delete vetoMuonPlotter;
//    allJetPlotter->Write(fout);
//    delete allJetPlotter;
//    selJetPlotter->Write(fout);
//    delete selJetPlotter;
//    allVertexPlotter->Write(fout);
//    delete allVertexPlotter;
//    selVertexPlotter->Write(fout);
//    delete selVertexPlotter;

    //important: free memory
    treeLoader.UnLoadDataset();

  }				//loop on datasets

  //Once everything is filled ...
  if (verbose > 0)
    cout << " We ran over all the data ;-)" << endl;

  ///////////////////
  // Writting
  //////////////////
  if (verbose > 1)
  	cout << " - Writting outputs on files ..." << endl;

  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
    MultiSamplePlot *temp = it->second;
    string name = it->first;
    temp->Draw(false, name);
    temp->Write(fout, name, true, "Plots/");
  }

  fout->cd();
  // Write plots
//  h_muonPt->Write();
//  h_muonPtWeighted->Write();

  //add configuration info
  fout->cd();
  configTree->Fill();
  configTree->Write();

  //
  if (verbose > 1)
    cout << " - Done with writing the module outputs in the ouput file ..." << endl;
    cout << " - Closing the output file now..." << endl;
//  fout->Write();
  fout->Close();

  //delete
  delete fout;
  delete event;
  delete tcdatasets;
  delete tcAnaEnv;
  delete configTree;

  cout << "**********************************************************************" << endl;
  cout << "           End of the program !!" << endl;
  cout << " 		      hasn't crashed yet ;-) " << endl;
  cout << "**********************************************************************" << endl;

  return 0;
}
