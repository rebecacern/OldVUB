
#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <stdlib.h>

#include <sys/time.h>
#include <sys/resource.h>

  //user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysis/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysis/Tools/interface/JetTools.h"
#include "TopTreeAnalysis/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysis/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysis/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysis/Tools/interface/AnalysisEnvironmentLoader.h"
#include "TopTreeAnalysis/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysis/Content/interface/Dataset.h"
#include "TopTreeAnalysis/MCInformation/interface/Lumi3DReWeighting.h"
#include "TopTreeAnalysis/BkgEstimationMethods/interface/VJetEstimation.h"
#include "Style.C"

using namespace std;
using namespace TopTree;

Double_t TCHE (TopTree::TRootJet* jet) { return jet->btag_trackCountingHighEffBJetTags(); }
Double_t TCHP (TopTree::TRootJet* jet) { return jet->btag_trackCountingHighPurBJetTags(); }
Double_t JP   (TopTree::TRootJet* jet) { return jet->btag_jetProbabilityBJetTags(); }
Double_t JBP  (TopTree::TRootJet* jet) { return jet->btag_jetBProbabilityBJetTags(); }
Double_t SSVHE(TopTree::TRootJet* jet) { return jet->btag_simpleSecondaryVertexHighEffBJetTags(); }
Double_t SSVHP(TopTree::TRootJet* jet) { return jet->btag_simpleSecondaryVertexHighPurBJetTags(); }
Double_t CSV  (TopTree::TRootJet* jet) { return jet->btag_combinedSecondaryVertexBJetTags(); }

int main(int argc, const char* argv[]) {
  
  clock_t start = clock();
  
  cout << "********************************" << endl;
  cout << " Backgrounds for stops searches " << endl;
  cout << "********************************" << endl;
  
  TH1::AddDirectory(kFALSE);
    //SetStyle if needed
    //setTDRStyle(); 
  setMyStyle();
  
    /////////////////////
    // Configuration
    /////////////////////
  
    //xml file
    //  string xmlfile ="config/myRefSelconfig.xml";
  string baseDirectory = "." ;
  if (argc >= 2) {
    baseDirectory = string(argv[1]);
  }
  string outputDirectory = "." ;
  if (argc >= 3) {
    outputDirectory = string(argv[2]);
  }
  string xmlfile = baseDirectory + "/config/StopBckg.xml";
  
    //Input ROOT file
  string inputRootFileName = outputDirectory + "/StopBckg.root";
    //Output ROOT file
  string outputRootFileName = outputDirectory + "/StopBckg_Output.root";
  
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
  
  Bool_t runOnEvents_ = kTRUE;
  Bool_t reloadFromPreviousRun_ = kFALSE;
  
  TTreeLoader treeLoader;
  cout << " - Load datasets ..." << endl;
  vector < Dataset* > datasets;
  treeLoader.LoadDatasets (datasets, xmlfile.c_str());
  for(unsigned int i=0;i<datasets.size();i++) {
    new ((*tcdatasets)[i]) Dataset(*datasets[i]);
  }
  
  float oldLuminosity = anaEnv.Luminosity;	// in 1/pb
  float Luminosity = 5000.;
  Double_t eventsPrescaleFactor = 500. ;
  Bool_t channel_Single_Muon = kTRUE;
  Bool_t channel_Single_Electron = !channel_Single_Muon;
  
  /**
   Uncertainty variables
   */
  int uncertainty_PU = 0;
  int uncertainty_JER = 0; //0:nominal, 1:minus, 2:plus
  int uncertainty_JES = 0;
  int nsigma = 1;
  
  TFile *fin = new TFile (inputRootFileName.c_str(), "READ");
  
  /*
   for (unsigned int d = 0; d < datasets.size (); d++)
   {
   //cout << "luminosity of dataset "<< d << " is " << datasets[d]->EquivalentLumi() << endl;
   //    if(Luminosity > datasets[d]->EquivalentLumi() ) Luminosity = datasets[d]->EquivalentLumi();
   string dataSetName = datasets[d]->Name();
   if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
   {
   Luminosity = datasets[d]->EquivalentLumi();
   }
   
   }
   */
  if(Luminosity != oldLuminosity) cout << "changed analysis environment luminosity to "<< Luminosity << endl;
  
  
    //vector of objects
  cout << " - Variable declaration ..." << endl;
  vector < TRootVertex* > vertex;
  vector < TRootMuon* > init_muons;
  vector < TRootElectron* > init_electrons;
  vector < TRootJet* > init_jets, init_jets_fromTopTree;
  vector < TRootMET* > mets;
  
    //Global variable
  TRootEvent* event = 0;
  
    ////////////////////////////////////
    /// Normal Plots (TH1F* and TH2F*)
    ////////////////////////////////////
  
  map<string,TH1*> histo1D;
  map<string,TH2*> histo2D;
  
    ////////////////////////////////////
    /// Normal Plots (TH1F* and TH2F*)
    ////////////////////////////////////
  
  cout << " - 3D-reweighting instantiation ..." << endl;
  
  Lumi3DReWeighting Lumi3DWeights = Lumi3DReWeighting(baseDirectory+"/macros/PileUpReweighting/pileup_MC_Flat10PlusTail.root",baseDirectory+"/macros/PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root", "pileup", "pileup");
  Lumi3DWeights.weight3D_init(1.0);
  
  switch (uncertainty_PU) {
    case 0: Lumi3DWeights.weight3D_init(1.0); break;
    case 1: Lumi3DWeights.weight3D_init(0.92); break;
    case 2: Lumi3DWeights.weight3D_init(1.08); break;
    default: std::cerr << "Wrong parameter for : " << "uncertainty_PU" << std::endl; break;
  }
  
  
    ////////////////////////////////////
    /// Selection table and Multisample control plots
    ////////////////////////////////////
  
  vector<string> CutsSelecTable;
  
  cout << " - Preparing the cutflow table" << endl;
  CutsSelecTable.push_back(string("No Cuts"));
  CutsSelecTable.push_back(string("Scrap"));
  CutsSelecTable.push_back(string("Trigger"));
  CutsSelecTable.push_back(string("Nmu == 1"));
  CutsSelecTable.push_back(string("Nel == 0"));
  CutsSelecTable.push_back(string("Npf $\\geq$ 1"));
  CutsSelecTable.push_back(string("Npf $\\geq$ 2"));
  CutsSelecTable.push_back(string("Npf $\\geq$ 3"));
  CutsSelecTable.push_back(string("Npf $\\geq$ 4"));
  CutsSelecTable.push_back(string("Npf $\\geq$ 5"));
  
  map<string,MultiSamplePlot*> cutCP;
  cutCP["Trigger"] = new MultiSamplePlot(datasets, string("Before Isolated Muon selection"), 10, 0., 10., string("Number of good muons"), string("#Events/bin"));
  cutCP["IsolatedMuon"] = new MultiSamplePlot(datasets, string("Before vetoing on electron"), 10, 0., 10., string("Number of vetoing electrons"), string("#Events/bin"));
  cutCP["VetoElectron"] = new MultiSamplePlot(datasets, string("Before asking for at least 1 good jet"), 150, 0., 150., string("1^{st} Good Jet p_{T}"), string("#Events/bin"));
  cutCP["AtLeast1Jet"] = new MultiSamplePlot(datasets, string("Before asking for at least 2 good jets"), 150, 0., 150., string("2^{nd} Good Jet p_{T}"), string("#Events/bin"));
  cutCP["AtLeast2Jets"] = new MultiSamplePlot(datasets, string("Before asking for at least 3 good jets"), 150, 0., 150., string("3^{rd} Good Jet p_{T}"), string("#Events/bin"));
  cutCP["AtLeast3Jets"] = new MultiSamplePlot(datasets, string("Before asking for at least 4 good jets"), 150, 0., 150., string("4^{th} Good Jet p_{T}"), string("#Events/bin"));
  cutCP["AtLeast4Jets"] = new MultiSamplePlot(datasets, string("Before asking for at least 5 good jets"), 150, 0., 150., string("5^{th} Good Jet p_{T}"), string("#Events/bin"));
    //  cutCP["AtLeast5Jets"] = new MultiSamplePlot(datasets, string(""), int Nbins, float Min, float Max, string XaxisLabel = string("Variable"), string YaxisLabel = string("#Events"));
  
  map<string,MultiSamplePlot*> btagDiscCP;
  {
    string plotname = "TCHE : ";
    int Nbins = 200 ;
    Float_t Min = -10. ;
    Float_t Max = 10. ;
    string btagDiscName = "TCHE discriminant";
    btagDiscCP["Trigger"]      = new MultiSamplePlot(datasets, plotname+"after trigger",       Nbins, Min, Max, btagDiscName, string("#Jets/bin"));
    btagDiscCP["IsolatedMuon"] = new MultiSamplePlot(datasets, plotname+"after muon",          Nbins, Min, Max, btagDiscName, string("#Jets/bin"));
    btagDiscCP["VetoElectron"] = new MultiSamplePlot(datasets, plotname+"after veto electron", Nbins, Min, Max, btagDiscName, string("#Jets/bin"));
    btagDiscCP["AtLeast1Jet"]  = new MultiSamplePlot(datasets, plotname+"after 1 jet",         Nbins, Min, Max, btagDiscName, string("#Jets/bin"));
    btagDiscCP["AtLeast2Jets"] = new MultiSamplePlot(datasets, plotname+"after 2 jets",        Nbins, Min, Max, btagDiscName, string("#Jets/bin"));
    btagDiscCP["AtLeast3Jets"] = new MultiSamplePlot(datasets, plotname+"after 3 jets",        Nbins, Min, Max, btagDiscName, string("#Jets/bin"));
    btagDiscCP["AtLeast4Jets"] = new MultiSamplePlot(datasets, plotname+"after 4 jets",        Nbins, Min, Max, btagDiscName, string("#Jets/bin"));
  }
  
  if (verbose > 0)
    cout << " - CutsSelectionTable instantiated ... size: " << CutsSelecTable.size() <<endl;
  SelectionTable selecTable(CutsSelecTable, datasets);
  selecTable.SetLuminosity(Luminosity);
  if (verbose > 0)
    cout << " - SelectionTable instantiated ..." << endl;
  
    ////////////////////////////////////
    /// V-jets estimation
    ////////////////////////////////////
  
  if (verbose > 0)
    cout << " - V-jets estimation instantiation ..." << endl;
  const UInt_t NofBtagWorkingPoint = 2;
  float* BtagWorkingPoint;
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
  BtagWorkingPoint = new Float_t[NofBtagWorkingPoint];
  BtagWorkingPoint[0] = 1.7;
  BtagWorkingPoint[1] = 3.3;
  const UInt_t NofJets = 3;       // Min. mult. of jets
  const UInt_t NofJetBins = 4;    // Nb of bins for mult. of jets
                                  //  double** EffXbq;   // fraction of b-tag jets, by multiplicity ([i][j] is b-mult "j" in "i+NofJets" jets) ?
  histo1D["JetMultiplicities_Vlike"]  = new TH1F("JetMultiplicities_Vlike", "JetMultiplicities_Vlike", 21, 0, 21) ;
  histo1D["JetMultiplicities_TTlike"] = new TH1F("JetMultiplicities_TTlike", "JetMultiplicities_TTlike" ,21, 0, 21) ;
  histo1D["b-JetMultiplicities_Vlike"]  = new TH1F("b-JetMultiplicities_Vlike", "b-JetMultiplicities_Vlike", 21, 0, 21) ;
  histo1D["b-JetMultiplicities_TTlike"] = new TH1F("b-JetMultiplicities_TTlike", "b-JetMultiplicities_TTlike" ,21, 0, 21) ;
  
  std::vector < Dataset > vDatasets;
  for (std::vector<Dataset*>::iterator it=datasets.begin(); it!=datasets.end(); it++) {
    vDatasets.push_back(*(*it));
  }
  int NofDatasets = vDatasets.size();
  vector<Int_t> iDTTLike; // Indices of TT-Like datasets
  vector<Int_t> iDWLike;  // 
  vector<Int_t> iDVLike;  // 
  vector<Int_t> iDVbLike; //
  vector<Int_t> iDTTSTLike; //
  
  std::vector<std::string> ttLikeNames;
  std::vector<std::string> ttStLikeNames;
  std::vector<std::string> vLikeNames;
  std::vector<std::string> wLikeNames;
  std::vector<std::string> vbLikeNames;
  
  UInt_t num_data_datasets=0;
  for (Int_t i=0; i<NofDatasets; i++) {
    if(datasets[i]->Name() == "Data" || datasets[i]->Name() == "data" || datasets[i]->Name() == "DATA") { // Data!
      num_data_datasets++;
      continue; //ie Data cannot be considered as any MC (define other iD*Like vectors !!!)
    }
    if ((datasets[i]->Name().find("TT_")==0) || (datasets[i]->Name().find("Stop_")==0)) {
      iDTTLike.push_back(i);
      ttLikeNames.push_back(datasets[i]->Name());
      iDTTSTLike.push_back(i);
      ttStLikeNames.push_back(datasets[i]->Name());
    } 
    if (datasets[i]->Name().find("W_")==0) {
      iDWLike.push_back(i);
      wLikeNames.push_back(datasets[i]->Name());
      iDVLike.push_back(i);
      vLikeNames.push_back(datasets[i]->Name());
    }
    if ((datasets[i]->Name().find("Z_")==0) || (datasets[i]->Name().find("DY_")==0)) {
      iDVLike.push_back(i);
      vLikeNames.push_back(datasets[i]->Name());
    }
  }
  printf("W-like datasets : ");
  for (vector<Int_t>::iterator iter=iDWLike.begin(); iter!=iDWLike.end(); iter++) {
    if (iter != iDWLike.begin()) {
      printf(", ");
    }
    printf("%s", datasets[(*iter)]->Name().c_str());
  }
  printf("\n");
  printf("V-like datasets : ");
  for (vector<Int_t>::iterator iter=iDVLike.begin(); iter!=iDVLike.end(); iter++) {
    if (iter != iDVLike.begin()) {
      printf(", ");
    }
    printf("%s", datasets[(*iter)]->Name().c_str());
  }
  printf("\n");
  printf("TT-like datasets : ");
  for (vector<Int_t>::iterator iter=iDTTLike.begin(); iter!=iDTTLike.end(); iter++) {
    if (iter != iDTTLike.begin()) {
      printf(", ");
    }
    printf("%s", datasets[(*iter)]->Name().c_str());
  }
  printf("\n");
  printf("Vb-like datasets : ");
  for (vector<Int_t>::iterator iter=iDVbLike.begin(); iter!=iDVbLike.end(); iter++) {
    if (iter != iDVbLike.begin()) {
      printf(", ");
    }
    printf("%s", datasets[(*iter)]->Name().c_str());
  }
  printf("\n");
  
  printf("\n\nBy dataset names :\n");
  printf("W-like datasets : ");
  for (std::vector<std::string>::iterator iter=wLikeNames.begin(); iter!=wLikeNames.end(); iter++) {
    if (iter != wLikeNames.begin()) {
      printf(", ");
    }
    printf("%s", iter->c_str());
  }
  printf("\n");
  printf("V-like datasets : ");
  for (std::vector<std::string>::iterator iter=vLikeNames.begin(); iter!=vLikeNames.end(); iter++) {
    if (iter != vLikeNames.begin()) {
      printf(", ");
    }
    printf("%s", iter->c_str());
  }
  printf("\n");
  printf("TT-like datasets : ");
  for (std::vector<std::string>::iterator iter=ttLikeNames.begin(); iter!=ttLikeNames.end(); iter++) {
    if (iter != ttLikeNames.begin()) {
      printf(", ");
    }
    printf("%s", iter->c_str());
  }
  printf("\n");
  printf("Vb-like datasets : ");
  for (std::vector<std::string>::iterator iter=vbLikeNames.begin(); iter!=vbLikeNames.end(); iter++) {
    if (iter != vbLikeNames.begin()) {
      printf(", ");
    }
    printf("%s", iter->c_str());
  }
  printf("\n");
  
  
  VJetEstimation *vJetEst__tt__W               = new VJetEstimation(NofBtagWorkingPoint, BtagWorkingPoint, NofJets, NofJetBins, vDatasets, ttLikeNames, wLikeNames, vbLikeNames);
  /**/

  VJetEstimation *vJetEst__tt__W__Z__SingleT__QCD = new VJetEstimation(NofBtagWorkingPoint, BtagWorkingPoint, NofJets, NofJetBins, vDatasets, ttLikeNames, wLikeNames, vbLikeNames);

  VJetEstimation *vJetEst__tt__W_Z             = new VJetEstimation(NofBtagWorkingPoint, BtagWorkingPoint, NofJets, NofJetBins, vDatasets, ttLikeNames, vLikeNames, vbLikeNames);
  VJetEstimation *vJetEst__tt__W_Z__QCD         = new VJetEstimation(NofBtagWorkingPoint, BtagWorkingPoint, NofJets, NofJetBins, vDatasets, ttLikeNames, vLikeNames, vbLikeNames);
  VJetEstimation *vJetEst__tt__W_Z__SingleT__QCD = new VJetEstimation(NofBtagWorkingPoint, BtagWorkingPoint, NofJets, NofJetBins, vDatasets, ttStLikeNames, vLikeNames, vbLikeNames);
  VJetEstimation *vJetEst__tt_SingleT__W_Z__QCD = new VJetEstimation(NofBtagWorkingPoint, BtagWorkingPoint, NofJets, NofJetBins, vDatasets, ttStLikeNames, vLikeNames, vbLikeNames);
  VJetEstimation *vJetEst__tt_SingleT__W_Z = new VJetEstimation(NofBtagWorkingPoint, BtagWorkingPoint, NofJets, NofJetBins, vDatasets, ttStLikeNames, vLikeNames, vbLikeNames);
  
  VJetEstimation *vJetEst__data = new VJetEstimation(NofBtagWorkingPoint, BtagWorkingPoint, NofJets, NofJetBins, vDatasets, std::vector<std::string>(), std::vector<std::string>(), std::vector<std::string>());
  /**/
    //  VJetEstimation vJetEstSTasV(NofBtagWorkingPoint, BtagWorkingPoint, NofJets, NofJetBins, NofDatasets, iDTTLike, iDVLike, iDVbLike);
    //  VJetEstimation vJetEstSTasTT(NofBtagWorkingPoint, BtagWorkingPoint, NofJets, NofJetBins, NofDatasets, iDTTLike, iDVLike, iDVbLike);
    //  VJetEstimation vJetEst(NofBtagWorkingPoint, BtagWorkingPoint, NofJets, NofJetBins, EffXbq, NofDatasets, iDTTLike, iDVLike, iDVbLike);
  
  if (reloadFromPreviousRun_) {
    char name[200];
    sprintf(name, "VJetEstimation_%djets", NofJets);
      //      vJetEst = (VJetEstimation*) fin->Get(name);
  }
  
  
  
    //<a type="ParamForVJetEstimation" BtagAlgo_vjEst="0" NofBtagWorkingPoint_vjEst="1" BtagWorkingPoint_vjEst="2.03,3.20" MinMethod="Minuit2" MinOption="Combined" useMJLE="0" useUnBinMLE="1" NVJetPE="500" TagEffInit="0.794,0.128,0.097-0.70,0.043,0.02-0.63,0.05,0.010/0.807,0.134,0.124-0.70,0.043,0.02-0.63,0.05,0.010" NVlikeInit="14./4." NTTlikeInit="6./8." EffEbsel="0.0515,0.4170,0.5281/0.0187,0.2604,0.7049" VJEstFixParam="0,1,2" NofIterationsVJestShapeEstim="40"/>
  
  Double_t initForMinVJetArray[NofJetBins][2+3*NofBtagWorkingPoint] = {
      //For each jet bin :
      // N-tt-like (=effectiveXsecTTLike*IntLumi), N-V-like (=effectiveXsecVLike*IntLumi), {e_b, e_udsc, e_uds} for each working point
    {(155293.*datasets[iDTTLike[0]]->NormFactor())*Luminosity ,(50659.*datasets[iDVLike[0]]->NormFactor())*Luminosity,
      0.77,0.15,0.128,
      0.69,0.064,0.039},
    {(122632.*datasets[iDTTLike[0]]->NormFactor())*Luminosity ,(9816.*datasets[iDVLike[0]]->NormFactor())*Luminosity,
      0.77,0.15,0.128,
      0.69,0.059,0.043},
    {(71963.3*datasets[iDTTLike[0]]->NormFactor())*Luminosity ,(1890.*datasets[iDVLike[0]]->NormFactor())*Luminosity,
      0.77,0.15,0.128,
      0.69,0.059,0.033}};
    //         *currentDatasetNormFactor*Luminosity
  /*  {6.*Luminosity ,14.*Luminosity,
   0.69,0.064,0.039},
   {8.*Luminosity ,4.*Luminosity,
   0.69,0.059,0.043},
   {10.*Luminosity ,2.*Luminosity,
   0.69,0.059,0.033}}; */
  
  
  
  Double_t** initForMinVJet = NULL;
  initForMinVJet = new Double_t*[NofJetBins];
  for (UInt_t nJet=0; nJet<NofJetBins; nJet++) {
    initForMinVJet[nJet] = new Double_t[3*NofBtagWorkingPoint+2];
    for (UInt_t iParam=0; iParam<3*NofBtagWorkingPoint+2; iParam++) {
      initForMinVJet[nJet][iParam] = initForMinVJetArray[nJet][iParam];
    }
  }
  
  /*
   // Initialisation of the inputs
   Double_t**** Njets;
   Double_t**** Njets_background; //multijets
   Njets[wp] = new Double_t***[NofBtagWorkingPoint];
   Njets_background[wp] = new Double_t***[NofBtagWorkingPoint];
   for (int wp=0; wp<NofBtagWorkingPoint; wp++) {
   Njets[wp] = new Double_t**[NofJetBins];
   Njets_background[wp] = new Double_t**[NofJetBins];
   for (int jetBin=0; jetBin<NofJetBins; jetBin++) {
   Njets[wp][jetBin] = new Double_t*[];
   Njets_background[wp][jetBin] = new Double_t*[];
   for (int bJetBin=0; bJetBin<; bJetBin++) {
   Njets[wp][jetBin][bJetBin] = new Double_t[datasets.size()];
   Njets_background[wp][jetBin][bJetBin] = new Double_t[datasets.size()];
   for (int iDataset=0; iDataset<datasets.size(); iDataset++) {
   Njets[wp][jetBin][bJetBin][iDataset] = 0.;
   Njets_background[wp][jetBin][bJetBin][iDataset] = 0.
   }
   }
   }
   }
   */ 
  
  
    //MC
  vector<string> triggerTestsListMC;
  triggerTestsListMC.push_back("HLT_Mu17_CentralJet30_v2");
  triggerTestsListMC.push_back("HLT_Mu17_DiCentralJet30_v2");
  triggerTestsListMC.push_back("HLT_Mu17_TriCentralJet30_v2");
  triggerTestsListMC.push_back("HLT_Mu17_QuadCentralJet30_v2");
    //DATA
  vector<string> triggerTestsListData;
    //MuHad/Run2011A-May10ReReco-v1
  triggerTestsListData.push_back("HLT_Mu17_CentralJet30_v1");
  triggerTestsListData.push_back("HLT_Mu17_DiCentralJet30_v1");
  triggerTestsListData.push_back("HLT_Mu17_TriCentralJet30_v1");
  triggerTestsListData.push_back("HLT_Mu17_QuadCentralJet30_v1");
  triggerTestsListData.push_back("HLT_Mu17_CentralJet30_v2");
  triggerTestsListData.push_back("HLT_Mu17_DiCentralJet30_v2");
  triggerTestsListData.push_back("HLT_Mu17_TriCentralJet30_v2");
  triggerTestsListData.push_back("HLT_Mu17_QuadCentralJet30_v2");
  triggerTestsListData.push_back("HLT_Mu17_CentralJet30_v4");
  triggerTestsListData.push_back("HLT_Mu17_DiCentralJet30_v4");
  triggerTestsListData.push_back("HLT_Mu17_TriCentralJet30_v4");
  triggerTestsListData.push_back("HLT_Mu17_QuadCentralJet30_v4");
    //MuHad/Run2011A-PromptReco-v4
  triggerTestsListData.push_back("HLT_Mu17_CentralJet30_v5");
  triggerTestsListData.push_back("HLT_Mu17_DiCentralJet30_v5");
  triggerTestsListData.push_back("HLT_Mu17_TriCentralJet30_v5");
  triggerTestsListData.push_back("HLT_Mu17_QuadCentralJet30_v5");
  triggerTestsListData.push_back("HLT_Mu17_CentralJet30_v6");
  triggerTestsListData.push_back("HLT_Mu17_DiCentralJet30_v6");
  triggerTestsListData.push_back("HLT_Mu17_TriCentralJet30_v6");
  triggerTestsListData.push_back("HLT_Mu17_QuadCentralJet30_v6");
  triggerTestsListData.push_back("HLT_Mu17_CentralJet30_v7");
  triggerTestsListData.push_back("HLT_Mu17_DiCentralJet30_v7");
  triggerTestsListData.push_back("HLT_Mu17_TriCentralJet30_v7");
  triggerTestsListData.push_back("HLT_Mu17_QuadCentralJet30_v7");
  triggerTestsListData.push_back("HLT_Mu17_CentralJet30_v8");
  triggerTestsListData.push_back("HLT_Mu17_DiCentralJet30_v8");
  triggerTestsListData.push_back("HLT_Mu17_TriCentralJet30_v8");
  triggerTestsListData.push_back("HLT_Mu17_QuadCentralJet30_v8");
    //MuHad/Run2011A-05Aug2011-v1 + //MuHad/Run2011B-PromptReco-v1
  triggerTestsListData.push_back("HLT_Mu17_CentralJet30_v10");
  triggerTestsListData.push_back("HLT_Mu17_DiCentralJet30_v10");
  triggerTestsListData.push_back("HLT_Mu17_TriCentralJet30_v10");
  triggerTestsListData.push_back("HLT_Mu17_QuadCentralJet30_v10");
    //MuHad/Run2011B-PromptReco-v1
  triggerTestsListData.push_back("HLT_IsoMu17_eta2p1_CentralJet30_v5");
  triggerTestsListData.push_back("HLT_IsoMu17_eta2p1_DiCentralJet30_v5");
  triggerTestsListData.push_back("HLT_IsoMu17_eta2p1_TriCentralJet30_v5");
  triggerTestsListData.push_back("HLT_IsoMu17_eta2p1_QuadCentralJet30_v5");
  triggerTestsListData.push_back("HLT_IsoMu17_eta2p1_CentralJet30_v6");
  triggerTestsListData.push_back("HLT_IsoMu17_eta2p1_DiCentralJet30_v6");
  triggerTestsListData.push_back("HLT_IsoMu17_eta2p1_TriCentralJet30_v6");
  triggerTestsListData.push_back("HLT_IsoMu17_eta2p1_QuadCentralJet30_v6");
  
  fin->Close();
  
  
  if (runOnEvents_) {
      ////////////////////////////////////
      //	Loop on datasets
      ////////////////////////////////////
    
    UInt_t n_data_datasets = 0;
    if (verbose > 0)
      cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
    for (unsigned int d = 0; d < datasets.size (); d++)
    {
      
      if (verbose > 1)
        cout << "   Dataset " << d << ": " << datasets[d]->Name () << "/ title : " << datasets[d]->Title () << endl;
      int iFile = -1;
      string previousFilename = "";
      if (verbose > 1)
        std::cout<<"      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
        //open files and load
      
      treeLoader.LoadDataset (datasets[d], anaEnv);
        //selecTable.Fill(d,0, datasets[d]->Xsection() * datasets[d]->EquivalentLumi() );
      
      string dataSetName = datasets[d]->Name();
      Double_t currentDatasetNormFactor = datasets[d]->NormFactor() ;
      
      if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") // Data!
        n_data_datasets++;
      UInt_t d_notData = d - n_data_datasets;
      
        // Trigger test
      vector<string> triggerTestsList;
      if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") // Data!
      {
        triggerTestsList = triggerTestsListData ;
      } else {
        triggerTestsList = triggerTestsListMC ;
      }
      
      ULong64_t num = (1 << triggerTestsList.size());
      printf("Trigger Size = %lu , 2**size = %llu", triggerTestsList.size(), num);
      histo1D["TriggerTests_"+datasets[d]->Name()]  = new TH1I(("TriggerTests_"+datasets[d]->Name()).c_str(), ("TriggerTests_"+datasets[d]->Name()).c_str(), 1, 0, 1) ;
      histo1D["TriggerTests_"+datasets[d]->Name()]->SetBit(TH1::kCanRebin);
      /*      for (ULong64_t i=0; i<num; i++) {
       string label = "";
       ULong64_t j = i;
       for (vector<string>::iterator iter=triggerTestsList.begin(); iter!=triggerTestsList.end(); iter++) {
       if((j&1)==1) {
       if (label!="") {
       label = label + " && ";
       }
       label = label + (*iter);
       }
       j = (j >> 1);
       }
       histo1D["TriggerTests_"+datasets[d]->Name()]->GetXaxis()->SetBinLabel((Int_t)(i+1), label.c_str());
       }
       */    
      
        /////////////////////////////////////
        /// Initialize JEC factors
        /////////////////////////////////////
      
      vector<JetCorrectorParameters> vCorrParam;
      
        // Create the JetCorrectorParameter objects, the order does not matter.
        // YYYY is the first part of the txt files: usually the global tag from which they are retrieved
      /*    JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters("macros/JECFiles/Jec11V2_db_AK5PFchs_L3Absolute.txt");
       JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters("macros/JECFiles/Jec11V2_db_AK5PFchs_L2Relative.txt");
       JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters("macros/JECFiles/Jec11V2_db_AK5PFchs_L1FastJet.txt");
       */
      JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters(baseDirectory+"/macros/JECFiles/START42_V17_AK5PFchs_L3Absolute.txt");
      JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters(baseDirectory+"/macros/JECFiles/START42_V17_AK5PFchs_L2Relative.txt");
      JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters(baseDirectory+"/macros/JECFiles/START42_V17_AK5PFchs_L1FastJet.txt");
      
      
        //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
      vCorrParam.push_back(*L1JetPar);
      vCorrParam.push_back(*L2JetPar);
      vCorrParam.push_back(*L3JetPar);
      
      if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") // Data!
      {
        JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters(baseDirectory+"/macros/JECFiles/Jec11V2_db_AK5PFchs_L2L3Residual.txt");
          //JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters(baseDirectory+"/macros/JECFiles/START42_V17_AK5PFchs_L2L3Residual.txt");
        vCorrParam.push_back(*ResJetCorPar);
      }
      
      
        //OLD
        //JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/JEC11_V10_AK5PF_UncertaintySources.txt");
        //NEW
      JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(baseDirectory+"/macros/JECFiles/START42_V17_AK5PFchs_Uncertainty.txt");
      
        // true means redo also the L1
      JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);
      
      
      
        ////////////////////////////////////
        //	Loop on events
        ////////////////////////////////////
      
        //nEvents[d] = 0;
      int itrigger = -1; //previousRun = -1;
      if (verbose > 1)
        cout << "	Loop over events " << endl;
        //          for (unsigned long ievt = 0ul; ievt < (unsigned long) datasets[d]->NofEvtsToRunOver()/eventsPrescaleFactor; ievt++)
      for ( ULong_t ievt = 0ul;
           (((dataSetName == "Data") && (ievt < ((ULong_t) datasets[d]->NofEvtsToRunOver())))
            || ((dataSetName != "Data") && (ievt < ((ULong_t) (datasets[d]->NofEvtsToRunOver()/eventsPrescaleFactor))))) ;
           ((dataSetName == "Data") ?  (ievt=ievt+((ULong_t) eventsPrescaleFactor)) : (ievt++) )
           )
      {
          //      printf("%lu\n", ievt);
          // Double_t eventWeight = ((Double_t)datasets[d]->NofEvtsToRunOver()) / nofEvtsRanOver ;
        Double_t eventWeight = 1. ;
        eventWeight *= eventsPrescaleFactor ;
        
          //nEvents[d]++;
        
          //cout << "anaEnv.JetType -> " << anaEnv.JetType << " --- " <<  endl;
          //cout << "anaEnv.METType -> " << anaEnv.METType << " --- " <<  endl;
        
        
        if(ievt%500 == 0)
        {
          struct rusage r_use;
          char line[400];
          if (getrusage(RUSAGE_SELF, &r_use) == 0 ) 
          {
            sprintf(line, "[ User time :%8ld.%03ds ; System time : %8ld.%03ds ; MaxResidentMem : %8ldkb ] >>> Processing the %15luth event\r", r_use.ru_utime.tv_sec, (int) r_use.ru_utime.tv_usec/1000, r_use.ru_stime.tv_sec, (int) r_use.ru_stime.tv_usec/1000, r_use.ru_maxrss/1024, ievt);
          }
          else
            sprintf(line, "Processing the %15luth event\r", ievt);
          std::cerr << line << std::flush <<"\r";
        }
        
          //load event
        
        event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets_fromTopTree, mets);
        
          /////////////
          // TRIGGER //
          /////////////
        
        bool trigged=false;
        
          //if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") {
        
        string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
        if(previousFilename != currentFilename)
        {
          previousFilename = currentFilename;
          iFile++;
          std::cerr << std::endl << "File changed!!! => iFile = "<< iFile << std::endl;
        }
        int currentRun = event->runId();
          //if(previousRun != currentRun) {
          //          previousRun = currentRun;
          //cout << currentRun << endl;
        
          //      itrigger = treeLoader.iTrigger ("HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30_v1", currentRun, iFile);
        
        if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") {
          /*
           if ( 160431<=currentRun && currentRun<=161176 ) //MuHad/Run2011A-May10ReReco-v1
           itrigger = treeLoader.iTrigger("HLT_Mu17_CentralJet30_v1", currentRun, iFile);
           else if ( 161217<=currentRun && currentRun<=163261 )
           itrigger = treeLoader.iTrigger("HLT_Mu17_CentralJet30_v2", currentRun, iFile);
           else if ( 163270<=currentRun && currentRun<=163869 )
           itrigger = treeLoader.iTrigger("HLT_Mu17_CentralJet30_v4", currentRun, iFile);
           else if ( 165088<=currentRun && currentRun<=165633 )//MuHad/Run2011A-PromptReco-v4
           itrigger = treeLoader.iTrigger("HLT_Mu17_CentralJet30_v5", currentRun, iFile);
           else if ( 165970<=currentRun && currentRun<=167043 )
           itrigger = treeLoader.iTrigger("HLT_Mu17_CentralJet30_v6", currentRun, iFile);
           else if ( 166346<=currentRun && currentRun<=166346 )
           itrigger = treeLoader.iTrigger("HLT_Mu17_CentralJet30_v7", currentRun, iFile);
           else if ( 167078<=currentRun && currentRun<=167913 )
           itrigger = treeLoader.iTrigger("HLT_Mu17_CentralJet30_v8", currentRun, iFile);
           else if ( 170722<=currentRun && currentRun<=172619 ) //MuHad/Run2011A-05Aug2011-v1
           itrigger = treeLoader.iTrigger("HLT_Mu17_CentralJet30_v10", currentRun, iFile);
           else if ( 172620<=currentRun && currentRun<=173198 ) //MuHad/Run2011A-PromptReco-v6
           itrigger = treeLoader.iTrigger("HLT_Mu17_CentralJet30_v10", currentRun, iFile);
           else if ( 173241<=currentRun && currentRun<=178380 ) //MuHad/Run2011B-PromptReco-v1
           //        else if ( 175860<=currentRun && currentRun<=178380 ) //MuHad/Run2011B-PromptReco-v1
           itrigger = treeLoader.iTrigger("HLT_Mu17_eta2p1_CentralJet30_v1", currentRun, iFile);
           else if ( 178420<=currentRun && currentRun<=179889 )
           itrigger = treeLoader.iTrigger("HLT_IsoMu17_eta2p1_CentralJet30_v5", currentRun, iFile);
           else if ( 179959<=currentRun && currentRun<=180252 )
           itrigger = treeLoader.iTrigger("HLT_IsoMu17_eta2p1_CentralJet30_v6", currentRun, iFile);
           */
          if ( 160431<=currentRun && currentRun<=161176 ) //MuHad/Run2011A-May10ReReco-v1
            itrigger = treeLoader.iTrigger("HLT_Mu17_TriCentralJet30_v1", currentRun, iFile);
          else if ( 161217<=currentRun && currentRun<=163261 )
            itrigger = treeLoader.iTrigger("HLT_Mu17_TriCentralJet30_v2", currentRun, iFile);
          else if ( 163270<=currentRun && currentRun<=163869 )
            itrigger = treeLoader.iTrigger("HLT_Mu17_TriCentralJet30_v4", currentRun, iFile);
          else if ( 165088<=currentRun && currentRun<=165633 )//MuHad/Run2011A-PromptReco-v4
            itrigger = treeLoader.iTrigger("HLT_Mu17_TriCentralJet30_v5", currentRun, iFile);
          else if ( 165970<=currentRun && currentRun<=167043 )
            itrigger = treeLoader.iTrigger("HLT_Mu17_TriCentralJet30_v6", currentRun, iFile);
          else if ( 166346<=currentRun && currentRun<=166346 )
            itrigger = treeLoader.iTrigger("HLT_Mu17_TriCentralJet30_v7", currentRun, iFile);
          else if ( 167078<=currentRun && currentRun<=167913 )
            itrigger = treeLoader.iTrigger("HLT_Mu17_TriCentralJet30_v8", currentRun, iFile);
          else if ( 170722<=currentRun && currentRun<=172619 ) //MuHad/Run2011A-05Aug2011-v1
            itrigger = treeLoader.iTrigger("HLT_Mu17_TriCentralJet30_v10", currentRun, iFile);
          else if ( 172620<=currentRun && currentRun<=173198 ) //MuHad/Run2011A-PromptReco-v6
            itrigger = treeLoader.iTrigger("HLT_Mu17_TriCentralJet30_v10", currentRun, iFile);
          else if ( 173241<=currentRun && currentRun<=178380 ) //MuHad/Run2011B-PromptReco-v1
                                                               //        else if ( 175860<=currentRun && currentRun<=178380 ) //MuHad/Run2011B-PromptReco-v1
            itrigger = treeLoader.iTrigger("HLT_Mu17_eta2p1_TriCentralJet30_v1", currentRun, iFile);
          else if ( 178420<=currentRun && currentRun<=179889 )
            itrigger = treeLoader.iTrigger("HLT_IsoMu17_eta2p1_TriCentralJet30_v5", currentRun, iFile);
          else if ( 179959<=currentRun && currentRun<=180252 )
            itrigger = treeLoader.iTrigger("HLT_IsoMu17_eta2p1_TriCentralJet30_v6", currentRun, iFile);
          
          else {
            printf("Trigger : run %6d not in the range and we have data for it !!!\n", currentRun);
              //exit(1);
          }
          
        }
        else {
          itrigger = treeLoader.iTrigger ("HLT_Mu17_TriCentralJet30_v2", currentRun, iFile);
        }
        
        
        /*
         TRootHLTInfo hltInfos = treeLoader.runInfos->getHLTinfo(currentRun);
         for(std::vector<std::string>::iterator hltName = hltInfos.hltNames.begin() ; hltName != hltInfos.hltNames.end() ; hltName++)
         printf("%s, ",hltName->c_str());
         printf("\n");
         exit 1;
         */
        if(itrigger == 9999)
        {
          cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (" << dataSetName << ") IN RUN " << event->runId() << endl;
            //            exit(1);
        }
          //}
        
        trigged = treeLoader.EventTrigged (itrigger);
        
          //} else 
          //trigged=true;
        
          // JES CORRECTION
        
          // Clone the init_jets vector, otherwise the corrections will be removed
        for(unsigned int i=0; i<init_jets.size(); i++)
          if(init_jets[i]) delete init_jets[i];
        init_jets.clear();
        
        for(unsigned int i=0; i<init_jets_fromTopTree.size(); i++)
          init_jets.push_back( (TRootJet*) init_jets_fromTopTree[i]->Clone() );
        
        
          // Apply Jet Corrections on-the-fly
        if( dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
        {
          jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho(),true); //last boolean: isData (needed for L2L3Residual...)
          jetTools->correctMETTypeOne(init_jets,mets[0],true);
          
        }
        else
        {
          jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho(),false); //last boolean: isData (needed for L2L3Residual...)
          jetTools->correctMETTypeOne(init_jets,mets[0],false);
          
          
            // Fall11 samples only ?
          Double_t lumiWeight3D = 1.;
          lumiWeight3D = Lumi3DWeights.weight3D(event->nPu(-1),event->nPu(0),event->nPu(+1));
          eventWeight *= lumiWeight3D;
          if (lumiWeight3D > 100 )
            printf("3D Reweighting : Event with very high weight : weight=%lf", lumiWeight3D);
            //lumiWeight3D = Lumi3DWeights.weight3D(event);
          
          
          vector<TRootGenJet*> genjets;
          genjets = treeLoader.LoadGenJet(ievt);
            // JER systematic
          switch (uncertainty_JER) {
            case 0: jetTools->correctJetJER(init_jets, genjets, mets[0], "nominal", nsigma); break;
            case 1: jetTools->correctJetJER(init_jets, genjets, mets[0], "minus", nsigma); break;
            case 2: jetTools->correctJetJER(init_jets, genjets, mets[0], "plus", nsigma); break;
            default: std::cerr << "Wrong parameter for : " << "uncertainty_JER" << std::endl; break;
          }
          
            //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JER correction:");	  
          
            // JES systematic 
          switch (uncertainty_JES) {
            case 0: ; break;
            case 1: jetTools->correctJetJESUnc(init_jets, mets[0], "minus", nsigma); break;
            case 2: jetTools->correctJetJESUnc(init_jets, mets[0], "plus", nsigma); break;
            default: std::cerr << "Wrong parameter for : " << "uncertainty_JES" << std::endl; break;
          }
            //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JER correction:");
          
        }
        
        selecTable.Fill(d,0,eventWeight);
        
        string labelFired = "";
        ULong64_t triggerSignature = 0;
        for (vector<string>::reverse_iterator iter=triggerTestsList.rbegin(); iter!=triggerTestsList.rend(); iter++) {
          Int_t iTrig= treeLoader.iTrigger(*iter, currentRun, iFile);
          if (iTrig != 9999) {
            triggerSignature = (triggerSignature<<1) | (treeLoader.EventTrigged(iTrig));
            if (treeLoader.EventTrigged(iTrig))
              labelFired = labelFired+ " + " +(*iter);
            
          } else {
            triggerSignature = triggerSignature<<1;
          }
        }
        
        string label = "" ;
        ULong64_t j = triggerSignature;
          //      printf("triggerSignature : %llu\n", triggerSignature);
        for (vector<string>::iterator iter=triggerTestsList.begin(); iter!=triggerTestsList.end(); iter++) {
          if((j&1)==1) {
            if (label!="") {
              label = label + " && ";
            }
            label = label + (*iter);
          }
          j = (j >> 1);
        }
          //printf("Fired : %s            \t Name : %s\n", labelFired.c_str(), label.c_str());
        histo1D["TriggerTests_"+datasets[d]->Name()]->Fill(label.c_str(), eventWeight);
        
        
        /* Beam scraping events veto and HBHE filter noise are implemented in the TopBrussels/AutoMaticTopTreeProducer/ConfigTemplates/PatTemplate_*_data_*_cfg.py */
        if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
        {
            // Apply the scraping veto. (Is it still needed?)
          
          
          bool isBeamBG = true;
          if(event->nTracks() > 10)
          {
            if( ( (float) event->nHighPurityTracks() ) / ( (float) event->nTracks() ) > 0.25 )
              isBeamBG = false;
          }
            //          selecTable.Fill(d,1,eventWeight);
          if(isBeamBG) {
            printf("Beam scraping event rejected !!!\n");
            continue;
          }
        }
        
        /*
         *
         *
         *
         *
         */
        
          /////////////////////////////
          //   Selection
          /////////////////////////////
        
          //Declare selection instance    
        Selection selection(init_jets, init_muons, init_electrons, mets);
        
          // Ne teste que le premier PV
        Bool_t isGoodPV = selection.isPVSelected(vertex, 4, 24, 2.);  
        
        if (channel_Single_Muon) {
            //selection.setMuonCuts(float Pt, float Eta, float RelIso, int NValidHits, float d0, float DRJets, int NMatchedStations, float DistVzPVz, int NPixelLayersWithMeas);
            // critères manquants :
            //selection.setMuonCuts(20,2.1,0.125,10,0.02,0.3,1,1,1);
            //        selection.setMuonCuts(20,2.1,FLT_MAX,10,0.02,0.3,1,1,1);    // Selection of the beginning (Synchronisation exercise with Gent)
          selection.setMuonCuts(20.,2.1,0.125, 10 ,0.02,0.3,1,1.,1); //relIso can be 0.1 ;
          
          
          
            //selection.setLooseMuonCuts(float Pt, float Eta, float RelIso)
            //selection.setLooseMuonCuts(10,2.5,0.2);
            //        // Selection of the beginning (Synchronisation exercise with Gent) : no loose muon veto
          selection.setLooseMuonCuts(10,2.5,0.2);
            //selection.setElectronCuts(float Et, float Eta, float RelIso, float d0, float DistVzPVz, float DRJets);
            //selection.setElectronCuts(20,2.5,0.125,0.02,1,0.3);
            //selection.setElectronCuts(20,2.5,FLT_MAX,0.02,1,0.3);// Selection of the beginning (Synchronisation exercise with Gent) : Hard electron veto
            //          selection.setElectronCuts(15,2.5,0.2,0.02,1,0.3); //ATTENTION DistVzPVz and DRJets not mentionned in https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel (2012/03/07)
          selection.setLooseElectronCuts(15,2.5,0.2);
            //selection.setJetCuts(float Pt, float Eta, float EMF, float n90Hits, float fHPD, float dRJetElectron, float dRJetMuon);
            //selection.setJetCuts(30.,2.4,0.01,1.,0.98,0.3,0.1);
            //          selection.setJetCuts(30.,2.4,(-1.)*FLT_MAX,(-1.)*FLT_MAX,(-1.)*FLT_MAX,0.3,0.1);//e+jets //RELIRE
            //          selection.setJetCuts(30.,2.4,(-1.)*FLT_MAX,(-1.)*FLT_MAX,(-1.)*FLT_MAX,-1.,0.3); //RELIRE : 0.3 ou 0.1 pour PF ???
          selection.setJetCuts(30.,2.4,(-1.)*FLT_MAX,(-1.)*FLT_MAX,(-1.)*FLT_MAX,0.3,0.1); // Partie en DeltaR : prise du code de Stijn
          /*
           for(std::vector<TopTree::Jet>::iterator ij = jets.begin(); ij!=jets.end(); ij++){    
           if ( ij->v4.Pt()>PtJetCut_ && fabs(ij->v4.Eta())<2.5 && ij->isLOOSE) {jets_out.push_back(*ij);}
           }
           */
          
        } else if (channel_Single_Electron) {
          printf("Electron channel analysis\n");
        } else {
          printf("Neither Single Muon nor Single Electron channel\n");
        }
        
        
          // FILL THE SELECTION TABLE //
        Bool_t selected = false;
        Bool_t selectedForVJetEstimation = false;
        
        vector<TRootJet*> selectedJets;
        vector<TRootMuon*> selectedMuons;
        vector<TRootMuon*> vetoMuons;
        vector<TRootElectron*> vetoElectrons;
        
        if(isGoodPV && trigged){
          selecTable.Fill(d,1,eventWeight);
          
          selectedJets = selection.GetSelectedJets(true);///////// Vient du code de Stijn
                                                         //GetSelectedJets(PtThr, EtaThr, electrons, dRElectronCut, applyJetID)
                                                         //          selection.GetSelectedJets(jet_pT_cut, jet_eta_cut, selectedMuonsSelection, dRMuonCut, kTRUE); ///
          
          sort(selectedJets.begin(), selectedJets.end(), HighestPt());
          selectedMuons = selection.GetSelectedMuons(vertex[0], selectedJets);
          vetoMuons = selection.GetSelectedLooseMuons();
          
          /* // Selection of the beginning (Synchronisation exercise with Gent)
           for(vector<TRootMuon*>::const_iterator muon=selectedMuons.begin() ; muon != selectedMuons.end() ; muon++)
           {
           if ((*muon)->relativeIso03() < 0.15)
           selectedOldIsoMuons.push_back((TRootMuon*) (*muon));
           } 
           */
          vetoElectrons = selection.GetSelectedLooseElectrons(false);
          vector<Double_t> btagDiscs_selectedJets;
          for (vector<TRootJet*>::iterator jet = selectedJets.begin(); jet != selectedJets.end(); jet++) {
            btagDiscs_selectedJets.push_back((*jet)->btag_trackCountingHighEffBJetTags());
          }
          cutCP["Trigger"]->Fill(selectedMuons.size(), datasets[d], true, Luminosity*eventWeight);/** if scale==true: histo are scaled by data.NormFactor()*Lumi */
          for (vector<Double_t>::iterator btagDisc = btagDiscs_selectedJets.begin(); btagDisc != btagDiscs_selectedJets.end(); btagDisc++)
            btagDiscCP["Trigger"]->Fill(*btagDisc, datasets[d], true, Luminosity*eventWeight);
          if(selectedMuons.size()==1){
            selecTable.Fill(d,2,eventWeight);
            cutCP["IsolatedMuon"]->Fill(vetoElectrons.size(), datasets[d], true, Luminosity*eventWeight);
            for (vector<Double_t>::iterator btagDisc = btagDiscs_selectedJets.begin(); btagDisc != btagDiscs_selectedJets.end(); btagDisc++)
              btagDiscCP["IsolatedMuon"]->Fill(*btagDisc, datasets[d], true, Luminosity*eventWeight);
            
            if(vetoMuons.size()==1){
              selecTable.Fill(d,3,eventWeight);
              cutCP["IsolatedMuon"]->Fill(vetoElectrons.size(), datasets[d], true, Luminosity*eventWeight);
              for (vector<Double_t>::iterator btagDisc = btagDiscs_selectedJets.begin(); btagDisc != btagDiscs_selectedJets.end(); btagDisc++)
                btagDiscCP["IsolatedMuon"]->Fill(*btagDisc, datasets[d], true, Luminosity*eventWeight);
              
              
              if(vetoElectrons.size()==0){
                selecTable.Fill(d,4,eventWeight);
                selectedForVJetEstimation = true;
                if (selectedJets.size()>0)
                  cutCP["VetoElectron"]->Fill(selectedJets[0]->Pt(), datasets[d], true, Luminosity*eventWeight);
                for (vector<Double_t>::iterator btagDisc = btagDiscs_selectedJets.begin(); btagDisc != btagDiscs_selectedJets.end(); btagDisc++)
                  btagDiscCP["VetoElectron"]->Fill(*btagDisc, datasets[d], true, Luminosity*eventWeight);
                if(selectedJets.size() >= 1) {
                  selecTable.Fill(d,5,eventWeight);
                  if (selectedJets.size()>1)
                    cutCP["AtLeast1Jet"]->Fill(selectedJets[1]->Pt(), datasets[d], true, Luminosity*eventWeight);
                  for (vector<Double_t>::iterator btagDisc = btagDiscs_selectedJets.begin(); btagDisc != btagDiscs_selectedJets.end(); btagDisc++)
                    btagDiscCP["AtLeast1Jet"]->Fill(*btagDisc, datasets[d], true, Luminosity*eventWeight);
                  if(selectedJets.size() >= 2) {
                    selecTable.Fill(d,6,eventWeight);
                    if (selectedJets.size()>2)
                      cutCP["AtLeast2Jets"]->Fill(selectedJets[2]->Pt(), datasets[d], true, Luminosity*eventWeight);
                    for (vector<Double_t>::iterator btagDisc = btagDiscs_selectedJets.begin(); btagDisc != btagDiscs_selectedJets.end(); btagDisc++)
                      btagDiscCP["AtLeast2Jets"]->Fill(*btagDisc, datasets[d], true, Luminosity*eventWeight);
                    if(selectedJets.size() >= 3) {
                      selecTable.Fill(d,7,eventWeight);
                      selected = true;
                      if (selectedJets.size()>3)
                        cutCP["AtLeast3Jets"]->Fill(selectedJets[3]->Pt(), datasets[d], true, Luminosity*eventWeight);
                      for (vector<Double_t>::iterator btagDisc = btagDiscs_selectedJets.begin(); btagDisc != btagDiscs_selectedJets.end(); btagDisc++)
                        btagDiscCP["AtLeast3Jets"]->Fill(*btagDisc, datasets[d], true, Luminosity*eventWeight);
                      if(selectedJets.size() >= 4) {
                        selecTable.Fill(d,8,eventWeight);
                        if (selectedJets.size()>4)
                          cutCP["AtLeast4Jets"]->Fill(selectedJets[4]->Pt(), datasets[d], true, Luminosity*eventWeight);
                        for (vector<Double_t>::iterator btagDisc = btagDiscs_selectedJets.begin(); btagDisc != btagDiscs_selectedJets.end(); btagDisc++)
                          btagDiscCP["AtLeast4Jets"]->Fill(*btagDisc, datasets[d], true, Luminosity*eventWeight);
                        if(selectedJets.size() >= 5) {
                          selecTable.Fill(d,9,eventWeight);
                          /*                      cutCP["AtLeast5Jets"]->Fill(, datasets[d], eventWeight, true, Lumi = Luminosity*eventWeight); */ /** if scale==true: histo are scaled by data.NormFactor()*Lumi */
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
        if (selectedForVJetEstimation)
        {
          if (find(iDVLike.begin(),iDVLike.end(),(int) d) != iDVLike.end())
            histo1D["JetMultiplicities_Vlike"]->Fill(selectedJets.size(), eventWeight*currentDatasetNormFactor*Luminosity);
          if (find(iDTTLike.begin(),iDTTLike.end(),(int) d) != iDTTLike.end())
            histo1D["JetMultiplicities_TTlike"]->Fill(selectedJets.size(), eventWeight*currentDatasetNormFactor*Luminosity);
          if (find(iDVLike.begin(),iDVLike.end(),(int) d) != iDVLike.end())
            histo1D["b-JetMultiplicities_Vlike"]->Fill(selectedJets.size(), eventWeight*currentDatasetNormFactor*Luminosity);
          if (find(iDTTLike.begin(),iDTTLike.end(),(int) d) != iDTTLike.end())
            histo1D["b-JetMultiplicities_TTlike"]->Fill(selectedJets.size(), eventWeight*currentDatasetNormFactor*Luminosity);
          /*
           printf("Vl decision : %d (index %d)   TTl decision : %d (index %d)\n",
           (find(iDVLike.begin(),iDVLike.end(),d) != iDVLike.end()),
           *find(iDVLike.begin(),iDVLike.end(),d),
           (find(iDTTLike.begin(),iDTTLike.end(),d) != iDTTLike.end()),
           *find(iDTTLike.begin(),iDTTLike.end(),d));
           */
          if (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") { // Data!
            /**/
              //            printf("Fill data dataset : %u", n_data_datasets-1);
            vJetEst__data->Fill(selectedJets, d/*n_data_datasets-1*/, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            /**/
              //         vJetEst->BckgdSubstraction(MCObsExpBjetMult, BckgdNames, Luminosity);
          }
          if ((dataSetName.find("TT_")==0) || (dataSetName.find("Stop_")==0)) {
              //            printf("Fill MC dataset (TT) : %u", d_notData);
            vJetEst__tt__W->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            /**/
            vJetEst__tt__W__Z__SingleT__QCD->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            vJetEst__tt__W_Z->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            vJetEst__tt__W_Z->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            vJetEst__tt__W_Z__QCD->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            vJetEst__tt__W_Z__SingleT__QCD->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            
            vJetEst__tt_SingleT__W_Z__QCD->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            vJetEst__tt_SingleT__W_Z->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            /**/
          } 
          if (dataSetName.find("W_")==0) {
              //            printf("Fill MC dataset (W) : %u", d);
            vJetEst__tt__W->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            /**/
            vJetEst__tt__W__Z__SingleT__QCD->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            vJetEst__tt__W_Z->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            vJetEst__tt__W_Z->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            vJetEst__tt__W_Z__QCD->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            vJetEst__tt__W_Z__SingleT__QCD->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);

            vJetEst__tt_SingleT__W_Z__QCD->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            vJetEst__tt_SingleT__W_Z->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            /**/
          }
          if ((dataSetName.find("Z_")==0) || (dataSetName.find("DY_")==0)) {
            ; /**/
              //            printf("Fill MC dataset (Z) : %u", d);
            vJetEst__tt__W__Z__SingleT__QCD->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            vJetEst__tt__W_Z->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            vJetEst__tt__W_Z->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            vJetEst__tt__W_Z__QCD->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            vJetEst__tt__W_Z__SingleT__QCD->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            
            vJetEst__tt_SingleT__W_Z__QCD->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            vJetEst__tt_SingleT__W_Z->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            /**/
          }
          if (dataSetName.find("ST_")==0) {
            ; /**/
            vJetEst__tt__W__Z__SingleT__QCD->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            vJetEst__tt__W_Z__SingleT__QCD->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            
            vJetEst__tt_SingleT__W_Z__QCD->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            vJetEst__tt_SingleT__W_Z->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);

            /**/
          }
          if (dataSetName.find("QCD_")==0) {
            ; /**/
            vJetEst__tt__W__Z__SingleT__QCD->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            vJetEst__tt__W_Z__QCD->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            vJetEst__tt__W_Z__SingleT__QCD->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
            vJetEst__tt_SingleT__W_Z__QCD->Fill(selectedJets, d, TCHE,(Float_t) eventWeight*currentDatasetNormFactor*Luminosity);
/**/
          }
        }
        if (! selected)
          continue;
        
        
      } // event loop
      
      
        // delete calls
      delete jetTools;
      
      std::cerr << std::flush <<"\n";
    } // dataset loop
  }  
  
  std::ofstream stopReport; 
  stopReport.open((outputDirectory + "/StopBckg.tex").c_str()) ;
  stopReport << "\\documentclass{article}" << std::endl;
  stopReport << "\\usepackage{supertabular}" << std::endl;
  stopReport << "\\usepackage{pdflscape}" << std::endl;
  stopReport << "\\begin{document}" << std::endl;
    //Selection tables
    //  selecTable.SetLuminosity(5000.);
    //  selecTable.Scale(5000.);
  selecTable.TableCalculator(false, true, true, true);
  selecTable.Write(stopReport);
  
    // V-jets estimation
  std::vector<VJetEstimation*> vJetEst_MC_methods;// crashe si Xlike pointe sur un dataset, et que ce dataset n'est pas rempli
  vJetEst_MC_methods.push_back(vJetEst__tt__W_Z);//OK with W_Z
  vJetEst_MC_methods.push_back(vJetEst__tt__W__Z__SingleT__QCD);//OK with W_Z
  vJetEst_MC_methods.push_back(vJetEst__tt__W);//KO with W_Z ; OK with W
  /**/
  vJetEst_MC_methods.push_back(vJetEst__tt__W_Z__QCD);
  vJetEst_MC_methods.push_back(vJetEst__tt__W_Z__SingleT__QCD);
  
  vJetEst_MC_methods.push_back(vJetEst__tt_SingleT__W_Z__QCD);
  vJetEst_MC_methods.push_back(vJetEst__tt_SingleT__W_Z);
  /**/
  
  for (std::vector<VJetEstimation*>::iterator iter=vJetEst_MC_methods.begin(); iter!=vJetEst_MC_methods.end(); iter++) {
    std::string dirName = "";
    if (*iter == vJetEst__tt__W) {
      dirName="vJetEst__tt__W";
    } else if (*iter == vJetEst__tt__W__Z__SingleT__QCD) {
      dirName="vJetEst__tt__W__Z__SingleT__QCD";
    } else if (*iter == vJetEst__tt__W_Z) {
      dirName="vJetEst__tt__W_Z";
    } else if (*iter == vJetEst__tt__W_Z__QCD) {
      dirName="vJetEst__tt__W_Z__QCD";
    } else if (*iter == vJetEst__tt__W_Z__SingleT__QCD) {
      dirName="vJetEst__tt__W_Z__SingleT__QCD";
    } else if (*iter == vJetEst__tt_SingleT__W_Z__QCD) {
      dirName="vJetEst__tt_SingleT__W_Z__QCD";
    } else if (*iter == vJetEst__tt_SingleT__W_Z) {
      dirName="vJetEst__tt_SingleT__W_Z";
    } else {
      printf("Houston, we have a problem !\n");
      exit(1);
    }
    printf("\n\n\n\n\n\n     FOR %s\n\n", dirName.c_str());
    std::vector<Bool_t> vMask;
    for (std::vector<Dataset>::const_iterator dataset=vDatasets.begin(); dataset!=vDatasets.end(); dataset++) {
      if(dataset->Name() == "Data" || dataset->Name() == "data" || dataset->Name() == "DATA") { // Data!
        vMask.push_back(kFALSE);
      }
      else {
        vMask.push_back(kTRUE);
      }
    }
    printf("vMask = { ");
    for (std::vector<Bool_t>::const_iterator it=vMask.begin(); it!=vMask.end(); it++) {
      if (*it ==kTRUE) {
        printf("kTRUE, ");
      } else {
        printf("kFALSE, ");
      }
    }
      printf("};\n");
    printf("SetProcesses();\n");
/*    if ( *iter == vJetEst__tt__W || *iter == vJetEst__tt__W__Z__SingleT__QCD)
      (*iter)->SetProcesses(vMask, ttLikeNames, vLikeNames, vbLikeNames);
      else*/
    (*iter)->SetProcesses(vMask, ttLikeNames, vLikeNames, vbLikeNames);
    printf("ComputeEffFromMC();\n");
    (*iter)->ComputeEffFromMC();
    printf("ComputeEffbqFromMC();\n");
    (*iter)->ComputeEffbqFromMC();  // From ttbar sample
  }
    //    vJetEst__data ->ComputeEffFromMC();
    //    vJetEst__data->ComputeEffbqFromMC(iDTTLike[0]);  // From ttbar sample
  
  printf("\n\nInitial values for the VJetEstimation : \n");
  for (UInt_t nJet=0; nJet<NofJetBins; nJet++) {
    for (UInt_t iParam=0; iParam<3*NofBtagWorkingPoint+2; iParam++) {
      if (iParam != 0)
        printf(" ; ");
      printf("%lf", initForMinVJet[nJet][iParam]);
    }
    printf("\n");
  }
  printf("\n\n\n");
  /**/ 
  vJetEst__data ->SumOverAllInputs();
  /**/
  for (std::vector<VJetEstimation*>::iterator iter=vJetEst_MC_methods.begin(); iter!=vJetEst_MC_methods.end(); iter++) {
    (*iter)->SumOverAllInputs();
  }
  /**/
  vJetEst__data ->PrintInputs();
  /**/
  for (std::vector<VJetEstimation*>::iterator iter=vJetEst_MC_methods.begin(); iter!=vJetEst_MC_methods.end(); iter++) {
    (*iter)->PrintInputs();
  }
  
  
  /*
   //OLD Method  
   vJetEst.SetInitialValues((Double_t**) initForMinVJet);
   */
  vector<Double_t> init;
    //  Double_t initForMinVJetArray[NofJetBins][2+3*NofBtagWorkingPoint]
    //init.clear();  for (UInt_t i=0; i<NofJetBins; i++) { init.push_back(initForMinVJet[i][0]); }
  init.clear();  for (UInt_t i=0; i<NofJetBins; i++) { init.push_back(vJetEst__tt__W->GetPredNtt(0,i)); }
  /**/
  vJetEst__data->SetInitialValues_Nttlike(init); //[njets][btagIdx]
  /**/
  for (std::vector<VJetEstimation*>::iterator iter=vJetEst_MC_methods.begin(); iter!=vJetEst_MC_methods.end(); iter++) {
    (*iter)->SetInitialValues_Nttlike(init); //[njets][btagIdx]
  }
    //  init.clear();  for (UInt_t i=0; i<NofJetBins; i++) { init.push_back(initForMinVJet[i][1]); }
  init.clear();  for (UInt_t i=0; i<NofJetBins; i++) { init.push_back(vJetEst__tt__W->GetPredNv(0,i)); }
  /**/
  vJetEst__data->SetInitialValues_Nvlike(init); //[njets][btagIdx]
  /**/
  for (std::vector<VJetEstimation*>::iterator iter=vJetEst_MC_methods.begin(); iter!=vJetEst_MC_methods.end(); iter++) {
    (*iter)->SetInitialValues_Nvlike(init); //[njets][btagIdx]
  }
  vector< vector<Double_t> > init2;
  init2.clear();
  for (UInt_t i=0; i<NofJetBins; i++) { 
    vector<Double_t> tmpInit ;
      //    for (UInt_t j=0; j<NofBtagWorkingPoint; j++) { tmpInit.push_back(initForMinVJet[i][2+3*j]); }
    for (UInt_t j=0; j<NofBtagWorkingPoint; j++) { tmpInit.push_back(vJetEst__tt__W->GetPredEb(j,i)); }
    init2.push_back(tmpInit);
  }
  /**/
  vJetEst__data->SetInitialValues_Eb(init2); //[njets][btagIdx]
  /**/
  for (std::vector<VJetEstimation*>::iterator iter=vJetEst_MC_methods.begin(); iter!=vJetEst_MC_methods.end(); iter++) {
    (*iter)->SetInitialValues_Eb(init2); //[njets][btagIdx]
  }
  
  init2.clear();
  for (UInt_t i=0; i<NofJetBins; i++) { 
    vector<Double_t> tmpInit ;
      //    for (UInt_t j=0; j<NofBtagWorkingPoint; j++) { tmpInit.push_back(initForMinVJet[i][2+1+3*j]); }
    for (UInt_t j=0; j<NofBtagWorkingPoint; j++) { tmpInit.push_back(vJetEst__tt__W->GetPredEudsc(j,i)); }
    init2.push_back(tmpInit);
  }
  /**/
  vJetEst__data->SetInitialValues_Eudsc(init2); //[njets][btagIdx]
  /**/
  for (std::vector<VJetEstimation*>::iterator iter=vJetEst_MC_methods.begin(); iter!=vJetEst_MC_methods.end(); iter++) {
    (*iter)->SetInitialValues_Eudsc(init2); //[njets][btagIdx]
  }
  init2.clear();
  for (UInt_t i=0; i<NofJetBins; i++) { 
    vector<Double_t> tmpInit ;
      //    for (UInt_t j=0; j<NofBtagWorkingPoint; j++) { tmpInit.push_back(initForMinVJet[i][2+2+3*j]); }
    for (UInt_t j=0; j<NofBtagWorkingPoint; j++) { tmpInit.push_back(vJetEst__tt__W->GetPredEuds(j,i)); }
    init2.push_back(tmpInit);
  }
  /**/
  vJetEst__data->SetInitialValues_Euds(init2); //[njets][btagIdx]
  /**/
  for (std::vector<VJetEstimation*>::iterator iter=vJetEst_MC_methods.begin(); iter!=vJetEst_MC_methods.end(); iter++) {
    (*iter)->SetInitialValues_Euds(init2); //[njets][btagIdx]
  }
  
    //  vector<int> fixedVar;    /*  0 : b-tag eff    ;    1 : mistag eff (udsc tagged as b)    ;    2 : mistag eff light (uds tagged as b) */
    //fixedVar.push_back(0);
    //  Bool_t doMinos = true;
  TFile *fout = new TFile (outputRootFileName.c_str(), "RECREATE");
    ////////////vJetEst->FillSummaryHistos();
    //  for (UInt_t index_wp=0; index_wp< NofBtagWorkingPoint ; index_wp++) {
    //  vJetEst__tt__W.UnBinnedMaximumNjetsLikelihoodEst(index_wp, fixedVar, doMinos, verbose);
    //vJetEst__tt__W->FillSummaryHistos();
  vJetEst__tt__W->PrintResults();
  vJetEst__tt__W->PrintResults_LatexFormat(stopReport, 5, 6);
  vJetEst__tt__W->PrintResults_LatexFormat(5, 6);
  
  /**/
   vJetEst__data->FillSummaryHistos();
  std::string dirName = "vJetEst__data";
  system( (std::string("mkdir -p ")+outputDirectory+"/"+dirName+" && mv "+outputDirectory+"/hNbjetsSummary_*.pdf "+outputDirectory+"/"+dirName+" ;").c_str());
   /**/
  for (std::vector<VJetEstimation*>::iterator iter=vJetEst_MC_methods.begin(); iter!=vJetEst_MC_methods.end(); iter++) {
    (*iter)->FillSummaryHistos();
    std::string dirName = "";
    if (*iter == vJetEst__tt__W) {
      dirName="vJetEst__tt__W";
    } else if (*iter == vJetEst__tt__W__Z__SingleT__QCD) {
      dirName="vJetEst__tt__W__Z__SingleT__QCD";
    } else if (*iter == vJetEst__tt__W_Z) {
      dirName="vJetEst__tt__W_Z";
    } else if (*iter == vJetEst__tt__W_Z__QCD) {
      dirName="vJetEst__tt__W_Z__QCD";
    } else if (*iter == vJetEst__tt__W_Z__SingleT__QCD) {
      dirName="vJetEst__tt__W_Z__SingleT__QCD";
    } else if (*iter == vJetEst__tt_SingleT__W_Z__QCD) {
      dirName="vJetEst__tt_SingleT__W_Z__QCD";
    } else if (*iter == vJetEst__tt_SingleT__W_Z) {
      dirName="vJetEst__tt_SingleT__W_Z";
    } else {
      printf("Houston, we have a problem !\n");
      exit(1);
    }
    system( (std::string("mkdir -p ")+outputDirectory+"/"+dirName+" && mv "+outputDirectory+"/hNbjetsSummary_*.pdf "+outputDirectory+"/"+dirName+" ;").c_str());
  }
  /*  vJetEst.Write(fout, "_allJetsInclusive_", (verbose>0));
   fout->cd();*/
  
    //  for (int njet=NofJets; njet<NofJets+NofJetBins; njet++) {
  /*    for (UInt_t njet=NofJetBins-1; njet<NofJetBins; njet++) {
   cout << std::endl << std::endl;
   stopReport << std::endl << std::endl;
   bool doPEwithRoofit = false;
   vJetEst.UnBinnedMaximumLikelihoodEst(index_wp, njet, fixedVar, doMinos, doPEwithRoofit, verbose);
   vJetEst.FillSummaryHistos();
   vJetEst.PrintResults();
   stopReport << "\\begin{landscape}" << std::endl;     
   vJetEst.PrintResults_LatexFormat(stopReport);
   stopReport << "\\end{landscape}" << std::endl;     
   stopReport << "\\pagebreak" << std::endl << std::endl;     
   char saveName[100]; sprintf(saveName, "VJetEstimation_%djets_TCHE_wp%f", njet+NofJets, BtagWorkingPoint[index_wp]);
   //      vJetEst.Write(fout, directorySuffix, (verbose>0));
   fout->cd();
   vJetEst.Write(saveName);
   
   }
   } */
  fout->cd();
  
  /**/
  vJetEst__data->Write("VJetEstimation-TCHE_LM--Data");
  /**/
  vJetEst__tt__W->Write("VJetEstimation-TCHE_LM--tt_W");
  /**/
  vJetEst__tt__W__Z__SingleT__QCD->Write("VJetEstimation-TCHE_LM--tt_W_Z_SingleT_QCD");
  vJetEst__tt__W_Z->Write("VJetEstimation-TCHE_LM--tt_W-Z");
  vJetEst__tt__W_Z__QCD->Write("VJetEstimation-TCHE_LM--tt_W-Z_QCD");
  vJetEst__tt__W_Z__SingleT__QCD->Write("VJetEstimation-TCHE_LM--tt_W-Z_SingleT_QCD");
  
  vJetEst__tt_SingleT__W_Z__QCD->Write("VJetEstimation-TCHE_LM--tt-SingleT_W-Z_QCD");
  vJetEst__tt_SingleT__W_Z->Write("VJetEstimation-TCHE_LM--tt-SingleT_W-Z");
  /**/
  
  stopReport << "\\end{document}" << std::endl;
  stopReport.close();
  
  fout->cd();
  
    //  TDirectory* th1dir = fout->mkdir("1D_histograms");
  
  fout->cd();
    //  th1dir->cd();
  
  if (runOnEvents_) {
    for(std::map<std::string,TH1*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++) {
      TH1 *temp = it->second;
        //int N = temp->GetNbinsX();
        //temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
        //temp->SetBinContent(N+1,0);
      temp->Write();
      TCanvas* tempCanvas = TCanvasCreator((TH1F*)temp, it->first);
      tempCanvas->SaveAs( (outputDirectory + "/"+it->first+".root").c_str() );
      if (it->first.find("TriggerTests_")==0) {
        UInt_t numBins = 0;
        for (Int_t bin = 1; bin <= temp->GetNbinsX(); bin++) {
          if (temp->GetBinContent(bin) != 0.) {
            numBins++;
          }
        }
        printf("\n%s (%u/%d)\n", it->first.c_str(), numBins, temp->GetNbinsX());
        for (Int_t bin = 1; bin <= temp->GetNbinsX(); bin++) {
          if (temp->GetBinContent(bin) != 0.) {
            printf("  %s \t : %10lf\n", temp->GetXaxis()->GetBinLabel(bin), temp->GetBinContent(bin));
          }
        }
      }
      
    }
    
    /** Saving control plots */
    fout->cd();
    for (map<string,MultiSamplePlot*>::iterator iter=cutCP.begin(); iter!=cutCP.end(); iter++) {
        //    iter->second->Draw(addRandomPseudoData = false, label = string("CMSPlot"), mergeTT = false, mergeQCD = false, mergeW = false, mergeZ = false, mergeST = false, scaleNPSignal = 1, addRatio = false, mergeVV = false, mergeTTV = false) ;
      iter->second->Draw(false, string("a plot"), false, false, false, false, true, 1, false, false, false) ;
      iter->second->Write(fout, iter->first, true, string("pdf")) ;
      delete iter->second ;
    }
    
    for (map<string,MultiSamplePlot*>::iterator iter=btagDiscCP.begin(); iter!=btagDiscCP.end(); iter++) {
      iter->second->Draw(false, string("a plot"), false, false, false, false, true, 1, false, false, false) ;
      iter->second->Write(fout, iter->first, true, string("pdf")) ;
      delete iter->second ;
    }
  }
    //fout->ls();
    //fout->Map();
  fout->Close();
  
  TFile *ftest = new TFile (outputRootFileName.c_str(), "READ");
    //ftest->ls();
    //ftest->Map();
  printf("Trying to read back\n");
  VJetEstimation *vje = NULL;
  vje = (VJetEstimation*) ftest->Get("VJetEstimation-TCHE_LM--Data");
  vje = (VJetEstimation*) ftest->Get("VJetEstimation-TCHE_LM--tt_W");
  vje = (VJetEstimation*) ftest->Get("VJetEstimation-TCHE_LM--tt_W_Z_SingleT_QCD");
  vje = (VJetEstimation*) ftest->Get("VJetEstimation-TCHE_LM--tt_W-Z");
  vje = (VJetEstimation*) ftest->Get("VJetEstimation-TCHE_LM--tt_W-Z_QCD");
  vje = (VJetEstimation*) ftest->Get("VJetEstimation-TCHE_LM--tt_W-Z_SingleT_QCD");
  vje = (VJetEstimation*) ftest->Get("VJetEstimation-TCHE_LM--tt-SingleT_W-Z_QCD");
  vje = (VJetEstimation*) ftest->Get("VJetEstimation-TCHE_LM--tt-SingleT_W-Z");

  printf("Could read back\n");
  ftest->Close();
  
  
  printf("It took us %7lgs to run the program\n", ((Double_t)clock() - start) / CLOCKS_PER_SEC);
  
}

