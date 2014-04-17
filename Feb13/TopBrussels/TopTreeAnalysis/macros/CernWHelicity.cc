///////////////////////////
///// TODO & COMMENTS /////
/////////////////////////// 

#include "TStyle.h"
#include "TGraphErrors.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <math.h>
#include "TLorentzVector.h" //Needs to be defined as " TLorentzVector V(Px,Py,Pz,E) " or " TLorentzVector V and V.SetPxPyPzE(Px,Py,Pz,E) "
#include "TF1.h"
#include "TH2.h"
#include "TString.h"
#include <stdio.h>
#include "TMath.h"

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "../Selection/interface/SelectionTable.h"
#include "../Tools/interface/PlottingTools.h"
#include "../Tools/interface/MultiSamplePlot.h"
#include "../Tools/interface/TTreeLoader.h"
#include "../Tools/interface/AnalysisEnvironmentLoader.h"
#include "../Content/interface/AnalysisEnvironment.h"
#include "../Content/interface/Dataset.h"
#include "../MCInformation/interface/JetPartonMatching.h"
#include "../Tools/interface/JetTools.h"
#include "Style.C"
#include "TVector3.h"
#include "TH1.h"
#include "Riostream.h"

//Class for reconstruction of topology
#include "../WHelicities/interface/WReconstruction.h"

//Includes necessary for kinFitter:
#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdio>
#include <algorithm>

#include "TMatrixD.h"

#include "../MCInformation/interface/ResolutionFit.h"
#include "../KinFitter/interface/TKinFitter.h"
#include "../KinFitter/interface/TFitConstraintM.h"
#include "../KinFitter/interface/TFitParticleEtThetaPhiEMomFix.h"

using namespace std;
using namespace TopTree;

int FirstEvent=0;

int main (int argc, char *argv[])
{
  clock_t start = clock();
  
  cout << "*****************************************" << endl;
  cout << " Beginning of the program for Cern WHelicity ! " << endl;
  cout << "*****************************************" << endl;
  
  //SetStyle if needed
  //setTDRStyle(); 
  setMyStyle();
   
  /////////////////////
  // Configuration
  /////////////////////
  
  bool doPF2PAT = false;
  
  //xml file
  string xmlFileName ="../config/CernMacroConfig.xml";   //Position of .xml file
  if (argc >= 2)
    xmlFileName=string(argv[1]);
  const char *xmlfile = xmlFileName.c_str();
  
  cout << "used config file: " << xmlfile << endl;
  
  //Output ROOT file
  string rootFileName ("MacroOutputPValue.root");         //Root output file
  const char *rootfile = rootFileName.c_str();
  
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
  float oldLuminosity = anaEnv.Luminosity;	// in 1/pb
  
  cout << "analysis environment luminosity for rescaling "<< oldLuminosity << endl;
  
  /////////////////////
  // Load Datasets
  /////////////////////
  bool SimulationSampleBoolean;
  
  TTreeLoader treeLoader;
  cout << " - Load datasets ..." << endl;
  vector < Dataset* > datasets;
  treeLoader.LoadDatasets (datasets, xmlfile);
  for(unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
  
  float Luminosity = oldLuminosity;
  bool isSemiMu = false;
  for (unsigned int d = 0; d < datasets.size (); d++){
    cout << "luminosity of dataset "<< d << " is " << datasets[d]->EquivalentLumi() << endl;
    if(Luminosity > datasets[d]->EquivalentLumi() ) Luminosity = datasets[d]->EquivalentLumi();
    string dataSetName = datasets[d]->Name();
    if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA"){
      Luminosity = datasets[d]->EquivalentLumi();
      break;
    }
    if(dataSetName.find("TTbarJets_SemiMu") == 0) isSemiMu = true;
  }
  if(Luminosity != oldLuminosity) cout << "changed analysis environment luminosity to "<< Luminosity << endl;
  
  //vector of objects
  cout << " - Variable declaration ..." << endl;          //All possible variables that can be used for analysis
  vector < TRootVertex* > vertex; 
  vector < TRootMuon* > init_muons;
  vector < TRootElectron* > init_electrons;
  vector < TRootJet* > init_jets;
  vector < TRootJet* > init_jets_corrected;
  vector < TRootMET* > mets;
  
  TFile *fout = new TFile (rootfile, "RECREATE");    //Making of a TFile where all information is stored
  //Global variable
  TRootEvent* event = 0;
  
  //nof selected events
  float NEvtsData = 0;
  Double_t *nEvents = new Double_t[datasets.size()];
  
  ////////////////////////////////////
  /// TGraphAsymmErrors
  ////////////////////////////////////
  
  map<string,TGraphAsymmErrors*> graphAsymmErr;         
  map<string,TGraphErrors*> graphErr;
  
  std::cout << " starting to define the histograms !! " << std::endl;
  ////////////////////////////////////
  /// Normal Plots (TH1F* and TH2F*)
  ////////////////////////////////////
  //All histograms can be defined as pre-programmed maps which makes definitions and looping easier
  map<string,TH1F*> histo1D;     
  map<string,TH2F*> histo2D;  
  
  // Histograms needed to calculate the sigma and Mc mass (from mean value) for W and top mass distribution
  //   --> Comment out after initializing most recent values ( also lines 1046 and 1356 )
  histo1D["WMass"]= new TH1F("WMass","WMass", 200,0,160);
  histo1D["TopMass"]= new TH1F("TopMass","TopMass", 200,0,350);

  // Give the theoretical distribution for SM values
  //  --> Always stay the same 
  //
  //histo1D["HelicityDistribution"] = new TH1F("HelicityDistribution","HelicityDistribution",200,-1,1);
  //histo1D["HelicityDistributionLong"] = new TH1F("HelicityDistributionLong","HelicityDistributionLong",200,-1,1);
  //histo1D["HelicityDistributionLeft"] = new TH1F("HelicityDistributionLeft","HelicityDistributionLeft",200,-1,1);
  //histo1D["HelicityDistributionRight"] = new TH1F("HelicityDistributionRight","HelicityDistributionRight",200,-1,1);  
  
  // Give the Cos theta distribution for Semi Mu sample before event selection
  // Fit histograms gives the helicity values that should be used as SM values from CMS data
  // --> Always stay the same 
  //
  histo1D["StandardCosTheta"]=new TH1F("StCosTheta","StCosTheta",200,-1,1);
  histo1D["StandardCosThetaFit"]=new TH1F("StCosThetaFit","StCosThetaFit",200,-1,1);   

  //Histograms which check the position in the descending order of the jets:
  //
  histo1D["Quark1Number"]=new TH1F("Quark1Number","Quark1Number",11,-0.5,10.5);
  histo1D["Quark2Number"]=new TH1F("Quark2Number","Quark2Number",11,-0.5,10.5);
  histo1D["BHadronicNumber"]=new TH1F("BHadronicNumber","BHadronicNumber",11,-0.5,10.5);
  histo1D["BLeptonicNumber"]=new TH1F("BLeptonicNumber","BLeptonicNumber",11,-0.5,10.5);  

  //Histograms to check the effect of the reconstruction on the helicity distribution:
  //
  histo1D["CosThetaRECO"]=new TH1F("CosThetaRECO","CosThetaRECO",200,-1,1);
  histo1D["CosThetaRECOBLeptCorrect"]=new TH1F("CosThetaRECOBLeptCorrect","CosThetaRECOBLeptCorrect",200,-1,1);
  histo1D["CosThetaBLept"]=new TH1F("CosThetaBLept","CosThetaBLept",200,-1,1);
  histo1D["CosThetaMuon"]=new TH1F("CosThetaMuon","CosThetaMuon",200,-1,1);
  histo1D["CosThetaNeutrino"]=new TH1F("CosThetaNeutrino","CosThetaNeutrino",200,-1,1);

  histo1D["TopMassForCorrect"]=new TH1F("TopMassForCorrect","TopMassForCorrect",200,0,400);
  histo1D["TopMassForWrong"]=new TH1F("TopMassForWrong","TopMassForWrong",500,0,1000);
  histo1D["TopMassReco"]=new TH1F("TopMassReco","TopMassReco",200,0,400);
  histo1D["TopMassNotWrong"]=new TH1F("TopMassNotWrong","TopMassNotWrong",500,0,1000);

  histo2D["BtagForTopWrong"]=new TH2F("BtagForTopWrong","BtagForTopWrong",200,0,400,100,0,1);
  histo2D["BtagForTopCorrect"]=new TH2F("BtagForTopCorrect","BtagForTopCorrect",200,0,400,100,0,1);
  histo2D["BtagForTop"]=new TH2F("BtagForTop","BtagForTop",200,0,400,100,0,1);
  histo2D["ChiSquaredForTop"]=new TH2F("ChiSquaredForTop","ChiSquaredForTop",200,0,400,100,0,150);  
  histo2D["MVACorrelationPlot"]=new TH2F("MVACorrelationPlot","MVACorrelationPlot",100,0,150,100,0,1);
  histo2D["MVACorrelationPlotWrong"]=new TH2F("MVACorrelationPlotWrong","MVACorrelationPlotWrong",100,0,150,100,0,1);

  histo1D["MVACombination"]=new TH1F("MVACombination","MVACombination",200,0,12);
  histo1D["MVACombinationBefore"]=new TH1F("MVACombinationBefore","MVACombinationBefore",200,0,12);
  histo1D["MVACombinationTwo"]=new TH1F("MVACombinationTwo","MVACombinationTwo",200,0,12);

  int CorrectZero=0;
  int CorrectTwo=0;
  int CorrectFour=0;
  int CorrectSix=0;
  int CorrectEight=0;
  int CorrectTen=0;
  histo1D["CorrectMVANumber"]=new TH1F("CorrectMVANumber","CorrectMVANumber",200,0,12);
  bool bTagSelected;
  bool ProblemEvent;

  histo1D["QuarkOnePtZero"]=new TH1F("QuarkOnePtZero","QuarkOnePtZero",200,0,400);
  histo1D["QuarkOnePtTwo"]=new TH1F("QuarkOnePtTwo","QuarkOnePtTwo",200,0,400);
  histo1D["QuarkOnePtFour"]=new TH1F("QuarkOnePtFour","QuarkOnePtFour",200,0,400);
  histo1D["QuarkOnePtSix"]=new TH1F("QuarkOnePtSix","QuarkOnePtSix",200,0,400);
  histo1D["QuarkOnePtEight"]=new TH1F("QuarkOnePtEight","QuarkOnePtEight",200,0,400);
  histo1D["QuarkOnePtTen"]=new TH1F("QuarkOnePtTen","QuarkOnePtTen",200,0,400);

  histo1D["QuarkTwoPtZero"]=new TH1F("QuarkTwoPtZero","QuarkTwoPtZero",200,0,400);
  histo1D["QuarkTwoPtTwo"]=new TH1F("QuarkTwoPtTwo","QuarkTwoPtTwo",200,0,400);
  histo1D["QuarkTwoPtFour"]=new TH1F("QuarkTwoPtFour","QuarkTwoPtFour",200,0,400);
  histo1D["QuarkTwoPtSix"]=new TH1F("QuarkTwoPtSix","QuarkTwoPtSix",200,0,400);
  histo1D["QuarkTwoPtEight"]=new TH1F("QuarkTwoPtEight","QuarkTwoPtEight",200,0,400);
  histo1D["QuarkTwoPtTen"]=new TH1F("QuarkTwoPtTen","QuarkTwoPtTen",200,0,400);

  histo1D["HadrBPtZero"]=new TH1F("HadrBPtZero","HadrBPtZero",200,0,400);
  histo1D["HadrBPtTwo"]=new TH1F("HadrBPtTwo","HadrBPtTwo",200,0,400);
  histo1D["HadrBPtFour"]=new TH1F("HadrBPtFour","HadrBPtFour",200,0,400);
  histo1D["HadrBPtSix"]=new TH1F("HadrBPtSix","HadrBPtSix",200,0,400);
  histo1D["HadrBPtEight"]=new TH1F("HadrBPtEight","HadrBPtEight",200,0,400);
  histo1D["HadrBPtTen"]=new TH1F("HadrBPtTen","HadrBPtTen",200,0,400);

  histo1D["LeptBPtZero"]=new TH1F("LeptBPtZero","LeptBPtZero",200,0,400);
  histo1D["LeptBPtTwo"]=new TH1F("LeptBPtTwo","LeptBPtTwo",200,0,400);
  histo1D["LeptBPtFour"]=new TH1F("LeptBPtFour","LeptBPtFour",200,0,400);
  histo1D["LeptBPtSix"]=new TH1F("LeptBPtSix","LeptBPtSix",200,0,400);
  histo1D["LeptBPtEight"]=new TH1F("LeptBPtEight","LeptBPtEight",200,0,400);
  histo1D["LeptBPtTen"]=new TH1F("LeptBPtTen","LeptBPtTen",200,0,400);

  int BTagCombinationEvents = 0;
  int BTagCorrectCombination = 0;

  ////////////////////////////////////
  /// MultiSamplePlot
  ////////////////////////////////////
  //Whenever there are more than one dataset this can be used to add up the different results of the datasets
  //
  map<string,MultiSamplePlot*> MSPlot;       
  
  //Histograms for chapter event selection: Compare Data vs MC for Cut values!
  //
  
  // Histograms that give the chi squared and BTagValueMVA values for the jet distribution selection method
  //
  MSPlot["ChiSquared"]= new MultiSamplePlot(datasets, "ChiSquared", 20,-1,150,"ChiSquared"); 
  MSPlot["NotSelectedChiSq"]=new MultiSamplePlot(datasets,"NotSelectedChiSq",100,-1,200,"NotSelectedChiSq");
  MSPlot["BTagValueMVA"]=new MultiSamplePlot(datasets,"BTagValueMVA",100,0,1,"BTagValueMVA");
  MSPlot["NotSelectedBTag"]=new MultiSamplePlot(datasets,"NotSelectedBTag",100,0,1,"NotSelectedBTag");
  //
  // For SemiMu MC sample can be checked how many times selected jet distribution is correct:
  //
  MSPlot["ChiSquaredCorrect"]= new MultiSamplePlot(datasets, "ChiSquaredCorrect", 200,-1,150,"ChiSquaredCorrect");
  MSPlot["ChiSquaredWrong"]= new MultiSamplePlot(datasets, "ChiSquaredWrong", 200,-1,150,"ChiSquaredWrong"); 
  MSPlot["BTagValueMVACorrect"]=new MultiSamplePlot(datasets,"BTagValueMVACorrect",100,0,1,"BTagValueMVACorrect");
  MSPlot["BTagValueMVAWrong"]=new MultiSamplePlot(datasets,"BTagValueMVAWrong",100,0,1,"BTagValueMVAWrong");

  MSPlot["RecoWMass"]= new MultiSamplePlot(datasets,"RecoWMass",100,20,250,"RecoWMass");
  MSPlot["RecoTopHadrMass"]= new MultiSamplePlot(datasets,"RecoTopHadrMass",400,60,600,"RecoTopHadrMass"); 
  MSPlot["RecoWMassDiff"]= new MultiSamplePlot(datasets,"RecoWMassDiff",40,-70,200,"RecoWMassDiff");
  MSPlot["RecoTopHadrMassDiff"]= new MultiSamplePlot(datasets,"RecoTopHadrMassDiff",80,-120,500,"RecoTopHadrMassDiff");
  
  int NumberSemiMuEvents = 0;
  int AllJetsCorrect = 0;
  int OneJetWrong = 0;
  int numberCorrectHadronicB = 0;
  int numberCorrectQuarks = 0;  

  //----------------------------------------------------------------
  //--   Definitions for comparing CosTheta between data and Mc   --
  //----------------------------------------------------------------
  //Luminosity should be rescaled  

  // Number of bins should be 10 since there are only about 200 data entries after all selections
  //
  int CosThetaBinNumber = 15;  
  int NumberOfHelicityBins=100; 
  MSPlot["CosTheta"]= new MultiSamplePlot(datasets, "CosTheta", CosThetaBinNumber,-1,1,"CosTheta");  
  histo1D["CosThetaMC"] = new TH1F("CosThetaMC","CosThetaMC",CosThetaBinNumber,-1,1);
  histo1D["CosThetaData"]=new TH1F("CosThetaData","CosThetaData",CosThetaBinNumber,-1,1);  

  std::string HelicityNumber[(NumberOfHelicityBins+1)];   
  for(int ii=0;ii<=NumberOfHelicityBins;ii++){
    std::stringstream out;
    out << ii;
    HelicityNumber[ii] = out.str();
  }
  
  histo1D["TheoreticalHelicityDistribution"]= new TH1F("TheoreticalHelicityDistribution","TheoreticalHelicityDistribution",10,-1,1);
  histo1D["ChiSquaredReweighting"]=new TH1F("ChiSquaredReweighting","ChiSquaredReweighting",50,0,100);

  for(int Longit =0;Longit <=NumberOfHelicityBins;Longit++){
    for(int Right =0; Right <=(NumberOfHelicityBins-Longit) ; Right++){
      int Left = NumberOfHelicityBins-Longit-Right;
      
      std::string HistoName = "CosThetaRight"+HelicityNumber[Right]+"Long"+HelicityNumber[Longit]+"Left"+HelicityNumber[Left];
      TString THistoName = "CosThetaRight"+HelicityNumber[Right]+"Long"+HelicityNumber[Longit]+"Left"+HelicityNumber[Left];
      std::string HelicityWeightHisto = "HelicityWeightRight"+HelicityNumber[Right]+"Long"+HelicityNumber[Longit]+"Left"+HelicityNumber[Left];
      TString THelicityWeightHisto = "HelicityWeightRight"+HelicityNumber[Right]+"Long"+HelicityNumber[Longit]+"Left"+HelicityNumber[Left];
      std::string ChiSquaredHelicity = "ChiSquaredRight"+HelicityNumber[Right]+"Long"+HelicityNumber[Longit]+"Left"+HelicityNumber[Left];
      
      histo1D[HistoName] = new TH1F(THistoName,THistoName,CosThetaBinNumber,-1,1);           
      histo1D[HelicityWeightHisto] =new TH1F(THelicityWeightHisto,THelicityWeightHisto,200,0,5);  

    }//end of second helicity loop
  }// end of first helicity loop
  
  // Histograms that give the p-value for the different helicity values
  //  
  float HighestPValue =0;
  histo1D["PValueAllBins"]=new TH1F("PValueAllBins","PValueAllBins",5155,0,5155);
  float HelicityPointsBin = (1./(NumberOfHelicityBins*2.));
  histo2D["HelicityPointsPlot"]=new TH2F("HelicityPointsPlot","HelicityPointsPlot",(NumberOfHelicityBins+1),(0.-HelicityPointsBin),(1.+HelicityPointsBin),(NumberOfHelicityBins+1),(0.-HelicityPointsBin),(1.+HelicityPointsBin)); 
  histo2D["ChiSqPlot"]=new TH2F("ChiSqPlot","ChiSqPlot",(NumberOfHelicityBins+1),(0.-HelicityPointsBin),(1.+HelicityPointsBin),(NumberOfHelicityBins+1),(0.-HelicityPointsBin),(1.+HelicityPointsBin)); 
  float LowestChiSqValue=1000000000000000;

  ///////////////////////////////////////////////////////
  //   Making the theoretical helicity distribution    // -> To know how the distribution should look like!! 
  ///////////////////////////////////////////////////////  --> Only Fit values are used in code!!
  //------   Standard Model   -----
  float LongitudinalFraction = 0.642324; // values obtained from fit
  float RightHandedFraction = 0.0327637;
  float LeftHandedFraction = 1 - LongitudinalFraction - RightHandedFraction;
  float SMDiffDistribution[3];//0:Longitudinal; 1: Righthanded; 2: Lefthanded
  float XVariable;
  // std::cout << " Used SM helicity variables: Longitudinal (f0) = " << LongitudinalFraction << " Righthanded (f+) = " << RightHandedFraction << " Lefthanded (f-) = " << LeftHandedFraction << std::endl;    

  // for(unsigned int jj=0;jj<200;jj++){
  //   XVariable=-1+0.01*jj;
    
  //   SMDiffDistribution[0] = (LongitudinalFraction*6*(1-pow(XVariable,2)))/8;
  //   SMDiffDistribution[2] = (pow((1-XVariable),2)*3*LeftHandedFraction)/8;
  //   SMDiffDistribution[1] = (RightHandedFraction*3*pow((1+XVariable),2))/8;
    
  //   histo1D["HelicityDistribution"]->SetBinContent(jj+1,(SMDiffDistribution[0]+SMDiffDistribution[1]+SMDiffDistribution[2]));
  //   histo1D["HelicityDistributionLong"]->SetBinContent(jj+1,SMDiffDistribution[0]);
  //   histo1D["HelicityDistributionLeft"]->SetBinContent(jj+1,SMDiffDistribution[2]);
  //   histo1D["HelicityDistributionRight"]->SetBinContent(jj+1,SMDiffDistribution[1]);
  // }
  
  /////////////////////////////
  /// ResolutionFit Stuff
  /////////////////////////////
  bool applyKinFit = false;

  // ResolutionFit *resFitLightJets_ = new ResolutionFit("LightJet");
  // resFitLightJets_->LoadResolutions("lightJetReso.root");
  // ResolutionFit *resFitBJets_ = new ResolutionFit("BJet");
  // resFitBJets_->LoadResolutions("bJetReso.root");  

  float TopMassKinFit;
  float WMassKinFit = 80.4;

  ////////////////////////////////
  // MVA Stuff
  ////////////////////////////////
  bool applyMVA = true;           //Remember: applyKinFit is defined above!!

  //Which event reconstruction applied:
  //
  if(applyKinFit==true){
    std::cout << "    Only Kinematic Fit applied as event selection " << std::endl;
  }
  else if(applyMVA==true){
    std::cout << "    MVA (btag) combined with Kinematic Fit as event selection " << std::endl;
  }
  else{
    std::cout << "    Minimal Chi squared method used as event selection " << std::endl;
  }
   
  ////////////////////////////////////
  /// Selection table
  ////////////////////////////////////
  
  vector<string> CutsSelecTable;//Writes information in a pdf file, is filled during the selection with the different layers
  CutsSelecTable.push_back(string("initial"));
  CutsSelecTable.push_back(string("preselected"));
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
  if (verbose > 0)
    cout << " - CutsSelectionTable instantiated ..." << endl;
  SelectionTable selecTable(CutsSelecTable, datasets);
  selecTable.SetLuminosity(Luminosity);
  if (verbose > 0)
    cout << " - SelectionTable instantiated ..." << endl; 

  std::cout << " starting to initialize Whelicities class !! " << std::endl;
  ////////////////////////////////////
  //Reconstruction class WHelicities
  ////////////////////////////////////
  WReconstruction wReconstruction = WReconstruction();
  
  std::cout << " WHelicities class initialized .. " << std::endl;
  
  ////////////////////////////////////
  //	Loop on datasets
  ////////////////////////////////////
  
  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
  for (unsigned int d = 0; d < datasets.size (); d++){
    if (verbose > 1)
      cout << "   Dataset " << d << ": " << datasets[d]->Name () << "/ title : " << datasets[d]->Title () << endl;
    if (verbose > 1)
      std::cout<<"      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
    
    //open files and load
    cout<<"LoadEvent"<<endl;
    treeLoader.LoadDataset (datasets[d], anaEnv);
    cout<<"LoadEvent"<<endl;
    
    string dataSetName = datasets[d]->Name();
    
    selecTable.Fill(d,0, datasets[d]->Xsection() * datasets[d]->EquivalentLumi() );
    
    string previousFilename = "";
    int iFile = -1;

    //Initializing top mass constraint for kinematic fit:
    //Is different for Data(measured by Tevatron) and MC(generated for this top mass)
    //
    if(dataSetName.find("Data") == 0){
      TopMassKinFit = 173.1;
      SimulationSampleBoolean = false;
    }
    else{
      TopMassKinFit = 172.5;
      SimulationSampleBoolean = true;
    }
    
    /////////////////////////////////////
    /// Initialize JEC factors
    /////////////////////////////////////
   	
    vector<JetCorrectorParameters> vCorrParam;    
    
    if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") // Data!
    {
      JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("JECFiles/Jec11V2_db_AK5PFchs_L2L3Residual.txt");
      vCorrParam.push_back(*ResJetCorPar);
    }
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/Jec11V2_db_AK5PFchs_Uncertainty.txt");       
    
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, false);
   
    ////////////////////////////////////
    //	Loop on events
    ////////////////////////////////////
    
    nEvents[d] = 0;
    int itrigger = -1, previousRun = -1;
    if (verbose > 1)
      cout << "	Loop over events " << endl;
    
    //for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++){     //In this loop plots before selection can be defined
    for (unsigned int ievt = 0; ievt < 200; ievt++){     
      
      nEvents[d]++;
      if(ievt%2000 == 0)
	std::cout<<"Processing the "<<ievt<<"th event" <<flush<<"\r";
	    
      //load event
      event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);
      vector<TRootGenJet*> genjets;
      if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
      {
        genjets = treeLoader.LoadGenJet(ievt);
        sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
      
      // scale factor for the event
      float scaleFactor = 1.;
      
      // Load the GenEvent and calculate the branching ratio correction
      if(dataSetName.find("TTbarJets") == 0)
      {
        TRootGenEvent* genEvt = treeLoader.LoadGenEvent(ievt);
        if( genEvt->isSemiLeptonic() )
	  scaleFactor *= (0.108*9.)*(0.676*1.5);
	else if( genEvt->isFullHadronic() )
	  scaleFactor *= (0.676*1.5)*(0.676*1.5);
	else if( genEvt->isFullLeptonic() )
	  scaleFactor *= (0.108*9.)*(0.108*9.);
      }
      
      // Clone the init_jets vector, otherwise the corrections will be removed
      for(unsigned int i=0; i<init_jets_corrected.size(); i++)
        if(init_jets_corrected[i]) delete init_jets_corrected[i];
      init_jets_corrected.clear();

      for(unsigned int i=0; i<init_jets.size(); i++)
        init_jets_corrected.push_back( (TRootJet*) init_jets[i]->Clone() );      
      
      // check which file in the dataset it is to have the HLTInfo right
      string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
      if(previousFilename != currentFilename)
      {
      	previousFilename = currentFilename;
	iFile++;
	cout<<"File changed!!! => iFile = "<<iFile<<endl;
      }
      
      int currentRun = event->runId();
      if(previousRun != currentRun){
	previousRun = currentRun;
	if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA"){
	  if( event->runId() <= 161016 )
	    itrigger = treeLoader.iTrigger (string ("HLT_Mu17_TriCentralJet30_v1"), currentRun, iFile);	  	  
	  else if( event->runId() >= 162762 && event->runId() <= 163261 )
	    itrigger = treeLoader.iTrigger (string ("HLT_Mu17_TriCentralJet30_v2"), currentRun, iFile);	  
	  else if( event->runId() >= 163270 && event->runId() <= 163869 )
	    itrigger = treeLoader.iTrigger (string ("HLT_Mu17_TriCentralJet30_v4"), currentRun, iFile);	  
	  else if( event->runId() >= 165088 && event->runId() <= 165633 )
	    itrigger = treeLoader.iTrigger (string ("HLT_Mu17_TriCentralJet30_v5"), currentRun, iFile);	  
	  else if( event->runId() == 166346 )
	    itrigger = treeLoader.iTrigger (string ("HLT_Mu17_QuadCentralJet30_v2"), currentRun, iFile);	 
	  else if( event->runId() >= 165970 && event->runId() <= 167043 )
	    itrigger = treeLoader.iTrigger (string ("HLT_Mu17_QuadCentralJet30_v1"), currentRun, iFile);	  
	  else if( event->runId() >= 167078 )
	    itrigger = treeLoader.iTrigger (string ("HLT_Mu17_QuadCentralJet30_v3"), currentRun, iFile);
	  else
	    cout << "Unknown run for HLTpath selection: " << event->runId() << endl;     
	}
	else{
	  itrigger = treeLoader.iTrigger (string ("HLT_Mu15_v2"), currentRun);
	  if (itrigger == 9999)
	    itrigger = treeLoader.iTrigger (string ("HLT_Mu15_v1"), currentRun); // Spring11: HLT_Mu15_v1
	}
      }
      
      if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA"){
        // Apply the scraping veto
        bool isBeamBG = true;
        if(event->nTracks() > 10){
          if( ( (float) event->nHighPurityTracks() ) / ( (float) event->nTracks() ) > 0.25 )
            isBeamBG = false;
        }
        if(isBeamBG) continue;
        
        // Apply the JSON
        // bool passedJSON = treeLoader.EventPassedJSON(datasets[d], event->runId(), event->lumiBlockId());
        // if(!passedJSON) continue;
      }
          
      // Apply Jet Corrections on-the-fly
      if( dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ){
        jetTools->correctJets(init_jets_corrected, vertex);
      }

      if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) ){
        // Correct for the difference in muon efficiency (HLT and Id) between Data and MC
        //scaleFactor *= 0.965;
        
        // Correct jets for JES uncertainy systematics
	//jetTools->correctJetJESUnc(init_jets_corrected, "minus");  //also with "plus"
        
        // Match RecoJets with GenJets
        vector< pair<size_t, size_t> > indexVector; //first index = RecoJet, second index = GenJet
        vector<bool> mLock(genjets.size(),false);   // when locked, genJet is already matched to a recoJet
        for(size_t i=0; i<init_jets_corrected.size(); i++){
          pair<size_t, size_t> tmpIndex;
          float minDR = 9999.;
          for(size_t j=0; j<genjets.size(); j++){
            if( ! mLock[j] ){
              if( init_jets_corrected[i]->DeltaR(*genjets[j]) < 0.4 && init_jets_corrected[i]->DeltaR(*genjets[j]) < minDR ){
                minDR = init_jets_corrected[i]->DeltaR(*genjets[j]);
                tmpIndex = pair<size_t, size_t>(i,j);
              }
            }
          }
          if(minDR < 999.){
            mLock[tmpIndex.second] = true;
            indexVector.push_back(tmpIndex);
          }
        }
	
        // Apply correction for jet energy resolution on-the-fly, only for recoJets matched with a genJet
        for(size_t i=0; i<indexVector.size(); i++){
          if( genjets[indexVector[i].second]->Pt() < 15 ) continue;
          float corrFactor = 0.1; // factor is either 0.1 for bias correction, 0.0 for JER_minus and 0.2 for JER_plus
          float deltapt = ( init_jets_corrected[indexVector[i].first]->Pt() - genjets[indexVector[i].second]->Pt() ) * corrFactor;
          float ptscale = max(0.0, ( init_jets_corrected[indexVector[i].first]->Pt() + deltapt) / init_jets_corrected[indexVector[i].first]->Pt() );
          if(ptscale > 0.0)
            init_jets_corrected[indexVector[i].first]->SetPxPyPzE(init_jets_corrected[indexVector[i].first]->Px()*ptscale, init_jets_corrected[indexVector[i].first]->Py()*ptscale,init_jets_corrected[indexVector[i].first]->Pz()*ptscale, init_jets_corrected[indexVector[i].first]->E()*ptscale);
        }
        
        //Scale jets with a certain factor
        //jetTools->scaleJets(init_jets_corrected, 1.); 
	//if(dataSetName.find("Data")==0){
	//jetTools->scaleJets(init_jets_corrected, 1.05); //1.05 for +5%
	//}
      }
      
      /////////////////////////////
      //   Selection
      /////////////////////////////
      
      //Declare selection instance    
      Selection selection(init_jets_corrected, init_muons, init_electrons, mets);
      
      bool trigged = treeLoader.EventTrigged (itrigger);
      bool isGoodPV = selection.isPVSelected(vertex, anaEnv.PVertexNdofCut, anaEnv.PVertexZCut, anaEnv.PVertexRhoCut);

      vector<TRootJet*> selectedJets;
      vector<TRootMuon*> selectedMuons;
	   
      if (init_jets_corrected.size() > 0) {
	if (init_jets_corrected[0]->jetType() == 1 || doPF2PAT) { // calojets
	  //cout << "Selecting for caloJets" << endl;
	  selectedJets = selection.GetSelectedJets(true);
	  selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
	}	
	else {
	  //cout << "Selecting for PF/JPT jets" << endl;
	  vector<TRootMuon*> overlapMuons = selection.GetSelectedMuons(vertex[0]);
	  //selection.setJetCuts(30.,2.4,0.01,1.,0.98,0.3,0.1); // refSelV4 values
	  //selection.setMuonCuts(20,2.1,0.12,10,0.02,0.3,1,1,1);  //WHelicity uses different PFMuon isolation
	  selectedJets = selection.GetSelectedJets(overlapMuons,true);
	  selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
	}
      }
      vector<TRootMuon*> vetoMuons = selection.GetSelectedLooseMuons();
      vector<TRootElectron*> vetoElectrons = selection.GetSelectedLooseElectrons(false);
      
      vector<TRootMCParticle*> mcParticles;
      if(dataSetName.find("TTbarJets_SemiMu") == 0){
        mcParticles = treeLoader.LoadMCPart(ievt);
        sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      }    
     
      ///////////////////////////////////////////////////////
      ///     Start of program: Defining variables        ///
      ///////////////////////////////////////////////////////
      int CorrectQuark1=999;  //Needed for Monte Carlo
      int CorrectQuark2=999;
      int CorrectBHadronic=999;
      int CorrectBLeptonic=999;
      	    
      float MassW=83.6103;
      float MassTop = 172.956;
      float SigmaW=11.1534;  //Obtained from gaussian fit on Top and W distribution with simulated information
      float SigmaTop=18.232;
      
      float ChiSquared[12];  //Needed for chi squared caclulation      
      int UsedCombination;
      int QuarkOneIndex[12];
      int QuarkTwoIndex[12];
      int BHadronicIndex[12];
      int BLeptonicIndex;      
      float ChiSquaredValue;
      double btag[12]; //Needed for MVA calculation
      for(int i =0; i<12;i++){
	btag[i]=0;
      }
      double BTagValueMVA;
      int BLeptonicIndexMVA[12];
      int MVACombination[2];

      float NeutrinoPx,NeutrinoPy;  //(No initialization needed since Px and Py can always be calculated)
      float NeutrinoPz=999;   //with this value it can be distinguished in plot! Has to be calculated with DCoefficient
      float NeutrinoE=999;
      
      TLorentzVector WLeptonic;
      TLorentzVector TopLeptonic;
      TLorentzVector Neutrino;
      TLorentzVector WLeptonicTZMF;
      TLorentzVector MuonWZMF;
      
      float CosTheta=0;
      float CosThetaData=0;
      float standardCosTheta=0;
      float TheoreticalDistributionValue = 0;  //Function value for SM helicities

      TRootMCParticle standardNeutrino, standardTop,standardMuon,standardWLeptonic;      

      vector<int> jetCombi;
      if(dataSetName.find("TTbarJets_SemiMu") == 0){      	
	
      	pair<unsigned int, unsigned int> leptonicBJet_, hadronicBJet_, hadronicWJet1_, hadronicWJet2_; //First index is the JET number, second one is the parton
      	leptonicBJet_ = hadronicBJet_ = hadronicWJet1_ = hadronicWJet2_ = pair<unsigned int, unsigned int>(9999,9999);
      	vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV;
      	vector<TRootMCParticle> mcParticlesMatching;
      	bool muPlusFromTop = false, muMinusFromTop = false;
      	for(unsigned int i=0; i<mcParticles.size(); i++){
      	  if( mcParticles[i]->status() != 3) continue;
	  
      	  if( mcParticles[i]->type() == 13 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -6 ){
      	    if(muMinusFromTop) cerr<<"muMinusFromTop was already true"<<endl;
      	    muMinusFromTop = true;
      	  }
      	  if( mcParticles[i]->type() == -13 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == 6 ){
      	    if(muPlusFromTop) cerr<<"muPlusFromTop was already true"<<endl;
      	    muPlusFromTop = true;
      	  }
	  
      	  if( abs(mcParticles[i]->type()) < 6 || abs(mcParticles[i]->type()) == 21 ){
      	    mcParticlesTLV.push_back(*mcParticles[i]);
      	    mcParticlesMatching.push_back(*mcParticles[i]);
      	  }
      	}
      	if(muPlusFromTop && muMinusFromTop)
      	  cerr<<"muPlusFromTop and muMinusFromTop are both true ?!\nCheck if you are using the right sample..."<<endl;
	
      	// take all the selectedJets_ to study the radiation stuff, selectedJets are already ordened in decreasing Pt()
      	for(unsigned int i=0; i<selectedJets.size(); i++)
      	  selectedJetsTLV.push_back(*selectedJets[i]);
	
      	JetPartonMatching matching = JetPartonMatching(mcParticlesTLV, selectedJetsTLV, 2, true, true, 0.3);
	
      	if(matching.getNumberOfAvailableCombinations() != 1)
      	  cerr << "matching.getNumberOfAvailableCombinations() = "<<matching.getNumberOfAvailableCombinations()<<"  This should be equal to 1 !!!"<<endl;
	
      	vector< pair<unsigned int, unsigned int> > JetPartonPair, ISRJetPartonPair; // First one is jet number, second one is mcParticle number
	
      	for(unsigned int i=0; i<mcParticlesTLV.size(); i++){
      	  int matchedJetNumber = matching.getMatchForParton(i, 0);
      	  if(matchedJetNumber != -1)
      	    JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
      	}
	
      	for(unsigned int i=0; i<JetPartonPair.size(); i++){
      	  unsigned int j = JetPartonPair[i].second;
	  
      	  if( fabs(mcParticlesMatching[j].type()) < 6 ){
      	    if( ( muPlusFromTop && mcParticlesMatching[j].motherType() == -24 && mcParticlesMatching[j].grannyType() == -6 ) || ( muMinusFromTop && mcParticlesMatching[j].motherType() == 24 && mcParticlesMatching[j].grannyType() == 6 ) ){
      	      if(hadronicWJet1_.first == 9999) 
      		hadronicWJet1_ = JetPartonPair[i];
      	      else if(hadronicWJet2_.first == 9999) 
      		hadronicWJet2_ = JetPartonPair[i];
      	      else cerr<<"Found a third jet coming from a W boson which comes from a top quark..."<<endl;
      	    }
      	  }
      	  if( fabs(mcParticlesMatching[j].type()) == 5 ){
      	    if( ( muPlusFromTop && mcParticlesMatching[j].motherType() == -6 ) || ( muMinusFromTop && mcParticlesMatching[j].motherType() == 6 ) )
      	      hadronicBJet_ = JetPartonPair[i];
      	    else if( ( muPlusFromTop && mcParticlesMatching[j].motherType() == 6 ) || ( muMinusFromTop && mcParticlesMatching[j].motherType() == -6 ) )
      	      leptonicBJet_ = JetPartonPair[i];
      	  }
      	}
	
      	jetCombi.push_back(hadronicWJet1_.first);
      	jetCombi.push_back(hadronicWJet2_.first);
      	jetCombi.push_back(hadronicBJet_.first);
      	jetCombi.push_back(leptonicBJet_.first);
	
      	CorrectQuark1=jetCombi[0];
      	CorrectQuark2=jetCombi[1];
      	CorrectBHadronic = jetCombi[2];
      	CorrectBLeptonic = jetCombi[3];
		
      	//------------------------------------------//
      	//    Identifying Monte Carlo particles     //
      	//------------------------------------------//
      	int EventParticleNumber[5]; //0:top; 1:b; 2: u,c,d,s; 3:W; 4:mu + neutrino
      	int MuonNeutrinoNumber[4]; //0:mu-; 1:mu+; 2:nu; 3: anti-nu 
      	int ChargeMuonW[4]; //0:mu-; 1:mu+; 2:W-; 3:W+
      	for(int ll = 0;ll<5;ll++){EventParticleNumber[ll]=0;}
      	for(int kk = 0;kk<4;kk++){MuonNeutrinoNumber[kk]=0;}
      	for(int mm = 0;mm<4;mm++){ChargeMuonW[mm]=0;}
      	int MuonNumber=0;
      	int EventChargeTop=0;  //1 for top, -1 for anti-top
      	int EventChargeAntiTop=0;
      	int EventChargeWPos=0; //1 for W+, -1 for W-
      	int EventChargeWNeg=0;
      	int EventChargeMuPos=0; //1 for mu+, -1 for mu-
      	int EventChargeMuNeg =0;
      	TRootMCParticle Top,AntiTop,WPos,WNeg,Muon,BQuark,BBarQuark,Neutrino,Bottom;
      	TLorentzVector standardMuonWZMF, standardWLeptonicTZMF;
      	TRootMCParticle hadrBottom,hadrTop,UpQuark,DownQuark,hadrW;	
      	TLorentzVector hadrWTZMF,UpQuarkWZMF;

      	for(unsigned int i=0; i<mcParticles.size(); i++){
	  if( mcParticles[i]->status() != 3) continue;

      	  if(fabs(mcParticles[i]->type()) == 6){//Identifying top quarks
      	    EventParticleNumber[0]++;	    
      	    if(mcParticles[i]->type()==6){
      	      Top=*mcParticles[i];
      	      EventChargeTop=1;
      	    }
      	    else if(mcParticles[i]->type()==-6){
      	      AntiTop=*mcParticles[i];
      	      EventChargeAntiTop=-1;
      	    }
      	  }

      	  if(fabs(mcParticles[i]->type()) == 5 && fabs(mcParticles[i]->motherType()) == 6){//Identifying bottom quarks
      	    EventParticleNumber[1]++;
      	    if(mcParticles[i]->type()==5){BQuark=*mcParticles[i];}
      	    else BBarQuark=*mcParticles[i];
      	  } 

      	  if(fabs(mcParticles[i]->type()) <= 4 && fabs(mcParticles[i]->motherType()) == 24 && fabs(mcParticles[i]->grannyType()) == 6){
      	    EventParticleNumber[2]++;
      	    if(fabs(mcParticles[i]->type())==1 ||fabs(mcParticles[i]->type())==3 ){
      	      UpQuark =*mcParticles[i];
      	    }
      	    else if(fabs(mcParticles[i]->type())==2 ||fabs(mcParticles[i]->type())==4){
      	      DownQuark=*mcParticles[i];
      	    }
      	  }

      	  if(fabs(mcParticles[i]->type()) == 24 && fabs(mcParticles[i]->motherType()) == 6){ //Identifying W bosons
      	    EventParticleNumber[3]++;
      	    if(mcParticles[i]->type()==24){
      	      WPos=*mcParticles[i];
      	      EventChargeWPos=1;
      	    }
      	    else if(mcParticles[i]->type()==-24){
      	      WNeg=*mcParticles[i];
      	      EventChargeWNeg=-1;
      	    }
      	  }

      	  if(fabs(mcParticles[i]->type()) == 14 && fabs(mcParticles[i]->motherType()) == 24 && fabs(mcParticles[i]->grannyType()) == 6){
      	    EventParticleNumber[4]++;  //Identifying neutrino's
      	    standardNeutrino=*mcParticles[i];
      	  }

      	  if(fabs(mcParticles[i]->type()) == 13 && fabs(mcParticles[i]->motherType()) == 24 && fabs(mcParticles[i]->grannyType()) == 6 && mcParticles[i]->status()==3){ //status: 1:stable; 2:shower; 3:hard scattering(coming from the studied hard proces)
      	    MuonNumber++; //Identifying muons
      	    EventParticleNumber[4]++;
      	    standardMuon=*mcParticles[i];
      	    if(mcParticles[i]->type()==13){EventChargeMuNeg=-1;}
      	    else if(mcParticles[i]->type()==-13){EventChargeMuPos=1;}
      	  }	  	  
      	}//  if 0 < i < mcParticles.size()
		
      	//////////////////////////////////////////////////
      	//   Selecting correct event (b b q q mu nu )   //
      	//////////////////////////////////////////////////
      	if(EventParticleNumber[0]==2 && EventParticleNumber[1]==2 && EventParticleNumber[2]==2 && EventParticleNumber[3]==2 && EventParticleNumber[4]==2){
	  
      	  //-----   Differentiating between proces from top and anti-top (choose leptonic):   -----
      	  if(EventChargeTop==1 && EventChargeWPos==1 && EventChargeMuPos==1){  // Proces: t -> b W+ -> mu+ nu
      	    standardWLeptonicTZMF=WPos;
      	    standardWLeptonic=WPos;
      	    standardTop=Top;
      	    Bottom = BQuark;
      	  }
      	  else if(EventChargeAntiTop==-1 && EventChargeWNeg==-1 && EventChargeMuNeg==-1){ //proces: anti-t -> anti-b W- -> mu- anti-nu	
      	    standardWLeptonicTZMF=WNeg;
      	    standardWLeptonic=WNeg;
      	    standardTop=AntiTop;
      	    Bottom = BBarQuark;
      	  }

      	  //-----   Applying boost on muon and W:   -----
      	  standardMuonWZMF=standardMuon;
      	  standardMuonWZMF.Boost(-standardWLeptonicTZMF.BoostVector());
      	  standardWLeptonicTZMF.Boost(-standardTop.BoostVector());

      	  //-----   Calculating cos theta:   -----
      	  standardCosTheta = ((standardWLeptonicTZMF.Vect()).Dot(standardMuonWZMF.Vect()))/(((standardWLeptonicTZMF.Vect()).Mag())*((standardMuonWZMF.Vect()).Mag()));
      	  TheoreticalDistributionValue = (LongitudinalFraction*6*(1-standardCosTheta*standardCosTheta) + (1-standardCosTheta)*(1-standardCosTheta)*3*LeftHandedFraction + RightHandedFraction*3*(1+standardCosTheta)*(1+standardCosTheta))/8;
  	  
      	  histo1D["StandardCosTheta"]->Fill(standardCosTheta);  // Histogram without fit
      	  histo1D["StandardCosThetaFit"]->Fill(standardCosTheta);  // Histogram with fit   

      	}
      }//if dataset Semi mu ttbar  
      
      /////////////////////////
      //   Event selection   //
      /////////////////////////
      bool eventSelected = false;      
      selecTable.Fill(d,1,1.);  
      if(trigged){   //selection steps start: https://twiki.cern.ch/twiki/bin/view/CMS/TopLeptonPlusJetsRefSel_mu#Selection_Version_SelV3_for_HCP2
      	selecTable.Fill(d,2,1.);	
      	if(isGoodPV){
      	  selecTable.Fill(d,3,1.);
      	  if(selectedMuons.size()==1){
      	    selecTable.Fill(d,4,1.);	
      	    if(vetoMuons.size()==1){
      	      selecTable.Fill(d,5,1.);	      
      	      if(vetoElectrons.size()==0){
      		selecTable.Fill(d,6,1.);
		if(init_jets.size()>=4){
		}
      		if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-3){
      		  selecTable.Fill(d,7,1.);
      		  if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-2){
      		    selecTable.Fill(d,8,1.);
      		    if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-1){
      		      selecTable.Fill(d,9,1.);
      		      if(selectedJets.size()>=(unsigned int)anaEnv.NofJets){   //we ask for 4 jets, but this way (with -3,-2,-1)it is possible to look at plots with other amount of jets
      			selecTable.Fill(d,10,1.);
      			eventSelected = true;
      			if(verbose>2) cout << event->runId() << ":" << event->eventId() << ":" << event->lumiBlockId()  << ":" << setprecision(8) << selectedMuons[0]->Pt() << endl;			
      		      }
      		    }
      		  }
      		}
      	      }
      	    }
      	  }
      	}
      }          
      
      if(eventSelected==true){  

	float CorrectRecMassW=0;
	float CorrectRecMassTop=0;
	vector<int> jetCombi;
	if(dataSetName.find("TTbarJets_SemiMu") == 0){
	  
	  pair<unsigned int, unsigned int> leptonicBJet_, hadronicBJet_, hadronicWJet1_, hadronicWJet2_; //First index is the JET number, second one is the parton
	  leptonicBJet_ = hadronicBJet_ = hadronicWJet1_ = hadronicWJet2_ = pair<unsigned int, unsigned int>(9999,9999);
	  vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV;
	  vector<TRootMCParticle> mcParticlesMatching;
	  bool muPlusFromTop = false, muMinusFromTop = false;
	  for(unsigned int i=0; i<mcParticles.size(); i++){
	    if( mcParticles[i]->status() != 3) continue;
	    
	    if( mcParticles[i]->type() == 13 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -6 ){
	      if(muMinusFromTop) cerr<<"muMinusFromTop was already true"<<endl;
	      muMinusFromTop = true;
	    }
	    if( mcParticles[i]->type() == -13 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == 6 ){
	      if(muPlusFromTop) cerr<<"muPlusFromTop was already true"<<endl;
	      muPlusFromTop = true;
	    }
	    
	    if( abs(mcParticles[i]->type()) < 6 || abs(mcParticles[i]->type()) == 21 ){
	      mcParticlesTLV.push_back(*mcParticles[i]);
	      mcParticlesMatching.push_back(*mcParticles[i]);
	    }
	  }
	  if(muPlusFromTop && muMinusFromTop)
	    cerr<<"muPlusFromTop and muMinusFromTop are both true ?!\nCheck if you are using the right sample..."<<endl;
	  
	  // take all the selectedJets_ to study the radiation stuff, selectedJets are already ordened in decreasing Pt()
	  for(unsigned int i=0; i<selectedJets.size(); i++)
	    selectedJetsTLV.push_back(*selectedJets[i]);
	  
	  JetPartonMatching matching = JetPartonMatching(mcParticlesTLV, selectedJetsTLV, 2, true, true, 0.3);
	  
	  if(matching.getNumberOfAvailableCombinations() != 1)
	    cerr << "matching.getNumberOfAvailableCombinations() = "<<matching.getNumberOfAvailableCombinations()<<"  This should be equal to 1 !!!"<<endl;
	  
	  vector< pair<unsigned int, unsigned int> > JetPartonPair, ISRJetPartonPair; //First one is jet number, second one is mcParticle number
	  
	  for(unsigned int i=0; i<mcParticlesTLV.size(); i++){
	    int matchedJetNumber = matching.getMatchForParton(i, 0);
	    if(matchedJetNumber != -1)
	      JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
	  }
	  
	  for(unsigned int i=0; i<JetPartonPair.size(); i++){
	    unsigned int j = JetPartonPair[i].second;
	    
	    if( fabs(mcParticlesMatching[j].type()) < 6 ){
	      if( ( muPlusFromTop && mcParticlesMatching[j].motherType() == -24 && mcParticlesMatching[j].grannyType() == -6 ) || ( muMinusFromTop && mcParticlesMatching[j].motherType() == 24 && mcParticlesMatching[j].grannyType() == 6 ) ){
		if(hadronicWJet1_.first == 9999) 
		  hadronicWJet1_ = JetPartonPair[i];
		else if(hadronicWJet2_.first == 9999) 
		  hadronicWJet2_ = JetPartonPair[i];
		else cerr<<"Found a third jet coming from a W boson which comes from a top quark..."<<endl;
	      }
	    }
	    if( fabs(mcParticlesMatching[j].type()) == 5 ){
	      if( ( muPlusFromTop && mcParticlesMatching[j].motherType() == -6 ) || ( muMinusFromTop && mcParticlesMatching[j].motherType() == 6 ) )
		hadronicBJet_ = JetPartonPair[i];
	      else if( ( muPlusFromTop && mcParticlesMatching[j].motherType() == 6 ) || ( muMinusFromTop && mcParticlesMatching[j].motherType() == -6 ) )
		leptonicBJet_ = JetPartonPair[i];
	    }
	  }
	  
	  jetCombi.push_back(hadronicWJet1_.first);
	  jetCombi.push_back(hadronicWJet2_.first);
	  jetCombi.push_back(hadronicBJet_.first);
	  jetCombi.push_back(leptonicBJet_.first);
	  	 
	  CorrectQuark1=jetCombi[0];
	  CorrectQuark2=jetCombi[1];
	  CorrectBHadronic = jetCombi[2];
	  CorrectBLeptonic = jetCombi[3];			  
	  
	  //Working on generator level (i.e. jets level):  
	  if(jetCombi[0]!=9999 && jetCombi[1]!=9999 && jetCombi[2]!=9999 && jetCombi[3]!=9999){    
	    CorrectRecMassW=(*selectedJets[jetCombi[0]]+*selectedJets[jetCombi[1]]).M();
	    CorrectRecMassTop=(*selectedJets[jetCombi[0]]+*selectedJets[jetCombi[1]]+*selectedJets[jetCombi[2]]).M();
	    
	    histo1D["WMass"]->Fill(CorrectRecMassW);
	    histo1D["TopMass"]->Fill(CorrectRecMassTop);
	  }	      	      	      	       	      
	}//if dataset Semi mu ttbar

	//Reconstruction class WHelicities
	//
	//std::cout << " Start of calling WHelicities class " << std::endl;
	//

	std::cout << " ************** Studied event : " << ievt << std::endl;
	
	std::vector<float> bTagValues;
	for(int ii = 0; ii<selectedJets.size();ii++){bTagValues.push_back(selectedJets[ii]->btag_trackCountingHighEffBJetTags());}

	std::vector<int> ChosenReconstructionJets;

	//vector of TLorentzVector should be given as input:
	std::vector<float> lorentzJetsPx;
	std::vector<float> lorentzJetsPy;
	std::vector<float> lorentzJetsPz;
	std::vector<float> lorentzJetsE;
	for(int ii=0; ii<selectedJets.size();ii++){
	  lorentzJetsPx.push_back((*selectedJets[ii]).Px());
	  lorentzJetsPy.push_back((*selectedJets[ii]).Py());
	  lorentzJetsPz.push_back((*selectedJets[ii]).Pz());
	  lorentzJetsE.push_back((*selectedJets[ii]).E());
	}

	std::vector<float> lorentzMuonsPx;
	std::vector<float> lorentzMuonsPy;
	std::vector<float> lorentzMuonsPz;
	std::vector<float> lorentzMuonsE;
	for(int ii=0;ii<selectedMuons.size();ii++){
	  lorentzMuonsPx.push_back((*selectedMuons[ii]).Px());
	  lorentzMuonsPy.push_back((*selectedMuons[ii]).Py());
	  lorentzMuonsPz.push_back((*selectedMuons[ii]).Pz());
	  lorentzMuonsE.push_back((*selectedMuons[ii]).E());
	}

	//wReconstruction.Initializing(SimulationSampleBoolean,lorentzJetsPx,lorentzJetsPy,lorentzJetsPz,lorentzJetsE,lorentzMuonsPx,lorentzMuonsPy,lorentzMuonsPz,lorentzMuonsE,bTagValues, applyKinFit, applyMVA);
	wReconstruction.Initializing(SimulationSampleBoolean,lorentzJetsPx,lorentzJetsPy,lorentzJetsPz,lorentzJetsE,bTagValues);

	std::cout << " --------------------------------------------- " << std::endl;

	ChosenReconstructionJets = wReconstruction.BTagAnalysis();
	std::cout << " value for boolean event : " << ChosenReconstructionJets[4] << std::endl;
	bool EventSelectedWithBTag;
	if(ChosenReconstructionJets[4]==1)EventSelectedWithBTag=true;
	else continue;

	std::cout << " hadronic b number : " << ChosenReconstructionJets[0] << std::endl;
	std::cout << " leptonic b number : " << ChosenReconstructionJets[1] << std::endl;
	std::cout << " quark1 number : " << ChosenReconstructionJets[2] << std::endl;
	std::cout << " quark2 number : " << ChosenReconstructionJets[3] << std::endl;
	std::cout << " value for boolean event : " << ChosenReconstructionJets[4] << std::endl;
	std::cout << " --------------------------------------------- " << std::endl;

	int bHadronicNumber = ChosenReconstructionJets[0];
	int bLeptonicNumber = ChosenReconstructionJets[1];
	int Quark1Number = ChosenReconstructionJets[2];
	int Quark2Number = ChosenReconstructionJets[3];
	
	//std::cout << " Clearing used vectors for next event : " << std::endl;
	lorentzMuonsPx.clear();
	lorentzMuonsPy.clear();
	lorentzMuonsPz.clear();
	lorentzMuonsE.clear();
	lorentzJetsPx.clear();
	lorentzJetsPy.clear();
	lorentzJetsPz.clear();
	lorentzJetsE.clear();
	ChosenReconstructionJets.clear();
	bTagValues.clear();
	
	/*//////////////////////////////////////////////////
	//     Calculating correct jet distribution     //
	//////////////////////////////////////////////////
	//
	std::cout << " Studied event : " << ievt << std::endl;
	int NumberCombinations=0;	
	for(int i=0;i<3;i++){
	  for(int j=i+1;j<4;j++){
	    for(int k=0;k<4;k++){
	      if(k!=i && k!=j){
		
		if(applyKinFit == true){
		  TLorentzVector lightJet1 = *selectedJets[i];
		  TLorentzVector lightJet2 = *selectedJets[j];
		  TLorentzVector bJet = *selectedJets[k];
		  
		  // prepare everything for the Kinematic Fit
		  TMatrixD Ml1(3,3), Ml2(3,3), Mb(3,3);
		  Ml1.Zero(); Ml2.Zero(); Mb.Zero();
		  Ml1(0,0) = pow(resFitLightJets_->EtResolution(&lightJet1), 2);
		  Ml1(1,1) = pow(resFitLightJets_->ThetaResolution(&lightJet1), 2);
		  Ml1(2,2) = pow(resFitLightJets_->PhiResolution(&lightJet1), 2);
		  Ml2(0,0) = pow(resFitLightJets_->EtResolution(&lightJet2), 2);
		  Ml2(1,1) = pow(resFitLightJets_->ThetaResolution(&lightJet2), 2);
		  Ml2(2,2) = pow(resFitLightJets_->PhiResolution(&lightJet2), 2);
		  Mb(0,0) = pow(resFitBJets_->EtResolution(&bJet), 2);
		  Mb(1,1) = pow(resFitBJets_->ThetaResolution(&bJet), 2);
		  Mb(2,2) = pow(resFitBJets_->PhiResolution(&bJet), 2);
		  
		  TKinFitter *theFitter = new TKinFitter("hadtopFit", "hadtopFit");
		  theFitter->setVerbosity(0);
		  
		  TFitParticleEtThetaPhiEMomFix *fitLight1 = new TFitParticleEtThetaPhiEMomFix("lightJet1", "lightJet1", &lightJet1, &Ml1);
		  TFitParticleEtThetaPhiEMomFix *fitLight2 = new TFitParticleEtThetaPhiEMomFix("lightJet2", "lightJet2", &lightJet2, &Ml2);
		  TFitParticleEtThetaPhiEMomFix *fitB = new TFitParticleEtThetaPhiEMomFix("bJet", "bJet", &bJet, &Mb);
		  theFitter->addMeasParticles(fitLight1,fitLight2,fitB);
		  
		  TFitConstraintM *consW = new TFitConstraintM("WBosonMass", "MassConstraint", 0, 0, WMassKinFit);
		  TFitConstraintM *consTop = new TFitConstraintM("TopQuarkMass", "MassConstraint", 0, 0, TopMassKinFit );//Different mass for MC and Data!!
		  consW->addParticles1(fitLight1,fitLight2);
		  consTop->addParticles1(fitB,fitLight1,fitLight2);
		  
		  theFitter->addConstraint(consW);
		  theFitter->addConstraint(consTop);
		  theFitter->setMaxNbIter(30);
		  theFitter->setMaxDeltaS(5e-5);
		  theFitter->setMaxF(1e-4);
		  
		  //do the fit!
		  theFitter->fit();
		  if (theFitter->getStatus() == 0) // if the fitter converged
		    ChiSquared[NumberCombinations]=theFitter->getS();
		  //else
		  //cout << "FIT NOT CONVERGED" << endl;
		  
		  delete theFitter;
		  delete fitLight1;
		  delete fitLight2;
		  delete fitB;
		  delete consW;
		  delete consTop;
		}//Kinematic fit applied
		else if(applyMVA==true){
		  for(int l = 0; l < 4 ; l++){ 
		    if(l!=i && l!=j && l!=k){
		      float btag_i = selectedJets[i]->btag_trackCountingHighEffBJetTags();
		      if(btag_i < -90) btag_i = 0;
		      float btag_j = selectedJets[j]->btag_trackCountingHighEffBJetTags();
		      if(btag_j < -90) btag_j = 0;
		      float btag_k = selectedJets[k]->btag_trackCountingHighEffBJetTags();
		      if(btag_k < -90) btag_k = 0;
		      float btag_l = selectedJets[l]->btag_trackCountingHighEffBJetTags();
		      if(btag_l < -90) btag_l = 0;		      

		      bTagSelected[NumberCombinations]=false;
		      if(btag_k > 1.7 || btag_l > 1.7){
			bTagSelected[NumberCombinations]=true;
		      }
		      
		      btag[NumberCombinations] = pow(btag_k,2) + pow(btag_l,2);
		      btag[NumberCombinations] = btag[NumberCombinations] / ( pow(btag_i,2) + pow(btag_j,2) + pow(btag_k,2) + pow(btag_l,2) );	

		      ProblemEvent = false;
		      if(btag_i ==0 && btag_j ==0 && btag_k ==0 && btag_k ==0){
			btag[NumberCombinations]=0;
			ProblemEvent=true;
		      }
		      
		      //Leptonic b quark is also obtained with this method:
		      BLeptonicIndexMVA[NumberCombinations]=l;
		    }
		  }		  
		}
		else{
		  float recMassW = (*selectedJets[i]+*selectedJets[j]).M();
		  float recMassTop=(*selectedJets[i]+*selectedJets[j]+*selectedJets[k]).M();
		  
		  ChiSquared[NumberCombinations]=pow(((recMassW-MassW)/SigmaW),2)+pow(((recMassTop-MassTop)/SigmaTop),2);
		}//No Kinematic Fit applied (minimal chi squared applied)
		QuarkOneIndex[NumberCombinations]=i;
		QuarkTwoIndex[NumberCombinations]=j;
		BHadronicIndex[NumberCombinations]=k;		
		
		NumberCombinations++;  //Always gives 12 as it should be!
	      }
	    }//end of k loop for jet combination selection
	  }//end of j loop for jet combination selection
	}//end of i loop for jet combination selection		
	
	//histo1D["MVACombinationBefore"]->Fill(MVACombination[0]);

	//
	//Select lowest chi squared value:
	//
	if(applyMVA ==  false){
	  ChiSquaredValue=ChiSquared[0];
	  for(int ii=0;ii<12;ii++){
	    if(ChiSquaredValue>ChiSquared[ii]){
	      ChiSquaredValue=ChiSquared[ii];
	      UsedCombination=ii;        
	    } 
	  }
	  //
	  //Jet not in Chisquared combination is the Leptonic B jet
	  //
	  for(int ll=0;ll<4;ll++){
	    if(ll!=QuarkOneIndex[UsedCombination] && ll!=QuarkTwoIndex[UsedCombination] && ll!=BHadronicIndex[UsedCombination]){
	      std::cout << " calculating leptonic b " << std::endl;
	      BLeptonicIndex=ll;
	    }
	  }
	}
	else if(applyMVA == true && bTagSelected==true){
	  std::cout << " ****************************************************** " << std::endl;

	  BTagValueMVA=0;
	  MVACombination[0]=12;
	  MVACombination[1]=12;
	  histo1D["MVACombinationBefore"]->Fill(MVACombination[0]);
	  for(int ii=0;ii<12;ii++){  //From the twelve combinations the combination with the highest BTag value is selected
	    if(BTagValueMVA<btag[ii] && !(fabs(BTagValueMVA-btag[ii])<0.000001)){
	      BTagValueMVA=btag[ii];
	      MVACombination[0]=ii;
	      MVACombination[1]=ii+1;	     
	      
	    }	
	    if(fabs(BTagValueMVA-btag[ii])<0.000001){
	      MVACombination[1]=ii;
	    }
	  }	 
	  if(ProblemEvent==true){//In this configuration the b-jets have the highest pT which is more reasonable
	    MVACombination[0]=10;
	    MVACombination[1]=11;
	    std::cout << " Filled MVA combination of problem event : " << ievt << std::endl;
	  }
	  std::cout << " two mva combinations : " << MVACombination[0] << " | " << MVACombination[1] << std::endl;

	  if(MVACombination[0]==CorrectBHadronic || MVACombination[0] == CorrectBLeptonic || MVACombination[1] == CorrectBHadronic || MVACombination[1] == CorrectBLeptonic){
	    BTagCorrectCombination++;
	  }
	  BTagCombinationEvents++;

	  histo1D["MVACombination"]->Fill(MVACombination[0]);
	  histo1D["MVACombinationTwo"]->Fill(MVACombination[1]);

	  if(MVACombination[0]==0){CorrectZero++;}
	    else if(MVACombination[0]==2){CorrectTwo++;}
	    else if(MVACombination[0]==4){CorrectFour++;}
	    else if(MVACombination[0]==6){CorrectSix++;}
	    else if(MVACombination[0]==8){CorrectEight++;}
	    else if(MVACombination[0]==10){CorrectTen++;}

	  for(int ii=0;ii<12;ii++){
	    if(ii!=MVACombination[0] && ii!=MVACombination[1]){
	      MSPlot["NotSelectedBTag"]->Fill(btag[ii], datasets[d], true, Luminosity*scaleFactor);
	    }
	  }
	  
	  //BLeptonicIndex = BLeptonicIndexMVA[UsedCombination];
	  //The two obtained values of the MVA have to be compared with a kinematic fit to select the correct event topology.
	  for(int jj=0;jj<2; jj++){
	    int MVAComb=MVACombination[jj];
	    
	    TLorentzVector lightJet1 = *selectedJets[QuarkOneIndex[MVAComb]];
	    TLorentzVector lightJet2 = *selectedJets[QuarkTwoIndex[MVAComb]];
	    TLorentzVector bJet = *selectedJets[BHadronicIndex[MVAComb]];
	    
	    // prepare everything for the Kinematic Fit
	    TMatrixD Ml1(3,3), Ml2(3,3), Mb(3,3);
	    Ml1.Zero(); Ml2.Zero(); Mb.Zero();
	    Ml1(0,0) = pow(resFitLightJets_->EtResolution(&lightJet1), 2);
	    Ml1(1,1) = pow(resFitLightJets_->ThetaResolution(&lightJet1), 2);
	    Ml1(2,2) = pow(resFitLightJets_->PhiResolution(&lightJet1), 2);
	    Ml2(0,0) = pow(resFitLightJets_->EtResolution(&lightJet2), 2);
	    Ml2(1,1) = pow(resFitLightJets_->ThetaResolution(&lightJet2), 2);
	    Ml2(2,2) = pow(resFitLightJets_->PhiResolution(&lightJet2), 2);
	    Mb(0,0) = pow(resFitBJets_->EtResolution(&bJet), 2);
	    Mb(1,1) = pow(resFitBJets_->ThetaResolution(&bJet), 2);
	    Mb(2,2) = pow(resFitBJets_->PhiResolution(&bJet), 2);
	    
	    TKinFitter *theFitter = new TKinFitter("hadtopFit", "hadtopFit");
	    theFitter->setVerbosity(0);
	    
	    TFitParticleEtThetaPhiEMomFix *fitLight1 = new TFitParticleEtThetaPhiEMomFix("lightJet1", "lightJet1", &lightJet1, &Ml1);
	    TFitParticleEtThetaPhiEMomFix *fitLight2 = new TFitParticleEtThetaPhiEMomFix("lightJet2", "lightJet2", &lightJet2, &Ml2);
	    TFitParticleEtThetaPhiEMomFix *fitB = new TFitParticleEtThetaPhiEMomFix("bJet", "bJet", &bJet, &Mb);
	    theFitter->addMeasParticles(fitLight1,fitLight2,fitB);
	    
	    TFitConstraintM *consW = new TFitConstraintM("WBosonMass", "MassConstraint", 0, 0, WMassKinFit);
	    TFitConstraintM *consTop = new TFitConstraintM("TopQuarkMass", "MassConstraint", 0, 0, TopMassKinFit );//Different mass for MC and Data!!
	    consW->addParticles1(fitLight1,fitLight2);
	    consTop->addParticles1(fitB,fitLight1,fitLight2);
	    
	    theFitter->addConstraint(consW);
	    theFitter->addConstraint(consTop);
	    theFitter->setMaxNbIter(30);
	    theFitter->setMaxDeltaS(5e-5);
	    theFitter->setMaxF(1e-4);
	    
	    //do the fit!
	    theFitter->fit();
	    if (theFitter->getStatus() == 0) // if the fitter converged
	      ChiSquared[jj]=theFitter->getS();	    
	    //else
	    //cout << "FIT NOT CONVERGED" << endl;
	    
	    delete theFitter;
	    delete fitLight1;
	    delete fitLight2;
	    delete fitB;
	    delete consW;
	    delete consTop;
	  }
	  ChiSquaredValue=ChiSquared[0];
	  int ChosenNumber;
	  UsedCombination=MVACombination[0];
	  for(int ii=0;ii<2;ii++){
	    if(ChiSquaredValue>ChiSquared[ii]){
	      ChiSquaredValue=ChiSquared[ii];
	      UsedCombination=MVACombination[ii];      
	      ChosenNumber=ii;
	    } 
	  }
	  BLeptonicIndex=BLeptonicIndexMVA[UsedCombination];
 
	  for(int ii=0;ii<2;ii++){
	    if(ii!=ChosenNumber){
	      MSPlot["NotSelectedChiSq"]->Fill(ChiSquared[ii], datasets[d], true, Luminosity*scaleFactor);
	      histo2D["MVACorrelationPlotWrong"]->Fill(ChiSquared[ii],BTagValueMVA);

	    }
	  }
	  //Make correlation plot for Chi2 and btag of MVA method:
	  //
	  histo2D["MVACorrelationPlot"]->Fill(ChiSquaredValue,BTagValueMVA);

	std::cout << " bhadronic number : " << BHadronicIndex[UsedCombination] << std::endl;
	std::cout << " bleptonic number : " << BLeptonicIndex << std::endl;
	std::cout << " quark1 number : " << QuarkOneIndex[UsedCombination] << std::endl;
	std::cout << " quark2 number : " << QuarkTwoIndex[UsedCombination] << std::endl;
	std::cout << " " << std::endl;
	std::cout << " ****************************************************** " << std::endl;
	std::cout << " " << std::endl;
	}*/
		
	// if(bTagSelected == true){
	//   //Only for Mc data, check if chosen combination is correct
	//   //
	//   if(dataSetName.find("TTbarJets_SemiMu") == 0){
	//     int EventWrong = 1;
	//     float RecoTopMassHadrTT = (*selectedJets[QuarkOneIndex[UsedCombination]] + *selectedJets[QuarkTwoIndex[UsedCombination]] + *selectedJets[BHadronicIndex[UsedCombination]]).M();
	    
	//     NumberSemiMuEvents++;
	//     if((QuarkOneIndex[UsedCombination]==CorrectQuark1 && QuarkTwoIndex[UsedCombination]==CorrectQuark2 && BHadronicIndex[UsedCombination]==CorrectBHadronic && BLeptonicIndex==CorrectBLeptonic) || (QuarkOneIndex[UsedCombination]==CorrectQuark2 && QuarkTwoIndex[UsedCombination]==CorrectQuark1 && BHadronicIndex[UsedCombination]==CorrectBHadronic && BLeptonicIndex==CorrectBLeptonic)){	        	    
	//       AllJetsCorrect++;
	//       histo1D["TopMassForCorrect"]->Fill(RecoTopMassHadrTT);
	//       MSPlot["ChiSquaredCorrect"]->Fill(ChiSquaredValue, datasets[d], true, Luminosity*scaleFactor);
	//       if(applyMVA==true && bTagSelected==true){
	// 	MSPlot["BTagValueMVACorrect"]->Fill(BTagValueMVA,datasets[d],true,Luminosity*scaleFactor);
	// 	histo2D["BtagForTopCorrect"]->Fill(RecoTopMassHadrTT,BTagValueMVA);
	//       }
	//     }
	//     if((QuarkOneIndex[UsedCombination]!=CorrectQuark1 && QuarkTwoIndex[UsedCombination]!=CorrectQuark2 && BHadronicIndex[UsedCombination]!=CorrectBHadronic && BLeptonicIndex!=CorrectBLeptonic) ||(QuarkOneIndex[UsedCombination]!=CorrectQuark2 && QuarkTwoIndex[UsedCombination]!=CorrectQuark1 && BHadronicIndex[UsedCombination]!=CorrectBHadronic && BLeptonicIndex!=CorrectBLeptonic)){	    
	//       OneJetWrong++;
	//       histo1D["TopMassForWrong"]->Fill(RecoTopMassHadrTT);
	//       MSPlot["ChiSquaredWrong"]->Fill(ChiSquaredValue, datasets[d], true, Luminosity*scaleFactor);
	//       if(applyMVA==true && bTagSelected==true){	      
	//       MSPlot["BTagValueMVAWrong"]->Fill(BTagValueMVA,datasets[d],true,Luminosity*scaleFactor);
	//       histo2D["BtagForTopWrong"]->Fill(RecoTopMassHadrTT,BTagValueMVA);
	//       }
	//       EventWrong=0;
	//     }
	//     if(QuarkOneIndex[UsedCombination]==CorrectQuark1 ||QuarkOneIndex[UsedCombination]==CorrectQuark2 || QuarkTwoIndex[UsedCombination]==CorrectQuark1 || QuarkTwoIndex[UsedCombination]==CorrectQuark2 || BHadronicIndex[UsedCombination]==CorrectBHadronic || BLeptonicIndex==CorrectBLeptonic ){
	//       histo1D["TopMassNotWrong"]->Fill(RecoTopMassHadrTT);
	      
	//     }	  
	//     //if(BLeptonicIndex==CorrectBLeptonic && bTagValue>1.7){
	//     //numberCorrectLeptonicBLowbTag++;
	//     //}
	//     if(BHadronicIndex[UsedCombination]==CorrectBHadronic){
	//       numberCorrectHadronicB++;
	      
	//     }
	//     if((QuarkOneIndex[UsedCombination]==CorrectQuark1 && QuarkTwoIndex[UsedCombination]==CorrectQuark2) ||(QuarkOneIndex[UsedCombination]==CorrectQuark2 && QuarkTwoIndex[UsedCombination]==CorrectQuark1)){
	//       numberCorrectQuarks++;
	      
	//     }
	//     //Compare reconstructed mass against corresponding chi squared value (and btag value in MVA case)
	//     histo2D["ChiSquaredForTop"]->Fill(RecoTopMassHadrTT,ChiSquaredValue);
	//     if(applyMVA==true && bTagSelected==true){
	//       histo2D["BtagForTop"]->Fill(RecoTopMassHadrTT,BTagValueMVA);
	//     }
	//     histo1D["TopMassReco"]->Fill(RecoTopMassHadrTT);
	//   }
	//   MSPlot["ChiSquared"]->Fill(ChiSquaredValue,datasets[d],true,Luminosity*scaleFactor);
	//   MSPlot["BTagValueMVA"]->Fill(BTagValueMVA,datasets[d],true,Luminosity*scaleFactor);
	  	  
	//   ////////////////////////////////////////////////////////////////
	//   //     Calculating MET_Pz() (equation ax² + bx + c = 0 ):     //
	//   ////////////////////////////////////////////////////////////////
	//   //
	//   // MET_Px() and MET_Py() is known since there is no Px() and Py() component before the collision.
	//   // The Pz() component before the collision is not known for proton-proton collisions.
	//   //
	//   NeutrinoPx = -(*selectedMuons[0]+*selectedJets[0]+*selectedJets[1]+*selectedJets[2]+*selectedJets[3]).Px();
	//   NeutrinoPy = -(*selectedMuons[0]+*selectedJets[0]+*selectedJets[1]+*selectedJets[2]+*selectedJets[3]).Py();	
	  
	//   //Calculating solutions for quadratic equation  
	//   //
	//   float aCoefficient = 4*pow(selectedMuons[0]->E(),2)-4*pow(selectedMuons[0]->Pz(),2);
	//   float bCoefficient = 4*(selectedMuons[0]->Pz())*(pow(selectedMuons[0]->M(),2)-pow(MassW,2)-2*(selectedMuons[0]->Px())*NeutrinoPx-2*(selectedMuons[0]->Py())*NeutrinoPy);
	//   float cCoefficient = -pow(selectedMuons[0]->M(),4)-pow(MassW,4)-4*pow(selectedMuons[0]->Px(),2)*pow(NeutrinoPx,2)-4*pow(selectedMuons[0]->Py(),2)*pow(NeutrinoPy,2)+4*(pow(selectedMuons[0]->M(),2)-pow(MassW,2))*((selectedMuons[0]->Px())*NeutrinoPx+(selectedMuons[0]->Py())*NeutrinoPy)-8*(selectedMuons[0]->Px())*NeutrinoPx*(selectedMuons[0]->Py())*NeutrinoPy+4*pow(selectedMuons[0]->E(),2)*pow(NeutrinoPx,2)+4*pow(selectedMuons[0]->E(),2)*pow(NeutrinoPy,2)+2*(pow(MassW,2))*(pow(selectedMuons[0]->M(),2));
	  
	//   float DCoefficient = pow(bCoefficient,2)-4*aCoefficient*cCoefficient; 		
	//   if(DCoefficient>=0){	  //Still need to find solution for D<0 !!!
	//     float NeutrinoPzOne = ((-bCoefficient + sqrt(DCoefficient))/(aCoefficient*2));
	//     float NeutrinoPzTwo = ((-bCoefficient - sqrt(DCoefficient))/(aCoefficient*2));	  
	    
	//     float NeutrinoEOne = sqrt(NeutrinoPx*NeutrinoPx+NeutrinoPy*NeutrinoPy+NeutrinoPzOne*NeutrinoPzOne);
	//     float NeutrinoETwo = sqrt(NeutrinoPx*NeutrinoPx+NeutrinoPy*NeutrinoPy+NeutrinoPzTwo*NeutrinoPzTwo);
	//     TLorentzVector NeutrinoOne;
	//     NeutrinoOne.SetPxPyPzE(NeutrinoPx,NeutrinoPy,NeutrinoPzOne,NeutrinoEOne);  //Has to be initialized like this !
	//     TLorentzVector NeutrinoTwo;
	//     NeutrinoTwo.SetPxPyPzE(NeutrinoPx,NeutrinoPy,NeutrinoPzTwo,NeutrinoETwo);
	    	  
	//     float TopMassDiffOne = (MassTop -((NeutrinoOne+*selectedMuons[0]+*selectedJets[BLeptonicIndex]).M()) );	  
	//     float TopMassDiffTwo = (MassTop -((NeutrinoTwo+*selectedMuons[0]+*selectedJets[BLeptonicIndex]).M()) );	   
	    
	//     //-----   Selected neutrino has the smallest derivation from the top mass   -----
	//     //
	//     if(fabs(TopMassDiffOne)<fabs(TopMassDiffTwo)){
	//       NeutrinoPz=NeutrinoPzOne;
	//       NeutrinoE = NeutrinoEOne;
	//     }
	//     else{
	//       NeutrinoPz=NeutrinoPzTwo;
	//       NeutrinoE = NeutrinoETwo;
	//     }	  	  
	//     Neutrino.SetPxPyPzE(NeutrinoPx,NeutrinoPy,NeutrinoPz,NeutrinoE);	  
	    
	//     WLeptonic = (Neutrino+*selectedMuons[0]);
	//     TopLeptonic = (Neutrino+*selectedMuons[0]+*selectedJets[BLeptonicIndex]);
	    
	//     //Reboost the particles to rest frames 
	//     //
	//     MuonWZMF = *selectedMuons[0]; // In W Zero Mass Frame (WZMF)
	//     WLeptonicTZMF = WLeptonic;  // In Top Zero Mass Frame (TZMF)	  
	    
	//     MuonWZMF.Boost(-WLeptonic.BoostVector());
	//     WLeptonicTZMF.Boost(-TopLeptonic.BoostVector());
	    
	//     //Calculating cos:	      
	//     CosTheta = ((WLeptonicTZMF.Vect()).Dot(MuonWZMF.Vect()))/(((WLeptonicTZMF.Vect()).Mag())*((MuonWZMF.Vect()).Mag()));
	//     MSPlot["CosTheta"]->Fill(CosTheta,datasets[d],true,Luminosity*scaleFactor);//No applied weight: CosTheta is the same for all helicities 	  
	    
	//     if(dataSetName.find("Data") == 0){ 	    //Make histogram with only CosTheta values for Data sample
	//       histo1D["CosThetaData"]->Fill(CosTheta); //No applied weight since this is considered as Data!!
	//       //histo1D["CosThetaData"]->Fill(CosTheta,scaleFactor*Luminosity*datasets[d]->NormFactor()); //Bckg systematics
	//     }
	//     //Systematics run:
	//     //	  	  
	//     histo1D["CosThetaData"]->Fill(CosTheta,scaleFactor*Luminosity*datasets[d]->NormFactor()); 	  
	    
	//     if(dataSetName.find("Data") != 0){ //Make histogram with CosTheta values for non-Data sample 
	//       histo1D["CosThetaMC"]->Fill(CosTheta,scaleFactor*Luminosity*datasets[d]->NormFactor());	    
	//     }	  	  
	    
	//     //////////////////////////////////////////////////
	//     //     Create loop with different helicities    //
	//     //////////////////////////////////////////////////
	//     float HelicityFraction[3];  //0:Longitudinal; 1: Righthanded; 2: Lefthanded
	//     float UsedDistributionValue = 0; 
	//     float HelicityWeight = 0;  
	    
	//     for(int Longit =0;Longit <=NumberOfHelicityBins;Longit++){
	//       for(int Right =0; Right <=(NumberOfHelicityBins-Longit) ; Right++){
	// 	int Left = NumberOfHelicityBins-Longit-Right;
		
	// 	double longit=Longit*(1.00/NumberOfHelicityBins);
	// 	double right=Right*(1.00/NumberOfHelicityBins);
	// 	double left=Left*(1.00/NumberOfHelicityBins);
	// 	std::string THistoName = "CosThetaRight"+HelicityNumber[Right]+"Long"+HelicityNumber[Longit]+"Left"+HelicityNumber[Left];     
	// 	std::string HelicityWeightHisto = "HelicityWeightRight"+HelicityNumber[Right]+"Long"+HelicityNumber[Longit]+"Left"+HelicityNumber[Left];	      	      
		
	// 	HelicityFraction[0]=longit;
	// 	HelicityFraction[1]=right;
	// 	HelicityFraction[2]=left;    
		
	// 	HelicityWeight=1;   //enkel voor ttMu via herweging, voor alle andere: HelicityWeight==1	       
	// 	if(dataSetName.find("TTbarJets_SemiMu") == 0){
	// 	  UsedDistributionValue = (HelicityFraction[0]*6*(1-standardCosTheta*standardCosTheta) + (1-standardCosTheta)*(1-standardCosTheta)*3*HelicityFraction[2] + HelicityFraction[1]*3*(1+standardCosTheta)*(1+standardCosTheta))/8;
	// 	  HelicityWeight = UsedDistributionValue/TheoreticalDistributionValue;	      	      
	// 	  histo1D[HelicityWeightHisto]->Fill(HelicityWeight);

	// 	  if(standardCosTheta == 0){
	// 	    std::cout << ievt << ") Still problem with Standard cos theta " << std::endl;
	// 	  }
	// 	}	       		       
		
	// 	if(dataSetName.find("Data") != 0){		  
	// 	  histo1D[THistoName]->Fill(CosTheta,(HelicityWeight*Luminosity*scaleFactor*datasets[d]->NormFactor()));  
	// 	}  // end of MC loop (not Data)	 	      

	//       } // end of second helicity loop
	//     } // end of first helicity loop  
	    
	//   }  //end of D>=0 loop
	// }// end of bTagSelected loop
      }  //delete selection;
    }//loop on events

    // /*//////////////////////////////////////////////////////////////////////////// --> No reweighting to PValue = 1 !!
    // //---    Compare Data values with simulated values for all helicities  ---//  
    // ////////////////////////////////////////////////////////////////////////////      
    // //  --> Should be outside eventloop since histograms should be filled
    // //    
    // std::string HighestPValueCosTheta;
    // std::string LowestChiSqValueCosTheta;
    // if(d==(datasets.size()-1)){ //Go in this loop when the last datasample is active
    //   int NumberHelicity = 0;
      
    //   //Loop to calculate the ChiSquaredValues:
    //   for(int Longit =0;Longit <=NumberOfHelicityBins;Longit++){
    // 	for(int Right =0; Right <=(NumberOfHelicityBins-Longit) ; Right++){	  
    // 	  int Left = NumberOfHelicityBins-Longit-Right;
	  
    // 	  double longit=Longit*(1.00/NumberOfHelicityBins);  
    // 	  double right=Right*(1.00/NumberOfHelicityBins);
    // 	  double left=Left*(1.00/NumberOfHelicityBins);	 
    // 	  std::string THistoName = "CosThetaRight"+HelicityNumber[Right]+"Long"+HelicityNumber[Longit]+"Left"+HelicityNumber[Left];
    // 	  std::string ChiSquaredHelicity = "ChiSquaredRight"+HelicityNumber[Right]+"Long"+HelicityNumber[Longit]+"Left"+HelicityNumber[Left];	       	  
	  
    // 	  //Loop over all bins of CosTheta (weighted) and CosThetaData to calculate ChiSquared for each helicity combination
    //  	  //
    // 	  //To compare MC and data scale them to equal number of entries:
    // 	  //
    // 	  float ScaleValue = histo1D["CosThetaData"]->Integral();
    // 	  histo1D[THistoName]->Scale(ScaleValue/histo1D[THistoName]->Integral());

    // 	  //std::cout << " Integral of Data : " << histo1D["CosThetaData"]->Integral() << " | Integral of MC (after scaling) : " << histo1D[THistoName]->Integral() << std::endl;

    // 	  float ChiSquaredAllBins =0;
    // 	  for(int ii=1;ii<=CosThetaBinNumber;ii++){  //Bin 0 = underflow bin  ==> Does not need to be considered!
    // 	    ChiSquaredAllBins=ChiSquaredAllBins+((((histo1D[THistoName]->GetBinContent(ii))-(histo1D["CosThetaData"]->GetBinContent(ii)))/(sqrt(histo1D["CosThetaData"]->GetBinContent(ii))))*(((histo1D[THistoName]->GetBinContent(ii))-(histo1D["CosThetaData"]->GetBinContent(ii)))/(sqrt(histo1D["CosThetaData"]->GetBinContent(ii)))));    //Total Chi Squared = Sum of chi squared for all bins  
    //  	  } // end of loop over all bins
	  
    // 	  NumberHelicity++;	  
    //    	  float PValueAllBins = TMath::Prob(ChiSquaredAllBins,CosThetaBinNumber);
    //    	  histo1D["PValueAllBins"]->SetBinContent(NumberHelicity,PValueAllBins);

    // 	  //Calculate highest p-value
    //    	  if(HighestPValue<PValueAllBins){
    //    	    HighestPValue = PValueAllBins;	    
    //    	    HighestPValueCosTheta = THistoName;
    //    	  }
	  
    //    	  float XBinValue = histo2D["HelicityPointsPlot"]->GetXaxis()->FindBin(right);	 
    //    	  float YBinValue = histo2D["HelicityPointsPlot"]->GetYaxis()->FindBin(longit);
    //    	  histo2D["HelicityPointsPlot"]->SetBinContent(XBinValue,YBinValue,PValueAllBins);
    //    	  histo2D["ChiSqPlot"]->SetBinContent(XBinValue,YBinValue,ChiSquaredAllBins);
    //    	  histo2D["HelicityPointsPlot"]->SetXTitle("f+");	  
    // 	  histo2D["HelicityPointsPlot"]->SetYTitle("f0");

    // 	} //end of second helicity loop
    //   } // end of first helicity loop 

    //  //Write out histogram with highest p-value
    //   TCanvas* HighestPValueCanvas = new TCanvas("HighestPValue","HighestPValue",1400,600);
    //   std::cout << " highest p value : " << HighestPValue << " corresponding to cos theta : " << HighestPValueCosTheta << std::endl;
    //   histo1D[HighestPValueCosTheta]->SetXTitle("Cos theta");
    //   histo1D[HighestPValueCosTheta]->SetYTitle("Entries");
    //   HighestPValueCanvas->cd();
    //   histo1D[HighestPValueCosTheta]->Draw();    
    //   HighestPValueCanvas->Modified();
    //   HighestPValueCanvas->SaveAs("PlotsMacro/HighestPValueCosTheta.png");    

    //   std::cout << "--------------- Effiency of b-tag analysis: -------------------- " << std::endl;
    //   std::cout << " Total number of studied events : " << BTagCombinationEvents << std::endl;
    //   std::cout << " Number of events with one of the b-tags correct : " << BTagCorrectCombination << std::endl;
    //   std::cout << " " << std::endl;
      
    //   std::cout << " Combination Zero  : " << CorrectZero << std::endl;
    //   std::cout << " Combination Two  : " << CorrectTwo << std::endl;
    //   std::cout << " Combination Four  : " << CorrectFour << std::endl;
    //   std::cout << " Combination Six  : " << CorrectSix << std::endl;
    //   std::cout << " Combination Eight  : " << CorrectEight << std::endl;
    //   std::cout << " Combination Ten  : " << CorrectTen << std::endl;
    //   std::cout << " " << std::endl;
      
    //   //---------------- Output for event selection efficiency ----------------
    //   std::cout << " total number of semi-muonic events : " << NumberSemiMuEvents << std::endl;
    //   std::cout << " All jets correctly matched : " << AllJetsCorrect << std::endl;
    //   std::cout << " One of the jets wrongly matched : " << OneJetWrong << std::endl;
    //   std::cout << " Hadronic B quark correct : " << numberCorrectHadronicB << std::endl;
    //   std::cout << " Light quarks correct : " << numberCorrectQuarks << std::endl;
      
    //   //--------------------  Sigma for W Mass and Top Mass  --------------------
    //   histo1D["WMass"]->Fit("gaus","Q");     
    //   histo1D["TopMass"]->Fit("gaus","Q");
    //   std::cout << " sigma values : " << histo1D["WMass"]->GetFunction("gaus")->GetParameter(2) << " " << histo1D["TopMass"]->GetFunction("gaus")->GetParameter(2) << std::endl;
    //   std::cout << " mass values : " << histo1D["WMass"]->GetFunction("gaus")->GetParameter(1) << " " << histo1D["TopMass"]->GetFunction("gaus")->GetParameter(1) << std::endl;
      
    //   //--------------------  Fit to collect used helicity values   --------------------
    //   // if(dataSetName.find("TTbarJets_SemiMu") == 0){
    //   //   TF1 *helicityFit = new TF1("helicityFit","((([0]*3*(1+x)*(1+x))+([1]*3*(1-x)*(1-x))+([2]*6*(1-x*x)))/8)",-1,1);
    //   //   histo1D["StandardCosThetaFit"]->Fit("helicityFit","Q");
    //   //   std::cout << " fit values (before event selection) : " << helicityFit->GetParameter(0) << " " << helicityFit->GetParameter(1) << " " << helicityFit->GetParameter(2) << std::endl;
    //   //   std::cout << " fit values error (before event selection) : " << helicityFit->GetParError(0) << " " << helicityFit->GetParError(1) << " " << helicityFit->GetParError(2) << std::endl;
    //   // }
    // }    // end of Data loop */ //--> Calculating pValues without putting the maximum on 1 !!
    
    // ////////////////////////////////////////////////////////////////////////////
    // //---    Compare Data values with simulated values for all helicities  ---//  
    // ////////////////////////////////////////////////////////////////////////////      
    // //  --> Should be outside eventloop since histograms should be filled
    // //    
    // std::string HighestPValueCosTheta;
    // std::string LowestChiSqValueCosTheta;
    // if(d==(datasets.size()-1)){ //Go in this loop when the last datasample is active
    //   int NumberHelicity = 0;
      
    //   //Loop to calculate the ChiSquaredValues:
    //   for(int Longit =0;Longit <=NumberOfHelicityBins;Longit++){
    // 	for(int Right =0; Right <=(NumberOfHelicityBins-Longit) ; Right++){	  
    // 	  int Left = NumberOfHelicityBins-Longit-Right;
	  
    // 	  double longit=Longit*(1.00/NumberOfHelicityBins);  
    // 	  double right=Right*(1.00/NumberOfHelicityBins);
    // 	  double left=Left*(1.00/NumberOfHelicityBins);	 
    // 	  std::string THistoName = "CosThetaRight"+HelicityNumber[Right]+"Long"+HelicityNumber[Longit]+"Left"+HelicityNumber[Left];
    // 	  std::string ChiSquaredHelicity = "ChiSquaredRight"+HelicityNumber[Right]+"Long"+HelicityNumber[Longit]+"Left"+HelicityNumber[Left];	       	  
	  
    // 	  //Loop over all bins of CosTheta (weighted) and CosThetaData to calculate ChiSquared for each helicity combination
    //  	  //
    // 	  float ChiSquaredAllBins =0;
    // 	  for(int ii=1;ii<=CosThetaBinNumber;ii++){  //Bin 0 = underflow bin  ==> Does not need to be considered!
    // 	    ChiSquaredAllBins=ChiSquaredAllBins+((((histo1D[THistoName]->GetBinContent(ii))-(histo1D["CosThetaData"]->GetBinContent(ii)))/(sqrt(histo1D["CosThetaData"]->GetBinContent(ii))))*(((histo1D[THistoName]->GetBinContent(ii))-(histo1D["CosThetaData"]->GetBinContent(ii)))/(sqrt(histo1D["CosThetaData"]->GetBinContent(ii)))));    //Total Chi Squared = Sum of chi squared for all bins  
    //  	  } // end of loop over all bins
	  
    // 	  //Calculate lowest ChiSquaredValue
    // 	  if(LowestChiSqValue>ChiSquaredAllBins){
    // 	    LowestChiSqValue = ChiSquaredAllBins;
    // 	    LowestChiSqValueCosTheta = THistoName;	    
    // 	  }
	  
    // 	} //end of second helicity loop
    //   } // end of first helicity loop 
      
    //   //Loop to calculate the ChiSquaredValue-LowestChiSquared  --> Makes it possible to compare the different results!!
    //   for(int Longit =0;Longit <=NumberOfHelicityBins;Longit++){
    //    	for(int Right =0; Right <=(NumberOfHelicityBins-Longit) ; Right++){
    //    	  int Left = NumberOfHelicityBins-Longit-Right;
	  
    //    	  double longit=Longit*(1.00/NumberOfHelicityBins);  
    //    	  double right=Right*(1.00/NumberOfHelicityBins);
    //   	  double left=Left*(1.00/NumberOfHelicityBins);	 
    //    	  std::string THistoName = "CosThetaRight"+HelicityNumber[Right]+"Long"+HelicityNumber[Longit]+"Left"+HelicityNumber[Left];
    //    	  std::string ChiSquaredHelicity = "ChiSquaredRight"+HelicityNumber[Right]+"Long"+HelicityNumber[Longit]+"Left"+HelicityNumber[Left];	       
	       	         	  
    //    	  //Loop over all bins of CosTheta (weighted) and CosThetaData to calculate ChiSquared for each helicity combination
    //    	  //
    //    	  float ChiSquaredAllBins =0;
    //    	  for(int ii=1;ii<=CosThetaBinNumber;ii++){  //Bin 0 = underflow bin  ==> Does not need to be considered!
    //    	    ChiSquaredAllBins=ChiSquaredAllBins+((((histo1D[THistoName]->GetBinContent(ii))-(histo1D["CosThetaData"]->GetBinContent(ii)))/(sqrt(histo1D["CosThetaData"]->GetBinContent(ii))))*(((histo1D[THistoName]->GetBinContent(ii))-(histo1D["CosThetaData"]->GetBinContent(ii)))/(sqrt(histo1D["CosThetaData"]->GetBinContent(ii)))));    //Total Chi Squared = Sum of chi squared for all bins  
    //    	  } // end of loop over all bins
	  
    //    	  ChiSquaredAllBins = ChiSquaredAllBins-LowestChiSqValue;

    //    	  NumberHelicity++;	  
    //    	  float PValueAllBins = TMath::Prob(ChiSquaredAllBins,CosThetaBinNumber);
    //    	  histo1D["PValueAllBins"]->SetBinContent(NumberHelicity,PValueAllBins);
	  
    //    	  //Calculate highest p-value
    //    	  if(HighestPValue<PValueAllBins){
    //    	    HighestPValue = PValueAllBins;	    
    //    	    HighestPValueCosTheta = THistoName;
    //    	  }
	  
    //    	  float XBinValue = histo2D["HelicityPointsPlot"]->GetXaxis()->FindBin(right);	 
    //    	  float YBinValue = histo2D["HelicityPointsPlot"]->GetYaxis()->FindBin(longit);
    //    	  histo2D["HelicityPointsPlot"]->SetBinContent(XBinValue,YBinValue,PValueAllBins);
    //    	  histo2D["ChiSqPlot"]->SetBinContent(XBinValue,YBinValue,ChiSquaredAllBins);
    //    	  histo2D["HelicityPointsPlot"]->SetXTitle("f+");	  
    // 	  histo2D["HelicityPointsPlot"]->SetYTitle("f0");	 
       	  
    //    	} //end of second helicity loop
    //   } // end of first helicity loop 
      
    //   //Write out histogram with highest p-value
    //   TCanvas* HighestPValueCanvas = new TCanvas("HighestPValue","HighestPValue",1400,600);
    //   std::cout << " highest p value : " << HighestPValue << " corresponding to cos theta : " << HighestPValueCosTheta << std::endl;
    //   histo1D[HighestPValueCosTheta]->SetXTitle("Cos theta");
    //   histo1D[HighestPValueCosTheta]->SetYTitle("Entries");
    //   HighestPValueCanvas->cd();
    //   histo1D[HighestPValueCosTheta]->Draw();    
    //   HighestPValueCanvas->Modified();
    //   HighestPValueCanvas->SaveAs("PlotsMacro/HighestPValueCosTheta.png");    

    //   std::cout << "--------------- Effiency of b-tag analysis: -------------------- " << std::endl;
    //   std::cout << " Total number of studied events : " << BTagCombinationEvents << std::endl;
    //   std::cout << " Number of events with one of the b-tags correct : " << BTagCorrectCombination << std::endl;
    //   std::cout << " " << std::endl;
      
    //   std::cout << " Combination Zero  : " << CorrectZero << std::endl;
    //   std::cout << " Combination Two  : " << CorrectTwo << std::endl;
    //   std::cout << " Combination Four  : " << CorrectFour << std::endl;
    //   std::cout << " Combination Six  : " << CorrectSix << std::endl;
    //   std::cout << " Combination Eight  : " << CorrectEight << std::endl;
    //   std::cout << " Combination Ten  : " << CorrectTen << std::endl;
    //   std::cout << " " << std::endl;
      
    //   //---------------- Output for event selection efficiency ----------------
    //   std::cout << " total number of semi-muonic events : " << NumberSemiMuEvents << std::endl;
    //   std::cout << " All jets correctly matched : " << AllJetsCorrect << std::endl;
    //   std::cout << " One of the jets wrongly matched : " << OneJetWrong << std::endl;
    //   std::cout << " Hadronic B quark correct : " << numberCorrectHadronicB << std::endl;
    //   std::cout << " Light quarks correct : " << numberCorrectQuarks << std::endl;
      
    //   //--------------------  Sigma for W Mass and Top Mass  --------------------
    //   histo1D["WMass"]->Fit("gaus","Q");     
    //   histo1D["TopMass"]->Fit("gaus","Q");
    //   std::cout << " sigma values : " << histo1D["WMass"]->GetFunction("gaus")->GetParameter(2) << " " << histo1D["TopMass"]->GetFunction("gaus")->GetParameter(2) << std::endl;
    //   std::cout << " mass values : " << histo1D["WMass"]->GetFunction("gaus")->GetParameter(1) << " " << histo1D["TopMass"]->GetFunction("gaus")->GetParameter(1) << std::endl;
      
    //   //--------------------  Fit to collect used helicity values   --------------------
    //   // if(dataSetName.find("TTbarJets_SemiMu") == 0){
    //   //   TF1 *helicityFit = new TF1("helicityFit","((([0]*3*(1+x)*(1+x))+([1]*3*(1-x)*(1-x))+([2]*6*(1-x*x)))/8)",-1,1);
    //   //   histo1D["StandardCosThetaFit"]->Fit("helicityFit","Q");
    //   //   std::cout << " fit values (before event selection) : " << helicityFit->GetParameter(0) << " " << helicityFit->GetParameter(1) << " " << helicityFit->GetParameter(2) << std::endl;
    //   //   std::cout << " fit values error (before event selection) : " << helicityFit->GetParError(0) << " " << helicityFit->GetParError(1) << " " << helicityFit->GetParError(2) << std::endl;
    //   // }
    // }    // end of Data loop  
    
    //important: free memory
    treeLoader.UnLoadDataset();
  }				//loop on datasets  -->Stops already here
  
  //Once everything is filled ...
  if (verbose > 0)
    cout << " We ran over all the data ;-)" << endl;
  
  //Selection tables
  selecTable.TableCalculator(false, true, true, true);
  string selectiontable = "SelectionTable_Macro";
  if (argc >= 3){
    string sample=string(argv[2]);
    selectiontable = selectiontable +"_"+sample;
  }
  selectiontable = selectiontable +".tex"; 	
  selecTable.Write(selectiontable.c_str());
  
  // Do some special things with certain plots (normalize, BayesDivide, ... )
  if (verbose > 0)
    cout << "Treating the special plots." << endl;
  
  ///////////////////
  // Writting
  //////////////////
  if (verbose > 1)
    cout << " - Writing outputs on files ..." << endl;  
  
  string pathPNG = "CernStudy/PValueMaximumInfluence/PlotsMacro";
  if (argc >= 3){
    string sample=string(argv[2]);
    pathPNG = pathPNG+"_"+sample;
  }
  pathPNG = pathPNG +"/"; 	
  
  mkdir(pathPNG.c_str(),0777);
  mkdir((pathPNG+"MSPlot/").c_str(),0777);
  
  // 1D 
  TDirectory* th1dir = fout->mkdir("1D_histograms");
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++){  //Because all the histograms are stored as maps it is quite easy to loop over them
    
    MultiSamplePlot *temp = it->second;
    string name = it->first;
    temp->Draw(false, name, true, true, true);
    temp->Write(fout, name, true, pathPNG+"MSPlot/");
  }
  
  //Write histograms
  fout->cd();
  th1dir->cd();
  
  fout->cd();
  
  for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++){
    TH1F *temp = it->second;
    string name = it->first;
    if( name.find("HelicityWeightRight") !=0 && name.find("CosThetaRight") != 0 ){//Does not want to write out these 5151 histograms (2X)
      int N = temp->GetNbinsX();
      temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
      temp->SetBinContent(N+1,0);
      temp->Write();
      TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
      tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
    }
  }
  
  // 2D
  TDirectory* th2dir = fout->mkdir("2D_histograms_graphs");
  th2dir->cd();
  for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++){
    TH2F *temp = it->second;
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
  }
  
  //Write TGraphAsymmErrors
  for(map<string,TGraphAsymmErrors*>::const_iterator it = graphAsymmErr.begin(); it != graphAsymmErr.end(); it++){
    TGraphAsymmErrors *temp = it->second;
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
  }  
  
  fout->cd();
  //add configuration info
  fout->cd();
  configTree->Fill();
  configTree->Write();
  
  //Write TGraphErrors
  fout->cd();
  for(map<string,TGraphErrors*>::const_iterator it = graphErr.begin(); it != graphErr.end(); it++){
    TGraphErrors *temp = it->second;
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
  }
  
  //
  if (verbose > 1)
    cout << " - Done with writing the module outputs in the ouput file ..." << endl;
  cout << " - Closing the output file now..." << endl;
  //  fout->Write();
  fout->Close();
  
  delete fout;
  delete tcdatasets;
  delete tcAnaEnv;
  delete configTree;
  
  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
  
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "           hasn't crashed yet ;-)           " << endl;
  cout << "********************************************" << endl;
  
  return 0;
}
