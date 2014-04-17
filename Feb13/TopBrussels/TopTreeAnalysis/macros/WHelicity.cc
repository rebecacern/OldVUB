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
#include "TLorentzVector.h"
#include "TF1.h"
#include "TString.h"
#include <stdio.h>

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
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1.h"

using namespace std;
using namespace TopTree;

int main (int argc, char *argv[])
{
  clock_t start = clock();
  
  cout << "*****************************************" << endl;
  cout << " Beginning of the program for WHelicity ! " << endl;
  cout << "*****************************************" << endl;
  
  //SetStyle if needed
  //setTDRStyle(); 
  setMyStyle();
  
  /////////////////////
  // Configuration
  /////////////////////
  
  bool doPF2PAT = false;
  
  //xml file
  string xmlFileName ="../config/myMacroConfig.xml";   //Position of .xml file
  if (argc >= 2)
    xmlFileName=string(argv[1]);
  const char *xmlfile = xmlFileName.c_str();
  
  cout << "used config file: " << xmlfile << endl;
  
  //Output ROOT file
  string rootFileName ("MacroOutput.root");         //Root output file
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
  
  TTreeLoader treeLoader;
  cout << " - Load datasets ..." << endl;
  vector < Dataset* > datasets;
  treeLoader.LoadDatasets (datasets, xmlfile);
  for(unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
  
  float Luminosity = oldLuminosity;
  bool isSemiMu = false;
  for (unsigned int d = 0; d < datasets.size (); d++){
    //cout << "luminosity of dataset "<< d << " is " << datasets[d]->EquivalentLumi() << endl;
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
  
  ////////////////////////////////////
  /// Normal Plots (TH1F* and TH2F*)
  ////////////////////////////////////
  
  map<string,TH1F*> histo1D;     //All histograms can be defined as pre-programmed maps which makes definitions and looping easier
  map<string,TH2F*> histo2D;
  
//  histo1D["WMass"]= new TH1F("WMass","WMass", 200,0,160);
//  histo1D["TopMass"]= new TH1F("TopMass","TopMass", 200,0,350);

  histo1D["HelicityDistribution"] = new TH1F("HelicityDistribution","HelicityDistribution",200,-1,1);
  histo1D["HelicityDistributionLong"] = new TH1F("HelicityDistributionLong","HelicityDistributionLong",200,-1,1);
  histo1D["HelicityDistributionLeft"] = new TH1F("HelicityDistributionLeft","HelicityDistributionLeft",200,-1,1);
  histo1D["HelicityDistributionRight"] = new TH1F("HelicityDistributionRight","HelicityDistributionRight",200,-1,1);
  histo2D["HelicityPoints"] = new TH2F("HelicityPoints","HelicityPoints",100,-0.2,1.2,100,-0.2,1.2);

  //histo1D["HelicityDistribution424"] = new TH1F("HelicityDistribution424","HelicityDistribution424",200,-1,1);
  //histo1D["HelicityDistributionLong424"] = new TH1F("HelicityDistributionLong424","HelicityDistributionLong424",200,-1,1);
  //histo1D["HelicityDistributionLeft424"] = new TH1F("HelicityDistributionLeft424","HelicityDistributionLeft424",200,-1,1);
  //histo1D["HelicityDistributionRight424"] = new TH1F("HelicityDistributionRight424","HelicityDistributionRight424",200,-1,1);
  
  histo1D["StandardCosTheta"]=new TH1F("StCosTheta","StCosTheta",200,-1,1);
  histo1D["StandardCosThetaFit"]=new TH1F("StCosThetaFit","StCosThetaFit",200,-1,1);
//  histo1D["McCosTheta"]= new TH1F("McCosTheta","McCosTheta",200,-1,1);
  
  //------  Control plots for Mc particle selection  ------
  histo1D["McTopMass"]=new TH1F("McTopMass","McTopMass",200,0,300);
  histo1D["McWMass"]=new TH1F("McWMass","McWMass",200,0,160);
  histo1D["McMuonMass"]=new TH1F("McMuonMass","McMuonMass",200,0,0.5);
  histo1D["McMuonPx"]=new TH1F("McMuonPx","McMuonPx",200,-100,100);  
  histo1D["McMuonPy"]=new TH1F("McMuonPy","McMuonPy",200,-100,100);
  histo1D["McMuonPz"]=new TH1F("McMuonPz","McMuonPz",200,-100,100);
  histo1D["McMuonE"]=new TH1F("McMuonE","McMuonE",200,-100,100);
  histo1D["McBottomMass"]=new TH1F("McBottomMass","McBottomMass",200,0,10);
  histo1D["McCorrectTopMass"]=new TH1F("McCorrectTopMass","McCorrectTopMass",200,0,300);
  histo1D["McTopMassMc"]=new TH1F("McTopMassMc","McTopMassMc",200,0,300);
  histo1D["TopEventTypes"]=new TH1F("TopEventTypes","TopEventTypes",200,-30,30);
  histo1D["AntiTopEventTypes"]=new TH1F("AntiTopEventTypes","AntiTopEventTypes",200,-30,30);

  histo1D["TopPt"]= new TH1F("TopPt","TopPt",200,0,600);
  histo1D["WPt"]= new TH1F("WPt","WPt",200,0,500);
  histo1D["MuonPt"]= new TH1F("MuonPt","MuonPt",200,0,300);
  //histo1D["BottomPt"]= new TH1F("BottomPt","BottomPt",200,-1000,1000);
  //histo1D["NeutrinoPt"]= new TH1F("NeutrinoPt","NeutrinoPt",200,-1000,1000);
  histo1D["TopEta"]= new TH1F("TopEta","TopEta",100,-7,7);
  histo1D["WEta"]= new TH1F("WEta","WEta",100,-7,7);
  histo1D["MuonEta"]= new TH1F("MuonEta","MuonEta",100,-7,7);
  //histo1D["BottomEta"]= new TH1F("BottomEta","BottomEta",100,-7,7);
  //histo1D["NeutrinoEta"]= new TH1F("NeutrinoEta","NeutrinoEta",100,-7,7);
  histo1D["TopPhi"]= new TH1F("TopPhi","TopPhi",50,-3.2,3.2);
  histo1D["WPhi"]= new TH1F("WPhi","WPhi",50,-3.2,3.2);
  histo1D["MuonPhi"]= new TH1F("MuonPhi","MuonPhi",50,-3.2,3.2);
  //histo1D["BottomPhi"]= new TH1F("BottomPhi","BottomPhi",50,-3.2,3.2);
  //histo1D["NeutrinoPhi"]= new TH1F("NeutrinoPhi","NeutrinoPhi",50,-3.2,3.2);
  histo1D["TopAndBottomAngle"]=new TH1F("TopAndBottomAngle","TopAndBottomAngle",100,-1,1.2);
  histo1D["TopAndWAngle"]=new TH1F("TopAndWAngle","TopAndWAngle",100,-1,1.2);
  histo1D["WAndBAngle"]= new TH1F("WAndBAngle","WAndBAngle",100,-1,1.2);

    
  ////////////////////////////////////
  /// MultiSamplePlot
  ////////////////////////////////////
  
  map<string,MultiSamplePlot*> MSPlot;        //Whenever there are more than one dataset this can be used to add up the different results of the datasets
  
  // Variables used during event selection
  
//  MSPlot["ChiSquared"]= new MultiSamplePlot(datasets, "ChiSquared", 200,0,50,"ChiSquared");
  //MSPlot["ChiSquaredCorrect"]= new MultiSamplePlot(datasets, "ChiSquaredCorrect", 200,0,50,"ChiSquaredCorrect");
  //MSPlot["ChiSquaredWrong"]= new MultiSamplePlot(datasets, "ChiSquaredWrong", 200,0,50,"ChiSquaredWrong");
  
  //MSPlot["NeutrinoPzOne"]= new MultiSamplePlot(datasets, "NeutrinoPzOne", 200,-1000,1000,"NeutrinoPzOne");
  //MSPlot["NeutrinoPzTwo"]= new MultiSamplePlot(datasets, "NeutrinoPzTwo", 1000,-1000,1000,"NeutrinoPzTwo"); 
//  MSPlot["McNeutrinoPz"] = new MultiSamplePlot(datasets, "McNeutrinoPz", 1000,-1000,1000,"McNeutrinoPz");
//  MSPlot["NeutrinoPzSelected"]=new MultiSamplePlot(datasets, "NeutrinoPzSelected", 1000,-1000,1000,"NeutrinoPzSelected");
  //MSPlot["NeutrinoPzCompared"]=new MultiSamplePlot(datasets, "NeutrinoPzCompared", 250,-20,20,"NeutrinoPzCompared");
  
  // MSPlot["WMassNeutrinoOne"]= new MultiSamplePlot(datasets, "WMassNeutrinoOne", 200,60,100,"WMassNeutrinoOne");
  //MSPlot["WMassNeutrinoTwo"]= new MultiSamplePlot(datasets, "WMassNeutrinoTwo", 200,60,100,"WMassNeutrinoTwo");
  //MSPlot["WMassNeutrinoMc"]= new MultiSamplePlot(datasets, "WMassNeutrinoMc", 200,0,110,"WMassNeutrinoMc");

  //  MSPlot["NeutrinoMassBoost"]= new MultiSamplePlot(datasets, "NeutrinoMassBoost", 200,0,150,"NeutrinoMassBoost");
  //  MSPlot["WMassVector"]= new MultiSamplePlot(datasets, "WMassVector", 200,0,150,"WMassVector");
  //  MSPlot["WMassBoost"]= new MultiSamplePlot(datasets, "WMassBoost", 200,0,110,"WMassBoost");
  //  MSPlot["CosTheta"]= new MultiSamplePlot(datasets, "CosTheta", 4,0,1,"CosTheta");
  //  MSPlot["CosThetaMcNeutrinoMcTop"]= new MultiSamplePlot(datasets, "CosThetaMcNeutrinoMcTop", 4,0,1,"CosThetaMcNeutrinoMcTop");
  //  //MSPlot["CosThetaMcTop"]= new MultiSamplePlot(datasets, "CosThetaMcTop", 4,0,1,"CosThetaMcTop");
  //  MSPlot["CosThetaMcNeutrino"]= new MultiSamplePlot(datasets, "CosThetaMcNeutrino", 4,0,1,"CosThetaMcNeutrino");
  //MSPlot["standardCosTheta"]= new MultiSamplePlot(datasets, "standardCosTheta", 200,-1,1,"standardCosTheta");
  //MSPlot["DifferentialDistribution"] = new MultiSamplePlot(datasets, "DifferentialDistribution",100,-1,1,"DifferentialDistribution");
  
  //MSPlot["TopMassNeutrinoOne"]=new MultiSamplePlot(datasets, "TopMassNeutrinoOne", 700,0,450,"TopMassNeutrinoOne");
  //MSPlot["TopMassNeutrinoTwo"]=new MultiSamplePlot(datasets, "TopMassNeutrinoTwo", 700,0,450,"TopMassNeutrinoTwo");
//  MSPlot["TopMassNeutrinoSelected"]=new MultiSamplePlot(datasets, "TopMassNeutrinoSelected", 700,0,450,"TopMassNeutrinoSelected");
  //MSPlot["TopMassGen"]=new MultiSamplePlot(datasets, "TopMassGen", 700,0,450,"TopMassGen");
//  MSPlot["TopMassMC"]=new MultiSamplePlot(datasets, "TopMassMc", 700,0,450,"TopMassMc");
  //MSPlot["TopAfterBoost"]=new MultiSamplePlot(datasets, "TopAfterBoost", 1200, 0, 650,"TopAfterBoost");
  //MSPlot["TopPzAfterBoost"]=new MultiSamplePlot(datasets, "TopPzAfterBoost", 100, 0, 4,"TopPzAfterBoost");
  
  //MSPlot["HighTopB"]=new MultiSamplePlot(datasets, "HighTopB", 250,0,20,"HighTopB");
  //MSPlot["HighTopMuon"]=new MultiSamplePlot(datasets, "HighTopMuon", 50,0,3,"HighTopMuon");
  //MSPlot["HighTopNeutrino"]=new MultiSamplePlot(datasets, "HighTopNeutrino", 500,0,250,"HighTopNeutrino");

//  MSPlot["PosMuonOneMass"]=new MultiSamplePlot(datasets, "PosMuonOneMass",500,-300,300,"PosMuonOneMass");
//  MSPlot["PosMuonTwoMass"]=new MultiSamplePlot(datasets, "PosMuonTwoMass",500,-300,300,"PosMuonTwoMass");
//  MSPlot["NegMuonOneMass"]=new MultiSamplePlot(datasets, "NegMuonOneMass",500,-300,300,"NegMuonOneMass");
//  MSPlot["NegMuonTwoMass"]=new MultiSamplePlot(datasets, "NegMuonTwoMass",500,-300,300,"NegMuonTwoMass");

//  MSPlot["McMuonOne"]= new MultiSamplePlot(datasets, "McMuonOne",500,-0.2,0.2,"McMuonOne");
//  MSPlot["McMuonTwo"]= new MultiSamplePlot(datasets, "McMuonTwo",500,-0.2,0.2,"McMuonTwo");

  ///////////////////////////////////////////////////////
  //   Making the theoretical helicity distribution    //
  ///////////////////////////////////////////////////////
  //------   Standard Model   -----
  float LongitudinalFraction = 0.642324; // values obtained from fit
  float RightHandedFraction = 0.0327637;
  float LeftHandedFraction = 1 - LongitudinalFraction - RightHandedFraction;
  float SMDiffDistribution[3];//0:Longitudinal; 1: Righthanded; 2: Lefthanded
  float XVariable;
  std::cout << " Used SM helicity variables: Longitudinal (f0) = " << LongitudinalFraction << " Righthanded (f+) = " << RightHandedFraction << " Lefthanded (f-) = " << LeftHandedFraction << std::endl;
  for(unsigned int jj=0;jj<200;jj++){
    XVariable=-1+0.01*jj;
    
    SMDiffDistribution[0] = (LongitudinalFraction*6*(1-pow(XVariable,2)))/8;
    SMDiffDistribution[2] = (pow((1-XVariable),2)*3*LeftHandedFraction)/8;
    SMDiffDistribution[1] = (RightHandedFraction*3*pow((1+XVariable),2))/8;
    
    histo1D["HelicityDistribution"]->SetBinContent(jj+1,(SMDiffDistribution[0]+SMDiffDistribution[1]+SMDiffDistribution[2]));
    histo1D["HelicityDistributionLong"]->SetBinContent(jj+1,SMDiffDistribution[0]);
    histo1D["HelicityDistributionLeft"]->SetBinContent(jj+1,SMDiffDistribution[2]);
    histo1D["HelicityDistributionRight"]->SetBinContent(jj+1,SMDiffDistribution[1]);
  }
  
  //-----  Non Standard Model values   -----
  float HelicityFraction[3];  //0:Longitudinal; 1: Righthanded; 2: Lefthanded
  float DiffDistribution[3];  //0:Longitudinal; 1: Righthanded; 2: Lefthanded
  std::string Number[11]={"0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"};
  for(double longit =0;longit <=1;longit=(longit+0.1)){
    for(double right =0; right <=(1-longit) ; right = right+0.1){
      double left = 1 - longit - right;
      HelicityFraction[0]=longit;
      HelicityFraction[1]=right;
      HelicityFraction[2]=left;

      // Name giving:
      int Longit=longit*10;
      int Right=right*10;
      int Left=left*10;
      TString Helicity = "Right"+Number[Right]+"Long"+Number[Longit]+"Left"+Number[Left];
      TString HistoName = "HelicityDistribution"+Helicity;
      TString WeightsName = "WeightsFor"+Helicity;
      TString CanvasName = "Canvas"+Helicity;
      TString WeightsCanvasName = "WeightsCanvas"+Helicity;

      TH1F * Histo = new TH1F(HistoName,HistoName,200,-1,1);
      TH1F * WeightsHisto = new TH1F(WeightsName,WeightsName,200,-1,1);
      TCanvas * Canvas = new TCanvas(CanvasName,CanvasName,1400,600);
      TCanvas * WeightsCanvas = new TCanvas(WeightsCanvasName,WeightsCanvasName,1400,600);

      for(unsigned int jj=0;jj<200;jj++){
	XVariable=-1+0.01*jj;
	
	DiffDistribution[0] = (HelicityFraction[0]*6*(1-pow(XVariable,2)))/8;
	DiffDistribution[2] = (pow((1-XVariable),2)*3*HelicityFraction[2])/8;
	DiffDistribution[1] = (HelicityFraction[1]*3*pow((1+XVariable),2))/8;
	  
	Histo->SetBinContent(jj+1,(DiffDistribution[0]+DiffDistribution[1]+DiffDistribution[2]));
      }
      Canvas->cd();
      Histo->Draw("hist");
      Histo->SetMaximum(1.6);
      Histo->SetMinimum(0);
      Histo->SetXTitle("Cos #theta");
      Histo->SetYTitle("Differential distribution");
      Canvas->Modified();
      TString PrintName = "PlotsMacro/HelicityDistribution/"+HistoName+".png";
      Canvas->Print(PrintName);
      
      WeightsCanvas->cd();
      WeightsHisto->Draw("hist");
      WeightsHisto->SetMinimum(0);
      WeightsHisto->Divide(Histo,histo1D["HelicityDistribution"],1,1);
      WeightsHisto->SetXTitle("Cos #theta");
      WeightsHisto->SetYTitle("Ratio");
      WeightsCanvas->Modified();
      TString WeightsPrintName = "PlotsMacro/Weights/"+WeightsName+".png";
      WeightsCanvas->Print(WeightsPrintName);
      
      histo2D["HelicityPoints"]->Fill(right,longit);
      histo2D["HelicityPoints"]->SetXTitle("f+");
      histo2D["HelicityPoints"]->SetYTitle("f0");      
    }
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
    
    /////////////////////////////////////
    /// Initialize JEC factors
    /////////////////////////////////////
   	
    vector<JetCorrectorParameters> vCorrParam;
    if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") // Data!
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("JECFiles/Fall10_L1Offset_AK5PF.txt");
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("JECFiles/GR_R_38X_V15_AK5PF_L2Relative.txt");
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("JECFiles/GR_R_38X_V15_AK5PF_L3Absolute.txt");
      JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("JECFiles/GR_R_38X_V15_AK5PF_L2L3Residual.txt");
      vCorrParam.push_back(*L1JetCorPar);
      vCorrParam.push_back(*L2JetCorPar);
      vCorrParam.push_back(*L3JetCorPar);
      vCorrParam.push_back(*ResJetCorPar);
    }
    else // MC!
    {
     	JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("JECFiles/START38_V14_AK5PF_L2Relative.txt");
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("JECFiles/START38_V14_AK5PF_L3Absolute.txt");
      vCorrParam.push_back(*L2JetCorPar);
      vCorrParam.push_back(*L3JetCorPar);
    }
    
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/START38_V14_AK5PF_Uncertainty.txt");
    
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);
    
    ////////////////////////////////////
    //	Loop on events
    ////////////////////////////////////
    
    nEvents[d] = 0;
    int itrigger = -1, previousRun = -1;
    if (verbose > 1)
      cout << "	Loop over events " << endl;

    int NumberOfStandardHits =0;
    int NumberOfAlternativeHits = 0;
    
    for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++){     //In this loop plots before selection can be defined
      //    for (unsigned int ievt = 0; ievt < 10000; ievt++)
      
      nEvents[d]++;
      if(ievt%1000 == 0)
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
      
      // Correct the JetId fractions for PFJets
      if (init_jets_corrected.size() > 0)
      {
        if(init_jets_corrected[0]->jetType() == 2 || doPF2PAT)
        {
          for(unsigned int j=0; j<init_jets_corrected.size(); j++)
          {
            TRootPFJet* pfJet = jetTools->convertToPFJets(init_jets_corrected[j]);
            pfJet->setChargedHadronEnergyFraction(pfJet->chargedHadronEnergyFraction()*pfJet->getJetCorrFactor("L1L2L3"));
            pfJet->setNeutralHadronEnergyFraction(pfJet->neutralHadronEnergyFraction()*pfJet->getJetCorrFactor("L1L2L3"));
            pfJet->setChargedEmEnergyFraction(pfJet->chargedEmEnergyFraction()*pfJet->getJetCorrFactor("L1L2L3"));
            pfJet->setChargedMuEnergyFraction(pfJet->chargedMuEnergyFraction()*pfJet->getJetCorrFactor("L1L2L3"));
            pfJet->setNeutralEmEnergyFraction(pfJet->neutralEmEnergyFraction()*pfJet->getJetCorrFactor("L1L2L3"));
          }
        }
      }
      
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
	  //cout << "data" << endl;
	  if (event->runId() < 147196)
	    itrigger = treeLoader.iTrigger (string ("HLT_Mu9"), currentRun, iFile);
	  else if (event->runId() >= 147196)
	    itrigger = treeLoader.iTrigger (string ("HLT_Mu15_v1"), currentRun, iFile);
	}
	else{
          itrigger = treeLoader.iTrigger (string ("HLT_Mu9"), currentRun);
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
      jetTools->correctJets(init_jets_corrected, vertex);
      
      if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) ){
        // Correct for the difference in muon efficiency (HLT and Id) between Data and MC
        scaleFactor *= 0.965;
        
        // Correct jets for JES uncertainy systematics
	//jetTools->correctJetJESUnc(init_jets_corrected, "minus");
        
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
            init_jets_corrected[indexVector[i].first]->SetPxPyPzE(init_jets_corrected[indexVector[i].first]->Px()*ptscale, init_jets_corrected[indexVector[i].first]->Py()*ptscale,
								  init_jets_corrected[indexVector[i].first]->Pz()*ptscale, init_jets_corrected[indexVector[i].first]->E()*ptscale);
        }
        
        //Scale jets with a certain factor
        jetTools->scaleJets(init_jets_corrected, 1.);
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
    
      //////////////////////////////////////////////////
      ///     Copy of Mc: Standard simulation        ///
      //////////////////////////////////////////////////

      float MassW=80.398;
      float MassTop = 172.5;
      float SigmaW=11.0354;  //Obtained from gaussian fit on Top and W distribution with simulated information
      float SigmaTop=18.0629;
      float CorrectRecMassW=0;
      float CorrectRecMassTop=0;
      int CorrectQuark1=999;
      int CorrectQuark2=999;
      int CorrectBHadronic=999;
      int CorrectBLeptonic=999;
      int neutrinoFound=0;
      int topFound=0;
      int muonFound=0;
      
      int NumberCombinations=0;
      float recMassW[12];
      float recMassTop[12];
      int UsedCombination;
      float ChiSquared[12];
      int QuarkOneIndex[12];
      int QuarkTwoIndex[12];
      int BHadronicIndex[12];
      int BLeptonicIndex;
      
      float ChiSquaredValue;
      float NeutrinoPx;
      float NeutrinoPy;
      float NeutrinoPz=999;//with this value it can be distinguished in plot!
      float NeutrinoPzOne;  //Pz belonging to first neutrino solution
      TLorentzVector NeutrinoOne;
      float NeutrinoEOne;
      TLorentzVector TopOne;
      float NeutrinoPzTwo;  //Pz belonging to second neutrino solution
      TLorentzVector NeutrinoTwo;
      float NeutrinoETwo;
      TLorentzVector TopTwo;
      
      TLorentzVector WLeptonic;
      TLorentzVector TopLeptonic;
      TLorentzVector Neutrino;
      TLorentzVector MuonInWFrame;
      TLorentzVector WInTopFrame;
      
      float CosTheta;
      float CosThetaMuonMc;  //Op deze manier effect van de drie deeltjes begrijpen (kunnen ook andere deeltjes nog zijn)
      float CosThetaWMc;
      float CosThetaTopMc;

      TRootMCParticle standardNeutrino, standardTop,standardMuon,standardWLeptonic;
      float standardCosTheta=0;
      int standardNeutrinoFound=0;
      int standardTopFound=0;
      int standardMuonFound=0;

      vector<int> jetCombi;
      if(dataSetName.find("TTbarJets_SemiMu") == 0){
	  
	// loop for helicities
	for(float left =0;left <=1;left=left+0.1){
	  for(float right =0; right <=1 ; right = right+0.1){
	    float longit = 1 - left - right;
	  }
	}
	
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
	int EventChargeTop=0;
	int EventChargeAntiTop=0;//1 for top, -1 for anti-top
	int EventChargeWPos=0;
	int EventChargeWNeg=0; //1 for W+, -1 for W-
	int EventChargeMuPos=0;
	int EventChargeMuNeg =0; //1 for mu+, -1 for mu-
	TRootMCParticle Top,AntiTop,WPos,WNeg,MuonNeg,MuonPos,BQuark,BBarQuark,Neutrino;
	//TLorentzVector standardMuonWZMF, standardWLeptonicTZMF;
	TRootMCParticle standardMuonWZMF, standardWLeptonicTZMF;

	for(unsigned int i=0; i<mcParticles.size(); i++){

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
	    histo1D["TopPt"]->Fill(mcParticles[i]->Pt());
	    histo1D["TopEta"]->Fill(mcParticles[i]->Eta());
	    histo1D["TopPhi"]->Fill(mcParticles[i]->Phi());
	  }

	  if(fabs(mcParticles[i]->type()) == 5 && fabs(mcParticles[i]->motherType()) == 6){//Identifying bottom quarks
	    EventParticleNumber[1]++;
	    if(mcParticles[i]->type()==5){BQuark=*mcParticles[i];}
	    else BBarQuark=*mcParticles[i];
	  } 
	  if(fabs(mcParticles[i]->type()) <= 4 && fabs(mcParticles[i]->motherType()) == 24 && fabs(mcParticles[i]->grannyType()) == 6){EventParticleNumber[2]++;}

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
	    histo1D["WPt"]->Fill(mcParticles[i]->Pt());
	    histo1D["WEta"]->Fill(mcParticles[i]->Eta());
	    histo1D["WPhi"]->Fill(mcParticles[i]->Phi());
	  }
	  if(fabs(mcParticles[i]->type()) == 14 && fabs(mcParticles[i]->motherType()) == 24 && fabs(mcParticles[i]->grannyType()) == 6){
	    EventParticleNumber[4]++;  //Identifying neutrino's
	    Neutrino=*mcParticles[i];
	  }
	  if(fabs(mcParticles[i]->type()) == 13 && fabs(mcParticles[i]->motherType()) == 24 && fabs(mcParticles[i]->grannyType()) == 6 && mcParticles[i]->status()==3){ //status: 1:stable; 2:shower; 3:hard scattering(coming from the studied hard proces)
	    MuonNumber++; //Identifying muons
	    EventParticleNumber[4]++;
	    if(mcParticles[i]->type()==13){
	      MuonNeg=*mcParticles[i];
	      EventChargeMuNeg=-1;
	    }
	    else if(mcParticles[i]->type()==-13){
	      MuonPos=*mcParticles[i];
	      EventChargeMuPos=1;
	    }
	    histo1D["MuonPt"]->Fill(mcParticles[i]->Pt());
	    histo1D["MuonEta"]->Fill(mcParticles[i]->Eta());
	    histo1D["MuonPhi"]->Fill(mcParticles[i]->Phi());
	  }	  	  
	}//  if 0 < i < mcParticles.size()
		
	//////////////////////////////////////////////////
	//   Selecting correct event (b b q q mu nu )   //
	//////////////////////////////////////////////////
	if(EventParticleNumber[0]==2 && EventParticleNumber[1]==2 && EventParticleNumber[2]==2 && EventParticleNumber[3]==2 && EventParticleNumber[4]==2){
	  TRootMCParticle Bottom;
	  
	  // Differentiating between proces from top and anti-top:
	  if(EventChargeTop==1 && EventChargeWPos==1 && EventChargeMuPos==1){  // Proces: t -> b W+ -> mu+ nu
	    standardMuonWZMF=MuonPos;
	    standardWLeptonicTZMF=WPos;
	    standardWLeptonic=WPos;
	    standardTop=Top;
	    Bottom = BQuark;
	  }
	  else if(EventChargeAntiTop==-1 && EventChargeWNeg==-1 && EventChargeMuNeg==-1){ //proces: anti-t -> anti-b W- -> mu- anti-nu	   
	    standardMuonWZMF=MuonNeg;
	    standardWLeptonicTZMF=WNeg;
	    standardWLeptonic=WNeg;
	    standardTop=AntiTop;
	    Bottom = BBarQuark;
	  }

	  if(standardTop.type()==6){  // Types are correct now!
	    histo1D["TopEventTypes"]->Fill(standardTop.type());
	    histo1D["TopEventTypes"]->Fill(standardMuonWZMF.type());
	    histo1D["TopEventTypes"]->Fill(standardWLeptonicTZMF.type());
	    histo1D["TopEventTypes"]->Fill(Bottom.type());
	    histo1D["TopEventTypes"]->Fill(Neutrino.type());
	  }
	  if(standardTop.type()==-6){
	    histo1D["AntiTopEventTypes"]->Fill(standardTop.type());
	    histo1D["AntiTopEventTypes"]->Fill(standardMuonWZMF.type());
	    histo1D["AntiTopEventTypes"]->Fill(standardWLeptonicTZMF.type());
	    histo1D["AntiTopEventTypes"]->Fill(Bottom.type());
	    histo1D["AntiTopEventTypes"]->Fill(Neutrino.type());
	  }

	  // Calculate W mass from muon and neutrino and top mass from muon, neutrino and b leptonic:
	  histo1D["McTopMass"]->Fill((Bottom+Neutrino+standardMuonWZMF).M());
	  histo1D["McBottomMass"]->Fill(Bottom.M());
	  histo1D["McTopMassMc"]->Fill(standardTop.M());
	  float TopMassMc=(Bottom+Neutrino+standardMuonWZMF).M(); //How can top mass differ so much when W and Bottom have the correct mass?
	  	    
	  // Applying boost on muon and W:
	  standardMuonWZMF.Boost(-standardWLeptonic.BoostVector());
	  standardWLeptonicTZMF.Boost(-standardTop.BoostVector());
	  
	  standardCosTheta = ((standardWLeptonicTZMF.Vect()).Dot(standardMuonWZMF.Vect()))/(((standardWLeptonicTZMF.Vect()).Mag())*((standardMuonWZMF.Vect()).Mag()));
	  histo1D["StandardCosTheta"]->Fill(standardCosTheta);  // Histogram without fit
	  histo1D["StandardCosThetaFit"]->Fill(standardCosTheta);  // Histogram with fit
	 	    
	
	  //-----  Check whether cuts are applied on MonteCarlo  -----
	  float TopAndBottomAngle = ((standardTop.Vect()).Dot(Bottom.Vect()))/(((standardTop.Vect()).Mag())*((Bottom.Vect()).Mag()));
	  histo1D["TopAndBottomAngle"]->Fill(TopAndBottomAngle);
	  float TopAndWAngle = ((standardTop.Vect()).Dot(standardWLeptonic.Vect()))/(((standardTop.Vect()).Mag())*((standardWLeptonic.Vect()).Mag()));
	  histo1D["TopAndWAngle"]->Fill(TopAndWAngle);
	  float WAndBAngle = (((Bottom.Vect()).Dot(standardWLeptonic.Vect()))/(((Bottom.Vect()).Mag())*((standardWLeptonic.Vect()).Mag())));
	  histo1D["WAndBAngle"]->Fill(WAndBAngle);	      								      	 	
	}
      }//if dataset Semi mu ttbar   
      
      // /////////////////////////
      // //   Event selection   //
      // /////////////////////////
      // bool eventSelected = false;      
      // selecTable.Fill(d,1,1.);   
      
      // if(trigged){   //selection steps start: https://twiki.cern.ch/twiki/bin/view/CMS/TopLeptonPlusJetsRefSel_mu#Selection_Version_SelV3_for_HCP2	
      // 	selecTable.Fill(d,2,1.);	
      // 	if(isGoodPV){
      // 	  selecTable.Fill(d,3,1.);	  
      // 	  if(selectedMuons.size()==1){
      // 	    selecTable.Fill(d,4,1.);	    
      // 	    if(vetoMuons.size()==1){
      // 	      selecTable.Fill(d,5,1.);	      
      // 	      if(vetoElectrons.size()==0){
      // 		selecTable.Fill(d,6,1.);		
      // 		if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-3){
      // 		  selecTable.Fill(d,7,1.);
      // 		  if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-2){
      // 		    selecTable.Fill(d,8,1.);
      // 		    if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-1){
      // 		      selecTable.Fill(d,9,1.);
      // 		      if(selectedJets.size()>=(unsigned int)anaEnv.NofJets){   //we ask for 4 jets, but this way (with -3,-2,-1)it is possible to look at plots with other amount of jets
      // 			selecTable.Fill(d,10,1.);
      // 			eventSelected = true;
      // 			if(verbose>2) cout << event->runId() << ":" << event->eventId() << ":" << event->lumiBlockId()  << ":" << setprecision(8) << selectedMuons[0]->Pt() << endl;			
      // 		      }
      // 		    }
      // 		  }
      // 		}
      // 	      }
      // 	    }
      // 	  }
      // 	}
      // }      
      
      // if(eventSelected==true){
      // 	vector<int> jetCombi;
      // 	TRootMCParticle neutrino;
      // 	TRootMCParticle top;          
      // 	TRootMCParticle muon;
      // 	TRootMCParticle WLeptonicMc;
      // 	if(dataSetName.find("TTbarJets_SemiMu") == 0){
	  
      // 	  NumberOfAlternativeHits++;
	  
      // 	  //if(dataSetName.find("TTbarJets") == 0){
      // 	  pair<unsigned int, unsigned int> leptonicBJet_, hadronicBJet_, hadronicWJet1_, hadronicWJet2_; //First index is the JET number, second one is the parton
      // 	  leptonicBJet_ = hadronicBJet_ = hadronicWJet1_ = hadronicWJet2_ = pair<unsigned int, unsigned int>(9999,9999);
      // 	  vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV;
      // 	  vector<TRootMCParticle> mcParticlesMatching;
      // 	  bool muPlusFromTop = false, muMinusFromTop = false;
      // 	  for(unsigned int i=0; i<mcParticles.size(); i++){
      // 	    if( mcParticles[i]->status() != 3) continue;
	    
      // 	    if( mcParticles[i]->type() == 13 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -6 ){
      // 	      if(muMinusFromTop) cerr<<"muMinusFromTop was already true"<<endl;
      // 	      muMinusFromTop = true;
      // 	    }
      // 	    if( mcParticles[i]->type() == -13 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == 6 ){
      // 	      if(muPlusFromTop) cerr<<"muPlusFromTop was already true"<<endl;
      // 	      muPlusFromTop = true;
      // 	    }
	    
      // 	    if( abs(mcParticles[i]->type()) < 6 || abs(mcParticles[i]->type()) == 21 ){
      // 	      mcParticlesTLV.push_back(*mcParticles[i]);
      // 	      mcParticlesMatching.push_back(*mcParticles[i]);
      // 	    }
      // 	  }
      // 	  if(muPlusFromTop && muMinusFromTop)
      // 	    cerr<<"muPlusFromTop and muMinusFromTop are both true ?!\nCheck if you are using the right sample..."<<endl;
	  
      // 	  // take all the selectedJets_ to study the radiation stuff, selectedJets are already ordened in decreasing Pt()
      // 	  for(unsigned int i=0; i<selectedJets.size(); i++)
      // 	    selectedJetsTLV.push_back(*selectedJets[i]);
	  
      // 	  JetPartonMatching matching = JetPartonMatching(mcParticlesTLV, selectedJetsTLV, 2, true, true, 0.3);
	  
      // 	  if(matching.getNumberOfAvailableCombinations() != 1)
      // 	    cerr << "matching.getNumberOfAvailableCombinations() = "<<matching.getNumberOfAvailableCombinations()<<"  This should be equal to 1 !!!"<<endl;
	  
      // 	  vector< pair<unsigned int, unsigned int> > JetPartonPair, ISRJetPartonPair; //First one is jet number, second one is mcParticle number
	  
      // 	  for(unsigned int i=0; i<mcParticlesTLV.size(); i++){
      // 	    int matchedJetNumber = matching.getMatchForParton(i, 0);
      // 	    if(matchedJetNumber != -1)
      // 	      JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
      // 	  }
	  
      // 	  for(unsigned int i=0; i<JetPartonPair.size(); i++){
      // 	    unsigned int j = JetPartonPair[i].second;
	    
      // 	    //cout<<"JetPartonPair["<<i<<"]: "<<endl;
      // 	    //cout<<"Jet Pt: "<<selectedJets[JetPartonPair[i].first].Pt()<<"  Eta: "<<selectedJets[JetPartonPair[i].first].Eta()<<"  Phi: "<<selectedJets[JetPartonPair[i].first].Phi()<<endl;
      // 	    //cout<<"MCParticle Pt: "<<mcParticlesMatching[JetPartonPair[i].second].Pt()<<"  Eta: "<<mcParticlesMatching[JetPartonPair[i].second].Eta()<<"  Phi: "<<mcParticlesMatching[JetPartonPair[i].second].Phi()<<endl;
      // 	    //cout<<"\ttype: "<<mcParticlesMatching[JetPartonPair[i].second].type()<<"  motherType: "<<mcParticlesMatching[JetPartonPair[i].second].motherType()<<"  grannyType: "<<mcParticlesMatching[JetPartonPair[i].second].grannyType()<<endl;
      // 	    //cout<<"DR(MCParticle, jet) = "<<selectedJets[JetPartonPair[i].first].DeltaR(mcParticlesMatching[JetPartonPair[i].second])<<endl;
	    
      // 	    if( fabs(mcParticlesMatching[j].type()) < 6 ){
      // 	      if( ( muPlusFromTop && mcParticlesMatching[j].motherType() == -24 && mcParticlesMatching[j].grannyType() == -6 ) || ( muMinusFromTop && mcParticlesMatching[j].motherType() == 24 && mcParticlesMatching[j].grannyType() == 6 ) ){
      // 		if(hadronicWJet1_.first == 9999) 
      // 		  hadronicWJet1_ = JetPartonPair[i];
      // 		else if(hadronicWJet2_.first == 9999) 
      // 		  hadronicWJet2_ = JetPartonPair[i];
      // 		else cerr<<"Found a third jet coming from a W boson which comes from a top quark..."<<endl;
      // 	      }
      // 	    }
      // 	    if( fabs(mcParticlesMatching[j].type()) == 5 ){
      // 	      if( ( muPlusFromTop && mcParticlesMatching[j].motherType() == -6 ) || ( muMinusFromTop && mcParticlesMatching[j].motherType() == 6 ) )
      // 		hadronicBJet_ = JetPartonPair[i];
      // 	      else if( ( muPlusFromTop && mcParticlesMatching[j].motherType() == 6 ) || ( muMinusFromTop && mcParticlesMatching[j].motherType() == -6 ) )
      // 		leptonicBJet_ = JetPartonPair[i];
      // 	    }
      // 	  }
	  
      // 	  jetCombi.push_back(hadronicWJet1_.first);
      // 	  jetCombi.push_back(hadronicWJet2_.first);
      // 	  jetCombi.push_back(hadronicBJet_.first);
      // 	  jetCombi.push_back(leptonicBJet_.first);
	  
      // 	  CorrectQuark1=jetCombi[0];
      // 	  CorrectQuark2=jetCombi[1];
      // 	  CorrectBHadronic = jetCombi[2];
      // 	  CorrectBLeptonic = jetCombi[3];
	  
      // 	  //Working on generator level:  //Waarom niet gewoon werken met de TRootMCParticles ???
      // 	  if(jetCombi[0]!=9999 && jetCombi[1]!=9999 && jetCombi[2]!=9999 && jetCombi[3]!=9999){    
      // 	    CorrectRecMassW=(*selectedJets[jetCombi[0]]+*selectedJets[jetCombi[1]]).M();
      // 	    CorrectRecMassTop=(*selectedJets[jetCombi[0]]+*selectedJets[jetCombi[1]]+*selectedJets[jetCombi[2]]).M();
	    
      // 	    histo1D["WMass"]->Fill(CorrectRecMassW);
      // 	    histo1D["TopMass"]->Fill(CorrectRecMassTop);
      // 	  }
	  
      // 	  for(unsigned int i=0; i<mcParticles.size(); i++){
      // 	    //cout << "mcParticles[" << i << "]  pdgId = " << mcParticles[i].type() << "  mother pdgId = " << mcParticles[i].motherType() << "  granny pdgId = " << mcParticles[i].grannyType() << endl;
	    
      // 	    if(fabs(mcParticles[i]->type()) == 14 && fabs(mcParticles[i]->motherType()) == 24 && fabs(mcParticles[i]->grannyType()) == 6){
      // 	      neutrino = *mcParticles[i];
      // 	      neutrinoFound=1;
      // 	    }
      // 	    if(fabs(mcParticles[i]->type()) == 6){
      // 	      top = *mcParticles[i];
      // 	      topFound=1;
      // 	    }
      // 	    if(fabs(mcParticles[i]->type()) == 13 && fabs(mcParticles[i]->motherType()) == 24 && fabs(mcParticles[i]->grannyType()) == 6){
      // 	      muon = *mcParticles[i];
      // 	      muonFound=1;
      // 	    }
      // 	    WLeptonicMc = (neutrino + muon);
      // 	    if(CorrectBLeptonic!=9999 && neutrinoFound==1){
      // 	      //MSPlot["TopMassGen"]->Fill((neutrino+selectedJets[jetCombi[3]]+muon).M(), datasets[d], true, Luminosity);
      // 	    }

      // 	    TLorentzVector McMuonWZMF = muon;
      // 	    TLorentzVector McWLeptonicTZMF = WLeptonicMc;
	    
      // 	    McMuonWZMF.Boost(-WLeptonicMc.BoostVector());
      // 	    McWLeptonicTZMF.Boost(-top.BoostVector());
	    
      // 	    float McCosTheta;
      // 	    McCosTheta = ((McWLeptonicTZMF.Vect()).Dot(McMuonWZMF.Vect()))/(((McWLeptonicTZMF.Vect()).Mag())*((McMuonWZMF.Vect()).Mag()));
      // 	    if(CorrectBLeptonic!=9999 && neutrinoFound==1 && muonFound==1 && topFound == 1){
      // 	      histo1D["McCosTheta"]->Fill(McCosTheta);
      // 	    }

      // 	  } //if 0 < i < mcParticles.size()    
      // 	}//if dataset Semi mu ttbar
	
      // 	//calculating the reconstructed mass of the top and the W:
      // 	for(int i=0;i<3;i++){
      // 	  for(int j=i+1;j<4;j++){
      // 	    for(int k=0;k<4;k++){
      // 	      if(k!=i && k!=j){
      // 		recMassW[NumberCombinations] = (*selectedJets[i]+*selectedJets[j]).M();
      // 		recMassTop[NumberCombinations]=(*selectedJets[i]+*selectedJets[j]+*selectedJets[k]).M();
		
      // 		ChiSquared[NumberCombinations]=pow(((recMassW[NumberCombinations]-MassW)/SigmaW),2)+pow(((recMassTop[NumberCombinations]-MassTop)/SigmaTop),2);
      // 		QuarkOneIndex[NumberCombinations]=i;
      // 		QuarkTwoIndex[NumberCombinations]=j;
      // 		BHadronicIndex[NumberCombinations]=k;
		
      // 		NumberCombinations++;  //Is steeds 12 zoals verwacht
      // 	      }
      // 	    }
      // 	  }
      // 	}
	
      // 	//select lowest chi squared value:
      // 	ChiSquaredValue=ChiSquared[0];
      // 	for(int ii=0;ii<12;ii++){
      // 	  if(ChiSquaredValue>ChiSquared[ii]){
      // 	    ChiSquaredValue=ChiSquared[ii];
      // 	    UsedCombination=ii;        
      // 	  } 
      // 	}
      // 	MSPlot["ChiSquared"]->Fill(ChiSquaredValue, datasets[d], true, Luminosity);
	
      // 	// Jet not in Chisquared combination is the Leptonic B jet
      // 	for(int ll=0;ll<4;ll++){
      // 	  if(ll!=QuarkOneIndex[UsedCombination] && ll!=QuarkTwoIndex[UsedCombination] && ll!=BHadronicIndex[UsedCombination]){
      // 	    BLeptonicIndex=ll;
      // 	  }
      // 	}
	
      // 	//Only for Mc data, check if chosen combination is correct
      // 	if(dataSetName.find("TTbarJets_SemiMu") == 0){
      // 	  if((QuarkOneIndex[UsedCombination]==CorrectQuark1 && QuarkTwoIndex[UsedCombination]==CorrectQuark2 && BHadronicIndex[UsedCombination]==CorrectBHadronic) ||(QuarkOneIndex[UsedCombination]==CorrectQuark2 && QuarkTwoIndex[UsedCombination]==CorrectQuark1 && BHadronicIndex[UsedCombination]==CorrectBHadronic)){
      // 	    //MSPlot["ChiSquaredCorrect"]->Fill(ChiSquaredValue, datasets[d], true, Luminosity);
      // 	  }
      // 	  else{
      // 	    //MSPlot["ChiSquaredWrong"]->Fill(ChiSquaredValue, datasets[d], true, Luminosity);
      // 	  }
      // 	}
	
      // 	//cout << "flavour: " << selectedJets[QuarkOneIndex[UsedCombination]].partonFlavour() << "  topJet: " << selectedJets[QuarkOneIndex[UsedCombination]].isTopJet() << endl;
	
      // 	//Calculating MET_Pz() (equation ax² + bx + c = 0 ):
      // 	NeutrinoPx = -(*selectedMuons[0]+*selectedJets[0]+*selectedJets[1]+*selectedJets[2]+*selectedJets[3]).Px();
      // 	NeutrinoPy = -(*selectedMuons[0]+*selectedJets[0]+*selectedJets[1]+*selectedJets[2]+*selectedJets[3]).Py();
	
      // 	float aCoefficient = 4*pow(selectedMuons[0]->E(),2)-4*pow(selectedMuons[0]->Pz(),2);
      // 	float bCoefficient = 4*(selectedMuons[0]->Pz())*(pow(selectedMuons[0]->M(),2)-pow(MassW,2)-2*(selectedMuons[0]->Px())*NeutrinoPx-2*(selectedMuons[0]->Py())*NeutrinoPy);
      // 	float cCoefficient = -pow(selectedMuons[0]->M(),4)-pow(MassW,4)-4*pow(selectedMuons[0]->Px(),2)*pow(NeutrinoPx,2)-4*pow(selectedMuons[0]->Py(),2)*pow(NeutrinoPy,2)+4*(pow(selectedMuons[0]->M(),2)-pow(MassW,2))*((selectedMuons[0]->Px())*NeutrinoPx+(selectedMuons[0]->Py())*NeutrinoPy)-8*(selectedMuons[0]->Px())*NeutrinoPx*(selectedMuons[0]->Py())*NeutrinoPy+4*pow(selectedMuons[0]->E(),2)*pow(NeutrinoPx,2)+4*pow(selectedMuons[0]->E(),2)*pow(NeutrinoPy,2);
	
      // 	float DCoefficient = pow(bCoefficient,2)-4*aCoefficient*cCoefficient;
	
      // 	int NeutrinoFound =0;
      // 	if(DCoefficient>0){
      // 	  NeutrinoFound =1;
      // 	  NeutrinoPzOne = ((-bCoefficient + sqrt(DCoefficient))/(aCoefficient*2));
      // 	  NeutrinoPzTwo = ((-bCoefficient - sqrt(DCoefficient))/(aCoefficient*2));
	  
      // 	  //MSPlot["NeutrinoPzOne"]->Fill(NeutrinoPzOne, datasets[d], true, Luminosity);
      // 	  //MSPlot["NeutrinoPzTwo"]->Fill(NeutrinoPzTwo, datasets[d], true, Luminosity);
	  
      // 	  NeutrinoEOne = sqrt(pow(NeutrinoPx,2)+pow(NeutrinoPy,2)+pow(NeutrinoPzOne,2));
      // 	  NeutrinoETwo = sqrt(pow(NeutrinoPx,2)+pow(NeutrinoPy,2)+pow(NeutrinoPzTwo,2));
      // 	  TLorentzVector NeutrinoOne(NeutrinoPx,NeutrinoPy,NeutrinoPzOne,NeutrinoEOne);
      // 	  TLorentzVector NeutrinoTwo(NeutrinoPx,NeutrinoPy,NeutrinoPzTwo,NeutrinoETwo);
      // 	  //MSPlot["WMassNeutrinoOne"]->Fill((NeutrinoOne+selectedMuons[0]).M(), datasets[d], true, Luminosity);
      // 	  //MSPlot["WMassNeutrinoTwo"]->Fill((NeutrinoTwo+selectedMuons[0]).M(), datasets[d], true, Luminosity);
      // 	  //MSPlot["WMassNeutrinoMc"]->Fill((neutrino+selectedMuons[0]).M(), datasets[d], true, Luminosity);
	  
      // 	  //Selecting which neutrino solution is the most correct one (neutrino + muon + bottom should give top mass)
      // 	  //MSPlot["TopMassNeutrinoOne"]->Fill((NeutrinoOne+selectedMuons[0]+selectedJets[BLeptonicIndex]).M(), datasets[d], true, Luminosity);
      // 	  //MSPlot["TopMassNeutrinoTwo"]->Fill((NeutrinoTwo+selectedMuons[0]+selectedJets[BLeptonicIndex]).M(), datasets[d], true, Luminosity);
	  
      // 	  float TopMassDiffOne = fabs(MassTop - (NeutrinoOne+*selectedMuons[0]+*selectedJets[BLeptonicIndex]).M());
      // 	  float TopMassDiffTwo = fabs(MassTop - (NeutrinoTwo+*selectedMuons[0]+*selectedJets[BLeptonicIndex]).M());
      // 	  if(TopMassDiffOne<TopMassDiffTwo){
      // 	    Neutrino = NeutrinoOne;
      // 	    NeutrinoPz=Neutrino.Pz();
      // 	  }
      // 	  else{
      // 	    Neutrino = NeutrinoTwo;
      // 	    NeutrinoPz=Neutrino.Pz();
      // 	  }
      // 	  //MSPlot["NeutrinoPzSelected"]->Fill(NeutrinoPz, datasets[d], true, Luminosity);
	  
      // 	  //Compare obtained Pz neutrino value with MC value:
      // 	  if(neutrinoFound==1 && NeutrinoPz!=999){
      // 	    float neutrinoPzDiff = (neutrino.Pz()-NeutrinoPz)/NeutrinoPz;
      // 	    //MSPlot["NeutrinoPzCompared"]->Fill(neutrinoPzDiff, datasets[d], true, Luminosity);
      // 	  }
      // 	  MSPlot["TopMassNeutrinoSelected"]->Fill((Neutrino+*selectedMuons[0]+*selectedJets[BLeptonicIndex]).M(), datasets[d], true, Luminosity);
      // 	}//end of D>0 loop
	
      // 	WLeptonic = (Neutrino+*selectedMuons[0]);
      // 	TopLeptonic = (Neutrino+*selectedMuons[0]+*selectedJets[BLeptonicIndex]);
	
      // 	//Reboost the particles to rest frames 
      // 	TLorentzVector LeptonWZMF = *selectedMuons[0]; // In W Zero Mass Frame (WZMF)
      // 	TLorentzVector WParticleTZMF = WLeptonic;  // In Top Zero Mass Frame (TZMF)
      // 	TLorentzVector NeutrinoWZMF = Neutrino;
      // 	TLorentzVector TopTZMF = TopLeptonic;
	
      // 	LeptonWZMF.Boost(-WLeptonic.BoostVector());
      // 	NeutrinoWZMF.Boost(-WLeptonic.BoostVector());  //Only gives correct result if a Neutrino is found, thus NeutrinoFound==1 !
      // 	WParticleTZMF.Boost(-TopLeptonic.BoostVector());
      // 	TopTZMF.Boost(-TopLeptonic.BoostVector());
	
      // 	if(NeutrinoFound==1){
      // 	  //MSPlot["TopAfterBoost"]->Fill(TopTZMF.M(), datasets[d], true, Luminosity);
      // 	  //MSPlot["TopPzAfterBoost"]->Fill(TopTZMF.Pz(), datasets[d], true, Luminosity);
      // 	}
      // 	//std::cout << " values of top after boost :" << TopTZMF.Px() << " " << TopTZMF.Py() << " " << TopTZMF.Pz() << " " << TopTZMF.E() << " " << TopTZMF.M() << std::endl;
	
      // 	//Calculating cos:
	
      // 	// CosTheta = fabs(WTopBoost.Px()*LeptonWBoost.Px() + WTopBoost.Py()*LeptonWBoost.Py() + WTopBoost.Pz()*LeptonWBoost.Pz())/(sqrt(pow(WTopBoost.Px(),2) + pow(WTopBoost.Py(),2) + pow(WTopBoost.Pz(),2))*sqrt(pow(LeptonWBoost.Px(),2) + pow(LeptonWBoost.Py(),2) + pow(LeptonWBoost.Pz(),2)));
      // 	// CosThetaMcNeutrino=fabs(WTopBoostMcNeutrino.Px()*LeptonWBoostMcNeutrino.Px() + WTopBoostMcNeutrino.Py()*LeptonWBoostMcNeutrino.Py() + WTopBoostMcNeutrino.Pz()*LeptonWBoostMcNeutrino.Pz())/(sqrt(pow(WTopBoostMcNeutrino.Px(),2) + pow(WTopBoostMcNeutrino.Py(),2) + pow(WTopBoostMcNeutrino.Pz(),2))*sqrt(pow(LeptonWBoostMcNeutrino.Px(),2) + pow(LeptonWBoostMcNeutrino.Py(),2) + pow(LeptonWBoostMcNeutrino.Pz(),2)));
      // 	// CosThetaMcNeutrino=fabs(WTopBoostMcNeutrinoMcTop.Px()*LeptonWBoostMcNeutrino.Px() + WTopBoostMcNeutrinoMcTop.Py()*LeptonWBoostMcNeutrino.Py() + WTopBoostMcNeutrinoMcTop.Pz()*LeptonWBoostMcNeutrino.Pz())/(sqrt(pow(WTopBoostMcNeutrinoMcTop.Px(),2) + pow(WTopBoostMcNeutrinoMcTop.Py(),2) + pow(WTopBoostMcNeutrinoMcTop.Pz(),2))*sqrt(pow(LeptonWBoostMcNeutrino.Px(),2) + pow(LeptonWBoostMcNeutrino.Py(),2) + pow(LeptonWBoostMcNeutrino.Pz(),2)));
      // 	// MSPlot["CosTheta"]->Fill(CosTheta, datasets[d], true, Luminosity);
      // 	// MSPlot["CosThetaMcNeutrino"]->Fill(CosThetaMcNeutrino, datasets[d], true, Luminosity);
      // 	// MSPlot["CosThetaMcNeutrinoMcTop"]->Fill(CosThetaMcNeutrinoMcTop, datasets[d], true, Luminosity);
	
      //      } 
      //delete selection;
    }//loop on events

    //--------------------  Sigma for W Mass and Top Mass  --------------------
    // histo1D["WMass"]->Fit("gaus","Q"); 
    // float SigmaWFit=histo1D["WMass"]->GetFunction("gaus")->GetParameter(2);
    // histo1D["TopMass"]->Fit("gaus","Q");
    // float SigmaTopFit=histo1D["TopMass"]->GetFunction("gaus")->GetParameter(2);
    // std::cout << SigmaWFit << " " << SigmaTopFit << std::endl;

    //--------------------  Scaling of Cos Theta histograms  --------------------
    histo1D["StandardCosTheta"]->Scale(100./(histo1D["StandardCosTheta"]->Integral()));
    histo1D["StandardCosTheta"]->SetMinimum(0);
    histo1D["StandardCosTheta"]->SetMaximum(0.8);
    histo1D["StandardCosThetaFit"]->Scale(100./(histo1D["StandardCosThetaFit"]->Integral()));
    histo1D["StandardCosThetaFit"]->SetMinimum(0);
    histo1D["StandardCosThetaFit"]->SetMaximum(0.8);
    //    histo1D["McCosTheta"]->Scale(100./(histo1D["McCosTheta"]->Integral()));
    //    histo1D["McCosTheta"]->SetMinimum(0);
    //    histo1D["McCosTheta"]->SetMaximum(0.8);

    histo1D["HelicityDistributionRight"]->SetMinimum(0);
    histo1D["HelicityDistributionRight"]->SetMaximum(0.8);
    histo1D["HelicityDistributionLeft"]->SetMinimum(0);
    histo1D["HelicityDistributionLeft"]->SetMaximum(0.8);
    histo1D["HelicityDistributionLong"]->SetMinimum(0);
    histo1D["HelicityDistributionLong"]->SetMaximum(0.8);
    histo1D["HelicityDistribution"]->SetMinimum(0);
    histo1D["HelicityDistribution"]->SetMaximum(0.8);
  
    //--------------------  Fit to collect used helicity values   --------------------
    TF1 *helicityFit = new TF1("helicityFit","((([0]*3*(1+x)*(1+x))+([1]*3*(1-x)*(1-x))+([2]*6*(1-x*x)))/8)",-1,1);
    histo1D["StandardCosThetaFit"]->Fit("helicityFit","Q");
    std::cout << " fit values : " << helicityFit->GetParameter(0) << " " << helicityFit->GetParameter(1) << " " << helicityFit->GetParameter(2) << std::endl;

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
  
  string pathPNG = "PlotsMacro";
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
    int N = temp->GetNbinsX();
    temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
    temp->SetBinContent(N+1,0);
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
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

