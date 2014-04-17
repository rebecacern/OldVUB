
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
#include "../Reconstruction/interface/JetCorrectorParameters.h"
#include "../Reconstruction/interface/JetCorrectionUncertainty.h"

//Include WTree information
#include "../WHelicities/interface/WTree.h"

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
#include "../KinFitter/interface/TFitParticleEtThetaPhi.h"

using namespace std;
using namespace TopTree;

int main (int argc, char *argv[])
{
  clock_t start = clock();
    
  cout << "***********************************************" << endl;
  cout << " Beginning of the program for Cern WHelicity ! " << endl;
  cout << "***********************************************" << endl;
    
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
       
  //Configuration output format
  TTree *configTree = new TTree("configTree","configuration Tree");
  TClonesArray* tcdatasets = new TClonesArray("Dataset",1000);
  configTree->Branch("Datasets","TClonesArray",&tcdatasets);
  TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000);
  configTree->Branch("AnaEnv","TClonesArray",&tcAnaEnv);

  string rootFileName;

  /////////////////////////
  // Which decay channel //
  /////////////////////////
  
  bool semiElectron = false; // use semiElectron channel,
  bool semiMuon = true; // use semiMuon channel?
  if(semiElectron && semiMuon) cout << "  --> Using semiMuon and semiElectron channel..." << endl;
  else
  {
    if(semiMuon) cout << " --> Using the semiMuon channel..." << endl;
    else{
      cout << " --> Using the semiElectron channel..." << endl;
      rootFileName = "MacroOutputElectronChannel.root";
    }
  }

  /////////////////////
  //  Which trigger  //
  /////////////////////
  bool TriCentralJet30Trigger = false;
  bool IsoMu172024Trigger = false;
  if(TriCentralJet30Trigger == true) rootFileName = "MacroOutputTriCentralJet30Trigger.root";         //Root output file
  else if(IsoMu172024Trigger == true) rootFileName = "MacroOutputIsoMu172024Trigger.root";
  else if(TriCentralJet30Trigger == false && IsoMu172024Trigger == false) rootFileName = "MacroOutputNoTrigger.root";
  
  ////////////////////////
  //  Which systematics //
  ////////////////////////

  int doJESShift = 0; // 0: off 1: minus 2: plus
  cout << "doJESShift: " << doJESShift << endl;

  int doJERShift = 0; // 0: off (except nominal scalefactor for jer) 1: minus 2: plus
  cout << "doJERShift: " << doJERShift << endl;

  int dobTagEffShift = 0; //0: off (except nominal scalefactor for btag eff) 1: minus 2: plus
  cout << "dobTagEffShift: " << dobTagEffShift << endl;

  int domisTagEffShift = 0; //0: off (except nominal scalefactor for mistag eff) 1: minus 2: plus
  cout << "domisTagEffShift: " << domisTagEffShift << endl;


  //Output ROOT file
  const char *rootfile = rootFileName.c_str();

  ////////////////////////////////////
  /// AnalysisEnvironment  
  ////////////////////////////////////
    
  AnalysisEnvironment anaEnv;
  cout<<"Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  cout << " 1 " << endl;
  new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);
  cout << " 2 " << endl;
  int verbose = anaEnv.Verbose;
  cout << " 3 " << endl;
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
  //All histograms can be defined as pre-programmed maps which makes definitions and looping easier
  map<string,TH1F*> histo1D;     
  map<string,TH2F*> histo2D;  
    
  // Histograms needed to calculate the sigma and Mc mass (from mean value) for W and top mass distribution
  //   --> Comment out after initializing most recent values ( also lines 1046 and 1356 )
  histo1D["WMass"]= new TH1F("WMass","Hadr W Mass distribution", 200,0,160);
  histo1D["TopMass"]= new TH1F("TopMass","Hadr Top mass distribution", 200,0,350);

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

  histo1D["ElectronPt"]=new TH1F("ElectronPt","ElectronPt",200,0,250);
  histo1D["ElectronEta"]=new TH1F("ElectronEta","ElectronEta",200,-5,5);

  ////////////////////////////////////
  /// MultiSamplePlot
  ////////////////////////////////////
  //Whenever there are more than one dataset this can be used to add up the different results of the datasets
  //
  map<string,MultiSamplePlot*> MSPlot;   

  //MSPlot["AllMuonsRelPFIso"] = new MultiSamplePlot(datasets, "AllMuonsRelPFIso", 100, 0, 5, "Muons PF Iso");
  //MSPlot["AllElectronsRelPFIso"] = new MultiSamplePlot(datasets, "AllElectronsRelPFIso", 100, 0, 5, "Electrons PF Iso");
  if(semiMuon == true){
    MSPlot["SelectedEventsMuonsRelPFIso"] = new MultiSamplePlot(datasets, "SelectedEventsMuonsRelPFIso", 100, 0, 5, "Muons PF Iso");

  }
  if(semiElectron == true){
    MSPlot["SelectedEventsElectronsRelPFIso"] = new MultiSamplePlot(datasets, "SelectedEventsElectronsRelPFIso", 100, 0, 5, "Electrons PF Iso");

  }
  MSPlot["nEventsAfterCutsSemiMu"] = new MultiSamplePlot(datasets, "nEventsAfterCutsSemiMu", 23, -0.5, 22.5, "Nr. of events after each cut, SemiMu");
  MSPlot["nEventsAfterCutsSemiEl"] = new MultiSamplePlot(datasets, "nEventsAfterCutsSemiEl", 23, -0.5, 22.5, "Nr. of events after each cut, SemiEl");
  MSPlot["nPrimaryVertices"] = new MultiSamplePlot(datasets, "nPrimaryVertices", 12, -0.5, 11.5, "Nr. of primary vertices");
  MSPlot["nGoodPrimaryVertices"] = new MultiSamplePlot(datasets, "nGoodPrimaryVertices", 12, -0.5, 11.5, "Nr. of good primary vertices");
  MSPlot["nSelectedMuons"] = new MultiSamplePlot(datasets, "nSelectedMuons", 5, -0.5, 4.5, "Nr. of selected Muons");
  MSPlot["nSelectedElectrons"] = new MultiSamplePlot(datasets, "nSelectedElectrons", 5, -0.5, 4.5, "Nr. of selected Electrons");
  MSPlot["nLooseOtherMuons"] = new MultiSamplePlot(datasets, "nLooseOtherMuons", 5, -0.5, 4.5, "Nr. of Loose Other Muons");
  MSPlot["nLooseElectronsSemiMu"] = new MultiSamplePlot(datasets, "nLooseElectronsSemiMu", 5, -0.5, 4.5, "Nr. of Loose Electrons, SemiMu");
  MSPlot["nLooseElectronsSemiEl"] = new MultiSamplePlot(datasets, "nLooseElectronsSemiEl", 5, -0.5, 4.5, "Nr. of Loose Electrons, SemiEl");
  MSPlot["nSelectedJets"] = new MultiSamplePlot(datasets, "nSelectedJets", 10, -0.5, 9.5, "Nr. of Selected Jets");  

  ///////////////////////////////////////////////////////
  //   Making the theoretical helicity distribution    // -> To know how the distribution should look like!! 
  ///////////////////////////////////////////////////////  --> Only Fit values are used in code!!
  //------   Standard Model   -----
  float LongitudinalFraction = 0.64491; // values obtained from fit (December 2011 -- RunA2011 result)
  float RightHandedFraction = 0.0334506;
  float LeftHandedFraction = 1 - LongitudinalFraction - RightHandedFraction;
  float SMDiffDistribution[3]; //0:Longitudinal; 1: Righthanded; 2: Lefthanded
  float XVariable;
  std::cout << " Used SM helicity variables: Longitudinal (f0) = " << LongitudinalFraction << " Righthanded (f+) = " << RightHandedFraction << " Lefthanded (f-) = " << LeftHandedFraction << std::endl;    

  /*for(unsigned int jj=0;jj<200;jj++){
    XVariable=-1+0.01*jj;
      
    SMDiffDistribution[0] = (LongitudinalFraction*6*(1-pow(XVariable,2)))/8;
    SMDiffDistribution[2] = (pow((1-XVariable),2)*3*LeftHandedFraction)/8;
    SMDiffDistribution[1] = (RightHandedFraction*3*pow((1+XVariable),2))/8;
     
    histo1D["HelicityDistribution"]->SetBinContent(jj+1,(SMDiffDistribution[0]+SMDiffDistribution[1]+SMDiffDistribution[2]));
    histo1D["HelicityDistributionLong"]->SetBinContent(jj+1,SMDiffDistribution[0]);
    histo1D["HelicityDistributionLeft"]->SetBinContent(jj+1,SMDiffDistribution[2]);
    histo1D["HelicityDistributionRight"]->SetBinContent(jj+1,SMDiffDistribution[1]);
    }*/

  /////////////////////////////
  /// ResolutionFit Stuff
  /////////////////////////////

  bool CalculateResolutions = false; // If false, the resolutions will be loaded from a previous calculation

  std::cout << " CalculateResolutions = " << CalculateResolutions << endl;

  ResolutionFit *resFitLightJets = 0, *resFitBJets = 0, *resFitMuon = 0, *resFitElectron = 0, *resFitNeutrino = 0;
    
  resFitLightJets = new ResolutionFit("LightJet");
  resFitBJets = new ResolutionFit("BJet");
  resFitMuon = new ResolutionFit("Muon");
  resFitElectron = new ResolutionFit("Electron");
  resFitNeutrino = new ResolutionFit("Neutrino");

  if(!CalculateResolutions){
    resFitLightJets->LoadResolutions("resolutions/lightJetReso.root");
    resFitBJets->LoadResolutions("resolutions/bJetReso.root");
    if(semiMuon == true){
      resFitMuon->LoadResolutions("resolutions/muonReso.root");
      resFitNeutrino->LoadResolutions("resolutions/neutrinoSemiMuReso.root");  //Once resolutions are newly created they will be split up for SemiMu and SemiEl for Neutrino !!
    }
    else if(semiElectron == true){
      resFitNeutrino->LoadResolutions("resolutions/neutrinoSemiElReso.root");
      resFitElectron->LoadResolutions("resolutions/electronReso.root");
    }
    std::cout << " Resolutions loaded " << std::endl;
  }

  if (verbose > 0)
    cout << " - ResolutionFit instantiated ..." << endl;

  float TopMassKinFit;
  float WMassKinFit = 80.4;

  ////////////////////////////////////
  /// Selection table
  ////////////////////////////////////
    
  vector<string> CutsSelecTableSemiMu;//Writes information in a pdf file, is filled during the selection with the different layers
  CutsSelecTableSemiMu.push_back(string("initial"));
  CutsSelecTableSemiMu.push_back(string("preselected"));
  CutsSelecTableSemiMu.push_back(string("trigged"));
  CutsSelecTableSemiMu.push_back(string("Good PV"));
  CutsSelecTableSemiMu.push_back(string("1 selected muon"));
  CutsSelecTableSemiMu.push_back(string("Veto 2nd muon"));
  CutsSelecTableSemiMu.push_back(string("Veto electron"));

  vector<string> CutsSelecTableSemiEl;
  CutsSelecTableSemiEl.push_back(string("preselected"));
  CutsSelecTableSemiEl.push_back(string("trigged"));
  CutsSelecTableSemiEl.push_back(string("Good PV"));
  CutsSelecTableSemiEl.push_back(string("1 selected electron"));
  CutsSelecTableSemiEl.push_back(string("Veto muon"));
  CutsSelecTableSemiEl.push_back(string("Veto 2nd electron from Z-decay"));
  CutsSelecTableSemiEl.push_back(string("Conversion veto"));
  

  char LabelNJets[100];
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-3);
  CutsSelecTableSemiMu.push_back(string(LabelNJets)); 
  CutsSelecTableSemiEl.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-2);
  CutsSelecTableSemiMu.push_back(string(LabelNJets));
  CutsSelecTableSemiEl.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-1);
  CutsSelecTableSemiMu.push_back(string(LabelNJets));
  CutsSelecTableSemiEl.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets);
  CutsSelecTableSemiMu.push_back(string(LabelNJets));
  CutsSelecTableSemiEl.push_back(string(LabelNJets));
  if (verbose > 0)
    cout << " - CutsSelectionTable instantiated ..." << endl;
  SelectionTable selecTableSemiMu(CutsSelecTableSemiMu, datasets);
  selecTableSemiMu.SetLuminosity(Luminosity);
  SelectionTable selecTableSemiEl(CutsSelecTableSemiEl, datasets);
  selecTableSemiEl.SetLuminosity(Luminosity);
  if (verbose > 0)
    cout << " - SelectionTable instantiated ..." << endl; 
    
  ////////////////////////////////////
  //	Loop on datasets
  ////////////////////////////////////
   
  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
  for (unsigned int d = 0; d < datasets.size (); d++){
      
    std::cout << " ********************** " << std::endl;
    std::cout << " value of data set : " << std::endl;
    std::cout << "    " << d << std::endl;
    std::cout << " ********************** " << std::endl;

    if (verbose > 1)
      cout << "   Dataset " << d << ": " << datasets[d]->Name () << "/ title : " << datasets[d]->Title () << endl;
    if (verbose > 1)
      std::cout<<"      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
      
    //open files and load
    cout<<"LoadEvent"<<endl;
    treeLoader.LoadDataset (datasets[d], anaEnv);
    cout<<"LoadEvent"<<endl;
      
    string dataSetName = datasets[d]->Name();
            
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
      
    // Create the JetCorrectorParameter objects, the order does not matter.
    // YYYY is the first part of the txt files: usually the global tag from which they are retrieved
    JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters("JECFiles/START42_V17_AK5PFchs_L3Absolute.txt");
    JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters("JECFiles/START42_V17_AK5PFchs_L2Relative.txt");
    JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters("JECFiles/START42_V17_AK5PFchs_L1FastJet.txt");

    //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
    vCorrParam.push_back(*L1JetPar);
    vCorrParam.push_back(*L2JetPar);
    vCorrParam.push_back(*L3JetPar);

    if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") // Data!
      {
	JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("JECFiles/START42_V17_AK5PFchs_L2L3Residual.txt");
	vCorrParam.push_back(*ResJetCorPar);
      }

    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/START42_V17_AK5PFchs_Uncertainty.txt");
    
    // true means redo also the L1      
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);
      
    ////////////////////////////////////////////////////////////
    // CREATE OUTPUT FILE AND TTREE FOR STORAGE OF THE NTUPLE //
    ////////////////////////////////////////////////////////////

    string decayChannel;
    if(semiElectron && semiMuon) decayChannel = "SemiLep";
    else{
      if(semiMuon) decayChannel = "SemiMu";
      if(semiElectron) decayChannel = "SemiEl";
    }

    string UsedTrigger;
    if(TriCentralJet30Trigger == true) UsedTrigger = "TriCentralJet30Trigger";
    else if(IsoMu172024Trigger == true) UsedTrigger = "IsoMu172024Trigger";
    else if(TriCentralJet30Trigger == false && IsoMu172024Trigger == false) UsedTrigger = "NoTrigger";

    //------------------------------------
    // Files for Nominal & JES up/down
    //------------------------------------
    string wTreeFileTitle;

    if(doJESShift == 0) wTreeFileTitle = "WTree/KinFit_WTree_"+UsedTrigger+"_"+dataSetName+"_"+decayChannel+".root";
    if(doJESShift == 1) wTreeFileTitle = "WTree/KinFit_WTree_"+UsedTrigger+"_JESMinus_1Sig_"+dataSetName+"_"+decayChannel+".root";  //JES systematics
    if(doJESShift == 2) wTreeFileTitle = "WTree/KinFit_WTree_"+UsedTrigger+"_JESPlus_1Sig_"+dataSetName+"_"+decayChannel+".root";  //JES systematics
        
    cout << "INFO: creating WTree file "+wTreeFileTitle << endl;
        
    TFile* WTreeFile = new TFile(wTreeFileTitle.c_str(),"RECREATE");
      
    WTree* wTree = 0;
    TTree* WTreeTree = new TTree("WTreeTree","Tree containing the W information");
    WTreeTree->Branch("TheWTree","WTree",&wTree);
     
    ////////////////////////////////////
    //	Loop on events
    ////////////////////////////////////
      
    nEvents[d] = 0;
    int itriggerSemiMu = -1, itriggerSemiEl = -1, previousRun = -1;
    if (verbose > 1)
      cout << "	Loop over events " << endl;
    
    for(unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++){     //In this loop plots before selection can be defined
    //for(unsigned int ievt = 0; ievt < 10000; ievt++){  

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

      //load mcParticles --> Put here to avoid overwriting of init_jets and mets after second treeLoader.Load !!
      vector<TRootMCParticle*> mcParticles;
      if(dataSetName.find("TTbarJets_SemiMu") == 0 || dataSetName.find("TTbarJets_SemiEl") == 0){
	mcParticles = treeLoader.LoadMCPart(ievt);
	sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      }         
      
      //Only keep init_jets with pt larger than 20 GeV (MC preselected at 20 and Data at 10 GeV --> Need to be consistent)
      for(int ii=0;ii<init_jets.size();ii++){
	if(init_jets[ii]->Pt()<20){
	  init_jets.erase(init_jets.begin()+ii);
	  ii=ii-1;
	}
      }
        
      // scale factor for the event
      float scaleFactor = 1.;
      float lumiWeight = 1;
      if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 )
	lumiWeight = 1;
      
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

	if(semiMuon==true){
	  if(TriCentralJet30Trigger == true){
	    //First bunch of triggers: TriCentralJet with (Iso)Mu17 resulting in an offline muon pt cut of 20 
	    //
	    if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA"){
	      if( event->runId() <= 161176 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_Mu17_TriCentralJet30_v1"), currentRun, iFile);
	      else if( event->runId() >= 161217 && event->runId() <= 163261 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_Mu17_TriCentralJet30_v2"), currentRun, iFile);
	      else if( event->runId() >= 163270 && event->runId() <= 163869 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_Mu17_TriCentralJet30_v4"), currentRun, iFile);
	      else if( event->runId() >= 165088 && event->runId() <= 165633 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_Mu17_TriCentralJet30_v5"), currentRun, iFile);
	      else if( event->runId() == 166346 )
		itriggerSemiMu= treeLoader.iTrigger (string ("HLT_IsoMu17_TriCentralJet30_v2"), currentRun, iFile);
	      else if( event->runId() >= 165970 && event->runId() <= 167043 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu17_TriCentralJet30_v1"), currentRun, iFile);
	      else if( event->runId() >= 167078 && event->runId() <= 167913 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu17_TriCentralJet30_v3"), currentRun, iFile);
	      else if( event->runId() >= 170826 && event->runId() <= 173198 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu17_TriCentralJet30_v5"), currentRun, iFile);
	      else if( event->runId() >= 173236 && event->runId() <= 178380 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralJet30_v1"), currentRun, iFile);
	      else if( event->runId() >= 178381 && event->runId() <= 179889 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFJet30_v2"), currentRun, iFile);
	      else if( event->runId() >= 179959 && event->runId() <= 180252 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFJet30_v3"), currentRun, iFile);
	      else
		cout << "Unknown run for HLTpath selection: " << event->runId() << endl;     
	    }
	    else{
	      itriggerSemiMu = treeLoader.iTrigger (string ("HLT_Mu17_TriCentralJet30_v2"), currentRun);
	    }
	    //End of first bunch of triggers!!
	  }
	  
	  else if(IsoMu172024Trigger == true){
	    
	    //Second bunch of triggers: IsoMu (17/20/24) resulting in an offline muon pt cut of minimum 27
	    //
	    if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") {
	      
	      // the first 1.1/fb part of RUN2011A (may10->promptv4)
	      if( event->runId() <= 163261 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu17_v5"), currentRun, iFile);
	      else if( event->runId() >= 163270 && event->runId() <= 163869 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu17_v6"), currentRun, iFile);
	      else if( event->runId() >= 165088 && event->runId() <= 165633 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu17_v8"), currentRun, iFile);
	      else if( event->runId() == 166346 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu17_v10"), currentRun, iFile);
	      else if( event->runId() >= 165970 && event->runId() <= 167043 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu17_v9"), currentRun, iFile);
	      else if( event->runId() >= 167078 && event->runId() <= 167913 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu17_v11"), currentRun, iFile);
	      
	      // the other part of RUN2011A (another 1/fb) (aug05,promptv6)
	      
	      else if( event->runId() >= 170249 && event->runId() <= 172619 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu20_v8"), currentRun, iFile);
	      else if( event->runId() >= 172620 && event->runId() <= 173198 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu20_v8"), currentRun, iFile);
	      else if (event->runId() >= 173236 && event->runId() <= 173692)
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_v9"), currentRun, iFile);	

	      // RUN2011B (promptv1)
	      
	      else if( event->runId() ==  176928 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v3"), currentRun, iFile);
	      else if( event->runId() == 176982 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v3"), currentRun, iFile);
	      
	      else if( event->runId() >= 175860 && event->runId() <= 176469 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v3"), currentRun, iFile);
	      else if( event->runId() >=  176548 && event->runId() <=  176702 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v3"), currentRun, iFile);
	      else if( event->runId() >=  176797 && event->runId() <=  176889 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v3"), currentRun, iFile);
	      else if( event->runId() >=  176929 && event->runId() <=  176959 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v3"), currentRun, iFile);
	      else if( event->runId() >=  177053 && event->runId() <=  177452 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v3"), currentRun, iFile);
	      
	      else if( event->runId() >=  176545 && event->runId() <=  176547 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v3"), currentRun, iFile);
	      else if( event->runId() >=  176765 && event->runId() <=  176796 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v3"), currentRun, iFile);
	      
	      else if( event->runId() >=  177718 && event->runId() <=  178380 ) // TopTree ID 804
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v3"), currentRun, iFile);
	      else if( event->runId() >=  178420 && event->runId() <=  178479 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v6"), currentRun, iFile);                                
	      else if( event->runId() >=  178703 && event->runId() <=  179889 ) // TopTree ID 816
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v6"), currentRun, iFile);
	      else if( event->runId() >=  179959 && event->runId() <=  180252 )
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v7"), currentRun, iFile);
	      	      
	      else
		cout << "Unknown run for HLTpath selection: " << event->runId() << endl;
	      
	      //  HLT_IsoMu17_v5        160431-163261                          32.88(/pb)    32.88(/pb)
	      //  HLT_IsoMu17_v6        163270-163869                         163.81(/pb)   163.81(/pb)
	      //  HLT_IsoMu17_v8        165088-165633                         137.51(/pb)   137.46(/pb)
	      //  HLT_IsoMu17_v9        165970-167043, except 166346          530.06(/pb)   529.50(/pb)
	      //  HLT_IsoMu17_v10       166346                                  4.26(/pb)     4.26(/pb)
	      //  HLT_IsoMu17_v11       167078-167151                          22.07(/pb)    22.00(/pb)
	      //  HLT_IsoMu17_v11       167281-167913                         222.84(/pb)   222.84(/pb)
	    } 
	    else {
	      itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu17_v5"), currentRun);
	      if (itriggerSemiMu == 9999)
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu17_v4"), currentRun); // Spring11: HLT_Mu15_v1
	      if (itriggerSemiMu == 9999)
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu9"), currentRun); // Fall10: HLT_Mu9
	      if (itriggerSemiMu == 9999)
		itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu17_v14"), currentRun); // Fall11

	      //  Summer11 MC:        HLT_IsoMu17_v5        HLT_Mu15_v2
	      //  Spring11 MC:        HLT_IsoMu17_v4        HLT_Mu15_v1 or HLT_Mu17_v1
	      //  Fall10 MC:          HLT_IsoMu9            HLT_Mu9 or HLT_Mu11
	    }
	    //End of second bunch of triggers!!
	  }
	  
	}  //End of triggers for semiMu case 
	
	else if(semiElectron == true){
	  if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA"){
	    
	    if( event->runId() <= 161176 )
	      itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30_v1"), currentRun, iFile);
	    else if( event->runId() >= 161177 && event->runId() <= 163261 )
	      itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30_v2"), currentRun, iFile);
	    else if( event->runId() >= 163262 && event->runId() <= 163869 )
	      itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30_v3"), currentRun, iFile);
	    else if( event->runId() >= 163870 && event->runId() <= 165633 )
	      itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_TrkIdT_TriCentralJet30_v3"), currentRun, iFile);
	    else if( event->runId() >= 165970 && event->runId() <= 166967 )
	      itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30_v1"), currentRun, iFile);
	    else if( event->runId() >= 167039 && event->runId() <= 167913 )
	      itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30_v2"), currentRun, iFile);
	    else if( event->runId() >= 170826 && event->runId() <= 173198 )
	      itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30_v4"), currentRun, iFile);
	    else if( event->runId() >= 173236 && event->runId() <= 178380 )
	      itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30_v5"), currentRun, iFile);
	    else if( event->runId() >= 178381 && event->runId() <= 178479 )
	      itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v2"), currentRun, iFile);
	    else
	      cout << "Unknown run for SemiEl HLTpath selection: " << event->runId() << endl;
	    if( itriggerSemiEl == 9999 ){
	      cout << "itriggerSemiEl == 9999 for SemiEl HLTpath selection: " << event->runId() << endl;
	      exit(-1);
	    }
	  }
	  else{
	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30_v2"), currentRun);
	  }
	}//End of triggers for semiEl case
	else{
	  cout << " ------------             ----------------                ------------------- " << endl;
	  cout << "                   Both electron channel and muon channel active !! " << endl;
 	  cout << " ------------             ----------------                ------------------- " << endl;
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
            
      //JES CORRECTION
      // Apply Jet Corrections on-the-fly
      if( dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
	jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho(),true); //last boolean: isData (needed for L2L3Residual...)
      else 
	jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho(),false); //last boolean: isData (needed for L2L3Residual...)
      
      //ordering is relevant; most probably 1) Type I MET correction, 2) JER where jet corrections are propagated to MET, 3) JES systematics where jet corrections are propagated to MET
      //----------------------------------------------------------
      // Apply type I MET corrections:  (Only for |eta| <= 4.7 )
      //---------------------------------------------------------
      if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
        jetTools->correctMETTypeOne(init_jets,mets[0],true);
      else
        jetTools->correctMETTypeOne(init_jets,mets[0],false);
      
      if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
      {	
	if(doJERShift == 1)
	  jetTools->correctJetJER(init_jets, genjets, mets[0], "minus",false);   //false means don't use old numbers but newer ones...
	else if(doJERShift == 2)
	  jetTools->correctJetJER(init_jets, genjets, mets[0], "plus",false);
	else
	  jetTools->correctJetJER(init_jets, genjets, mets[0], "nominal",false);
	
	// JES systematic! 
	if (doJESShift == 1)
	  jetTools->correctJetJESUnc(init_jets, mets[0], "minus");
	else if (doJESShift == 2)
	  jetTools->correctJetJESUnc(init_jets, mets[0], "plus");	       
      }      
              
      ///////////////////////////////////////////////////////
      ///     Start of program: Defining variables        ///
      ///////////////////////////////////////////////////////  
    	           
      vector<float> ChiSquared[2];  //Needed for chi squared caclulation      
      vector<float> ChiSquaredFull[2];
      
      //Vectors which contain the particles with kinematics changed due to KinFit:
      vector<TLorentzVector> fittedLepton[2];
      vector<TLorentzVector> fittedNeutrino[2];
      vector<TLorentzVector> fittedBLept[2];
      vector<TLorentzVector> fittedBHadr[2];
      vector<TLorentzVector> fittedLight1[2];
      vector<TLorentzVector> fittedLight2[2];
      
      vector<TLorentzVector> fittedFullLepton[2];
      vector<TLorentzVector> fittedFullNeutrino[2];
      vector<TLorentzVector> fittedFullBLept[2];
      vector<TLorentzVector> fittedFullBHadr[2];
      vector<TLorentzVector> fittedFullLight1[2];
      vector<TLorentzVector> fittedFullLight2[2];

      float MassW= 83.3924;
      float MassTop = 172.452;
      float MassTopLept = 180.349;
           
      float standardCosTheta=0; 
      TRootMCParticle standardNeutrino, standardTop,standardLepton,standardWLeptonic;      

      if((dataSetName.find("TTbarJets_SemiMu") == 0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") == 0 && semiElectron == true)){     

	//------------------------------------------//
	//    Identifying Monte Carlo particles     //
	//------------------------------------------//
	int EventParticleNumber[5]; //0:top; 1:b; 2: u,c,d,s; 3:W; 4:mu + neutrino
	for(int ll = 0;ll<5;ll++){EventParticleNumber[ll]=0;}
	int EventChargeTop=0;  //1 for top, -1 for anti-top
	int EventChargeAntiTop=0;
	int EventChargeWPos=0; //1 for W+, -1 for W-
	int EventChargeWNeg=0;
	int EventChargeLepPos=0; //1 for mu/el+, -1 for mu/el-
	int EventChargeLepNeg =0;
	TRootMCParticle Top,AntiTop,WPos,WNeg;
	TLorentzVector standardLeptonWZMF, standardWLeptonicTZMF;

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
	  } 

	  if(fabs(mcParticles[i]->type()) <= 4 && fabs(mcParticles[i]->motherType()) == 24 && fabs(mcParticles[i]->grannyType()) == 6){
	    EventParticleNumber[2]++;
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

	  if(semiMuon == true){
	    if(fabs(mcParticles[i]->type()) == 14 && fabs(mcParticles[i]->motherType()) == 24 && fabs(mcParticles[i]->grannyType()) == 6){
	      EventParticleNumber[4]++;  //Identifying neutrino's
	      standardNeutrino=*mcParticles[i];
	    }
	    
	    if(fabs(mcParticles[i]->type()) == 13 && fabs(mcParticles[i]->motherType()) == 24 && fabs(mcParticles[i]->grannyType()) == 6 && mcParticles[i]->status()==3){ //status: 1:stable; 2:shower; 3:hard scattering(coming from the studied hard proces)
	      EventParticleNumber[4]++;
	      standardLepton=*mcParticles[i];
	      if(mcParticles[i]->type()==13){EventChargeLepNeg=-1;}
	      else if(mcParticles[i]->type()==-13){EventChargeLepPos=1;}
	    }	 
	  }
	  else if(semiElectron == true){
	    if(fabs(mcParticles[i]->type()) == 12 && fabs(mcParticles[i]->motherType()) == 24 && fabs(mcParticles[i]->grannyType()) == 6){
	      EventParticleNumber[4]++;  //Identifying neutrino's
	      standardNeutrino=*mcParticles[i];
	    }

	    if(fabs(mcParticles[i]->type()) == 11 && fabs(mcParticles[i]->motherType()) == 24 && fabs(mcParticles[i]->grannyType()) == 6 && mcParticles[i]->status()==3){ //status: 1:stable; 2:shower; 3:hard scattering(coming from the studied hard proces)
	      EventParticleNumber[4]++;
	      standardLepton=*mcParticles[i];
	      if(mcParticles[i]->type()==11){EventChargeLepNeg=-1;}
	      else if(mcParticles[i]->type()==-11){EventChargeLepPos=1;}
	    }	  	  
	  }
	  
	}//  if 0 < i < mcParticles.size()
	
	//////////////////////////////////////////////////
	//   Selecting correct event (b b q q mu nu )   //
	//////////////////////////////////////////////////
	if(EventParticleNumber[0]==2 && EventParticleNumber[1]==2 && EventParticleNumber[2]==2 && EventParticleNumber[3]==2 && EventParticleNumber[4]==2){
	    
	  //-----   Differentiating between proces from top and anti-top (choose leptonic):   -----
	  if(EventChargeTop==1 && EventChargeWPos==1 && EventChargeLepPos==1){  // Proces: t -> b W+ -> mu/el+ nu
	    standardWLeptonicTZMF=WPos;
	    standardWLeptonic=WPos;
	    standardTop=Top;
	  }
	  else if(EventChargeAntiTop==-1 && EventChargeWNeg==-1 && EventChargeLepNeg==-1){ //proces: anti-t -> anti-b W- -> mu/el- anti-nu	
	    standardWLeptonicTZMF=WNeg;
	    standardWLeptonic=WNeg;
	    standardTop=AntiTop;
	  }

	  //-----   Applying boost on muon and W:   -----
	  standardLeptonWZMF=standardLepton;
	  standardLeptonWZMF.Boost(-standardWLeptonicTZMF.BoostVector());
	  standardWLeptonicTZMF.Boost(-standardTop.BoostVector());

	  //-----   Calculating cos theta:   -----
	  standardCosTheta = ((standardWLeptonicTZMF.Vect()).Dot(standardLeptonWZMF.Vect()))/(((standardWLeptonicTZMF.Vect()).Mag())*((standardLeptonWZMF.Vect()).Mag()));

	  histo1D["StandardCosTheta"]->Fill(standardCosTheta);  // Histogram without fit
	  histo1D["StandardCosThetaFit"]->Fill(standardCosTheta);  // Histogram with fit   	  

	}
      }//if dataset Semi mu ttbar  

      /////////////////////////////
      //   Selection
      /////////////////////////////        
      //Declare selection instance    
      Selection selection(init_jets, init_muons, init_electrons, mets);
      selection.setJetCuts(30.,2.4,0.01,1.,0.98,0.3,0.1);   //CIEMAT values, not refSel values !!!!
      if(TriCentralJet30Trigger == true) selection.setMuonCuts(20,2.1,0.15,10,0.02,0.3,1,1,1); //Values for TriCentralJet trigger
      if(IsoMu172024Trigger == true) selection.setMuonCuts(35,2.1,0.15,10,0.02,0.3,1,1,1); //Values for IsoMu(17/20/24) trigger -- Should be 27, but put on 25 to match CIEMAT constraints
      if(TriCentralJet30Trigger == false && IsoMu172024Trigger == false) selection.setMuonCuts(10,2.1,0.15,10,0.02,0.3,1,1,1);
      selection.setLooseMuonCuts(10,2.1,0.15);
      selection.setElectronCuts(30,2.5,0.15,0.02,1,0.3);
      selection.setLooseElectronCuts(15,2.5,0.2); // semiMu looseMuon cuts
        
      bool triggedSemiMu = false;
      bool triggedSemiEl = false;
      if(TriCentralJet30Trigger != false && IsoMu172024Trigger != false){
	if(semiMuon == true){ triggedSemiMu = treeLoader.EventTrigged (itriggerSemiMu);}
	else if(semiElectron == true){triggedSemiEl = treeLoader.EventTrigged (itriggerSemiEl);}
      }
      else if(TriCentralJet30Trigger == false && IsoMu172024Trigger == false){
	//Look at semi-mu sample without trigger influence:
	if(semiMuon == true) triggedSemiMu = true;
	else if(semiElectron == true) triggedSemiEl = true;
      }

      bool isGoodPV = selection.isPVSelected(vertex, anaEnv.PVertexNdofCut, anaEnv.PVertexZCut, anaEnv.PVertexRhoCut);

      vector<TRootElectron*> vetoElectronsSemiMu = selection.GetSelectedLooseElectrons(false);
      vector<TRootElectron*> vetoElectronsSemiEl = selection.GetSelectedLooseElectrons(20,2.5,1.0,true);      

      vector<TRootJet*> selectedJets = selection.GetSelectedJets(true);
      vector<TRootMuon*> selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
      vector<TRootElectron*> selectedElectrons = selection.GetSelectedElectrons(vertex[0],selectedJets);
      vector<TRootMuon*> vetoMuons = selection.GetSelectedLooseMuons();
        
      //      vector<TRootMuon*> looseMuonCIEMAT = selection.GetSelectedDiMuons();
      //      vector<TRootMuon*> tightMuonCIEMAT = selection.GetSelectedDiMuonsTight();
      //      vector<TRootElectron*> goodElectronCIEMAT = selection.GetSelectedLooseElectrons();   
        
      /////////////////////////
      //   Event selection   //
      /////////////////////////
 
      bool eventSelectedSemiMu = false;
      bool eventSelectedSemiEl = false;
      
      //selecTable.Fill(d,0, datasets[d]->Xsection() * datasets[d]->EquivalentLumi() );
      
      MSPlot["nPrimaryVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nGoodPrimaryVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nSelectedMuons"]->Fill(selectedMuons.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nSelectedElectrons"]->Fill(selectedElectrons.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nLooseOtherMuons"]->Fill(vetoMuons.size()-1, datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nLooseElectronsSemiMu"]->Fill(vetoElectronsSemiMu.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nLooseElectronsSemiEl"]->Fill(vetoElectronsSemiEl.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nSelectedJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
      
      MSPlot["nEventsAfterCutsSemiMu"]->Fill(0, datasets[d], true, Luminosity*scaleFactor);
      selecTableSemiMu.Fill(d,0,scaleFactor*lumiWeight);
      if( triggedSemiMu && semiMuon )
	{
        MSPlot["nEventsAfterCutsSemiMu"]->Fill(1, datasets[d], true, Luminosity*scaleFactor);
        selecTableSemiMu.Fill(d,1,scaleFactor*lumiWeight);
        if( isGoodPV )
	  {
	    MSPlot["nEventsAfterCutsSemiMu"]->Fill(2, datasets[d], true, Luminosity*scaleFactor);
	    selecTableSemiMu.Fill(d,2,scaleFactor*lumiWeight);
	    if( selectedMuons.size() == 1 )
	      {
		MSPlot["nEventsAfterCutsSemiMu"]->Fill(3, datasets[d], true, Luminosity*scaleFactor);
		selecTableSemiMu.Fill(d,3,scaleFactor*lumiWeight);
		if( vetoMuons.size() == 1 ) 
		  {
		    MSPlot["nEventsAfterCutsSemiMu"]->Fill(4, datasets[d], true, Luminosity*scaleFactor);
		    selecTableSemiMu.Fill(d,4,scaleFactor*lumiWeight);
		    if( vetoElectronsSemiMu.size() == 0 )
		      {
			MSPlot["nEventsAfterCutsSemiMu"]->Fill(5, datasets[d], true, Luminosity*scaleFactor);
			selecTableSemiMu.Fill(d,5,scaleFactor*lumiWeight);
			if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-3)
			  {
			    MSPlot["nEventsAfterCutsSemiMu"]->Fill(6, datasets[d], true, Luminosity*scaleFactor);
			    selecTableSemiMu.Fill(d,6,scaleFactor*lumiWeight);
			    if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-2)
			      {
				MSPlot["nEventsAfterCutsSemiMu"]->Fill(7, datasets[d], true, Luminosity*scaleFactor);
				selecTableSemiMu.Fill(d,7,scaleFactor*lumiWeight);
				if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-1)
				  {
				    MSPlot["nEventsAfterCutsSemiMu"]->Fill(8, datasets[d], true, Luminosity*scaleFactor);
				    selecTableSemiMu.Fill(d,8,scaleFactor*lumiWeight);
				    if(selectedJets.size()>=(unsigned int)anaEnv.NofJets)
				      {
					MSPlot["nEventsAfterCutsSemiMu"]->Fill(9, datasets[d], true, Luminosity*scaleFactor);
					selecTableSemiMu.Fill(d,9,scaleFactor*lumiWeight);
					eventSelectedSemiMu = true;
					float reliso = (selectedMuons[0]->chargedHadronIso()+selectedMuons[0]->neutralHadronIso()+selectedMuons[0]->photonIso())/selectedMuons[0]->Pt();
					MSPlot["SelectedEventsMuonsRelPFIso"]->Fill(reliso, datasets[d], true, Luminosity*scaleFactor);               
				      }
				  }
			      }
			  }
		      }
		  }
	      }
	  }
	}
      
      MSPlot["nEventsAfterCutsSemiEl"]->Fill(0, datasets[d], true, Luminosity*scaleFactor);
      selecTableSemiEl.Fill(d,0,scaleFactor*lumiWeight);
      if( semiElectron && triggedSemiEl )
	{
	  MSPlot["nEventsAfterCutsSemiEl"]->Fill(1, datasets[d], true, Luminosity*scaleFactor);
	  selecTableSemiEl.Fill(d,1,scaleFactor*lumiWeight);
	  if( isGoodPV )
	    {
	      MSPlot["nEventsAfterCutsSemiEl"]->Fill(2, datasets[d], true, Luminosity*scaleFactor);
	      selecTableSemiEl.Fill(d,2,scaleFactor*lumiWeight);
	      if( selectedElectrons.size() == 1 )
		{
		  MSPlot["nEventsAfterCutsSemiEl"]->Fill(3, datasets[d], true, Luminosity*scaleFactor);
		  selecTableSemiEl.Fill(d,3,scaleFactor*lumiWeight);
		  if( vetoMuons.size() == 0 )
		    {
		      MSPlot["nEventsAfterCutsSemiEl"]->Fill(4, datasets[d], true, Luminosity*scaleFactor);
		      selecTableSemiEl.Fill(d,4,scaleFactor*lumiWeight);
		      if( !selection.foundZCandidate(selectedElectrons[0], vetoElectronsSemiEl) )
			{
			  MSPlot["nEventsAfterCutsSemiEl"]->Fill(5, datasets[d], true, Luminosity*scaleFactor);
			  selecTableSemiEl.Fill(d,5,scaleFactor*lumiWeight);
			  if( selection.passConversionRejection(selectedElectrons[0]) )
			    {
			      MSPlot["nEventsAfterCutsSemiEl"]->Fill(6, datasets[d], true, Luminosity*scaleFactor);
			      selecTableSemiEl.Fill(d,6,scaleFactor*lumiWeight);
			      if( selectedJets.size()>=(unsigned int)anaEnv.NofJets-3 )
				{
				  MSPlot["nEventsAfterCutsSemiEl"]->Fill(7, datasets[d], true, Luminosity*scaleFactor);
				  selecTableSemiEl.Fill(d,7,scaleFactor*lumiWeight);
				  if( selectedJets.size()>=(unsigned int)anaEnv.NofJets-2 )
				    {
				      MSPlot["nEventsAfterCutsSemiEl"]->Fill(8, datasets[d], true, Luminosity*scaleFactor);
				      selecTableSemiEl.Fill(d,8,scaleFactor*lumiWeight);
				      if( selectedJets.size()>=(unsigned int)anaEnv.NofJets-1 )
					{
					  MSPlot["nEventsAfterCutsSemiEl"]->Fill(9, datasets[d], true, Luminosity*scaleFactor);
					  selecTableSemiEl.Fill(d,9,scaleFactor*lumiWeight);
					  if( selectedJets.size()>=(unsigned int)anaEnv.NofJets )
					    {
					      MSPlot["nEventsAfterCutsSemiEl"]->Fill(10, datasets[d], true, Luminosity*scaleFactor);
					      selecTableSemiEl.Fill(d,10,scaleFactor*lumiWeight);
					      eventSelectedSemiEl = true;
					      float reliso = (selectedElectrons[0]->chargedHadronIso()+selectedElectrons[0]->neutralHadronIso()+selectedElectrons[0]->photonIso())/selectedElectrons[0]->Pt();
					      MSPlot["SelectedEventsElectronsRelPFIso"]->Fill(reliso, datasets[d], true, Luminosity*scaleFactor);
					      
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
      
      if(eventSelectedSemiEl==true || eventSelectedSemiMu==true){  
	
	float CorrectRecMassW=0;
	float CorrectRecMassTop=0;
	vector<int> jetCombi;
	TLorentzVector hadrBQuark,hadrLQuark1,hadrLQuark2,leptBQuark;
	pair<unsigned int, unsigned int> leptonicBJet_, hadronicBJet_, hadronicWJet1_, hadronicWJet2_; //First index is the JET number, second the parton

	if(dataSetName.find("TTbarJets_SemiMu") == 0 || dataSetName.find("TTbarJets_SemiEl") == 0){
	  
	  leptonicBJet_ = hadronicBJet_ = hadronicWJet1_ = hadronicWJet2_ = pair<unsigned int, unsigned int>(9999,9999);
	  
	  vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV;
	  vector<TRootMCParticle*> mcParticlesMatching_;
	  mcParticlesMatching_.clear();	  
	  TLorentzVector topQuark, antiTopQuark;
	  
	  bool muPlusFromTop = false, muMinusFromTop = false, elPlusFromTop = false, elMinusFromTop = false;
	  int nTTbarQuarks = 0;
	  
	  for(unsigned int i=0; i<mcParticles.size(); i++){
	    if( mcParticles[i]->status() != 3) continue;
	    
	    if( mcParticles[i]->type() == 6 )
	      topQuark = *mcParticles[i];
	    else if( mcParticles[i]->type() == -6 )
	      antiTopQuark = *mcParticles[i];
	    
	    if( mcParticles[i]->type() == 13 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -6 )
	      muMinusFromTop = true;
	    if( mcParticles[i]->type() == -13 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == 6 )
	      muPlusFromTop = true;
	    if( mcParticles[i]->type() == 11 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -6 )
	      elMinusFromTop = true;
	    if( mcParticles[i]->type() == -11 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == 6 )
	      elPlusFromTop = true;
	    
	    if( abs(mcParticles[i]->type()) < 6 || abs(mcParticles[i]->type()) == 21 ){
	      mcParticlesTLV.push_back(*mcParticles[i]);
	      mcParticlesMatching_.push_back(mcParticles[i]);
	      
	      if( fabs(mcParticles[i]->motherType()) == 6 || fabs(mcParticles[i]->grannyType()) == 6 ){
		nTTbarQuarks++;
	      }
	    }
	  }
	  
	  // take all the selectedJets_ to study the radiation stuff, selectedJets_ are already ordened in decreasing Pt()
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
	    
	    if( fabs(mcParticlesMatching_[j]->type()) < 6 ){
	      if( ( ( muPlusFromTop || elPlusFromTop ) && mcParticlesMatching_[j]->motherType() == -24 && mcParticlesMatching_[j]->grannyType() == -6 ) || 
		  ( ( muMinusFromTop || elMinusFromTop ) && mcParticlesMatching_[j]->motherType() == 24 && mcParticlesMatching_[j]->grannyType() == 6 ) ){
		if(hadronicWJet1_.first == 9999) 
		  hadronicWJet1_ = JetPartonPair[i];
		else if(hadronicWJet2_.first == 9999) 
		  hadronicWJet2_ = JetPartonPair[i];
		else cerr<<"Found a third jet coming from a W boson which comes from a top quark..."<<endl;
	      }
	    }
	    if( fabs(mcParticlesMatching_[j]->type()) == 5 ){
	      if( ( ( muPlusFromTop || elPlusFromTop ) && mcParticlesMatching_[j]->motherType() == -6 )
		  || ( ( muMinusFromTop || elMinusFromTop ) && mcParticlesMatching_[j]->motherType() == 6 ) )
		hadronicBJet_ = JetPartonPair[i];
	      else if( ( ( muPlusFromTop || elPlusFromTop ) && mcParticlesMatching_[j]->motherType() == 6 )
		       || ( ( muMinusFromTop || elMinusFromTop ) && mcParticlesMatching_[j]->motherType() == -6 ) )
		leptonicBJet_ = JetPartonPair[i];
	    }
	    
	    // look for ISR stuff
	    if( fabs(mcParticlesMatching_[j]->type()) != 6 && fabs(mcParticlesMatching_[j]->motherType()) != 24 && fabs(mcParticlesMatching_[j]->motherType()) != 6 &&
		fabs(mcParticlesMatching_[j]->grannyType()) != 24 && fabs(mcParticlesMatching_[j]->grannyType()) != 6 ){
	      ISRJetPartonPair.push_back(JetPartonPair[i]);
	    }
	  }
	  
	  /*if(hadronicWJet1_.first != 9999 && hadronicWJet2_.first != 9999 && hadronicBJet_.first != 9999 && leptonicBJet_.first != 9999){
	    
	    all4PartonsMatched = true;
	    if(hadronicWJet1_.first < 4 && hadronicWJet2_.first < 4 && hadronicBJet_.first < 4 && leptonicBJet_.first < 4)
	      all4JetsMatched_MCdef_ = true;
	  }
	  if(hadronicWJet1_.first < 4 && hadronicWJet2_.first < 4 && hadronicBJet_.first < 4)
	    hadronictopJetsMatched_MCdef_ = true;
	  */

	  /*for(int i=1; i<=histo1D_["PtJetCut_nEventsBefore"]->GetNbinsX(); i++){
	    float binCenter = histo1D_["PtJetCut_nEventsBefore"]->GetBinCenter(i);
	    histo1D_["PtJetCut_nEventsBefore"]->Fill(binCenter);
	    if(selectedJets_.size() > 3 && selectedJets_[3]->Pt() > binCenter){
	      histo1D_["PtJetCut_nEventsAfter"]->Fill(binCenter);
	      if( hadronicWJet1_.first > 3 || hadronicWJet2_.first > 3 || hadronicBJet_.first > 3 || leptonicBJet_.first > 3 ){
		bool foundISR = false;
		for(unsigned int j=0; j<ISRJetPartonPair.size(); j++)
		  if(ISRJetPartonPair[j].first < 4)
		    foundISR = true;
		
	      }
	    }
	    }*/
	  
	  /*if(all4PartonsMatched){
	    for(int i=1; i<=histo1D_["PtJetCut_nEventsBefore"]->GetNbinsX(); i++){
	      float binCenter = histo1D_["PtJetCut_nEventsBefore"]->GetBinCenter(i);
	      histo1D_["PtJetCut_nEvents4PartonsMatchedBefore"]->Fill(binCenter);
	      if(selectedJets_.size() > 3 && selectedJets_[leptonicBJet_.first]->Pt() > binCenter && selectedJets_[hadronicBJet_.first]->Pt() > binCenter &&
		 selectedJets_[hadronicWJet1_.first]->Pt() > binCenter && selectedJets_[hadronicWJet2_.first]->Pt() > binCenter)
		histo1D_["PtJetCut_nEvents4PartonsMatchedAfter"]->Fill(binCenter);
	    }
	    
	    if(all4JetsMatched_MCdef_){
	      // all 4 jets found and matched, now do something with them!
	      
	      for(int i=1; i<=histo1D_["PtJetCut_nEventsBefore"]->GetNbinsX(); i++){
		float binCenter = histo1D_["PtJetCut_nEventsBefore"]->GetBinCenter(i);
		histo1D_["PtJetCut_nEvents4JetsMatchedBefore"]->Fill(binCenter);
		if(selectedJets_.size() > 3 && selectedJets_[3]->Pt() > binCenter)
		  histo1D_["PtJetCut_nEvents4JetsMatchedAfter"]->Fill(binCenter);
	      }
	    }
	    }*/

	  jetCombi.push_back(hadronicWJet1_.first);
	  jetCombi.push_back(hadronicWJet2_.first);
	  jetCombi.push_back(hadronicBJet_.first);
	  jetCombi.push_back(leptonicBJet_.first);
	  
	  if(hadronicBJet_.second !=9999){ hadrBQuark = *mcParticles[hadronicBJet_.second];}
	  if(hadronicWJet1_.second != 9999){hadrLQuark1 = *mcParticles[hadronicWJet1_.second];}
	  if(hadronicWJet2_.second != 9999){hadrLQuark2 = *mcParticles[hadronicWJet2_.second];}
	  if(leptonicBJet_.second != 9999){leptBQuark = *mcParticles[leptonicBJet_.second];}	  
	  
	  //Working on generator level (i.e. jets level):  
	  if(jetCombi[0]!=9999 && jetCombi[1]!=9999 && jetCombi[2]!=9999 && jetCombi[3]!=9999){    
	    CorrectRecMassW=(*selectedJets[jetCombi[0]]+*selectedJets[jetCombi[1]]).M();
	    CorrectRecMassTop=(*selectedJets[jetCombi[0]]+*selectedJets[jetCombi[1]]+*selectedJets[jetCombi[2]]).M();
	    
	    histo1D["WMass"]->Fill(CorrectRecMassW);
	    histo1D["TopMass"]->Fill(CorrectRecMassTop);
	  }	      	      	      	       	      	  
	  
	}//if Semi mu or semi el ttbar
	
	if(CalculateResolutions){
	  // Fill the resolution-stuff for events where the 4 ttbar semi-lep partons are all matched to jets
	  if(hadronicWJet1_.first < 9999 && hadronicWJet2_.first < 9999 && hadronicBJet_.first < 9999 && leptonicBJet_.first < 9999 ){
	    resFitLightJets->Fill(selectedJets[hadronicWJet1_.first], mcParticles[hadronicWJet1_.second]);
	    resFitLightJets->Fill(selectedJets[hadronicWJet2_.first], mcParticles[hadronicWJet2_.second]);
	    resFitBJets->Fill(selectedJets[hadronicBJet_.first], mcParticles[hadronicBJet_.second]);
	    resFitBJets->Fill(selectedJets[leptonicBJet_.first], mcParticles[leptonicBJet_.second]);
	    if(semiMuon == true) resFitMuon->Fill(selectedMuons[0], &standardLepton);
	    else if(semiElectron == true) resFitElectron->Fill(selectedElectrons[0], &standardLepton);
	    //histo1D["ElectronPt"]->Fill(selectedElectrons[0]->Pt());
	    //histo1D["ElectronEta"]->Fill(selectedElectrons[0]->Eta());
	    resFitNeutrino->Fill(mets[0], &standardNeutrino);
	  }
	}

	if(!CalculateResolutions){
	  //oooooooooOOOOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOOOOoooooooooooooOOOOOOOOOOOOOOOOooooooooooooOOOOOOOOOO
	  //                  Perform Kinematic Fit --> Need ChiSquared value for all 12 combinations
	  //oooooooooOOOOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOOOOoooooooooooooOOOOOOOOOOOOOOOOooooooooooooOOOOOOOOOO  
	
	  int NumberCombinations =0;
	  int BLeptIndex[12];
	  int BHadrIndex[12];
	  int Light1Index[12];
	  int Light2Index[12];

	  for(int i=0;i<3;i++){
	    for(int j=i+1;j<4;j++){
	      for(int k=0;k<4;k++){
	        if(k!=i && k!=j){
		
		  TLorentzVector lightJet1 = *selectedJets[i];
		  TLorentzVector lightJet2 = *selectedJets[j];
		  TLorentzVector bHadrJet = *selectedJets[k];
		  TLorentzVector leptonKinFit;
		  if(semiMuon == true) leptonKinFit = *selectedMuons[0];
		  else if(semiElectron == true) leptonKinFit = *selectedElectrons[0];
		  else cout << " Both false " << endl;
		  TLorentzVector neutrinoKinFit = *mets[0];
		  TLorentzVector bLeptJet;
		  int LeptNumber;
		  for(int ii = 0; ii < 4; ii++){
		    if( ii != i && ii != j && ii != k){
		      LeptNumber = ii;
		    }
		  }
		  bLeptJet = *selectedJets[LeptNumber];

		  BHadrIndex[NumberCombinations] = k;
		  Light1Index[NumberCombinations] = i;
		  Light2Index[NumberCombinations] = j;
		  BLeptIndex[NumberCombinations] = LeptNumber;
		  		
		  // prepare everything for the Kinematic Fit
		  TMatrixD Ml1(3,3), Ml2(3,3), Mbh(3,3), Mbl(3,3), Mlep(3,3), Mne(3,3);
		  Ml1.Zero(); Ml2.Zero(); Mbh.Zero(); Mlep.Zero(); Mne.Zero();
		  Ml1(0,0) = pow(resFitLightJets->EtResolution(&lightJet1), 2);
		  Ml1(1,1) = pow(resFitLightJets->ThetaResolution(&lightJet1), 2);
		  Ml1(2,2) = pow(resFitLightJets->PhiResolution(&lightJet1), 2);
		  Ml2(0,0) = pow(resFitLightJets->EtResolution(&lightJet2), 2);
		  Ml2(1,1) = pow(resFitLightJets->ThetaResolution(&lightJet2), 2);
		  Ml2(2,2) = pow(resFitLightJets->PhiResolution(&lightJet2), 2);
		  Mbh(0,0) = pow(resFitBJets->EtResolution(&bHadrJet), 2);
		  Mbh(1,1) = pow(resFitBJets->ThetaResolution(&bHadrJet), 2);
		  Mbh(2,2) = pow(resFitBJets->PhiResolution(&bHadrJet), 2);
		  Mbl(0,0) = pow(resFitBJets->EtResolution(&bLeptJet), 2);
		  Mbl(1,1) = pow(resFitBJets->ThetaResolution(&bLeptJet), 2);
		  Mbl(2,2) = pow(resFitBJets->PhiResolution(&bLeptJet), 2);
		  if(semiMuon == true){
		    Mlep(0,0) = pow(resFitMuon->EtResolution(&leptonKinFit), 2);
		    Mlep(1,1) = pow(resFitMuon->ThetaResolution(&leptonKinFit), 2);
		    Mlep(2,2) = pow(resFitMuon->PhiResolution(&leptonKinFit),2);
		  }
		  else if(semiElectron == true){
		    Mlep(0,0) = pow(resFitElectron->EtResolution(&leptonKinFit), 2);
		    Mlep(1,1) = pow(resFitElectron->ThetaResolution(&leptonKinFit), 2);
		    Mlep(2,2) = pow(resFitElectron->PhiResolution(&leptonKinFit),2);
		  }
		  Mne(0,0) = pow(resFitNeutrino->EtResolution(&neutrinoKinFit),2);
		  Mne(1,1) = pow(9999.,2);
		  Mne(2,2) = pow(resFitNeutrino->PhiResolution(&neutrinoKinFit),2);

		  float WMassConstraintFit;
		  float TopMassConstraintFit;
		  float LeptTopMassConstraintFit;
		  for(int ii = 0; ii<2; ii++){
		    if(ii == 0){
		      WMassConstraintFit = WMassKinFit;
		      TopMassConstraintFit = TopMassKinFit; //Different mass for MC and Data!!
		    }
		    else{
		      WMassConstraintFit = MassW;
		      TopMassConstraintFit = MassTop;
		    }		    
            		
		    TKinFitter *theFitter = new TKinFitter("hadtopFit", "hadtopFit");
		    TKinFitter *theFullFitter = new TKinFitter("hadAndLepTopFit", "hadAndLepTopFit");
		    theFitter->setVerbosity(0);
		    theFullFitter->setVerbosity(0);
		
		    TFitParticleEtThetaPhiEMomFix *fitLight1 = new TFitParticleEtThetaPhiEMomFix("lightJet1", "lightJet1", &lightJet1, &Ml1);
		    TFitParticleEtThetaPhiEMomFix *fitLight2 = new TFitParticleEtThetaPhiEMomFix("lightJet2", "lightJet2", &lightJet2, &Ml2);
		    TFitParticleEtThetaPhiEMomFix *fitBHadr = new TFitParticleEtThetaPhiEMomFix("bHadrJet", "bHadrJet", &bHadrJet, &Mbh);
		    TFitParticleEtThetaPhiEMomFix *fitBLept = new TFitParticleEtThetaPhiEMomFix("bLeptJet", "bLeptJet", &bLeptJet, &Mbl);
		    TFitParticleEtThetaPhiEMomFix *fitLepton = new TFitParticleEtThetaPhiEMomFix("leptonKinFit", "leptonKinFit", &leptonKinFit, &Mlep);
		    TFitParticleEtThetaPhi *fitNeutrino = new TFitParticleEtThetaPhi("neutrinoKinFit", "neutrinoKinFit", &neutrinoKinFit, &Mne);
		    theFitter->addMeasParticles(fitLight1,fitLight2,fitBHadr);
		    theFullFitter->addMeasParticles(fitLight1,fitLight2,fitBHadr,fitBLept,fitLepton,fitNeutrino);
		  
		    TFitConstraintM *consWHadr = new TFitConstraintM("WBosonMassHadr", "MassConstraintHadr", 0, 0, WMassConstraintFit);
		    TFitConstraintM *consTopHadr = new TFitConstraintM("TopQuarkMassHadr", "MassConstraintHadr", 0, 0, TopMassConstraintFit );
		    TFitConstraintM *consWLept = new TFitConstraintM("WBosonMassLept", "MassConstraintLept", 0, 0, WMassConstraintFit);
		    TFitConstraintM *consTopLept = new TFitConstraintM("TopQuarkMassLept", "MassConstraintLept", 0, 0, TopMassConstraintFit );
		    consWHadr->addParticles1(fitLight1,fitLight2);
		    consTopHadr->addParticles1(fitBHadr,fitLight1,fitLight2);
		    consWLept->addParticles1(fitLepton,fitNeutrino);
		    consTopLept->addParticles1(fitBLept,fitLepton,fitNeutrino);	
		
		    theFitter->addConstraint(consWHadr);
		    theFitter->addConstraint(consTopHadr);
		    theFullFitter->addConstraint(consWHadr);
		    theFullFitter->addConstraint(consTopHadr);
		    theFullFitter->addConstraint(consWLept);
		    //Only apply leptonic top mass constraint in the theoretical value case!!!
		    if(ii==0){theFullFitter->addConstraint(consTopLept);}
		    
		    theFitter->setMaxNbIter(200);
		    theFitter->setMaxDeltaS(5e-5);
		    theFitter->setMaxF(1e-4);
		    theFullFitter->setMaxNbIter(200);
		    theFullFitter->setMaxDeltaS(5e-5);
		    theFullFitter->setMaxF(1e-4);
		
		    //-----------------------------------
		    // Execution of hadronic fitter:
		    //-----------------------------------
		    theFitter->fit();
		    if(theFitter->getStatus() == 0){ChiSquared[ii].push_back(theFitter->getS());}// if the fitter converged
		    else{ChiSquared[ii].push_back(9999);}	//Need to push_back 9999 in stead of -9999 to avoid that this configuration is selected as the one with the lowest ChiSquare value

		    fittedLepton[ii].push_back( *(fitLepton->getCurr4Vec()) );
		    fittedNeutrino[ii].push_back( *(fitNeutrino->getCurr4Vec()) );
		    fittedBLept[ii].push_back( *(fitBLept->getCurr4Vec()) );
		    fittedBHadr[ii].push_back( *(fitBHadr->getCurr4Vec()) );
		    fittedLight1[ii].push_back( *(fitLight1->getCurr4Vec()) );
		    fittedLight2[ii].push_back( *(fitLight2->getCurr4Vec()) );		    
		    
		    //-------------------------------------------
		    // Execution of hadronic+leptonic fitter:
		    //-------------------------------------------
		    theFullFitter->fit();
		    if( theFullFitter->getStatus() == 0){ChiSquaredFull[ii].push_back(theFitter->getS());}
		    else{ChiSquaredFull[ii].push_back(9999);}
		    
		    fittedFullLepton[ii].push_back( *(fitLepton->getCurr4Vec()) );
		    fittedFullNeutrino[ii].push_back( *(fitNeutrino->getCurr4Vec()) );
		    fittedFullBLept[ii].push_back( *(fitBLept->getCurr4Vec()) );
		    fittedFullBHadr[ii].push_back( *(fitBHadr->getCurr4Vec()) );
		    fittedFullLight1[ii].push_back( *(fitLight1->getCurr4Vec()) );
		    fittedFullLight2[ii].push_back( *(fitLight2->getCurr4Vec()) );
		    		    		
		    delete theFitter;
		    delete theFullFitter;
		    delete fitLight1;
		    delete fitLight2;
		    delete fitBHadr;
		    delete fitBLept;
		    delete fitLepton;
		    delete fitNeutrino;
		    delete consWHadr;
		    delete consTopHadr;
		    delete consWLept;
		    delete consTopLept;

		  }//end of loop for mass constraint fit for two different mass options (fixed value or value obtained from fit)
		
		  NumberCombinations++;
	        }
	      }//end of k loop for jet combination selection
	    }//end of j loop for jet combination selection
	  }//end of i loop for jet combination selection		
	  
	  //oooooooooOOOOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOOOOoooooooooooooOOOOOOOOOOOOOOOO
	  //           Initialize and define all variables necessary for WTree
	  //oooooooooOOOOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOOOOoooooooooooooOOOOOOOOOOOOOOOO	
	
	  std::vector<float> TCHEbTagValues, TCHPbTagValues, SSVHEbTagValues, SSVHPbTagValues, CSVbTagValues;
	  for(int ii = 0; ii<selectedJets.size();ii++){
	    TCHEbTagValues.push_back(selectedJets[ii]->btag_trackCountingHighEffBJetTags());
	    TCHPbTagValues.push_back(selectedJets[ii]->btag_trackCountingHighPurBJetTags());
	    SSVHEbTagValues.push_back(selectedJets[ii]->btag_simpleSecondaryVertexHighEffBJetTags());
	    SSVHPbTagValues.push_back(selectedJets[ii]->btag_simpleSecondaryVertexHighPurBJetTags());
	    CSVbTagValues.push_back(selectedJets[ii]->btag_combinedSecondaryVertexBJetTags());
	  }
	
	  vector<TLorentzVector> SelectedJets;
	  for(int ii=0; ii<selectedJets.size();ii++){	
	    SelectedJets.push_back(*selectedJets[ii]);
	  }

	  //Initialize jetCombi values for not SemiMu sample:
	  if(dataSetName.find("TTbarJets_SemiMu") != 0){
	    jetCombi.push_back(-9999);
	    jetCombi.push_back(-9999);
	    jetCombi.push_back(-9999);
	    jetCombi.push_back(-9999);
	  }

	  //Set all variables defined in WTree
	  wTree = new WTree();
	  wTree->setEventID( event->eventId() );
	  wTree->setRunID( event->runId() );
	  wTree->setLumiBlockID( event->lumiBlockId() );
	  wTree->setNPV(vertex.size());
	  wTree->setNPUBXm1(event->nPu(-1));
	  wTree->setNPU(event->nPu(0));
	  wTree->setNPUBXp1(event->nPu(1));
	  wTree->setKinFitResults(ChiSquared[0],fittedLepton[0],fittedNeutrino[0],fittedBLept[0],fittedBHadr[0],fittedLight1[0],fittedLight2[0]);
	  wTree->setFullKinFitResults(ChiSquaredFull[0],fittedFullLepton[0],fittedFullNeutrino[0],fittedFullBLept[0],fittedFullBHadr[0],fittedFullLight1[0],fittedFullLight2[0]);
	  wTree->setKinFitResultsMassFit(ChiSquared[1],fittedLepton[1],fittedNeutrino[1],fittedBLept[1],fittedBHadr[1],fittedLight1[1],fittedLight2[1]);
	  wTree->setFullKinFitResultsMassFit(ChiSquaredFull[1],fittedFullLepton[1],fittedFullNeutrino[1],fittedFullBLept[1],fittedFullBHadr[1],fittedFullLight1[1],fittedFullLight2[1]);
	  wTree->setHadrBJet( jetCombi[2] );
	  wTree->setHadrLJet1( jetCombi[0] );
	  wTree->setHadrLJet2( jetCombi[1] );
	  wTree->setLeptBJet( jetCombi[3] );
	  wTree->setMET( *mets[0] );
	  //cout << "met: " << mets[0]->Pt() << endl;
	  wTree->setSelectedJets( SelectedJets );
	  wTree->setBTagTCHE(TCHEbTagValues);  
	  wTree->setBTagTCHP(TCHPbTagValues);
	  wTree->setBTagSSVHE(SSVHEbTagValues);
	  wTree->setBTagSSVHP(SSVHPbTagValues);
	  wTree->setBTagCSV(CSVbTagValues);
	  if(semiMuon == true) wTree->setMuon( *selectedMuons[0] );
	  else if(semiElectron == true) wTree->setMuon(*selectedElectrons[0]);
	  wTree->setHadrBQuark( hadrBQuark );  
	  wTree->setHadrLQuark1( hadrLQuark1 );
	  wTree->setHadrLQuark2( hadrLQuark2 );
	  wTree->setLeptBQuark( leptBQuark );
	  wTree->setStandardCosTheta( standardCosTheta );
	  wTree->setStandardNeutrino( standardNeutrino);
	  wTree->setStandardLepton( standardLepton);
          
	  WTreeTree->Fill();
	  delete wTree;
	}  // end of !CalculateResolutions
      }  //delete selection;    
    }//loop on events

    //////////////////////////////
    //  Executing fits          //
    //////////////////////////////

    if(!CalculateResolutions && ((semiMuon == true && dataSetName.find("TTbarJets_SemiMu") == 0 ) || (semiElectron == true && dataSetName.find("TTbarJets_SemiEl") == 0)) ){
      histo1D["WMass"]->Fit("gaus","Q");     
      histo1D["TopMass"]->Fit("gaus","Q");
      histo1D["WMass"]->Fit("gaus","Q","",histo1D["WMass"]->GetFunction("gaus")->GetParameter(1)-histo1D["WMass"]->GetFunction("gaus")->GetParameter(2),histo1D["WMass"]->GetFunction("gaus")->GetParameter(1)+histo1D["WMass"]->GetFunction("gaus")->GetParameter(2));
      histo1D["TopMass"]->Fit("gaus","Q","",histo1D["TopMass"]->GetFunction("gaus")->GetParameter(1)-histo1D["TopMass"]->GetFunction("gaus")->GetParameter(2),histo1D["TopMass"]->GetFunction("gaus")->GetParameter(1)+histo1D["TopMass"]->GetFunction("gaus")->GetParameter(2));
      
      std::cout << " sigma values : Top = " << histo1D["WMass"]->GetFunction("gaus")->GetParameter(2) << " , W = " << histo1D["TopMass"]->GetFunction("gaus")->GetParameter(2) << std::endl;
      std::cout << " mass values : Top = " << histo1D["WMass"]->GetFunction("gaus")->GetParameter(1) << " , W = " << histo1D["TopMass"]->GetFunction("gaus")->GetParameter(1) << std::endl;

      cout << " -----------------------------------------------------------------------------------------------------------------------" << endl;
      cout << " Performing helicity Generator fit : " << endl;
      cout << " ------------------------------------" << endl;
      TF1 *helicityFit = new TF1("helicityFit","[0]*((((1-[1]-[2])*3*(1+x)*(1+x))+([1]*3*(1-x)*(1-x))+([2]*6*(1-x*x)))/8)",-1,1);
      histo1D["StandardCosThetaFit"]->Fit("helicityFit","Q");
      std::cout << " fit values (before event selection) : Norm =" <<helicityFit->GetParameter(0) << " , Left = " << helicityFit->GetParameter(1) << " Long = " << helicityFit->GetParameter(2) << " ==> Right = " << 1-(helicityFit->GetParameter(1))-(helicityFit->GetParameter(2))<< std::endl;
      std::cout << " fit values error (before event selection) : " << helicityFit->GetParError(0) << " " << helicityFit->GetParError(1) << " " << helicityFit->GetParError(2) << std::endl;
      cout << "                      ------------------------------------" << endl;
      histo1D["StandardCosThetaFit"]->Scale(100./(histo1D["StandardCosThetaFit"]->Integral()));
      histo1D["StandardCosThetaFit"]->SetMinimum(0);
      histo1D["StandardCosThetaFit"]->SetMaximum(0.8);
      TF1 *helicityFit2 = new TF1("helicityFit2","((([0]*3*(1+x)*(1+x))+([1]*3*(1-x)*(1-x))+([2]*6*(1-x*x)))/8)",-1,1);
      histo1D["StandardCosThetaFit"]->Fit("helicityFit2","Q");
      std::cout << " fit values 2 (before event selection) : " << helicityFit2->GetParameter(0) << " " << helicityFit2->GetParameter(1) << " " << helicityFit2->GetParameter(2) << std::endl;
      std::cout << " fit values error 2 (before event selection) : " << helicityFit2->GetParError(0) << " " << helicityFit2->GetParError(1) << " " << helicityFit2->GetParError(2) << std::endl;
      cout << " -----------------------------------------------------------------------------------------------------------------------" << endl;    
    }
    
    WTreeFile->cd();
      
    TTree *configTreeWTreeFile = new TTree("configTreeWTreeFile","configuration Tree in WTree File");
    TClonesArray* tcdatasetwTreefile = new TClonesArray("Dataset",1);  //If filename is lost, this makes it possible to retrieve the corresponding data sample which is used to create the TFile.
    configTreeWTreeFile->Branch("Dataset","TClonesArray",&tcdatasetwTreefile);
    TClonesArray* tcAnaEnvWTreeFile = new TClonesArray("AnalysisEnvironment",1);
    configTreeWTreeFile->Branch("AnaEnv","TClonesArray",&tcAnaEnvWTreeFile);
    new ((*tcAnaEnvWTreeFile)[0]) AnalysisEnvironment(anaEnv);
    new ((*tcdatasetwTreefile)[0]) Dataset(*datasets[d]);
      
    configTreeWTreeFile->Fill();
    configTreeWTreeFile->Write();
    WTreeTree->Write();
    WTreeFile->Close();
    delete WTreeFile;      
      
    //important: free memory
    treeLoader.UnLoadDataset();
  }				//loop on datasets  -->Stops already here
    
  //Once everything is filled ...
  if (verbose > 0)
    cout << " We ran over all the data ;-)" << endl;
    
  //Selection tables
  selecTableSemiMu.TableCalculator(false, true, true, true);
  selecTableSemiEl.TableCalculator(false, true, true, true);
  string selectiontableSemiMu = "SelectionTable_SemiMu_Macro";
  string selectiontableSemiEl = "SelectionTable_SemiEl_Macro";
  if (argc >= 3){
    string sample=string(argv[2]);
    selectiontableSemiMu = selectiontableSemiMu +"_"+sample;
    selectiontableSemiEl = selectiontableSemiEl +"_"+sample;
  }
  selectiontableSemiMu = selectiontableSemiMu +".tex"; 	
  selectiontableSemiEl = selectiontableSemiEl +".tex"; 	
  selecTableSemiMu.Write(selectiontableSemiMu.c_str());
  selecTableSemiEl.Write(selectiontableSemiEl.c_str());
    
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

  // Fill the resolution histograms and calculate the resolutions
  if(CalculateResolutions)
    {
      mkdir((pathPNG+"resFit_LightJet/").c_str(),0777);
      mkdir((pathPNG+"resFit_BJet/").c_str(),0777);
      if(semiMuon == true){
	mkdir((pathPNG+"resFit_Muon/").c_str(),0777);
	mkdir((pathPNG+"resFit_NeutrinoSemiMu/").c_str(),0777);     
	resFitMuon->WritePlots(fout, true, pathPNG+"resFit_Muon/");
	resFitMuon->WriteResolutions("muonReso.root");
	resFitNeutrino->WritePlots(fout, true, pathPNG+"resFit_NeutrinoSemiMu/");
	resFitNeutrino->WriteResolutions("neutrinoSemiMuReso.root");
      }
      else if(semiElectron == true){
	mkdir((pathPNG+"resFit_NeutrinoSemiEl/").c_str(),0777);      
	mkdir((pathPNG+"resFit_Electron/").c_str(),0777);
	resFitNeutrino->WritePlots(fout, true, pathPNG+"resFit_Neutrino/");
	resFitNeutrino->WriteResolutions("neutrinoSemiElReso.root");
	resFitElectron->WritePlots(fout, true, pathPNG+"resFit_Electron/");
	resFitElectron->WriteResolutions("electronReso.root");      
      }

      resFitLightJets->WritePlots(fout, true, pathPNG+"resFit_LightJet/");
      resFitLightJets->WriteResolutions("lightJetReso.root");
      resFitBJets->WritePlots(fout, true, pathPNG+"resFit_BJet/");
      resFitBJets->WriteResolutions("bJetReso.root");
    }

  fout->cd();  
  
  // MSPlots 
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++){  //Because all the histograms are stored as maps it is quite easy to loop over them      
    MultiSamplePlot *temp = it->second;
    string name = it->first;
    cout << " name = " << name << endl;
    temp->Draw(false, name, false, true, true, true, true);
    temp->Write(fout, name, true, pathPNG+"MSPlot/");
  }
    
  //Write 1D histograms
  TDirectory* th1dir = fout->mkdir("1D_histograms");
  fout->cd();
  th1dir->cd();
        
  for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++){
    TH1F *temp = it->second;
    string name = it->first;
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

