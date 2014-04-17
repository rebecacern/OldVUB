#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

// Root stuff
#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TF2.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TBranch.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TDirectory.h"

#include "../Tools/interface/PlottingTools.h"
#include "../Tools/interface/MultiSamplePlot.h"
#include "../Content/interface/AnalysisEnvironment.h"
#include "../Content/interface/Dataset.h"
#include "../JESMeasurement/interface/Monster.h"
#include "../JESMeasurement/interface/MonsterCombination.h"
#include "../JESMeasurement/interface/MonsterTools.h"
#include "../JESMeasurement/interface/BinnedMonsterCombination.h"
#include "../JESMeasurement/interface/MVABinnedMonsterCombination.h"
#include "../JESMeasurement/interface/LightJetCombiner.h"
#include "../Reconstruction/interface/JetCorrectorParameters.h"
#include "../Reconstruction/interface/FactorizedJetCorrector.h"

#include "Style.C"

using namespace std;

/// TGraphAsymmErrors
map<string,TGraphAsymmErrors*> graphAsymmErr;
map<string,TGraphErrors*> graphErr;

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

int main (int argc, char *argv[])
{
  clock_t start = clock();

  cout << "***************************************************" << endl;
  cout << " Beginning of the program for TTbar JES Analysis ! " << endl;
  cout << "***************************************************" << endl;
 
//  setTDRStyle();
  setMyStyle();

  string pathPNG = "PlotsJESAnalysis/";
  mkdir(pathPNG.c_str(),0777);
  mkdir((pathPNG+"MSPlot/").c_str(),0777);
  
  // configuration
  bool measureTopMass = true;
  bool measureTopMassDifference = false;
  
//  float Luminosity = 36.1389;
  float Luminosity = 5215.88; // lumi of full ttbar-semimu-analysis sample
//  float Luminosity = 350; // lumi of the smallest sample (WJets), ignoring QCD
  
  if(measureTopMass)
    cout << "Executing the Top Mass analysis ";
  else if(measureTopMassDifference)
    cout << "Executing the Top Mass difference analysis";
  else
    cout << "Executing the JES analysis ";
  cout << "for an integrated luminosity of " << Luminosity << " pb^-1" << endl;
  
  vector<string> inputMonsters;
  if(measureTopMass)
  {
//    inputMonsters.push_back("Monsters/KinFit_Monsters_TopMass_TTbarJets_SemiMuon_Mass_178_5.root");
//    inputMonsters.push_back("Monsters/KinFit_Monsters_TopMass_TTbarJets_Other_Mass_178_5.root");
//    inputMonsters.push_back("Monsters/KinFit_Monsters_TopMass_TTbarJets_SemiMuon_Mass_166_5.root");
//    inputMonsters.push_back("Monsters/KinFit_Monsters_TopMass_TTbarJets_Other_Mass_166_5.root");
    inputMonsters.push_back("Monsters/KinFit_Monsters_TopMass_TTbarJets_SemiMuon_Analysis.root");
//    inputMonsters.push_back("Monsters/KinFit_Monsters_TopMass_TTbarJets_SemiMuon_Training.root");
//    inputMonsters.push_back("Monsters/KinFit_Monsters_TopMass_TTbarJets_Other.root");
//    inputMonsters.push_back("Monsters/KinFit_Monsters_TopMass_ST_tWChannel.root");
//    inputMonsters.push_back("Monsters/KinFit_Monsters_TopMass_ST_tChannel.root");
//    inputMonsters.push_back("Monsters/KinFit_Monsters_TopMass_WJets.root");
//    inputMonsters.push_back("Monsters/KinFit_Monsters_TopMass_ZJets.root");
//    inputMonsters.push_back("Monsters/KinFit_Monsters_TopMass_QCD_Mu15.root");
//    inputMonsters.push_back("Monsters/KinFit_Monsters_TopMass_Data.root");
  }
  else
  {
//    inputMonsters.push_back("Monsters/KinFit_Monsters_TTbarJets_SemiMuon_Analysis_PtBinnedIter8.root");
//    inputMonsters.push_back("Monsters/KinFit_Monsters_TTbarJets_SemiMuon_Analysis_JES_m5.root");
//    inputMonsters.push_back("Monsters/KinFit_Monsters_TTbarJets_SemiMuon_Analysis_JES_p5.root");
    inputMonsters.push_back("Monsters/KinFit_Monsters_TTbarJets_SemiMuon_Analysis.root");
//    inputMonsters.push_back("Monsters/KinFit_Monsters_TTbarJets_SemiMuon_Analysis_Iter6.root");
//    inputMonsters.push_back("Monsters/KinFit_Monsters_TTbarJets_Other.root");
//    inputMonsters.push_back("Monsters/KinFit_Monsters_ST_tWChannel.root");
//    inputMonsters.push_back("Monsters/KinFit_Monsters_ST_tChannel.root");
//    inputMonsters.push_back("Monsters/KinFit_Monsters_WJets.root");
//    inputMonsters.push_back("Monsters/KinFit_Monsters_ZJets.root");
//    inputMonsters.push_back("Monsters/KinFit_Monsters_QCD_Mu15.root");
//    inputMonsters.push_back("Monsters/KinFit_Monsters_Data.root");
  }
  
  TFile *fout = new TFile ("TTbarJES_Analysis.root", "RECREATE");
//  TFile *fout = new TFile ("TTbarJES_Analysis_Mass_PtBinnedIter8.root", "RECREATE");
  fout->cd();
  
  FullKinFit* kinFit = new FullKinFit(0, 0, 0, measureTopMass, measureTopMassDifference, true);
  TH2F* dummyMonster = kinFit->DummyMonster();
  delete kinFit;
  
  TH2F* expDEl2D = 0;
  TH2F* expDEb2D = 0;
  TFile *inJECFile = new TFile("ExpectedDEdiff.root","READ");
  if( inJECFile->Get("ExpectedDEl_nocut_PtBjetBins") != 0 && inJECFile->Get("ExpectedDEb_nocut_PtBjetBins") != 0 )
	{
	  cout << "Loading the expected corrections" << endl;
    expDEl2D = (TH2F*) inJECFile->Get("ExpectedDEl_nocut_2DPtBins");
    expDEb2D = (TH2F*) inJECFile->Get("ExpectedDEb_nocut_2DPtBins");
	  expDEl2D->SetDirectory(NULL);
	  expDEb2D->SetDirectory(NULL);
  }
  else
    cout << "inJECFile:  Problems while loading the expected corrections!" << endl;
	  
  inJECFile->Close();
  delete inJECFile;
  
  // Official L7 JEC, applied after L2L3
  vector<JetCorrectorParameters> vCorrParamLight, vCorrParamB;
  JetCorrectorParameters *LightJetCorPar = new JetCorrectorParameters("JECFiles/START38_V14_AK5PF_L7Parton_qT.txt");
  vCorrParamLight.push_back(*LightJetCorPar);
  FactorizedJetCorrector *JEClight = new FactorizedJetCorrector(vCorrParamLight);
  
  JetCorrectorParameters *BJetCorPar = new JetCorrectorParameters("JECFiles/START38_V14_AK5PF_L7Parton_bT.txt");
  vCorrParamB.push_back(*BJetCorPar);
  FactorizedJetCorrector *JECb = new FactorizedJetCorrector(vCorrParamB);
  
  // initialize histograms
  cout << "Initializing histograms" << endl;
  
 	//control plots
 	histo1D["mW_gen"] = new TH1F("mW_gen","mW_gen;mW;#events",50,50,100);
	histo1D["mTop_gen"] = new TH1F("mTop_gen","mTop_gen;mTop;#events",50,150,200);
 	histo1D["mW_MCcomb"] = new TH1F("mW_MCcomb","mW_MCcomb;mW;#events",75,00,150);
	histo1D["mTop_MCcomb"] = new TH1F("mTop_MCcomb","mTop_MCcomb;mTop;#events",75,100,250);
 	histo1D["mTop_MCcomb_bCorr"] = new TH1F("mTop_MCcomb_bCorr","mTop_MCcomb_bCorr;mTop;#events",75,100,250);
 	histo1D["mTop_MCcomb_AllCorr"] = new TH1F("mTop_MCcomb_AllCorr","mTop_MCcomb_AllCorr;mTop;#events",75,100,250);
 	histo1D["mTop_MCcomb_JECbCorr"] = new TH1F("mTop_MCcomb_JECbCorr","mTop_MCcomb_JECbCorr;mTop;#events",75,100,250);
 	histo1D["mTop_MCcomb_JECAllCorr"] = new TH1F("mTop_MCcomb_JECAllCorr","mTop_MCcomb_JECAllCorr;mTop;#events",75,100,250);
	histo1D["mTop_MCcomb_MaxProb"] = new TH1F("mTop_MCcomb_MaxProb","mTop_MCcomb_MaxProb;mTop;#events",75,100,250);
 	histo1D["mTop_MCcomb_MaxProbProbNoCorr"] = new TH1F("mTop_MCcomb_MaxProbProbNoCorr","mTop_MCcomb_MaxProbProbNoCorr;mTop;#events",75,100,250);
 	histo1D["mTop_MCcomb_MaxProbProbNoCorr_bCorr"] = new TH1F("mTop_MCcomb_MaxProbProbNoCorr_bCorr","mTop_MCcomb_MaxProbProbNoCorr_bCorr;mTop;#events",75,100,250);
 	histo1D["mTop_MCcomb_MaxProbProbNoCorr_AllCorr"] = new TH1F("mTop_MCcomb_MaxProbProbNoCorr_AllCorr","mTop_MCcomb_MaxProbProbNoCorr_AllCorr;mTop;#events",75,100,250);
 	histo1D["mTop_MCcomb_MaxProbProbNoCorr_JECbCorr"] = new TH1F("mTop_MCcomb_MaxProbProbNoCorr_JECbCorr","mTop_MCcomb_MaxProbProbNoCorr_JECbCorr;mTop;#events",75,100,250);
 	histo1D["mTop_MCcomb_MaxProbProbNoCorr_JECAllCorr"] = new TH1F("mTop_MCcomb_MaxProbProbNoCorr_JECAllCorr","mTop_MCcomb_MaxProbProbNoCorr_JECAllCorr;mTop;#events",75,100,250);
  
	histo1D["pT_bjet_aftermva_goodcomb"] = new TH1F("pT_bjet_aftermva_goodcomb","pT bjets good comb;p_{T}^b;#events",150,0,300);
	histo1D["pT_ljet_aftermva_goodcomb"] = new TH1F("pT_ljet_aftermva_goodcomb","pT ljets good comb;p_{T}^l;#events",150,0,300);
	histo1D["mW_aftermva_goodcomb"] = new TH1F("mW_aftermva_goodcomb","W boson mass good comb;m_{W};#events",100,0,200);
	histo1D["mtop_aftermva_goodcomb"] = new TH1F("mtop_aftermva_goodcomb","top quark mass good comb;m_{t};#events",75,100,250);
	histo1D["mtop_afterMaxProb_goodcomb"] = new TH1F("mtop_afterMaxProb_goodcomb","top quark mass good comb;m_{t};#events",75,100,250);
	histo1D["mvavalue_aftermva_goodcomb"] = new TH1F("mvavalue_aftermva_goodcomb","mva value good comb;mva-value;#events",100,0,1);
	histo1D["mvavalue_afterMaxProb_goodcomb"] = new TH1F("mvavalue_afterMaxProb_goodcomb","mva value good comb after MaxProb;mva-value;#events",100,0,1);	
	histo1D["mvavalue_afterprobcuts_goodcomb"] = new TH1F("mvavalue_afterprobcuts_goodcomb","mva value good comb after probcuts;mva-value;#events",100,0,1);
	histo1D["njets_aftermva_goodcomb"] = new TH1F("njets_aftermva_goodcomb","njets_aftermva_goodcomb;# selected jets;#events",8,2.5,10.5);		
	histo2D["mWVSmtop_aftermva_goodcomb"] = new TH2F("mWVSmtop_aftermva_goodcomb","mWVSmtop_aftermva_goodcomb",100,0,200,200,0,400);
	histo1D["Mlb_aftermva_goodcomb"] = new TH1F("Mlb_aftermva_goodcomb","Mlb_aftermva_goodcomb;m_{l,b};#events",100,0,400);
	histo1D["Mlb_afterMaxProb_goodcomb"] = new TH1F("Mlb_afterMaxProb_goodcomb","Mlb_afterMaxProb_goodcomb;m_{l,b};#events",100,0,400);
	histo1D["ProbNoCorr_afterMaxProb_goodcomb"] = new TH1F("ProbNoCorr_afterMaxProb_goodcomb","ProbNoCorr_afterMaxProb_goodcomb;ProbNoCorr;#events",100,0,1);
	histo1D["Mlb_afterprobcuts_goodcomb"] = new TH1F("Mlb_afterprobcuts_goodcomb","Mlb_afterprobcuts_goodcomb;m_{l,b};#events",100,0,400);
	histo1D["MaxProb_afterprobcuts_goodcomb"] = new TH1F("MaxProb_afterprobcuts_goodcomb","MaxProb_afterprobcuts_goodcomb;MaxProb;#events",100,0.98,1.);
	histo1D["ProbNoCorr_afterprobcuts_goodcomb"] = new TH1F("ProbNoCorr_afterprobcuts_goodcomb","ProbNoCorr_afterprobcuts_goodcomb;ProbNoCorr;#events",100,0,1);
	histo1D["nHoles_afterMaxProb_goodcomb"] = new TH1F("nHoles_afterMaxProb_goodcomb","nHoles_afterMaxProb_goodcomb;nHoles;#events",80,-0.5,79.5);
	histo1D["nHoles_afterProbNoCorr_goodcomb"] = new TH1F("nHoles_afterProbNoCorr_goodcomb","nHoles_afterProbNoCorr_goodcomb;nHoles;#events",80,-0.5,79.5);
	histo1D["mtop_afterProbNoCorr_goodcomb"] = new TH1F("mtop_afterProbNoCorr_goodcomb","mtop_afterProbNoCorr_goodcomb;mTop;#events",75,100,250);
  histo1D["mtop_afterProbNoCorr_TCHPT_goodcomb"] = new TH1F("mtop_afterProbNoCorr_TCHPT_goodcomb","mtop_afterProbNoCorr_TCHPT_goodcomb;mTop;#events",75,100,250);
  histo1D["mtop_afterProbNoCorr_SSVHPT_goodcomb"] = new TH1F("mtop_afterProbNoCorr_SSVHPT_goodcomb","mtop_afterProbNoCorr_SSVHPT_goodcomb;mTop;#events",75,100,250);
  histo1D["MET_afterProbNoCorr_goodcomb"] = new TH1F("MET_afterProbNoCorr_goodcomb","MET_afterProbNoCorr_goodcomb;MET;#events",75,0,300);
  histo1D["MHT_afterProbNoCorr_goodcomb"] = new TH1F("MHT_afterProbNoCorr_goodcomb","MHT_afterProbNoCorr_goodcomb;MHT;#events",75,0,300);
  histo1D["TCHP_hadrBJet_afterProbNoCorr_goodcomb"] = new TH1F("TCHP_hadrBJet_afterProbNoCorr_goodcomb","TCHP_hadrBJet_afterProbNoCorr_goodcomb;TCHP hadr b-jet;#events",50,-5,20);
  histo1D["TCHP_lepBJet_afterProbNoCorr_goodcomb"] = new TH1F("TCHP_lepBJet_afterProbNoCorr_goodcomb","TCHP_lepBJet_afterProbNoCorr_goodcomb;TCHP lep b-jet;#events",50,-5,20);
  histo1D["nJets_afterProbNoCorr_goodcomb"] = new TH1F("nJets_afterProbNoCorr_goodcomb","nJets_afterProbNoCorr_goodcomb;#jets,#events",6,3.5,9.5);
  histo1D["DPhi_MET_MHT_afterProbNoCorr_goodcomb"] = new TH1F("DPhi_MET_MHT_afterProbNoCorr_goodcomb","DPhi_MET_MHT_afterProbNoCorr_goodcomb;#Delta#Phi;#events",60,0,6);
  
	histo1D["pT_bjet_aftermva_badcomb"] = new TH1F("pT_bjet_aftermva_badcomb","pT bjets bad comb;p_{T}^b;#events",150,0,300);
	histo1D["pT_ljet_aftermva_badcomb"] = new TH1F("pT_ljet_aftermva_badcomb","pT ljets bad comb;p_{T}^l;#events",150,0,300);
	histo1D["mW_aftermva_badcomb"] = new TH1F("mW_aftermva_badcomb","W boson mass bad comb;m_{W};#events",100,0,200);
	histo1D["mtop_aftermva_badcomb"] = new TH1F("mtop_aftermva_badcomb","top quark mass bad comb;m_{t};#events",75,100,250);
	histo1D["mtop_afterMaxProb_badcomb"] = new TH1F("mtop_afterMaxProb_badcomb","top quark mass bad comb;m_{t};#events",75,100,250);
	histo1D["mvavalue_aftermva_badcomb"] = new TH1F("mvavalue_aftermva_badcomb","mva value bad comb;mva-value;#events",100,0,1);
	histo1D["mvavalue_afterMaxProb_badcomb"] = new TH1F("mvavalue_afterMaxProb_badcomb","mva value bad comb after MaxProb;mva-value;#events",100,0,1);	
	histo1D["mvavalue_afterprobcuts_badcomb"] = new TH1F("mvavalue_afterprobcuts_badcomb","mva value bad comb after probcuts;mva-value;#events",100,0,1);
	histo1D["njets_aftermva_badcomb"] = new TH1F("njets_aftermva_badcomb","njets_aftermva_badcomb;# selected jets;#events",8,2.5,10.5);		
	histo2D["mWVSmtop_aftermva_badcomb"] = new TH2F("mWVSmtop_aftermva_badcomb","mWVSmtop_aftermva_badcomb",100,0,200,200,0,400);	
	histo1D["Mlb_aftermva_badcomb"] = new TH1F("Mlb_aftermva_badcomb","Mlb_aftermva_badcomb;m_{l,b};#events",100,0,400);
	histo1D["Mlb_afterMaxProb_badcomb"] = new TH1F("Mlb_afterMaxProb_badcomb","Mlb_afterMaxProb_badcomb;m_{l,b};#events",100,0,400);
	histo1D["ProbNoCorr_afterMaxProb_badcomb"] = new TH1F("ProbNoCorr_afterMaxProb_badcomb","ProbNoCorr_afterMaxProb_badcomb;ProbNoCorr;#events",100,0,1);
	histo1D["Mlb_afterprobcuts_badcomb"] = new TH1F("Mlb_afterprobcuts_badcomb","Mlb_afterprobcuts_badcomb;m_{l,b};#events",100,0,400);
	histo1D["MaxProb_afterprobcuts_badcomb"] = new TH1F("MaxProb_afterprobcuts_badcomb","MaxProb_afterprobcuts_badcomb;MaxProb;#events",100,0.98,1.);
	histo1D["ProbNoCorr_afterprobcuts_badcomb"] = new TH1F("ProbNoCorr_afterprobcuts_badcomb","ProbNoCorr_afterprobcuts_badcomb;ProbNoCorr;#events",100,0,1);
	histo1D["nHoles_afterMaxProb_badcomb"] = new TH1F("nHoles_afterMaxProb_badcomb","nHoles_afterMaxProb_badcomb;nHoles;#events",80,-0.5,79.5);
	histo1D["nHoles_afterProbNoCorr_badcomb"] = new TH1F("nHoles_afterProbNoCorr_badcomb","nHoles_afterProbNoCorr_badcomb;nHoles;#events",80,-0.5,79.5);
	histo1D["mtop_afterProbNoCorr_badcomb"] = new TH1F("mtop_afterProbNoCorr_badcomb","mtop_afterProbNoCorr_badcomb;mTop;#events",75,100,250);
  histo1D["mtop_afterProbNoCorr_TCHPT_badcomb"] = new TH1F("mtop_afterProbNoCorr_TCHPT_badcomb","mtop_afterProbNoCorr_TCHPT_badcomb;mTop;#events",75,100,250);
  histo1D["mtop_afterProbNoCorr_SSVHPT_badcomb"] = new TH1F("mtop_afterProbNoCorr_SSVHPT_badcomb","mtop_afterProbNoCorr_SSVHPT_badcomb;mTop;#events",75,100,250);
  histo1D["MET_afterProbNoCorr_badcomb"] = new TH1F("MET_afterProbNoCorr_badcomb","MET_afterProbNoCorr_badcomb;MET;#events",75,0,300);
  histo1D["MHT_afterProbNoCorr_badcomb"] = new TH1F("MHT_afterProbNoCorr_badcomb","MHT_afterProbNoCorr_badcomb;MHT;#events",75,0,300);  
  histo1D["TCHP_hadrBJet_afterProbNoCorr_badcomb"] = new TH1F("TCHP_hadrBJet_afterProbNoCorr_badcomb","TCHP_hadrBJet_afterProbNoCorr_badcomb;TCHP hadr b-jet;#events",50,-5,20);
  histo1D["TCHP_lepBJet_afterProbNoCorr_badcomb"] = new TH1F("TCHP_lepBJet_afterProbNoCorr_badcomb","TCHP_lepBJet_afterProbNoCorr_badcomb;TCHP lep b-jet;#events",50,-5,20);
  histo1D["nJets_afterProbNoCorr_badcomb"] = new TH1F("nJets_afterProbNoCorr_badcomb","nJets_afterProbNoCorr_badcomb;#jets,#events",6,3.5,9.5);
  histo1D["DPhi_MET_MHT_afterProbNoCorr_badcomb"] = new TH1F("DPhi_MET_MHT_afterProbNoCorr_badcomb","DPhi_MET_MHT_afterProbNoCorr_badcomb;#Delta#Phi;#events",60,0,6);
  
	histo1D["mvavalue_afterprobcuts_goodcomb_template"] = new TH1F("mvavalue_afterprobcuts_goodcomb_template","mva value good comb after probcuts;mva-value;#events",50,0,1);
	histo1D["mvavalue_afterprobcuts_badcomb_template"] = new TH1F("mvavalue_afterprobcuts_badcomb_template","mva value bad comb after probcuts;mva-value;#events",50,0,1);	
	histo1D["mvavalue_afterprobcuts_allcomb_MC"] = new TH1F("mvavalue_afterprobcuts_allcomb_MC","mva value all comb after probcuts;mva-value;#events",50,0,1);
	
	histo1D["mvavalue_afterMaxProb_all4JetsMCMatched_goodcomb"] = new TH1F("mvavalue_afterMaxProb_all4JetsMCMatched_goodcomb","mvavalue_afterMaxProb_all4JetsMCMatched_goodcomb;maxMVA;#events",50,0,1);
	histo1D["mvavalue_afterMaxProb_allHadrJetsMCMatched_goodcomb"] = new TH1F("mvavalue_afterMaxProb_allHadrJetsMCMatched_goodcomb","mvavalue_afterMaxProb_allHadrJetsMCMatched_goodcomb;maxMVA;#events",50,0,1);
	histo1D["mvavalue_afterMaxProb_NotMCMatched_goodcomb"] = new TH1F("mvavalue_afterMaxProb_NotMCMatched_goodcomb","mvavalue_afterMaxProb_NotMCMatched_goodcomb;maxMVA;#events",50,0,1);
	
	histo1D["mvavalue_afterMaxProb_all4JetsMCMatched_badcomb"] = new TH1F("mvavalue_afterMaxProb_all4JetsMCMatched_badcomb","mvavalue_afterMaxProb_all4JetsMCMatched_badcomb;maxMVA;#events",50,0,1);
	histo1D["mvavalue_afterMaxProb_allHadrJetsMCMatched_badcomb"] = new TH1F("mvavalue_afterMaxProb_allHadrJetsMCMatched_badcomb","mvavalue_afterMaxProb_allHadrJetsMCMatched_badcomb;maxMVA;#events",50,0,1);
	histo1D["mvavalue_afterMaxProb_NotMCMatched_badcomb"] = new TH1F("mvavalue_afterMaxProb_NotMCMatched_badcomb","mvavalue_afterMaxProb_NotMCMatched_badcomb;maxMVA;#events",50,0,1);
	
	for(int iCombi=0; iCombi<12; iCombi++)
	{
	  stringstream ss1;
    ss1 << iCombi;
    string tmpName = "mvavalue_afterMaxProb_comb" + ss1.str();
    histo1D[tmpName] = new TH1F(tmpName.c_str(),(tmpName+";mva-value;#events").c_str(),50,0,1);
    histo1D[tmpName+"_goodcomb"] = new TH1F((tmpName+"_goodcomb").c_str(),(tmpName+"_goodcomb;mva-value;#events").c_str(),50,0,1);
    histo1D[tmpName+"_badcomb"] = new TH1F((tmpName+"_badcomb").c_str(),(tmpName+"_badcomb;mva-value;#events").c_str(),50,0,1);
	}
	
	histo1D["nGoodCombi_jetCombiBins"] = new TH1F("nGoodCombi_jetCombiBins","nGoodCombi_jetCombiBins;jetCombi;#goodCombi",12,-0.5,11.5);
	histo1D["nBadCombi_jetCombiBins"] = new TH1F("nBadCombi_jetCombiBins","nBadCombi_jetCombiBins;jetCombi;#badCombi",12,-0.5,11.5);
	
  // Holes stuff
  histo1D["minDistHoleMax"] = new TH1F("minDistHoleMax","minDistHoleMax",80,0,80);
  histo1D["nHoles"] = new TH1F("nHoles","nHoles",80,0,80);
  histo1D["mTopMax"] = new TH1F("mTopMax","mTopMax",200,0,400);
  histo2D["mTopMaxVSminDistHoleMax"] = new TH2F("mTopMaxVSminDistHoleMax","mTopMax VS minDistHoleMax",200,0,400,80,0,80);
  histo2D["nHolesVSminDistHoleMax"] = new TH2F("nHolesVSminDistHoleMax","nHoles VS minDistHoleMax",80,0,80,80,0,80);
  histo2D["mTopMaxVSnHoles"] = new TH2F("mTopMaxVSnHoles","mTopMax VS nHoles",200,0,400,80,0,80);
  
  histo1D["holes_PtJets"] = new TH1F("holes_PtJets","holes_PtJets",100,0,200);
  histo1D["holes_mW"] = new TH1F("holes_mW","holes_mW",100,0,200);
  histo1D["holes_mTop"] = new TH1F("holes_mTop","holes_mTop",200,0,400);
  histo1D["holes_EtaJets"] = new TH1F("holes_EtaJets","holes_EtaJets",100,-2.5,2.5);
  histo1D["holes_MVAvalue"] = new TH1F("holes_MVAvalue","holes_MVAvalue",100,0,1);
  histo1D["holes_DRlightjets"] = new TH1F("holes_DRlightjets","holes_DRlightjets",70,0,7);
  histo1D["holes_DRlightbjet"] = new TH1F("holes_DRlightbjet","holes_DRlightbjet",70,0,7);
  histo1D["noholes_PtJets"] = new TH1F("noholes_PtJets","noholes_PtJets",100,0,200);
  histo1D["noholes_mW"] = new TH1F("noholes_mW","noholes_mW",100,0,200);
  histo1D["noholes_mTop"] = new TH1F("noholes_mTop","noholes_mTop",200,0,400);
  histo1D["noholes_EtaJets"] = new TH1F("noholes_EtaJets","noholes_EtaJets",100,-2.5,2.5);
  histo1D["noholes_MVAvalue"] = new TH1F("noholes_MVAvalue","noholes_MVAvalue",100,0,1);
  histo1D["noholes_DRlightjets"] = new TH1F("noholes_DRlightjets","noholes_DRlightjets",70,0,7);
  histo1D["noholes_DRlightbjet"] = new TH1F("noholes_DRlightbjet","noholes_DRlightbjet",70,0,7);

  vector<float> topMass, topMassUnc, DEl, DElUnc; // to store the results of each Pseudo-experiment
  
  cout << "Initializing MSPlots" << endl;
  vector<Dataset*> datasets; // needed for MSPlots
  for(unsigned int iDataSet=0; iDataSet<inputMonsters.size(); iDataSet++)
  {
    TFile* inFile = new TFile(inputMonsters[iDataSet].c_str(),"READ");
    TTree* inConfigTree = (TTree*) inFile->Get("configTreeMonsterFile");
    TBranch* d_br = (TBranch*) inConfigTree->GetBranch("Dataset");
    TClonesArray* tc_dataset = new TClonesArray("Dataset",0);
    d_br->SetAddress(&tc_dataset);
    inConfigTree->GetEvent(0);
    Dataset* dataSet = (Dataset*) tc_dataset->At(0);
    datasets.push_back( (Dataset*) dataSet->Clone() );
  }
  
  //MSPlots
	MSPlot["pT_bjet_aftermva"] = new MultiSamplePlot(datasets,"pT_bjet_aftermva",150,0,300,"b Jet P_{T} after MVA");
	MSPlot["pT_ljet_aftermva"] = new MultiSamplePlot(datasets,"pT_ljet_aftermva",150,0,300,"l Jets P_{T} after MVA");
	MSPlot["mW_aftermva"] = new MultiSamplePlot(datasets,"mW_aftermva",100,0,200,"m_{W} after MVA");
	MSPlot["mtop_aftermva"] = new MultiSamplePlot(datasets,"mtop_aftermva",200,0,400,"m_{top} after MVA");
	MSPlot["mvavalue_aftermva"] = new MultiSamplePlot(datasets,"mvavalue_aftermva",100,0,1,"maxMVA value");
	MSPlot["massmub_aftermva"] = new MultiSamplePlot(datasets,"massmub_aftermva",200,0,400,"m_{b#mu} after MVA");
	MSPlot["njets_aftermva"] = new MultiSamplePlot(datasets,"njets_aftermva",8,2.5,10.5,"Number of selected jets after MVA");
	MSPlot["DR_ljets_aftermva"] = new MultiSamplePlot(datasets,"DR_ljets_aftermva",15,0,5,"#DeltaR(light jets)");

  // additional cuts or weights
  MSPlot["maxProbKinFit"] = new MultiSamplePlot(datasets, "maxProbKinFit", 100, -0.0001, 1, "max Probability of the KinFit");
  MSPlot["maxMVA"] = new MultiSamplePlot(datasets, "maxMVA",100,0,1,"maxMVA");
  MSPlot["ProbNoCorr_MaxProbCut"] = new MultiSamplePlot(datasets, "ProbNoCorr_MaxProbCut",100,0,1,"ProbNoCorr");

  MSPlot["maxMVA_HadrWMass_MaxProbCut"] = new MultiSamplePlot(datasets, "maxMVA_HadrWMass_MaxProbCut", 40, -0.0001, 200, "Hadronic W Mass, after MaxProb");
  MSPlot["maxMVA_HadrTopMass_MaxProbCut"] = new MultiSamplePlot(datasets, "maxMVA_HadrTopMass_MaxProbCut", 40, -0.0001, 400, "Hadronic Top Mass, after MaxProb");
  MSPlot["maxMVA_Mlb_MaxProbCut"] = new MultiSamplePlot(datasets, "maxMVA_Mlb_MaxProbCut", 40, -0.0001, 400, "M_{l,b}, after MaxProb");
  MSPlot["maxMVA_HadrWMass_probcuts"] = new MultiSamplePlot(datasets, "maxMVA_HadrWMass_probcuts", 40, -0.0001, 200, "Hadronic W Mass, after probcuts");
  MSPlot["maxMVA_HadrTopMass_probcuts"] = new MultiSamplePlot(datasets, "maxMVA_HadrTopMass_probcuts", 40, -0.0001, 400, "Hadronic Top Mass, after probcuts");
  MSPlot["maxMVA_Mlb_probcuts"] = new MultiSamplePlot(datasets, "maxMVA_Mlb_probcuts", 40, -0.0001, 400, "M_{l,b}, after probcuts");
  
  // mu+ and mu- differences
  MSPlot["mWHadr_MuPlus"] = new MultiSamplePlot(datasets, "mWHadr_MuPlus", 40, 0, 400, "m_{W^{-}}");
  MSPlot["mWHadr_MuMinus"] = new MultiSamplePlot(datasets, "mWHadr_MuMinus", 40, 0, 400, "m_{W^{+}}");
  MSPlot["mTopHadr_MuPlus"] = new MultiSamplePlot(datasets, "mTopHadr_MuPlus", 50, 0, 800, "m_{Top}");
  MSPlot["mTopHadr_MuMinus"] = new MultiSamplePlot(datasets, "mTopHadr_MuMinus", 50, 0, 800, "m_{Top}");
  MSPlot["mWHadr_MuPlus_MaxProb"] = new MultiSamplePlot(datasets, "mWHadr_MuPlus_MaxProb", 35, 0, 280, "m_{W^{-}}");
  MSPlot["mWHadr_MuMinus_MaxProb"] = new MultiSamplePlot(datasets, "mWHadr_MuMinus_MaxProb", 35, 0, 280, "m_{W^{+}}");
  MSPlot["mTopHadr_MuPlus_MaxProb"] = new MultiSamplePlot(datasets, "mTopHadr_MuPlus_MaxProb", 25, 0, 400, "m_{Top}");
  MSPlot["mTopHadr_MuMinus_MaxProb"] = new MultiSamplePlot(datasets, "mTopHadr_MuMinus_MaxProb", 25, 0, 400, "m_{Top}");
  MSPlot["mWHadr_MuPlus_ProbNoCorr"] = new MultiSamplePlot(datasets, "mWHadr_MuPlus_ProbNoCorr", 35, 0, 280, "m_{W^{-}}");
  MSPlot["mWHadr_MuMinus_ProbNoCorr"] = new MultiSamplePlot(datasets, "mWHadr_MuMinus_ProbNoCorr", 35, 0, 280, "m_{W^{+}}");
  MSPlot["mTopHadr_MuPlus_ProbNoCorr"] = new MultiSamplePlot(datasets, "mTopHadr_MuPlus_ProbNoCorr", 25, 0, 400, "m_{Top}");
  MSPlot["mTopHadr_MuMinus_ProbNoCorr"] = new MultiSamplePlot(datasets, "mTopHadr_MuMinus_ProbNoCorr", 25, 0, 400, "m_{Top}");
  
  cout << "Initializing BinnedMonsterCombinations" << endl;
  // initialize the monstercombination stuff
  vector<float> cutBins;
  cutBins.push_back(0);
  cutBins.push_back(1);
  cutBins.push_back(2);
  cutBins.push_back(3);
  cutBins.push_back(4);
  BinnedMonsterCombination* superMonsterBinnedCuts = new BinnedMonsterCombination("cutBins",measureTopMass,cutBins);
  BinnedMonsterCombination* superMonsterBinnedCuts_SemiMuMatched = new BinnedMonsterCombination("cutBins_SemiMuMatched",measureTopMass,cutBins);
  BinnedMonsterCombination* superMonsterBinnedCuts_Data = new BinnedMonsterCombination("cutBins_Data",measureTopMass,cutBins);
  
  vector<float> jetCombiBins;
  jetCombiBins.push_back(-0.5);
  jetCombiBins.push_back(0.5);
  jetCombiBins.push_back(1.5);
  jetCombiBins.push_back(2.5);
  jetCombiBins.push_back(3.5);
  jetCombiBins.push_back(4.5);
  jetCombiBins.push_back(5.5);
  jetCombiBins.push_back(6.5);
  jetCombiBins.push_back(7.5);
  jetCombiBins.push_back(8.5);
  jetCombiBins.push_back(9.5);
  jetCombiBins.push_back(10.5);
  jetCombiBins.push_back(11.5);
  BinnedMonsterCombination* superMonsterBinnedJetCombi = new BinnedMonsterCombination("jetCombiBins",measureTopMass,jetCombiBins);
  BinnedMonsterCombination* superMonsterBinnedJetCombi_SemiMuMatched = new BinnedMonsterCombination("jetCombiBins_SemiMuMatched",measureTopMass,jetCombiBins);
  BinnedMonsterCombination* superMonsterBinnedJetCombi_BadCombi_MaxProb = new BinnedMonsterCombination("jetCombiBins_BadCombi_MaxProb",measureTopMass,jetCombiBins);
  BinnedMonsterCombination* superMonsterBinnedJetCombi_BadCombi_MaxProbProbNoCorr = new BinnedMonsterCombination("jetCombiBins_BadCombi_MaxProbProbNoCorr",measureTopMass,jetCombiBins);
  BinnedMonsterCombination* superMonsterBinnedJetCombi_BadCombi_MaxProbProbNoCorr0p05 = new BinnedMonsterCombination("jetCombiBins_BadCombi_MaxProbProbNoCorr0p05",measureTopMass,jetCombiBins);
  BinnedMonsterCombination* superMonsterBinnedJetCombi_BadCombi_MaxProbProbNoCorr0p1 = new BinnedMonsterCombination("jetCombiBins_BadCombi_MaxProbProbNoCorr0p1",measureTopMass,jetCombiBins);
  BinnedMonsterCombination* superMonsterBinnedJetCombi_BadCombi_MaxProbProbNoCorr0p2 = new BinnedMonsterCombination("jetCombiBins_BadCombi_MaxProbProbNoCorr0p2",measureTopMass,jetCombiBins);
  BinnedMonsterCombination* superMonsterBinnedJetCombi_BadCombi_MaxProbProbNoCorrNoHoles = new BinnedMonsterCombination("jetCombiBins_BadCombi_MaxProbProbNoCorrNoHoles",measureTopMass,jetCombiBins);
  BinnedMonsterCombination* superMonsterBinnedJetCombi_BadCombi_MaxProbProbNoCorrNJets = new BinnedMonsterCombination("jetCombiBins_BadCombi_MaxProbProbNoCorrNJets",measureTopMass,jetCombiBins);
  BinnedMonsterCombination* superMonsterBinnedJetCombi_BadCombi_MaxProbProbNoCorrBtagHadrJet = new BinnedMonsterCombination("jetCombiBins_BadCombi_MaxProbProbNoCorrBtagHadrJet",measureTopMass,jetCombiBins);
  BinnedMonsterCombination* superMonsterBinnedJetCombi_BadCombi_MaxProbProbNoCorrBtagHadrLeptJet = new BinnedMonsterCombination("jetCombiBins_BadCombi_MaxProbProbNoCorrBtagHadrLeptJet",measureTopMass,jetCombiBins);
  
  BinnedMonsterCombination* superMonsterBinnedJetCombi_0W0B_MaxProb = new BinnedMonsterCombination("jetCombiBins_0W0B_MaxProb",measureTopMass,jetCombiBins);
  BinnedMonsterCombination* superMonsterBinnedJetCombi_0W1B_MaxProb = new BinnedMonsterCombination("jetCombiBins_0W1B_MaxProb",measureTopMass,jetCombiBins);
  BinnedMonsterCombination* superMonsterBinnedJetCombi_1W0B_MaxProb = new BinnedMonsterCombination("jetCombiBins_1W0B_MaxProb",measureTopMass,jetCombiBins);
  BinnedMonsterCombination* superMonsterBinnedJetCombi_1W1B_MaxProb = new BinnedMonsterCombination("jetCombiBins_1W1B_MaxProb",measureTopMass,jetCombiBins);
  BinnedMonsterCombination* superMonsterBinnedJetCombi_2W0B_MaxProb = new BinnedMonsterCombination("jetCombiBins_2W0B_MaxProb",measureTopMass,jetCombiBins);
  BinnedMonsterCombination* superMonsterBinnedJetCombi_2W1B_MaxProb = new BinnedMonsterCombination("jetCombiBins_2W1B_MaxProb",measureTopMass,jetCombiBins);
 
  BinnedMonsterCombination* superMonsterBinnedJetCombi_0W0B_MaxProbProbNoCorr = new BinnedMonsterCombination("jetCombiBins_0W0B_MaxProbProbNoCorr",measureTopMass,jetCombiBins);
  BinnedMonsterCombination* superMonsterBinnedJetCombi_0W1B_MaxProbProbNoCorr = new BinnedMonsterCombination("jetCombiBins_0W1B_MaxProbProbNoCorr",measureTopMass,jetCombiBins);
  BinnedMonsterCombination* superMonsterBinnedJetCombi_1W0B_MaxProbProbNoCorr = new BinnedMonsterCombination("jetCombiBins_1W0B_MaxProbProbNoCorr",measureTopMass,jetCombiBins);
  BinnedMonsterCombination* superMonsterBinnedJetCombi_1W1B_MaxProbProbNoCorr = new BinnedMonsterCombination("jetCombiBins_1W1B_MaxProbProbNoCorr",measureTopMass,jetCombiBins);
  BinnedMonsterCombination* superMonsterBinnedJetCombi_2W0B_MaxProbProbNoCorr = new BinnedMonsterCombination("jetCombiBins_2W0B_MaxProbProbNoCorr",measureTopMass,jetCombiBins);
  BinnedMonsterCombination* superMonsterBinnedJetCombi_2W1B_MaxProbProbNoCorr = new BinnedMonsterCombination("jetCombiBins_2W1B_MaxProbProbNoCorr",measureTopMass,jetCombiBins);
  
  vector<float> inclJetCombiBins;
  inclJetCombiBins.push_back(-0.5);
  inclJetCombiBins.push_back(5.5);
  inclJetCombiBins.push_back(6.5);
  inclJetCombiBins.push_back(11.5);
  BinnedMonsterCombination* superMonsterBinnedInclJetCombi = new BinnedMonsterCombination("inclJetCombiBins",measureTopMass,inclJetCombiBins);
  BinnedMonsterCombination* superMonsterBinnedInclJetCombi_SemiMuMatched = new BinnedMonsterCombination("inclJetCombiBins_SemiMuMatched",measureTopMass,inclJetCombiBins);
  BinnedMonsterCombination* superMonsterBinnedInclJetCombi_BadCombi = new BinnedMonsterCombination("inclJetCombiBins_BadCombi",measureTopMass,inclJetCombiBins);
  
  BinnedMonsterCombination* superMonsterBinnedInclJetCombi_BadCombiSubtr = new BinnedMonsterCombination("inclJetCombiBins_BadCombiSubtr",measureTopMass,inclJetCombiBins);
  BinnedMonsterCombination* superMonsterBinnedInclJetCombi_AllCombiSubtr = new BinnedMonsterCombination("inclJetCombiBins_AllCombiSubtr",measureTopMass,inclJetCombiBins);
  
  vector<float> pTbins;
  pTbins.push_back(30);
  pTbins.push_back(50);
  pTbins.push_back(70);
  pTbins.push_back(90);
  pTbins.push_back(120);
  pTbins.push_back(250);
  BinnedMonsterCombination* superMonsterBinnedPtL = new BinnedMonsterCombination("PtLbins",measureTopMass,pTbins);
  BinnedMonsterCombination* superMonsterBinnedPtB = new BinnedMonsterCombination("PtBbins",measureTopMass,pTbins);
  BinnedMonsterCombination* superMonsterBinnedPtL_Data = new BinnedMonsterCombination("PtLbins_Data",measureTopMass,pTbins);
  BinnedMonsterCombination* superMonsterBinnedPtB_Data = new BinnedMonsterCombination("PtBbins_Data",measureTopMass,pTbins);
  BinnedMonsterCombination* superMonsterBinnedPtL_SemiMuMatched = new BinnedMonsterCombination("PtLbins_SemiMuMatched",measureTopMass,pTbins);
  BinnedMonsterCombination* superMonsterBinnedPtB_SemiMuMatched = new BinnedMonsterCombination("PtBbins_SemiMuMatched",measureTopMass,pTbins);
  BinnedMonsterCombination* superMonsterBinnedPtL_BadCombi = new BinnedMonsterCombination("PtLbins_BadCombi",measureTopMass,pTbins);
  BinnedMonsterCombination* superMonsterBinnedPtB_BadCombi = new BinnedMonsterCombination("PtBbins_BadCombi",measureTopMass,pTbins);
  
  vector<float> chargeBins;
  chargeBins.push_back(-2);
  chargeBins.push_back(0);
  chargeBins.push_back(2);
  BinnedMonsterCombination* superMonsterBinnedCharge = new BinnedMonsterCombination("chargeBins",measureTopMass,chargeBins);
  BinnedMonsterCombination* superMonsterBinnedCharge_SemiMuMatched = new BinnedMonsterCombination("chargeBins_SemiMuMatched",measureTopMass,chargeBins);
  BinnedMonsterCombination* superMonsterBinnedCharge_Data = new BinnedMonsterCombination("chargeBins_Data",measureTopMass,chargeBins);
  
  vector<float> maxMVAbins;
  maxMVAbins.push_back(0);
  maxMVAbins.push_back(0.6);
  maxMVAbins.push_back(0.8);
  maxMVAbins.push_back(0.9);
  maxMVAbins.push_back(0.93);
  maxMVAbins.push_back(0.96);
  maxMVAbins.push_back(1.);
  BinnedMonsterCombination* superMonsterBinnedMaxMVA = new BinnedMonsterCombination("maxMVAbins",measureTopMass,maxMVAbins);
  BinnedMonsterCombination* superMonsterBinnedMaxMVA_SemiMuMatched = new BinnedMonsterCombination("maxMVAbins_SemiMuMatched",measureTopMass,maxMVAbins);
  BinnedMonsterCombination* superMonsterBinnedMaxMVA_Data = new BinnedMonsterCombination("maxMVAbins_Data",measureTopMass,maxMVAbins);
  
  MVABinnedMonsterCombination* mvaBinSuperMonsterMaxProb = new MVABinnedMonsterCombination("mvaBinned_MaxProb",measureTopMass,maxMVAbins);
  MVABinnedMonsterCombination* mvaBinSuperMonsterProbNoCorr = new MVABinnedMonsterCombination("mvaBinned_ProbNoCorr",measureTopMass,maxMVAbins);
  MVABinnedMonsterCombination* mvaBinSuperMonsterMaxProb_SemiMuMatched = new MVABinnedMonsterCombination("mvaBinned_MaxProb_SemiMuMatched",measureTopMass,maxMVAbins);
  MVABinnedMonsterCombination* mvaBinSuperMonsterProbNoCorr_SemiMuMatched = new MVABinnedMonsterCombination("mvaBinned_ProbNoCorr_SemiMuMatched",measureTopMass,maxMVAbins);
  MVABinnedMonsterCombination* mvaBinSuperMonsterMaxProb_BadCombi = new MVABinnedMonsterCombination("mvaBinned_MaxProb_BadCombi",measureTopMass,maxMVAbins);
  MVABinnedMonsterCombination* mvaBinSuperMonsterProbNoCorr_BadCombi = new MVABinnedMonsterCombination("mvaBinned_ProbNoCorr_BadCombi",measureTopMass,maxMVAbins);
  
  // check input/output of MVA
  vector<float> MVAbins;
  MVAbins.push_back(0);
  MVAbins.push_back(0.1);
  MVAbins.push_back(0.2);
  MVAbins.push_back(0.3);
  MVAbins.push_back(0.4);
  MVAbins.push_back(0.5);
  MVAbins.push_back(0.6);
  MVAbins.push_back(0.7);
  MVAbins.push_back(0.8);
  MVAbins.push_back(0.9);
  MVAbins.push_back(1.);
  BinnedMonsterCombination* superMonsterBinnedMVA_AllGoodCombi = new BinnedMonsterCombination("MVAbins_AllGoodCombi",measureTopMass,MVAbins);
  BinnedMonsterCombination* superMonsterBinnedMVA_AllBadCombi = new BinnedMonsterCombination("MVAbins_AllBadCombi",measureTopMass,MVAbins);
  
  BinnedMonsterCombination* superMonsterBinnedBtag_AllGoodCombi = new BinnedMonsterCombination("btagBins_AllGoodCombi",measureTopMass,MVAbins);
  BinnedMonsterCombination* superMonsterBinnedBtag_AllBadCombi = new BinnedMonsterCombination("btagBins_AllBadCombi",measureTopMass,MVAbins);
  
  BinnedMonsterCombination* superMonsterBinnedProbNoCorr_AllGoodCombi = new BinnedMonsterCombination("ProbNoCorrbins_AllGoodCombi",measureTopMass,MVAbins);
  BinnedMonsterCombination* superMonsterBinnedProbNoCorr_AllBadCombi = new BinnedMonsterCombination("ProbNoCorrbins_AllBadCombi",measureTopMass,MVAbins);
  
  vector<float> MaxProbbins;
  MaxProbbins.push_back(0.98);
  MaxProbbins.push_back(0.9825);
  MaxProbbins.push_back(0.985);
  MaxProbbins.push_back(0.9875);
  MaxProbbins.push_back(0.99);
  MaxProbbins.push_back(0.9925);
  MaxProbbins.push_back(0.995);
  MaxProbbins.push_back(0.9975);
  MaxProbbins.push_back(1.);
  BinnedMonsterCombination* superMonsterBinnedMaxProb_AllGoodCombi = new BinnedMonsterCombination("MaxProbbins_AllGoodCombi",measureTopMass,MaxProbbins);
  BinnedMonsterCombination* superMonsterBinnedMaxProb_AllBadCombi = new BinnedMonsterCombination("MaxProbbins_AllBadCombi",measureTopMass,MaxProbbins);
  
  vector<float> ThPtSumPtBins;
  ThPtSumPtBins.push_back(0);
  ThPtSumPtBins.push_back(0.08);
  ThPtSumPtBins.push_back(0.16);
  ThPtSumPtBins.push_back(0.24);
  ThPtSumPtBins.push_back(0.32);
  ThPtSumPtBins.push_back(0.4);
  ThPtSumPtBins.push_back(0.48);
  ThPtSumPtBins.push_back(0.56);
  ThPtSumPtBins.push_back(0.64);
  ThPtSumPtBins.push_back(0.72);
  ThPtSumPtBins.push_back(0.8);
  BinnedMonsterCombination* superMonsterBinnedThPtSumPt_AllGoodCombi = new BinnedMonsterCombination("ThPtSumPtBins_AllGoodCombi",measureTopMass,ThPtSumPtBins);
  BinnedMonsterCombination* superMonsterBinnedThPtSumPt_AllBadCombi = new BinnedMonsterCombination("ThPtSumPtBins_AllBadCombi",measureTopMass,ThPtSumPtBins);

  vector<float> angleBins;
  angleBins.push_back(0);
  angleBins.push_back(0.32);
  angleBins.push_back(0.64);
  angleBins.push_back(0.96);
  angleBins.push_back(1.28);
  angleBins.push_back(1.6);
  angleBins.push_back(1.92);
  angleBins.push_back(2.24);
  angleBins.push_back(2.56);
  angleBins.push_back(2.88);
  angleBins.push_back(3.2);
  BinnedMonsterCombination* superMonsterBinnedAngleBlMu_AllGoodCombi = new BinnedMonsterCombination("AngleBlMuBins_AllGoodCombi",measureTopMass,angleBins);
  BinnedMonsterCombination* superMonsterBinnedAngleBlMu_AllBadCombi = new BinnedMonsterCombination("AngleBlMuBins_AllBadCombi",measureTopMass,angleBins);

  BinnedMonsterCombination* superMonsterBinnedAngleThBl_AllGoodCombi = new BinnedMonsterCombination("AngleThBlBins_AllGoodCombi",measureTopMass,angleBins);
  BinnedMonsterCombination* superMonsterBinnedAngleThBl_AllBadCombi = new BinnedMonsterCombination("AngleThBlBins_AllBadCombi",measureTopMass,angleBins);
  
  BinnedMonsterCombination* superMonsterBinnedAngleThMu_AllGoodCombi = new BinnedMonsterCombination("AngleThMuBins_AllGoodCombi",measureTopMass,angleBins);
  BinnedMonsterCombination* superMonsterBinnedAngleThMu_AllBadCombi = new BinnedMonsterCombination("AngleThMuBins_AllBadCombi",measureTopMass,angleBins);
  
  vector<float> mLBbins;
  mLBbins.push_back(0);
  mLBbins.push_back(50);
  mLBbins.push_back(70);
  mLBbins.push_back(90);
  mLBbins.push_back(110);
  mLBbins.push_back(130);
  mLBbins.push_back(150);
  mLBbins.push_back(200);
  BinnedMonsterCombination* superMonsterBinnedMlb_AllGoodCombi = new BinnedMonsterCombination("mLBbins_AllGoodCombi",measureTopMass,mLBbins);
  BinnedMonsterCombination* superMonsterBinnedMlb_AllBadCombi = new BinnedMonsterCombination("mLBbins_AllBadCombi",measureTopMass,mLBbins);
  
  string MVAmethod = "Likelihood"; // MVAmethod to be used to get the good jet combi calculation
  LightJetCombiner* jetCombiner = new LightJetCombiner(false, Luminosity, datasets, measureTopMass, MVAmethod);
  cout << " - LightJetCombiner instantiated ..." << endl;
  
  for(unsigned int iDataSet=0; iDataSet<inputMonsters.size(); iDataSet++)
  {
    string dataSetName = datasets[iDataSet]->Name();
    cout << "Processing DataSet: " << dataSetName << endl;

    // Plots which need to be made for each dataset
    histo2D[dataSetName+"_maxProbVSmaxMVA"] = new TH2F((dataSetName+"_maxProbVSmaxMVA").c_str(),(dataSetName+"_maxProbVSmaxMVA").c_str(),100,-0.0001,1,100,-0.0001,1);
    
    TFile* inFile = new TFile(inputMonsters[iDataSet].c_str(),"READ");
    
    TTree* inMonstersTree = (TTree*) inFile->Get("MonsterTree");
    TBranch* m_br = (TBranch*) inMonstersTree->GetBranch("TheMonster");
    
    Monster* monster = 0;
    m_br->SetAddress(&monster);
    
    int nEvent = inMonstersTree->GetEntries();
    
    TTree* inConfigTree = (TTree*) inFile->Get("configTreeMonsterFile");
    TBranch* d_br = (TBranch*) inConfigTree->GetBranch("Dataset");
    TClonesArray* tc_dataset = new TClonesArray("Dataset",0);
    d_br->SetAddress(&tc_dataset);
    
    inConfigTree->GetEvent(0);
    Dataset* dataSet = (Dataset*) tc_dataset->At(0);
    
    for(unsigned int iEvt=0; iEvt<nEvent; iEvt++)
//    for(unsigned int iEvt=0; iEvt<3000; iEvt++)
    {
      inMonstersTree->GetEvent(iEvt);
//      cout << "event: " << iEvt << endl;
      if(iEvt%1000 == 0)
        std::cout<<"Processing the "<<iEvt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC <<endl;
      
//      if( ! ( monster->eventID() == 3150806 && monster->runID() == 1 && monster->lumiBlockID() == 7 ) )
//        continue;
      
      jetCombiner->ProcessEvent(dataSet, monster);
      
      float MaxProb = maxProb(monster);
      float mvaVal = monster->mvaVal(0);
      unsigned int* mvaResult = monster->mvaResult(0);
      bool SemiMuMatched = hadrJetsMVAMatched(monster);
      
      histo2D[dataSetName+"_maxProbVSmaxMVA"]->Fill(MaxProb, mvaVal);
      MSPlot["maxMVA"]->Fill(mvaVal, dataSet, true, Luminosity*monster->eventWeight());
      
/*      superMonsterBinnedCuts->FillExpected(monster, 0.5, dataSet->NormFactor()*Luminosity*monster->eventWeight());
      superMonsterBinnedCuts->FillExpected(monster, 1.5, dataSet->NormFactor()*Luminosity*monster->eventWeight());
      superMonsterBinnedCuts->FillExpected(monster, 2.5, dataSet->NormFactor()*Luminosity*monster->eventWeight());
      superMonsterBinnedCuts->FillExpected(monster, 3.5, dataSet->NormFactor()*Luminosity*monster->eventWeight());
      superMonsterBinnedPtL->FillExpected(monster, monster->hadrLJet1().Pt(), dataSet->NormFactor()*Luminosity*monster->eventWeight());
      superMonsterBinnedPtL->FillExpected(monster, monster->hadrLJet2().Pt(), dataSet->NormFactor()*Luminosity*monster->eventWeight());
      superMonsterBinnedPtB->FillExpected(monster, monster->hadrBJet().Pt(), dataSet->NormFactor()*Luminosity*monster->eventWeight());
      superMonsterBinnedCharge->FillExpected(monster, monster->muonCharge(), dataSet->NormFactor()*Luminosity*monster->eventWeight());
      superMonsterBinnedMaxMVA->FillExpected(monster, mvaVal.first, dataSet->NormFactor()*Luminosity*monster->eventWeight());
      
      if(monster->allHadronicJetsMVAMatched())
      {
        superMonsterBinnedCuts_SemiMuMatched->FillExpected(monster, 0.5, dataSet->NormFactor()*Luminosity*monster->eventWeight());
        superMonsterBinnedCuts_SemiMuMatched->FillExpected(monster, 1.5, dataSet->NormFactor()*Luminosity*monster->eventWeight());
        superMonsterBinnedCuts_SemiMuMatched->FillExpected(monster, 2.5, dataSet->NormFactor()*Luminosity*monster->eventWeight());
        superMonsterBinnedCuts_SemiMuMatched->FillExpected(monster, 3.5, dataSet->NormFactor()*Luminosity*monster->eventWeight());
        superMonsterBinnedPtL_SemiMuMatched->FillExpected(monster, monster->hadrLJet1().Pt(), dataSet->NormFactor()*Luminosity*monster->eventWeight());
        superMonsterBinnedPtL_SemiMuMatched->FillExpected(monster, monster->hadrLJet2().Pt(), dataSet->NormFactor()*Luminosity*monster->eventWeight());
        superMonsterBinnedPtB_SemiMuMatched->FillExpected(monster, monster->hadrBJet().Pt(), dataSet->NormFactor()*Luminosity*monster->eventWeight());
        superMonsterBinnedMaxMVA_SemiMuMatched->FillExpected(monster, mvaVal.first, dataSet->NormFactor()*Luminosity*monster->eventWeight());
        superMonsterBinnedCharge_SemiMuMatched->FillExpected(monster, monster->muonCharge(), dataSet->NormFactor()*Luminosity*monster->eventWeight());
      }
      else
      {
        superMonsterBinnedPtL_BadCombi->FillExpected(monster, monster->hadrLJet1().Pt(), dataSet->NormFactor()*Luminosity*monster->eventWeight());
        superMonsterBinnedPtL_BadCombi->FillExpected(monster, monster->hadrLJet2().Pt(), dataSet->NormFactor()*Luminosity*monster->eventWeight());
        superMonsterBinnedPtB_BadCombi->FillExpected(monster, monster->hadrBJet().Pt(), dataSet->NormFactor()*Luminosity*monster->eventWeight());
      }*/
      
      if( mvaVal >= 0. )
      {
        if( monster->hadrBQuark().Pt() > 0.1 && monster->hadrLQuark1().Pt() > 0.1 && monster->hadrLQuark2().Pt() > 0.1 )
        {
          TLorentzVector lightJet1 = monster->selectedJet(monster->hadrLJet1());
          TLorentzVector lightJet2 = monster->selectedJet(monster->hadrLJet2());
          TLorentzVector bJet = monster->selectedJet(monster->hadrBJet());
          float lightCorr = expDEl2D->GetBinContent( expDEl2D->FindBin( ( lightJet1.Pt() + lightJet2.Pt() )/2, bJet.Pt() ) );
          float bCorr = expDEb2D->GetBinContent( expDEb2D->FindBin( ( lightJet1.Pt() + lightJet2.Pt() )/2, bJet.Pt() ) );

          float mTopAllCorr = ( lightJet1*(1+lightCorr) + lightJet2*(1+lightCorr) + bJet*(1+bCorr) ).M();
          float mTopbCorr = ( lightJet1 + lightJet2 + bJet*(1+bCorr) ).M();
        	
       	  JEClight->setJetEta(lightJet1.Eta());
          JEClight->setJetPt(lightJet1.Pt());
          float corrJEClight1 = JEClight->getCorrection();
          JEClight->setJetEta(lightJet2.Eta());
          JEClight->setJetPt(lightJet2.Pt());
          float corrJEClight2 = JEClight->getCorrection();
          
          JECb->setJetEta(bJet.Eta());
          JECb->setJetPt(bJet.Pt());
          float corrJECb = JECb->getCorrection();
          
          float mTopJECAllCorr = ( lightJet1*corrJEClight1 + lightJet2*corrJEClight2 + bJet*corrJECb ).M();
          float mTopJECbCorr = ( lightJet1 + lightJet2 + bJet*corrJECb ).M();
          
          histo1D["mW_gen"]->Fill( ( monster->hadrLQuark1() + monster->hadrLQuark2() ).M() );
	        histo1D["mTop_gen"]->Fill( ( monster->hadrLQuark1() + monster->hadrLQuark2() + monster->hadrBQuark() ).M() );
	        
         	histo1D["mW_MCcomb"]->Fill( ( lightJet1 + lightJet2 ).M() );
        	histo1D["mTop_MCcomb"]->Fill( ( lightJet1 + lightJet2 + bJet ).M() );
        	
       	 	histo1D["mTop_MCcomb_bCorr"]->Fill( mTopbCorr );
         	histo1D["mTop_MCcomb_AllCorr"]->Fill( mTopAllCorr );
         	histo1D["mTop_MCcomb_JECbCorr"]->Fill( mTopJECbCorr );
         	histo1D["mTop_MCcomb_JECAllCorr"]->Fill( mTopJECAllCorr );
        }
        
        TLorentzVector Wmass = monster->selectedJet(mvaResult[0])+monster->selectedJet(mvaResult[1]);
				TLorentzVector topmass = Wmass+monster->selectedJet(mvaResult[2]);
				TLorentzVector mub = monster->selectedJet(mvaResult[3]) + monster->muon();
				
				if(monster->muonCharge() == 1)
				{
			    MSPlot["mWHadr_MuPlus"]->Fill(Wmass.M(), dataSet, true, Luminosity*monster->eventWeight());
			    MSPlot["mTopHadr_MuPlus"]->Fill(topmass.M(), dataSet, true, Luminosity*monster->eventWeight());
			  }
				else if(monster->muonCharge() == -1)
				{
          MSPlot["mWHadr_MuMinus"]->Fill(Wmass.M(), dataSet, true, Luminosity*monster->eventWeight());
          MSPlot["mTopHadr_MuMinus"]->Fill(topmass.M(), dataSet, true, Luminosity*monster->eventWeight());
				}
        
				if(dataSetName.find("TTbarJets_SemiMu") == 0 && SemiMuMatched)
				{
      		histo1D["pT_bjet_aftermva_goodcomb"]->Fill(monster->selectedJet(mvaResult[2]).Pt());
      		histo1D["pT_ljet_aftermva_goodcomb"]->Fill(monster->selectedJet(mvaResult[0]).Pt());
     			histo1D["pT_ljet_aftermva_goodcomb"]->Fill(monster->selectedJet(mvaResult[1]).Pt());
					histo1D["mW_aftermva_goodcomb"]->Fill(Wmass.M());
					histo1D["mtop_aftermva_goodcomb"]->Fill(topmass.M());
					histo1D["mvavalue_aftermva_goodcomb"]->Fill(mvaVal);
					histo1D["njets_aftermva_goodcomb"]->Fill( (monster->selectedJets()).size() );
					histo1D["Mlb_aftermva_goodcomb"]->Fill( mub.M() );
					
      		histo2D["mWVSmtop_aftermva_goodcomb"]->Fill(Wmass.M(),topmass.M());
				}
				else if(dataSetName.find("TTbarJets_SemiMu") == 0 && SemiMuMatched == false)
				{
      		histo1D["pT_bjet_aftermva_badcomb"]->Fill(monster->selectedJet(mvaResult[2]).Pt());
      		histo1D["pT_ljet_aftermva_badcomb"]->Fill(monster->selectedJet(mvaResult[0]).Pt());
     			histo1D["pT_ljet_aftermva_badcomb"]->Fill(monster->selectedJet(mvaResult[1]).Pt());
					histo1D["mW_aftermva_badcomb"]->Fill(Wmass.M());
					histo1D["mtop_aftermva_badcomb"]->Fill(topmass.M());
					histo1D["mvavalue_aftermva_badcomb"]->Fill(mvaVal);
					histo1D["njets_aftermva_badcomb"]->Fill( (monster->selectedJets()).size() );
					histo1D["Mlb_aftermva_badcomb"]->Fill( mub.M() );
					
      		histo2D["mWVSmtop_aftermva_badcomb"]->Fill(Wmass.M(),topmass.M());
				}
				MSPlot["pT_bjet_aftermva"]->Fill(monster->selectedJet(mvaResult[2]).Pt(), dataSet, true, Luminosity*monster->eventWeight());
				MSPlot["pT_ljet_aftermva"]->Fill(monster->selectedJet(mvaResult[0]).Pt(), dataSet, true, Luminosity*monster->eventWeight());
				MSPlot["pT_ljet_aftermva"]->Fill(monster->selectedJet(mvaResult[1]).Pt(), dataSet, true, Luminosity*monster->eventWeight());
				MSPlot["mW_aftermva"]->Fill(Wmass.M(), dataSet, true, Luminosity*monster->eventWeight());
				MSPlot["mtop_aftermva"]->Fill(topmass.M(), dataSet, true, Luminosity*monster->eventWeight());
				MSPlot["massmub_aftermva"]->Fill(mub.M(), dataSet, true, Luminosity*monster->eventWeight());
				MSPlot["mvavalue_aftermva"]->Fill(mvaVal, dataSet, true, Luminosity*monster->eventWeight());
				MSPlot["njets_aftermva"]->Fill( (monster->selectedJets()).size(), dataSet, true, Luminosity*monster->eventWeight());
				MSPlot["DR_ljets_aftermva"]->Fill( (monster->selectedJet(mvaResult[0])).DeltaR(monster->selectedJet(mvaResult[1]).Pt()), dataSet, true, Luminosity*monster->eventWeight());
      	
				MSPlot["maxProbKinFit"]->Fill(MaxProb, dataSet, true, Luminosity*monster->eventWeight());
        
        for(int iCombi=0; iCombi<12; iCombi++)
        {
          if( maxProb(monster, iCombi) > 0.98 && probNoCorr(monster, dummyMonster, measureTopMass, iCombi) > 0.02 )
          {
            unsigned int* MVAres = monster->mvaResult(iCombi);
            
            TLorentzVector lightJet1 = monster->selectedJet( MVAres[0] );
            TLorentzVector lightJet2 = monster->selectedJet( MVAres[1] );
            TLorentzVector hadrBJet = monster->selectedJet( MVAres[2] );
            TLorentzVector leptBJet = monster->selectedJet( MVAres[3] );
            TLorentzVector muon = monster->muon();
            
            float btag_ljet1 = monster->bTagSSVHE()[ MVAres[0] ];
            if(btag_ljet1 < -0.90) btag_ljet1 = 0;
            float btag_ljet2 = monster->bTagSSVHE()[ MVAres[1] ];
            if(btag_ljet2 < -0.90) btag_ljet2 = 0;
            float btag_hadrbjet = monster->bTagSSVHE()[ MVAres[2] ];
            if(btag_hadrbjet < -0.90) btag_hadrbjet = 0;
            float btag_leptbjet = monster->bTagSSVHE()[ MVAres[3] ];
            if(btag_leptbjet < -0.90) btag_leptbjet = 0;
            
            float btag = pow(btag_hadrbjet,2) + pow(btag_leptbjet,2);
            btag = btag / ( pow(btag_ljet1,2) + pow(btag_ljet1,2) + pow(btag_hadrbjet,2) + pow(btag_leptbjet,2) );
            
            float sumPt = (lightJet1+lightJet2+hadrBJet).Pt() + (lightJet1+lightJet2+leptBJet).Pt() + (lightJet1+hadrBJet+leptBJet).Pt() + (lightJet2+hadrBJet+leptBJet).Pt();
            float ThPtOverSumPt = (lightJet1+lightJet2+hadrBJet).Pt()/sumPt;
            
            float AngleBlMu = leptBJet.Angle(muon.Vect());
            float AngleThBl = (lightJet1+lightJet2+hadrBJet).Angle(leptBJet.Vect());
            float AngleThMu = (lightJet1+lightJet2+hadrBJet).Angle(muon.Vect());
            float Mlb = (leptBJet+muon).M();
            
            if( hadrJetsMVAMatched(monster, iCombi) )
            {
              superMonsterBinnedMVA_AllGoodCombi->Fill(monster, monster->mvaVal(iCombi), dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              superMonsterBinnedBtag_AllGoodCombi->Fill(monster, btag, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              superMonsterBinnedThPtSumPt_AllGoodCombi->Fill(monster, ThPtOverSumPt, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              superMonsterBinnedAngleBlMu_AllGoodCombi->Fill(monster, AngleBlMu, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              superMonsterBinnedAngleThBl_AllGoodCombi->Fill(monster, AngleThBl, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              superMonsterBinnedAngleThMu_AllGoodCombi->Fill(monster, AngleThMu, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              superMonsterBinnedMlb_AllGoodCombi->Fill(monster, Mlb, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              superMonsterBinnedMaxProb_AllGoodCombi->Fill(monster, maxProb(monster, iCombi), dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              superMonsterBinnedProbNoCorr_AllGoodCombi->Fill(monster, probNoCorr(monster, dummyMonster, measureTopMass, iCombi), dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
            }
            else
            {
              superMonsterBinnedMVA_AllBadCombi->Fill(monster, monster->mvaVal(iCombi), dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              superMonsterBinnedBtag_AllBadCombi->Fill(monster, btag, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              superMonsterBinnedThPtSumPt_AllBadCombi->Fill(monster, ThPtOverSumPt, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              superMonsterBinnedAngleBlMu_AllBadCombi->Fill(monster, AngleBlMu, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              superMonsterBinnedAngleThBl_AllBadCombi->Fill(monster, AngleThBl, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              superMonsterBinnedAngleThMu_AllBadCombi->Fill(monster, AngleThMu, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              superMonsterBinnedMlb_AllBadCombi->Fill(monster, Mlb, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              superMonsterBinnedMaxProb_AllBadCombi->Fill(monster, maxProb(monster, iCombi), dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              superMonsterBinnedProbNoCorr_AllBadCombi->Fill(monster, probNoCorr(monster, dummyMonster, measureTopMass, iCombi), dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
            }
          }
          
          if( maxProb(monster, iCombi) > 0.98 )
          {
            stringstream ss1;
            ss1 << iCombi;
            string tmpName = "mvavalue_afterMaxProb_comb" + ss1.str();
            histo1D[tmpName]->Fill( monster->mvaVal(iCombi) );
            
            if( ! hadrJetsMVAMatched(monster, iCombi) )
              superMonsterBinnedJetCombi_BadCombi_MaxProb->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
            
            int nWJets = nWJetsMVAMatched(monster, iCombi);
            int nBJets = nHadrBJetsMVAMatched(monster, iCombi);

            if(nBJets == 0)
            {
              if(nWJets == 0) superMonsterBinnedJetCombi_0W0B_MaxProb->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              else if(nWJets == 1) superMonsterBinnedJetCombi_1W0B_MaxProb->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              else if(nWJets == 2) superMonsterBinnedJetCombi_2W0B_MaxProb->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              else cout << "nWJets = " << nWJets << endl;
            }
            else if(nBJets == 1)
            {
              if(nWJets == 0) superMonsterBinnedJetCombi_0W1B_MaxProb->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              else if(nWJets == 1) superMonsterBinnedJetCombi_1W1B_MaxProb->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              else if(nWJets == 2) superMonsterBinnedJetCombi_2W1B_MaxProb->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              else cout << "nWJets = " << nWJets << endl;
            }
            else cout << "nBJets = " << nBJets << endl;
            
            if( probNoCorr(monster, dummyMonster, measureTopMass, iCombi) > 0.02 /*&& nHoles(monster, measureTopMass, iCombi) < 1  && monster->selectedJets().size() < 5*/ )
            {
              superMonsterBinnedJetCombi->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              superMonsterBinnedInclJetCombi->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              superMonsterBinnedInclJetCombi_BadCombiSubtr->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              superMonsterBinnedInclJetCombi_AllCombiSubtr->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
              
              if(nBJets == 0)
              {
                if(nWJets == 0) superMonsterBinnedJetCombi_0W0B_MaxProbProbNoCorr->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
                else if(nWJets == 1) superMonsterBinnedJetCombi_1W0B_MaxProbProbNoCorr->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
                else if(nWJets == 2) superMonsterBinnedJetCombi_2W0B_MaxProbProbNoCorr->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
                else cout << "nWJets = " << nWJets << endl;
              }
              else if(nBJets == 1)
              {
                if(nWJets == 0) superMonsterBinnedJetCombi_0W1B_MaxProbProbNoCorr->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
                else if(nWJets == 1) superMonsterBinnedJetCombi_1W1B_MaxProbProbNoCorr->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
                else if(nWJets == 2) superMonsterBinnedJetCombi_2W1B_MaxProbProbNoCorr->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
                else cout << "nWJets = " << nWJets << endl;
              }
              else cout << "nBJets = " << nBJets << endl;
              
              if( hadrJetsMVAMatched(monster, iCombi) )
              {
                superMonsterBinnedJetCombi_SemiMuMatched->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
                superMonsterBinnedInclJetCombi_SemiMuMatched->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
                histo1D[tmpName+"_goodcomb"]->Fill( monster->mvaVal(iCombi), dataSet->NormFactor()*Luminosity*monster->eventWeight() );
                histo1D["nGoodCombi_jetCombiBins"]->Fill(iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight());
              }
              else
              {
                if( probNoCorr(monster, dummyMonster, measureTopMass, iCombi) > 0.02 )
                {
                  superMonsterBinnedJetCombi_BadCombi_MaxProbProbNoCorr->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
                  if( probNoCorr(monster, dummyMonster, measureTopMass, iCombi) > 0.05 )
                    superMonsterBinnedJetCombi_BadCombi_MaxProbProbNoCorr0p05->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
                  if( probNoCorr(monster, dummyMonster, measureTopMass, iCombi) > 0.1 )
                    superMonsterBinnedJetCombi_BadCombi_MaxProbProbNoCorr0p1->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
                  if( probNoCorr(monster, dummyMonster, measureTopMass, iCombi) > 0.2 )
                    superMonsterBinnedJetCombi_BadCombi_MaxProbProbNoCorr0p2->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
                  if( nHoles(monster, measureTopMass) < 1 )
                    superMonsterBinnedJetCombi_BadCombi_MaxProbProbNoCorrNoHoles->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
                  if( monster->selectedJets().size() < 5 )
                    superMonsterBinnedJetCombi_BadCombi_MaxProbProbNoCorrNJets->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
                  if( monster->bTagTCHP()[ monster->mvaResult(iCombi)[2] ] > 3.41 )
                  {
                    superMonsterBinnedJetCombi_BadCombi_MaxProbProbNoCorrBtagHadrJet->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
                    if( monster->bTagTCHP()[ monster->mvaResult(iCombi)[3] ] > 3.41 )
                      superMonsterBinnedJetCombi_BadCombi_MaxProbProbNoCorrBtagHadrLeptJet->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
                  }
                }
                superMonsterBinnedInclJetCombi_BadCombi->Fill(monster, iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight(), iCombi);
                histo1D[tmpName+"_badcomb"]->Fill( monster->mvaVal(iCombi), dataSet->NormFactor()*Luminosity*monster->eventWeight() );
                histo1D["nBadCombi_jetCombiBins"]->Fill(iCombi, dataSet->NormFactor()*Luminosity*monster->eventWeight());
              }
            }
          }
        }
        
        if( MaxProb > 0.98 )
        {
          if( monster->hadrBQuark().Pt() > 0.1 && monster->hadrLQuark1().Pt() > 0.1 && monster->hadrLQuark2().Pt() > 0.1 )
          	histo1D["mTop_MCcomb_MaxProb"]->Fill( ( monster->selectedJet(monster->hadrLJet1()) + monster->selectedJet(monster->hadrLJet2()) + monster->selectedJet(monster->hadrBJet()) ).M() );
          
				  if(monster->muonCharge() == 1)
				  {
			      MSPlot["mWHadr_MuPlus_MaxProb"]->Fill(Wmass.M(), dataSet, true, Luminosity*monster->eventWeight());
            MSPlot["mTopHadr_MuPlus_MaxProb"]->Fill(topmass.M(), dataSet, true, Luminosity*monster->eventWeight());
				  }
				  else if(monster->muonCharge() == -1)
				  {
            MSPlot["mWHadr_MuMinus_MaxProb"]->Fill(Wmass.M(), dataSet, true, Luminosity*monster->eventWeight());
            MSPlot["mTopHadr_MuMinus_MaxProb"]->Fill(topmass.M(), dataSet, true, Luminosity*monster->eventWeight());
				  }
          
          if( dataSetName == "Data" || dataSetName == "DATA" || dataSetName == "data" )
            superMonsterBinnedCuts_Data->Fill(monster, 0.5, dataSet->NormFactor()*Luminosity*monster->eventWeight());
          else
          {
            superMonsterBinnedCuts->Fill(monster, 0.5, dataSet->NormFactor()*Luminosity*monster->eventWeight());
//            mvaBinSuperMonsterMaxProb->Fill(monster, dataSet->NormFactor()*Luminosity*monster->eventWeight());
            
            if( SemiMuMatched )
            {
//              mvaBinSuperMonsterMaxProb_SemiMuMatched->Fill(monster, dataSet->NormFactor()*Luminosity*monster->eventWeight());
              superMonsterBinnedCuts_SemiMuMatched->Fill(monster, 0.5, dataSet->NormFactor()*Luminosity*monster->eventWeight());
            }
//            else //if( monster->hadrLQuark1().E() > 1 && monster->hadrLQuark2().E() > 1 && monster->allHadronicJetsMCMatched() && (
//              ( monster->hadrLJet1().DeltaR(monster->hadrLQuark1()) > 0.6 && monster->hadrLJet1().DeltaR(monster->hadrLQuark2()) > 0.6 ) ||
//              ( monster->hadrLJet2().DeltaR(monster->hadrLQuark1()) > 0.6 && monster->hadrLJet2().DeltaR(monster->hadrLQuark2()) > 0.6 ) ) )
//            {
//              mvaBinSuperMonsterMaxPrmeasureTopMassob_BadCombi->Fill(monster, dataSet->NormFactor()*Luminosity*monster->eventWeight());
//            }
          }
          
          float ProbNoCorr = probNoCorr(monster, dummyMonster, measureTopMass);
          int NHoles = nHoles(monster, measureTopMass);
          
          MSPlot["maxMVA_HadrWMass_MaxProbCut"]->Fill(Wmass.M(), dataSet, true, Luminosity*monster->eventWeight());
          MSPlot["maxMVA_HadrTopMass_MaxProbCut"]->Fill(topmass.M(), dataSet, true, Luminosity*monster->eventWeight());
          MSPlot["maxMVA_Mlb_MaxProbCut"]->Fill(mub.M(), dataSet, true, Luminosity*monster->eventWeight());
          MSPlot["ProbNoCorr_MaxProbCut"]->Fill(ProbNoCorr, dataSet, true, Luminosity*monster->eventWeight());
          
          if( dataSetName.find("TTbarJets_SemiMu") == 0 )
          {
            if( SemiMuMatched )
            {
            	histo1D["Mlb_afterMaxProb_goodcomb"]->Fill( mub.M() );
            	histo1D["mvavalue_afterMaxProb_goodcomb"]->Fill( mvaVal );
             	histo1D["nHoles_afterMaxProb_goodcomb"]->Fill(NHoles);
             	histo1D["ProbNoCorr_afterMaxProb_goodcomb"]->Fill( ProbNoCorr );
             	histo1D["mtop_afterMaxProb_goodcomb"]->Fill( topmass.M() );
             	
             	if( monster->all4JetsMCMatched() )
               	histo1D["mvavalue_afterMaxProb_all4JetsMCMatched_goodcomb"]->Fill( mvaVal );
	            if( monster->allHadronicJetsMCMatched() )  
	              histo1D["mvavalue_afterMaxProb_allHadrJetsMCMatched_goodcomb"]->Fill( mvaVal );
	            else
	              histo1D["mvavalue_afterMaxProb_NotMCMatched_goodcomb"]->Fill( mvaVal );
            }
            else
            {
              histo1D["Mlb_afterMaxProb_badcomb"]->Fill( mub.M() );
              histo1D["mvavalue_afterMaxProb_badcomb"]->Fill( mvaVal );
              histo1D["nHoles_afterMaxProb_badcomb"]->Fill(NHoles);
              histo1D["ProbNoCorr_afterMaxProb_badcomb"]->Fill( ProbNoCorr );
              histo1D["mtop_afterMaxProb_badcomb"]->Fill( topmass.M() );
              
              if( monster->all4JetsMCMatched() )
              	histo1D["mvavalue_afterMaxProb_all4JetsMCMatched_badcomb"]->Fill( mvaVal );
              if( monster->allHadronicJetsMCMatched() )
              	histo1D["mvavalue_afterMaxProb_allHadrJetsMCMatched_badcomb"]->Fill( mvaVal );
              else
              	histo1D["mvavalue_afterMaxProb_NotMCMatched_badcomb"]->Fill( mvaVal );
            }
          }
          
          if( ProbNoCorr > 0.02 )
//          if(NHoles < 1)
          {
            if( monster->hadrBQuark().Pt() > 0.1 && monster->hadrLQuark1().Pt() > 0.1 && monster->hadrLQuark2().Pt() > 0.1 )
            {
              TLorentzVector lightJet1 = monster->selectedJet(monster->hadrLJet1());
              TLorentzVector lightJet2 = monster->selectedJet(monster->hadrLJet2());
              TLorentzVector bJet = monster->selectedJet(monster->hadrBJet());
              float lightCorr = expDEl2D->GetBinContent( expDEl2D->FindBin( ( lightJet1.Pt() + lightJet2.Pt() )/2, bJet.Pt() ) );
              float bCorr = expDEb2D->GetBinContent( expDEb2D->FindBin( ( lightJet1.Pt() + lightJet2.Pt() )/2, bJet.Pt() ) );

              float mTopAllCorr = ( lightJet1*(1+lightCorr) + lightJet2*(1+lightCorr) + bJet*(1+bCorr) ).M();
              float mTopbCorr = ( lightJet1 + lightJet2 + bJet*(1+bCorr) ).M();
            	
           	  JEClight->setJetEta(lightJet1.Eta());
              JEClight->setJetPt(lightJet1.Pt());
              float corrJEClight1 = JEClight->getCorrection();
              JEClight->setJetEta(lightJet2.Eta());
              JEClight->setJetPt(lightJet2.Pt());
              float corrJEClight2 = JEClight->getCorrection();
              
              JECb->setJetEta(bJet.Eta());
              JECb->setJetPt(bJet.Pt());
              float corrJECb = JECb->getCorrection();
              
              float mTopJECAllCorr = ( lightJet1*corrJEClight1 + lightJet2*corrJEClight2 + bJet*corrJECb ).M();
              float mTopJECbCorr = ( lightJet1 + lightJet2 + bJet*corrJECb ).M();
              
            	histo1D["mTop_MCcomb_MaxProbProbNoCorr"]->Fill( ( lightJet1 + lightJet2 + bJet ).M() );
           	 	histo1D["mTop_MCcomb_MaxProbProbNoCorr_bCorr"]->Fill( mTopbCorr );
             	histo1D["mTop_MCcomb_MaxProbProbNoCorr_AllCorr"]->Fill( mTopAllCorr );
             	histo1D["mTop_MCcomb_MaxProbProbNoCorr_JECbCorr"]->Fill( mTopJECbCorr );
             	histo1D["mTop_MCcomb_MaxProbProbNoCorr_JECAllCorr"]->Fill( mTopJECAllCorr );
            }
            
            float minDistHoleMax = -1; // minDisHoleMax(&monsterHisto);
            float MTopMax = -1; // mTopMax(&monsterHisto);
            
            if(monster->muonCharge() == 1)
				    {
			        MSPlot["mWHadr_MuPlus_ProbNoCorr"]->Fill(Wmass.M(), dataSet, true, Luminosity*monster->eventWeight());
              MSPlot["mTopHadr_MuPlus_ProbNoCorr"]->Fill(topmass.M(), dataSet, true, Luminosity*monster->eventWeight());
				    }
				    else if(monster->muonCharge() == -1)
				    {
              MSPlot["mWHadr_MuMinus_ProbNoCorr"]->Fill(Wmass.M(), dataSet, true, Luminosity*monster->eventWeight());
              MSPlot["mTopHadr_MuMinus_ProbNoCorr"]->Fill(topmass.M(), dataSet, true, Luminosity*monster->eventWeight());
				    }
            
            if( dataSetName == "Data" || dataSetName == "DATA" || dataSetName == "data" )
            {
              superMonsterBinnedCuts_Data->Fill(monster, 1.5, dataSet->NormFactor()*Luminosity*monster->eventWeight());
//              superMonsterBinnedPtL_Data->Fill(monster, monster->selectedJet(mvaResult[0]).Pt(), 1);
//              superMonsterBinnedPtL_Data->Fill(monster, monster->selectedJet(mvaResult[1]).Pt(), 1);
//              superMonsterBinnedPtB_Data->Fill(monster, monster->selectedJet(mvaResult[2]).Pt(), 1);
//              superMonsterBinnedMaxMVA_Data->Fill(monster, mvaVal, 1);
//              superMonsterBinnedCharge_Data->Fill(monster, monster->muonCharge(), 1);
            }
            else
            {
              if(dataSetName.find("TTbarJets_SemiMu") == 0) //semi-mu
              {
                if(minDistHoleMax < 999)
              	{
              	  histo1D["minDistHoleMax"]->Fill(minDistHoleMax);
                  histo1D["nHoles"]->Fill(NHoles);
                  histo1D["mTopMax"]->Fill(MTopMax);
                  histo2D["mTopMaxVSminDistHoleMax"]->Fill(MTopMax,minDistHoleMax);
                  histo2D["nHolesVSminDistHoleMax"]->Fill(NHoles,minDistHoleMax);
                  histo2D["mTopMaxVSnHoles"]->Fill(MTopMax,NHoles);
                  
                  histo1D["holes_PtJets"]->Fill(monster->selectedJet(mvaResult[1]).Pt());
                  histo1D["holes_PtJets"]->Fill(monster->selectedJet(mvaResult[2]).Pt());
                  histo1D["holes_PtJets"]->Fill(monster->selectedJet(mvaResult[0]).Pt());
                  histo1D["holes_mW"]->Fill( Wmass.M() );
                  histo1D["holes_mTop"]->Fill( topmass.M() );
                  histo1D["holes_EtaJets"]->Fill(monster->selectedJet(mvaResult[0]).Eta());
                  histo1D["holes_EtaJets"]->Fill(monster->selectedJet(mvaResult[1]).Eta());
                  histo1D["holes_EtaJets"]->Fill(monster->selectedJet(mvaResult[2]).Eta());
                  histo1D["holes_MVAvalue"]->Fill(mvaVal);
                  histo1D["holes_DRlightjets"]->Fill(monster->selectedJet(mvaResult[0]).DeltaR(monster->selectedJet(mvaResult[1])));
                  histo1D["holes_DRlightbjet"]->Fill(monster->selectedJet(mvaResult[2]).DeltaR(monster->selectedJet(mvaResult[0])));
                  histo1D["holes_DRlightbjet"]->Fill(monster->selectedJet(mvaResult[2]).DeltaR(monster->selectedJet(mvaResult[1])));
              	}
              	else
              	{
                  histo1D["noholes_PtJets"]->Fill(monster->selectedJet(mvaResult[0]).Pt());
                  histo1D["noholes_PtJets"]->Fill(monster->selectedJet(mvaResult[1]).Pt());
                  histo1D["noholes_PtJets"]->Fill(monster->selectedJet(mvaResult[2]).Pt());
                  histo1D["noholes_mW"]->Fill( Wmass.M() );
                  histo1D["noholes_mTop"]->Fill( topmass.M() );
                  histo1D["noholes_EtaJets"]->Fill(monster->selectedJet(mvaResult[0]).Eta());
                  histo1D["noholes_EtaJets"]->Fill(monster->selectedJet(mvaResult[1]).Eta());
                  histo1D["noholes_EtaJets"]->Fill(monster->selectedJet(mvaResult[2]).Eta());
                  histo1D["noholes_MVAvalue"]->Fill(mvaVal);
                  histo1D["noholes_DRlightjets"]->Fill(monster->selectedJet(mvaResult[0]).DeltaR(monster->selectedJet(mvaResult[1])));
                  histo1D["noholes_DRlightbjet"]->Fill(monster->selectedJet(mvaResult[2]).DeltaR(monster->selectedJet(mvaResult[0])));
                  histo1D["noholes_DRlightbjet"]->Fill(monster->selectedJet(mvaResult[2]).DeltaR(monster->selectedJet(mvaResult[1])));
              	}
                
                if( SemiMuMatched )
                {
                 	histo1D["Mlb_afterprobcuts_goodcomb"]->Fill( mub.M() );
               		histo1D["mvavalue_afterprobcuts_goodcomb"]->Fill( mvaVal );
               		histo1D["MaxProb_afterprobcuts_goodcomb"]->Fill(MaxProb);
             			histo1D["ProbNoCorr_afterprobcuts_goodcomb"]->Fill(ProbNoCorr);
             			histo1D["nHoles_afterProbNoCorr_goodcomb"]->Fill(NHoles);
             			histo1D["mtop_afterProbNoCorr_goodcomb"]->Fill(topmass.M());
             			if(monster->bTagTCHP()[0] > 3.41)
             			  histo1D["mtop_afterProbNoCorr_TCHPT_goodcomb"]->Fill(topmass.M());
             			if(monster->bTagSSVHP()[0] > 2.)
             			  histo1D["mtop_afterProbNoCorr_SSVHPT_goodcomb"]->Fill(topmass.M());
             			histo1D["MET_afterProbNoCorr_goodcomb"]->Fill(monster->met().Pt());
             			histo1D["MHT_afterProbNoCorr_goodcomb"]->Fill( ( topmass + mub ).Pt() );
                  histo1D["TCHP_hadrBJet_afterProbNoCorr_goodcomb"]->Fill(monster->bTagTCHP()[mvaResult[2]]);
                  histo1D["TCHP_lepBJet_afterProbNoCorr_goodcomb"]->Fill(monster->bTagTCHP()[mvaResult[3]]);
                  histo1D["nJets_afterProbNoCorr_goodcomb"]->Fill( monster->selectedJets().size() );
                  histo1D["DPhi_MET_MHT_afterProbNoCorr_goodcomb"]->Fill( ( topmass + mub ).DeltaPhi(monster->met()) );
                }
                else
                {
                  histo1D["Mlb_afterprobcuts_badcomb"]->Fill( mub.M() );
                  histo1D["mvavalue_afterprobcuts_badcomb"]->Fill( mvaVal );
                  histo1D["MaxProb_afterprobcuts_badcomb"]->Fill(MaxProb);
             			histo1D["ProbNoCorr_afterprobcuts_badcomb"]->Fill(ProbNoCorr);
             			histo1D["nHoles_afterProbNoCorr_badcomb"]->Fill(NHoles);
             			histo1D["mtop_afterProbNoCorr_badcomb"]->Fill(topmass.M());
             			if(monster->bTagTCHP()[0] > 3.41)
             			  histo1D["mtop_afterProbNoCorr_TCHPT_badcomb"]->Fill(topmass.M());
             			if(monster->bTagSSVHP()[0] > 2.)
             			  histo1D["mtop_afterProbNoCorr_SSVHPT_badcomb"]->Fill(topmass.M());
                  histo1D["MET_afterProbNoCorr_badcomb"]->Fill(monster->met().Pt());
             			histo1D["MHT_afterProbNoCorr_badcomb"]->Fill( ( topmass + mub ).Pt() );
                  histo1D["TCHP_hadrBJet_afterProbNoCorr_badcomb"]->Fill(monster->bTagTCHP()[mvaResult[2]]);
                  histo1D["TCHP_lepBJet_afterProbNoCorr_badcomb"]->Fill(monster->bTagTCHP()[mvaResult[3]]);
                  histo1D["nJets_afterProbNoCorr_badcomb"]->Fill( monster->selectedJets().size() );
                  histo1D["DPhi_MET_MHT_afterProbNoCorr_badcomb"]->Fill( ( topmass + mub ).DeltaPhi(monster->met()) );
                }
                
                if( SemiMuMatched )
                 	histo1D["mvavalue_afterprobcuts_goodcomb_template"]->Fill( mvaVal, dataSet->NormFactor()*Luminosity*monster->eventWeight());
                else
                	histo1D["mvavalue_afterprobcuts_badcomb_template"]->Fill( mvaVal, dataSet->NormFactor()*Luminosity*monster->eventWeight());
              }
              histo1D["mvavalue_afterprobcuts_allcomb_MC"]->Fill( mvaVal, dataSet->NormFactor()*Luminosity*monster->eventWeight());
              
//              mvaBinSuperMonsterProbNoCorr->Fill(monster, dataSet->NormFactor()*Luminosity*monster->eventWeight());
              
              superMonsterBinnedCuts->Fill(monster, 1.5, dataSet->NormFactor()*Luminosity*monster->eventWeight());
              
//              superMonsterBinnedPtL->Fill(monster, monster->selectedJet(mvaResult[0]).Pt(), dataSet->NormFactor()*Luminosity*monster->eventWeight());
//              superMonsterBinnedPtL->Fill(monster, monster->selectedJet(mvaResult[1]).Pt(), dataSet->NormFactor()*Luminosity*monster->eventWeight());
//              superMonsterBinnedPtB->Fill(monster, monster->selectedJet(mvaResult[2]).Pt(), dataSet->NormFactor()*Luminosity*monster->eventWeight());
//              superMonsterBinnedCharge->Fill(monster, monster->muonCharge(), dataSet->NormFactor()*Luminosity*monster->eventWeight());
//              superMonsterBinnedMaxMVA->Fill(monster, mvaVal, dataSet->NormFactor()*Luminosity*monster->eventWeight());
              
              if(SemiMuMatched)
              {
//                mvaBinSuperMonsterProbNoCorr_SemiMuMatched->Fill(monster, dataSet->NormFactor()*Luminosity*monster->eventWeight());
                superMonsterBinnedCuts_SemiMuMatched->Fill(monster, 1.5, dataSet->NormFactor()*Luminosity*monster->eventWeight());
                
//                superMonsterBinnedPtL_SemiMuMatched->Fill(monster, monster->selectedJet(mvaResult[0]).Pt(), dataSet->NormFactor()*Luminosity*monster->eventWeight());
//                superMonsterBinnedPtL_SemiMuMatched->Fill(monster, monster->selectedJet(mvaResult[1]).Pt(), dataSet->NormFactor()*Luminosity*monster->eventWeight());
//                superMonsterBinnedPtB_SemiMuMatched->Fill(monster, monster->selectedJet(mvaResult[2]).Pt(), dataSet->NormFactor()*Luminosity*monster->eventWeight());
//                superMonsterBinnedMaxMVA_SemiMuMatched->Fill(monster, mvaVal, dataSet->NormFactor()*Luminosity*monster->eventWeight());
//                superMonsterBinnedCharge_SemiMuMatched->Fill(monster, monster->muonCharge(), dataSet->NormFactor()*Luminosity*monster->eventWeight());
              }
//              else //if(monster->hadrBQuark().E() < 0.1 && monster->hadrLQuark1().E() < 0.1 && monster->hadrLQuark2().E() < 0.1)
//              {
//                mvaBinSuperMonsterProbNoCorr_BadCombi->Fill(monster, dataSet->NormFactor()*Luminosity*monster->eventWeight());
//                superMonsterBinnedPtL_BadCombi->Fill(monster, monster->selectedJet(mvaResult[0]).Pt(), dataSet->NormFactor()*Luminosity*monster->eventWeight());
//                superMonsterBinnedPtL_BadCombi->Fill(monster, monster->selectedJet(mvaResult[1]).Pt(), dataSet->NormFactor()*Luminosity*monster->eventWeight());
//                superMonsterBinnedPtB_BadCombi->Fill(monster, monster->selectedJet(mvaResult[2]).Pt(), dataSet->NormFactor()*Luminosity*monster->eventWeight());
//              }
            }
            
            MSPlot["maxMVA_HadrWMass_probcuts"]->Fill(Wmass.M(), dataSet, true, Luminosity*monster->eventWeight());
            MSPlot["maxMVA_HadrTopMass_probcuts"]->Fill(topmass.M(), dataSet, true, Luminosity*monster->eventWeight());
            MSPlot["maxMVA_Mlb_probcuts"]->Fill(mub.M(), dataSet, true, Luminosity*monster->eventWeight());
            
            if( monster->bTagTCHP()[mvaResult[2]] > 3.41 )
            {
              if( dataSetName == "Data" || dataSetName == "DATA" || dataSetName == "data" )
              {
                superMonsterBinnedCuts_Data->Fill(monster, 2.5, dataSet->NormFactor()*Luminosity*monster->eventWeight());
              }
              else
              {
                superMonsterBinnedCuts->Fill(monster, 2.5, dataSet->NormFactor()*Luminosity*monster->eventWeight());
                if(SemiMuMatched)
                {
                  superMonsterBinnedCuts_SemiMuMatched->Fill(monster, 2.5, dataSet->NormFactor()*Luminosity*monster->eventWeight());
                }
              }
              
              if( monster->bTagTCHP()[mvaResult[3]] > 3.41 )
              {
                if( dataSetName == "Data" || dataSetName == "DATA" || dataSetName == "data" )
                {
                  superMonsterBinnedCuts_Data->Fill(monster, 3.5, dataSet->NormFactor()*Luminosity*monster->eventWeight());
                }
                else
                {
                  superMonsterBinnedCuts->Fill(monster, 3.5, dataSet->NormFactor()*Luminosity*monster->eventWeight());
                  if(SemiMuMatched)
                  {
                    superMonsterBinnedCuts_SemiMuMatched->Fill(monster, 3.5, dataSet->NormFactor()*Luminosity*monster->eventWeight());
                  }
                }
              }
            } // nHoles cut
          } // probNoCorr cut
        } // MaxProb cut
      } // maxMva cut
    } // end loop over monsters
    cout << "Not closing and not deleting!" << endl;
//    inFile->Close();
//    delete inFile;
  } // end loop over datasets
  
  string pathPNGsupermonsters = pathPNG + "SuperMonsters/";
  mkdir(pathPNGsupermonsters.c_str(),0777);
  
  cout << "Finished running over all datasets..." << endl;
  
  cout << "Subtract effect of bad comi's" << endl;
  
  float nBadCombis = histo1D["nBadCombi_jetCombiBins"]->Integral(1,1);
  
  cout << " --> nBadCombis = " << nBadCombis << endl;
  
//  void SubTractAvMonster(TH2F* avMonster, float binValue, float nBadCombis);
//  TH2F* AverageMonster(float binValue);
  
  superMonsterBinnedInclJetCombi_BadCombiSubtr->SubTractAvMonster( superMonsterBinnedInclJetCombi_BadCombi->AverageMonster(11.), 0., nBadCombis);
  superMonsterBinnedInclJetCombi_AllCombiSubtr->SubTractAvMonster( superMonsterBinnedInclJetCombi->AverageMonster(11.), 0., nBadCombis);
  
  histo1D["mTop_MCcomb"]->Fit("gaus","Q","",155,190);
  histo1D["mTop_MCcomb_bCorr"]->Fit("gaus","Q","",155,190);
 	histo1D["mTop_MCcomb_AllCorr"]->Fit("gaus","Q","",155,190);
 	histo1D["mTop_MCcomb_JECbCorr"]->Fit("gaus","Q","",155,190);
 	histo1D["mTop_MCcomb_JECAllCorr"]->Fit("gaus","Q","",155,190);
  histo1D["mTop_MCcomb_MaxProb"]->Fit("gaus","Q","",155,190);
  histo1D["mTop_MCcomb_MaxProbProbNoCorr"]->Fit("gaus","Q","",155,190);
  histo1D["mTop_MCcomb_MaxProbProbNoCorr_bCorr"]->Fit("gaus","Q","",155,190);
  histo1D["mTop_MCcomb_MaxProbProbNoCorr_AllCorr"]->Fit("gaus","Q","",155,190);
 	histo1D["mtop_afterProbNoCorr_goodcomb"]->Fit("gaus","Q","",155,190);
 	histo1D["mtop_afterMaxProb_goodcomb"]->Fit("gaus","Q","",155,190);
  histo1D["mtop_aftermva_goodcomb"]->Fit("gaus","Q","",155,190);
  histo1D["mTop_MCcomb_MaxProbProbNoCorr_JECbCorr"]->Fit("gaus","Q","",155,190);
  histo1D["mTop_MCcomb_MaxProbProbNoCorr_JECAllCorr"]->Fit("gaus","Q","",155,190);
  
  cout << "Writing out..." << endl;
  fout->cd();
  
  // write out binned results
  string pathPNGbinnedResults = pathPNG + "BinnedResults/";
  mkdir(pathPNGbinnedResults.c_str(),0777);
  
//  superMonsterBinnedPtL->Write(fout, pathPNGbinnedResults, "Pt light jets", false);
//  superMonsterBinnedPtB->Write(fout, pathPNGbinnedResults, "Pt b-jet", false);
//  superMonsterBinnedPtL_Data->Write(fout, pathPNGbinnedResults, "Pt light jets", true);
//  superMonsterBinnedPtB_Data->Write(fout, pathPNGbinnedResults, "Pt b-jet", true);
//  superMonsterBinnedPtL_SemiMuMatched->Write(fout, pathPNGbinnedResults, "Pt light jets", true);
//  superMonsterBinnedPtB_SemiMuMatched->Write(fout, pathPNGbinnedResults, "Pt b-jet", true);
//  superMonsterBinnedPtL_BadCombi->Write(fout, pathPNGbinnedResults, "Pt light jets", true);
//  superMonsterBinnedPtB_BadCombi->Write(fout, pathPNGbinnedResults, "Pt b-jet", true);

//  superMonsterBinnedCuts->Write(fout, pathPNGbinnedResults, "Cut Number", true);
//  superMonsterBinnedCuts_SemiMuMatched->Write(fout, pathPNGbinnedResults, "Cut Number", true);
//  superMonsterBinnedCuts_Data->Write(fout, pathPNGbinnedResults, "Cut Number", true);
  
  superMonsterBinnedJetCombi->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
  superMonsterBinnedJetCombi_SemiMuMatched->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
  superMonsterBinnedJetCombi_BadCombi_MaxProb->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
  superMonsterBinnedJetCombi_BadCombi_MaxProbProbNoCorr->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
//  superMonsterBinnedJetCombi_BadCombi_MaxProbProbNoCorr0p05->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
//  superMonsterBinnedJetCombi_BadCombi_MaxProbProbNoCorr0p1->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
//  superMonsterBinnedJetCombi_BadCombi_MaxProbProbNoCorr0p2->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
//  superMonsterBinnedJetCombi_BadCombi_MaxProbProbNoCorrNoHoles->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
//  superMonsterBinnedJetCombi_BadCombi_MaxProbProbNoCorrNJets->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
//  superMonsterBinnedJetCombi_BadCombi_MaxProbProbNoCorrBtagHadrJet->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
//  superMonsterBinnedJetCombi_BadCombi_MaxProbProbNoCorrBtagHadrLeptJet->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
  
//  superMonsterBinnedJetCombi_0W0B_MaxProb->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
//  superMonsterBinnedJetCombi_0W1B_MaxProb->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
//  superMonsterBinnedJetCombi_1W0B_MaxProb->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
//  superMonsterBinnedJetCombi_1W1B_MaxProb->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
//  superMonsterBinnedJetCombi_2W0B_MaxProb->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
//  superMonsterBinnedJetCombi_2W1B_MaxProb->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
 
//  superMonsterBinnedJetCombi_0W0B_MaxProbProbNoCorr->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
//  superMonsterBinnedJetCombi_0W1B_MaxProbProbNoCorr->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
//  superMonsterBinnedJetCombi_1W0B_MaxProbProbNoCorr->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
//  superMonsterBinnedJetCombi_1W1B_MaxProbProbNoCorr->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
//  superMonsterBinnedJetCombi_2W0B_MaxProbProbNoCorr->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
//  superMonsterBinnedJetCombi_2W1B_MaxProbProbNoCorr->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
  
  superMonsterBinnedInclJetCombi->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
  superMonsterBinnedInclJetCombi_SemiMuMatched->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
  superMonsterBinnedInclJetCombi_BadCombi->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
  
  superMonsterBinnedInclJetCombi_BadCombiSubtr->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
  superMonsterBinnedInclJetCombi_AllCombiSubtr->Write(fout, pathPNGbinnedResults, "Jet Combi Number", true);
  
//  superMonsterBinnedCharge->Write(fout, pathPNGbinnedResults, "Muon Charge", true);
//  superMonsterBinnedCharge_SemiMuMatched->Write(fout, pathPNGbinnedResults, "Muon Charge", true);
//  superMonsterBinnedCharge_Data->Write(fout, pathPNGbinnedResults, "Muon Charge", true);
  
//  superMonsterBinnedMaxMVA->Write(fout, pathPNGbinnedResults, "Max MVA", true);
//  superMonsterBinnedMaxMVA_SemiMuMatched->Write(fout, pathPNGbinnedResults, "Max MVA", true);
//  superMonsterBinnedMaxMVA_Data->Write(fout, pathPNGbinnedResults, "Max MVA", true);
  
  superMonsterBinnedMVA_AllGoodCombi->Write(fout, pathPNGbinnedResults, "MVA value", true);
  superMonsterBinnedMVA_AllBadCombi->Write(fout, pathPNGbinnedResults, "MVA value", true);
  superMonsterBinnedBtag_AllGoodCombi->Write(fout, pathPNGbinnedResults, "btag", true);
  superMonsterBinnedBtag_AllBadCombi->Write(fout, pathPNGbinnedResults, "btag", true);
  superMonsterBinnedThPtSumPt_AllGoodCombi->Write(fout, pathPNGbinnedResults, "ThPtOverSumPt", true);
  superMonsterBinnedThPtSumPt_AllBadCombi->Write(fout, pathPNGbinnedResults, "ThPtOverSumPt", true);
  superMonsterBinnedAngleBlMu_AllGoodCombi->Write(fout, pathPNGbinnedResults, "AngleBlMu", true);
  superMonsterBinnedAngleBlMu_AllBadCombi->Write(fout, pathPNGbinnedResults, "AngleBlMu", true);
  superMonsterBinnedAngleThBl_AllGoodCombi->Write(fout, pathPNGbinnedResults, "AngleThBl", true);
  superMonsterBinnedAngleThBl_AllBadCombi->Write(fout, pathPNGbinnedResults, "AngleThBl", true);
  superMonsterBinnedAngleThMu_AllGoodCombi->Write(fout, pathPNGbinnedResults, "AngleThMu", true);
  superMonsterBinnedAngleThMu_AllBadCombi->Write(fout, pathPNGbinnedResults, "AngleThMu", true);
  superMonsterBinnedMlb_AllGoodCombi->Write(fout, pathPNGbinnedResults, "Mlb", true);
  superMonsterBinnedMlb_AllBadCombi->Write(fout, pathPNGbinnedResults, "Mlb", true);
//  superMonsterBinnedMaxProb_AllGoodCombi->Write(fout, pathPNGbinnedResults, "MaxProb", true);
//  superMonsterBinnedMaxProb_AllBadCombi->Write(fout, pathPNGbinnedResults, "MaxProb", true);
//  superMonsterBinnedProbNoCorr_AllGoodCombi->Write(fout, pathPNGbinnedResults, "ProbNoCorr", true);
//  superMonsterBinnedProbNoCorr_AllBadCombi->Write(fout, pathPNGbinnedResults, "ProbNoCorr", true);
  
  string pathPNGMVABinnedResults = pathPNG + "MVABinnedResults/";
  mkdir(pathPNGMVABinnedResults.c_str(),0777);
  
//  mvaBinSuperMonsterMaxProb->Write(fout, pathPNGMVABinnedResults, true, histo1D["mvavalue_afterprobcuts_goodcomb_template"], histo1D["mvavalue_afterprobcuts_badcomb_template"]);
//  mvaBinSuperMonsterProbNoCorr->Write(fout, pathPNGMVABinnedResults, true, histo1D["mvavalue_afterprobcuts_goodcomb_template"], histo1D["mvavalue_afterprobcuts_badcomb_template"]);

//  mvaBinSuperMonsterMaxProb_SemiMuMatched->Write(fout, pathPNGMVABinnedResults, true, histo1D["mvavalue_afterprobcuts_goodcomb_template"], histo1D["mvavalue_afterprobcuts_badcomb_template"]);
//  mvaBinSuperMonsterProbNoCorr_SemiMuMatched->Write(fout, pathPNGMVABinnedResults, true, histo1D["mvavalue_afterprobcuts_goodcomb_template"], histo1D["mvavalue_afterprobcuts_badcomb_template"]);
  
//  mvaBinSuperMonsterMaxProb_BadCombi->Write(fout, pathPNGMVABinnedResults, false, histo1D["mvavalue_afterprobcuts_goodcomb_template"], histo1D["mvavalue_afterprobcuts_badcomb_template"]);
//  mvaBinSuperMonsterProbNoCorr_BadCombi->Write(fout, pathPNGMVABinnedResults, false, histo1D["mvavalue_afterprobcuts_goodcomb_template"], histo1D["mvavalue_afterprobcuts_badcomb_template"]);
  
//  float** resultsTemplateFit = estimateGoodCombisMVA(histo1D["mvavalue_afterprobcuts_goodcomb_template"], histo1D["mvavalue_afterprobcuts_badcomb_template"], histo1D["mvavalue_afterprobcuts_allcomb_MC"], pathPNG, fout, "maxMVA_templateFit");
  
  // special 1D histo's
  histo1D["nCombi_jetCombiBins"] = (TH1F*) histo1D["nGoodCombi_jetCombiBins"]->Clone();
  histo1D["nCombi_jetCombiBins"]->Add( histo1D["nBadCombi_jetCombiBins"] );
  histo1D["nCombi_jetCombiBins"]->SetNameTitle("nCombi_jetCombiBins","nCombi_jetCombiBins");
  histo1D["fracGoodCombi_jetCombiBins"] = (TH1F*) histo1D["nGoodCombi_jetCombiBins"]->Clone();
  histo1D["fracGoodCombi_jetCombiBins"]->Divide( histo1D["nCombi_jetCombiBins"] );
  histo1D["fracGoodCombi_jetCombiBins"]->SetNameTitle("fracGoodCombi_jetCombiBins","fracGoodCombi_jetCombiBins");
  
  string pathPNGMVAstuff = pathPNG + "MVA/";
  mkdir(pathPNGMVAstuff.c_str(),0777);
  jetCombiner->Write(fout, true, pathPNGMVAstuff);
  
  TDirectory* th1dir = fout->mkdir("1D_histograms");
  th1dir->cd();
  // Write 1D histo's
  for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
	{
		TH1F *temp = it->second;
		temp->Write();
		TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
		tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
	}
	
	fout->cd();
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
    MultiSamplePlot *temp = it->second;
    string name = it->first;
    temp->Draw(false, name, true, true, true, true, true);
    temp->Write(fout, name, true, pathPNG+"MSPlot/");
  }
	
	// 2D
  TDirectory* th2dir = fout->mkdir("2D_histograms_graphs");
  th2dir->cd();
	for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
	{
		TH2F *temp = it->second;
		temp->Write();
		TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
		tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
	}

	//Write TGraphAsymmErrors
	fout->cd();
	for(map<string,TGraphAsymmErrors*>::const_iterator it = graphAsymmErr.begin(); it != graphAsymmErr.end(); it++)
	{
	  TGraphAsymmErrors *temp = it->second;
	  temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
		tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
	}
	
  //Write TGraphErrors
  fout->cd();
  for(map<string,TGraphErrors*>::const_iterator it = graphErr.begin(); it != graphErr.end(); it++)
  {
    TGraphErrors *temp = it->second;
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
  }
  
//  fout->Close();
  
  delete dummyMonster;
//  delete JEClight;
//  delete JECb;
  
  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;

  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "           hasn't crashed yet ;-)           " << endl;
  cout << "********************************************" << endl;
  
  return 0;
}
