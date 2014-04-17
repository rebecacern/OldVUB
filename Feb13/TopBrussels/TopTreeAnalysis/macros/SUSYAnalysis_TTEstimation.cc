#include "iomanip.h"
#include <cmath>
#include "TF1.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TBranch.h"
#include "TLegend.h"
#include "TStyle.h"

//user code
#include "TopTreeProducer/interface/TRootGenEvent.h"
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "../config/Datasets.cc"
#include "../Selection/interface/SelectionTable.h"
#include "../BkgEstimationMethods/interface/WJetEstimation.h"
#include "../BkgEstimationMethods/interface/WJetEstPseudoExp.h"
#include "../BkgEstimationMethods/interface/TtJetEstimation.h"
#include "../BkgEstimationMethods/interface/ABCDEstimation.h"
#include "../BkgEstimationMethods/interface/QCDShapeEstimation.h"
#include "../BkgEstimationMethods/interface/WJetShapeEstimation.h"
#include "../BkgEstimationMethods/interface/BkgEstimationSummary.h"
#include "../Tools/interface/PlottingTools.h"
#include "../Tools/interface/MultiSamplePlot.h"
#include "../Tools/interface/MultiCutPlot.h"
#include "../Tools/interface/CutImpactEvaluation.h"
#include "../Reconstruction/interface/MEzCalculator.h"

#include "Style.C"

using namespace std;
using namespace TopTree;


  //to be added

  //Choose cut for selection !!
  //Efficiency - Rejection - S/(S+B)
  // WARNING: code to be developed
  
  //Expected Estimation w/o SUSY contamination: WARNING: had possibility to choose !!

vector<TRootJet> JES(vector<TRootJet> jets, float factor){
	vector<TRootJet> output;
	for(unsigned int i=0;i<jets.size();i++){
		TRootJet jet(TLorentzVector(jets[i].Px()*factor,jets[i].Py()*factor,jets[i].Pz()*factor,jets[i].Energy()*factor));
		jet.setHcalEnergyFraction(jets[i].hcalEnergyFraction());
		jet.setBtag_trackCountingHighEffBJetTags(jets[i].btag_trackCountingHighEffBJetTags());
		jet.setBtag_trackCountingHighPurBJetTags(jets[i].btag_trackCountingHighPurBJetTags());
		jet.setBtag_jetProbabilityBJetTags(jets[i].btag_jetProbabilityBJetTags());
		jet.setBtag_jetBProbabilityBJetTags(jets[i].btag_jetBProbabilityBJetTags());
		jet.setBtag_simpleSecondaryVertexBJetTags(jets[i].btag_simpleSecondaryVertexBJetTags());
		output.push_back(jet);
	}
	return output;
}


void Scale(vector<Dataset> datasets, float*& numbers, float Luminosity ){
	for(unsigned int i=0;i<datasets.size();i++){
		numbers[i] = numbers[i]*datasets[i].NormFactor()*Luminosity;
	}
}

void Scale(vector<Dataset> datasets, double*& numbers, float Luminosity ){
	for(unsigned int i=0;i<datasets.size();i++){
		numbers[i] = numbers[i]*datasets[i].NormFactor()*Luminosity;
	}
}

void Scale( Dataset dataset, TH1F* histo, float Luminosity ){
	if(histo!=0) histo->Scale(dataset.NormFactor()*Luminosity);
	cout<<"factor = "<<dataset.NormFactor()*Luminosity<<endl;
}


float Mttbar(TRootMuon muon, std::vector< TRootJet > jets, TRootMET met)
{
        if(jets.size()<4) return 0;
        TLorentzVector TTBar = (TLorentzVector) (muon);
        TTBar += (TLorentzVector)(jets[0]);
        TTBar += (TLorentzVector)(jets[1]);
        TTBar += (TLorentzVector)(jets[2]);
        TTBar += (TLorentzVector)(jets[3]);
        TTBar += (TLorentzVector)(met);
        MEzCalculator NuPz;
        NuPz.SetMET(met);
        NuPz.SetMuon(muon);
        TTBar.SetPz(TTBar.Pz()+NuPz.Calculate());
        float Mtt = TTBar.M();
        return Mtt;
}


int main(int argc, char*argv[]){


  cout<<"**********************************************************************"<<endl;
  cout<<"Begining of the program for bkg estimation in lepton+jets channels !" <<endl;
  cout<<"**********************************************************************"<<endl;

  //setTDRStyle(); 
  setMyStyle();

  /////////////////////
  // Configuration
  /////////////////////

  int verbose = 2;
  float Luminosity = 100; // in 1/fb
  bool addLMX = false;
  cout<<"The results will be obtained for a luminosity of "<<Luminosity<<" x 1/pb"<<endl;
  
  //choice of variable
  bool isMET = false;
  bool isHT = true;
  bool isMttbar = false;
 
  //Binning for the plots of the variable (ex: MET)
  //int NbinsX = 10;
  //float binsX[] = {0,25,50,75,100,150,200,300,500,700,1000};
  //int NbinsX = 40;
  //float binsX[] = {0,25,50,75,100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500,525,550,575,600,625,650,675,700,725,750,775,800,825,850,875,900,925,950,975,1000};
  int NbinsX = 20;
  float binsX[] = {0,25,50,75,100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500};
  
  //int NbinsX = 21;
  //float binsX[] = {120,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200};
  string VarLabel;
  string XaxisLabel;
  if(isMET){
  	VarLabel = string("MET");
  	XaxisLabel = string("MET [GeV]");
  }
  if(isHT){
  	VarLabel = string("H_{T}");
  	XaxisLabel = string("H_{T} [GeV]");
  }
  if(isMttbar){
  	VarLabel = string("Mttbar");
  	XaxisLabel = string("Mttbar [GeV]");
  }

  //Output files
  string selecTableFileName("SUSY_SelectionTable.log"); 
  string selecTableCRFileName("SUSY_SelectionTableCR.log"); 
  string bkgEstimPerJetsBJetsFileName("Estim_Jets_BJets.log");
  string rootFileName("TtJetEstimation.root");


  /////////////////////
  // Load Datasets
  /////////////////////
  if(verbose>0) cout<<" - Load datasets ..."<<endl;
  vector<Dataset> datasets;
  LoadDatasets(datasets);
  vector<Dataset> LMX;
  if(addLMX){
  	LoadDatasetsLMX(datasets);
  	LoadDatasetsLMX(LMX);
  }
  /////////////////////
 
  ////////////////////////////
  //Selection table
  ////////////////////////////
  
  if(verbose>0) cout<<" - Variable decleration ..."<<endl;

  //vector of objects
  vector<TRootMuon> init_muons;
  vector<TRootElectron> init_electrons;
  vector<TRootJet> init_jets;
  vector<TRootMET> mets;
 
  TBranch *jets_br = 0;
  TClonesArray *tcjets;
  TBranch *muons_br = 0;
  TClonesArray *tcmuons;
  TBranch *electrons_br = 0;
  TClonesArray *tcelectrons;
  TBranch *mets_br = 0;
  TClonesArray *tcmets;
  TBranch *genEvt_br = 0;
  TClonesArray *tcgenEvt;
  TBranch* run_br = 0;
  TRootRun* runInfos = 0;
  TBranch* event_br = 0;
  TRootEvent* event = 0;



 ////////////////////////////////////
 /// Conditions 
 ////////////////////////////////////

  bool doJES = false;
  float JESFactor = 1.1;
  for(int i=1;i<argc;i++){
  	if(!strcmp(argv[i],"-JES") && i<argc-1) JESFactor = atof(argv[i+1]); 
  }
  
  bool doWJetEstimation = true;
	

  //b-tagging
  int btaggingWP = 0;
  int btaggingAlgo = 1;
  for(int i=1;i<argc;i++){
  	if(!strcmp(argv[i],"-btag") && i<argc-1) btaggingAlgo = atoi(argv[i+1]); 
  	if(!strcmp(argv[i],"-WP") && i<argc-1) btaggingWP = atoi(argv[i+1]); 
  }
  //0: loose - 1: medium - 2: tight
  int NofBtagAlgo = 5;
  //0: btag_trackCountingHighEffBJetTags
  //1: btag_trackCountingHighPurBJetTags
  //2: btag_jetProbabilityBJetTags
  //3: btag_jetBProbabilityBJetTags
  //4: btag_simpleSecondaryVertexBJetTags 
  float** WorkingPoint = new float*[NofBtagAlgo];
  for(int i=0;i<NofBtagAlgo;i++) WorkingPoint[i] = new float[3];
  WorkingPoint[0][0] = 2.03 ; WorkingPoint[0][1] = 4.38; WorkingPoint[0][2] = 5.36;
  WorkingPoint[1][0] = 2.03 ; WorkingPoint[1][1] = 4.38; WorkingPoint[1][2] = 5.36;
  WorkingPoint[2][0] = 0.241 ; WorkingPoint[2][1] = 0.49; WorkingPoint[2][2] = 0.795;
  WorkingPoint[3][0] = 1.1 ; WorkingPoint[3][1] = 1.37; WorkingPoint[3][2] = 1.39;
  WorkingPoint[4][0] = 1.25 ; WorkingPoint[4][1] = 2.05; WorkingPoint[4][2] = 4.07;
  //WorkingPoint[5][0] = 0.387 ; WorkingPoint[5][1] = 0.838; WorkingPoint[5][2] = 0.94;
  //////////////////////


 ////////////////////////////////////
 /// Selection 
 ////////////////////////////////////
  
   float ElectronPtCut = 10;
   float ElectronEtaCut = 2.0;
   float ElectronRelIsoCut = 0.1;
   //SR
   float MuonPtCutSR = 20;//20
   float MuonEtaCutSR = 2.1;
   float MuonRelIsoCutSR = 0.1;
   float JetsPtCut3SR = 30;//75
   float JetsPtCut4SR = 30;//50
   float JetsEtaCutSR = 2.4;
   //CR: ttjetEstimation
   float MuonPtCutCR = 20;
   float MuonEtaCutCR = 2.1;
   float MuonRelIsoCutCR = 0.1;
   float JetsPtCut3CR = 30;//75
   float JetsPtCut4CR = 30;
   float JetsEtaCutCR = 2.4;
  for(int i=1;i<argc;i++){
  	if(!strcmp(argv[i],"-PtJets") && i<argc-1){
		JetsPtCut3SR = atoi(argv[i+1]); 
		JetsPtCut3CR = atoi(argv[i+1]); 
  	}
	if(!strcmp(argv[i],"-PtMuon") && i<argc-1){
		MuonPtCutSR = atoi(argv[i+1]); 
		MuonPtCutCR = atoi(argv[i+1]); 
  	}
  }

 ////////////////////////////////////
 /// Selection Table
 ////////////////////////////////////
  if(verbose>0) cout<<" - Configuration of SelectionTable ..."<<endl;
 //Signal Region: selection table
  vector<string> CutsSelecTable;
  CutsSelecTable.push_back(string("all"));
  CutsSelecTable.push_back(string("trigged"));
  CutsSelecTable.push_back(string("1 muon"));
  CutsSelecTable.push_back(string("0 electron"));
  CutsSelecTable.push_back(string("4 jets"));
  CutsSelecTable.push_back(string("MET> 100 GeV"));
  SelectionTable selecTable(CutsSelecTable, datasets);
  selecTable.SetLuminosity(Luminosity);
  //Control Region: selection table
  vector<string> CutsSelecTableCR;
  CutsSelecTableCR.push_back(string("all"));
  CutsSelecTableCR.push_back(string("trigged"));
  CutsSelecTableCR.push_back(string("1 muon"));
  CutsSelecTableCR.push_back(string("0 electron"));
  CutsSelecTableCR.push_back(string("4 jets"));
  CutsSelecTableCR.push_back(string("2 b-jets"));
  CutsSelecTableCR.push_back(string("Mbl>"));
  CutsSelecTableCR.push_back(string("HLTBB>"));
  CutsSelecTableCR.push_back(string("DRBB>"));
  CutsSelecTableCR.push_back(string("MET>100 GeV"));
  SelectionTable selecTableCR(CutsSelecTableCR, datasets);
  selecTableCR.SetLuminosity(Luminosity);
 
 ////////////////////////////////////
 // Estimation per #jets & #b-jets
 ///////////////////////////////////
 EstimationPerJets_BJetsMultiplicity EstimJetsBJets;
 
 ////////////////////////////////////
 /// TtJetEstimation
 ////////////////////////////////////
  if(verbose>0) cout<<" - Configuration of TTJetEstimation ..."<<endl;
  //Declare TtJetEstimation
  TtJetEstimation ttjEstimation;
  // This values can be tuned !!
  //pair<float,float> NRegion(0.,100.); 
  pair<float,float> NRegion(120.,600.); 
  float METCut = 100; 
  //Varibles used to define the Control Region
  //float BtagCut = 4.38;
  float MblCut = 160;//150
  float DRBBCut = 2.3;//2.0
  float HTBBCut = 500.;//250
  //////////////////////////////
  vector<Dataset> datasetsBkg;
  vector<Dataset> datasetsTTJets;
  vector<Dataset> datasetsNP;
  for(unsigned int d=0;d<datasets.size();d++){
	if(datasets[d].Name()=="QCD" || datasets[d].Name()=="WJets" || datasets[d].Name()=="SingleTop") datasetsBkg.push_back(datasets[d]);
	if(datasets[d].Name()=="TTJets") datasetsTTJets.push_back(datasets[d]); 
	//if(datasets[d].Name()=="SUSY") datasetsNP.push_back(datasets[d]);
	if(datasets[d].Name().find("LM")<datasets[d].Name().size()) datasetsNP.push_back(datasets[d]);
  }
  //follow that order !!
  ttjEstimation.SetListOfDatasets(datasetsBkg, datasetsTTJets, datasetsNP);
  ttjEstimation.ConfigHistos(NbinsX, binsX, VarLabel, XaxisLabel, true);


  /////////////////////
  ///  Plots for ttbar
  /////////////////////

  MultiSamplePlot msplotMET(datasets, string("Met"), 50, 0., 500., string("MET [GeV]"));
  MultiSamplePlot msplotMbl(datasets, string("Mbl"), 50, 0., 500., string("M(b,l)"));
  MultiSamplePlot msplotDRBB(datasets, string("DRBB"), 50, 0., 5., string("#DeltaR(b,b)"));
  MultiSamplePlot msplotHTBB(datasets, string("HTBB"), 120, 0., 1200., string("H_{T}(b,b)"));
  MultiSamplePlot msplotDRbl(datasets, string("DRbl"), 50, 0., 5., string("#DeltaR(b,l)"));
  MultiSamplePlot msplotDPhiMETlepton(datasets, string("DPhiMETlepton"), 64, 0., 3.2, string("#Phi(MET,l)"));
  MultiSamplePlot msplotDEtaBB(datasets, string("DEtaBB"), 50, 0.,5. , string("#Delta#Eta(b,b)"));
  MultiSamplePlot msplotDPhiBB(datasets, string("DPhiBB"), 64, 0., 3.2, string("#Delta#phi(b,b)"));

  vector<pair<string,int> > METCutsList;
  METCutsList.push_back(pair<string,int>(string("SR"),1));
  METCutsList.push_back(pair<string,int>(string("2 b-jets"),2));
  METCutsList.push_back(pair<string,int>(string("Mbl"),3));
  METCutsList.push_back(pair<string,int>(string("DRBB"),4));
  METCutsList.push_back(pair<string,int>(string("HTBB"),5));
  //MultiCutPlot mcplotMET(METCutsList, string("MET"), 20, 0., 500., string("MET [GeV]"));
  MultiCutPlot mcplotMET(METCutsList, string("MET_Cuts"), NbinsX, binsX, string("MET [GeV]"));
  MultiCutPlot mcplotMET2(METCutsList, string("MET_Cuts_OneByOne"), NbinsX, binsX, string("MET [GeV]"));

  //semi - vs  di
  vector<pair<string,int> > ttbarComponents;
  ttbarComponents.push_back(pair<string,int>(string("Inclusif - SR"),1));
  ttbarComponents.push_back(pair<string,int>(string("Inclusif - CR"),2));
  ttbarComponents.push_back(pair<string,int>(string("Semi-leptonic - SR"),3));
  ttbarComponents.push_back(pair<string,int>(string("Semi-leptonic - CR"),4));
  ttbarComponents.push_back(pair<string,int>(string("Di-leptonic - SR"),5));
  ttbarComponents.push_back(pair<string,int>(string("Di-leptonic - CR"),6));
  //ttbarComponents.push_back(pair<string,int>(string("Hadronic - SR"),7));
  //ttbarComponents.push_back(pair<string,int>(string("Hadronic - CR"),8));
  MultiCutPlot mcplotMETttbarComponents(ttbarComponents,string("MET"), NbinsX, binsX, string("MET [GeV]"));

  //effect of the cuts in CR on shape
  CutImpactEvaluation cimpeMbl(string("cimpeMbl"), NbinsX, binsX, 50, 0., 500., LMX, string("Mass(b,l) [GeV]"), false);
  CutImpactEvaluation cimpeHTBB(string("cimpeHTBB"), NbinsX, binsX, 50, 0., 1000., LMX, string("H_{T}(b,b) [GeV]"), false);
  CutImpactEvaluation cimpeDRBB(string("cimpeDRBB"), NbinsX, binsX, 38, 0., 3.8, LMX, string("#DeltaR(b,b)"), true);


 ////////////////////////////////////
 /// WEstimation
 ////////////////////////////////////
  if(verbose>0) cout<<" - Configuration of WJetEstimation ..."<<endl;
  //Declare TtJetEstimation
  //parameters
  //Can be tuned !!
  int Iterations = 40;
  //
  vector<int> iDTTLike;
  vector<int> iDWLike;
  for(unsigned int d=0;d<datasets.size();d++){
	if(datasets[d].Name()=="TTJets") iDTTLike.push_back(d); 
	if(datasets[d].Name()=="WJets") iDWLike.push_back(d); 
  }
  WJetEstimation wjEstimation_SR(Iterations, (int) datasets.size(), iDTTLike, iDWLike);

  ////////////////////////////////
  int NofJetBins = 3;
  double** MultiJets_Estimated_Nj_SR = new double*[NofJetBins];
  for(int i=0;i<NofJetBins;i++){
    	MultiJets_Estimated_Nj_SR[i] = new double[4];
	for(unsigned int j=0;j<4;j++){ MultiJets_Estimated_Nj_SR[i][j] = 0 ;}
  }
  //nof b-tagged jets per jet multiplicity and per datasets
  double** N0_SR = new double*[NofJetBins];
  double** N1_SR = new double*[NofJetBins];
  double** N2_SR = new double*[NofJetBins];
  double** N3_SR = new double*[NofJetBins];
  for(int i=0;i<NofJetBins;i++){
  	N0_SR[i] = new double[datasets.size()];
  	N1_SR[i] = new double[datasets.size()];
  	N2_SR[i] = new double[datasets.size()];
  	N3_SR[i] = new double[datasets.size()];
  	for(unsigned int j=0;j<datasets.size();j++){
  		N0_SR[i][j] = 0; N1_SR[i][j] = 0; N2_SR[i][j] = 0; N3_SR[i][j] = 0;
  	}
  }
  ///////////////////////////////////////////
  ///  Iterative Shape Estimator for WJet
  ///////////////////////////////////////////
  //histo for the variable needed - here: MET
  //CR
  TH1F* histo_inclusif_CR = new TH1F("histo_inclusif_CR","histo_inclusif", NbinsX, binsX);
  TH1F* histo_0bj_CR = new TH1F("histo_0bj_CR","histo_0bj", NbinsX, binsX);
  TH1F* histo_1bj_CR = new TH1F("histo_1bj_CR","histo_1bj", NbinsX, binsX);
  TH1F* histo_Wjets_CR = new TH1F("histo_Wjets_CR","histo_Wjets", NbinsX, binsX);
  TH1F* histo_ttjets_CR = new TH1F("histo_ttjets_CR","histo_ttjets", NbinsX, binsX);
  //SR
  TH1F* histo_inclusif_SR = new TH1F("histo_inclusif_SR","histo_inclusif", NbinsX, binsX);
  TH1F* histo_0bj_SR = new TH1F("histo_0bj_SR","histo_0bj", NbinsX, binsX);
  TH1F* histo_1bj_SR = new TH1F("histo_1bj_SR","histo_1bj", NbinsX, binsX);
  TH1F* histo_Wjets_SR = new TH1F("histo_Wjets_SR","histo_Wjets", NbinsX, binsX);
  TH1F* histo_ttjets_SR = new TH1F("histo_ttjets_SR","histo_ttjets", NbinsX, binsX);
 
  //to get proper errors from MC expected stat
  //CR
  histo_inclusif_CR->Sumw2();
  histo_0bj_CR->Sumw2();
  histo_1bj_CR->Sumw2();
  histo_Wjets_CR->Sumw2();
  histo_ttjets_CR->Sumw2();
  //SR
  histo_inclusif_SR->Sumw2();
  histo_0bj_SR->Sumw2();
  histo_1bj_SR->Sumw2();
  histo_Wjets_SR->Sumw2();
  histo_ttjets_SR->Sumw2();
 
 ////////////////////////////////////
 /// ABCD Estimation for QCD
 ////////////////////////////////////
  if(verbose>0) cout<<" - Configuration of ABCDEstimation ..."<<endl;
  int NXbinsABCD = 200;
  int NYbinsABCD = 200;
  float XbinMinABCD = 0.;
  float XbinMaxABCD = 20.;
  float YbinMinABCD = 0.;
  float YbinMaxABCD = 20.;//20
  ABCDEstimation abcdEstim(TString("ABCD"), TString("ABCD"), TString("RelIso"), TString("IP Significance"), NXbinsABCD, XbinMinABCD, XbinMaxABCD, NYbinsABCD, YbinMinABCD,YbinMaxABCD);
  //ABCDEstimation abcdEstim(TString("ABCD"), TString("ABCD"), TString("RelIso"), TString("d0"), NXbinsABCD, XbinMinABCD, XbinMaxABCD, NYbinsABCD, YbinMinABCD,YbinMaxABCD);

  TH2F* h_ABCD_QCD = new TH2F("h_ABCD_QCD","",NXbinsABCD,XbinMinABCD,XbinMaxABCD,NYbinsABCD,YbinMinABCD,YbinMaxABCD);	
  TH2F* h_ABCD_WJets = new TH2F("h_ABCD_WJets","",NXbinsABCD,XbinMinABCD,XbinMaxABCD,NYbinsABCD,YbinMinABCD,YbinMaxABCD);	
  TH2F* h_ABCD_TTJets = new TH2F("h_ABCD_TTJets","",NXbinsABCD,XbinMinABCD,XbinMaxABCD,NYbinsABCD,YbinMinABCD,YbinMaxABCD);	
	
  int QCDNbinsX = 50;
  float*  QCDbinsX = new float[QCDNbinsX+1];
  float QCDxMax = 1.;
  for(int i=0;i<QCDNbinsX+1;i++) QCDbinsX[i] = i*QCDxMax/QCDNbinsX;
  int QCDNbinsY = 50;
  float*  QCDbinsY = new float[QCDNbinsX+1];
  float QCDyMax = 10.;
  for(int i=0;i<QCDNbinsY+1;i++) QCDbinsY[i] = i*QCDyMax/QCDNbinsY;
  //parameters to compute estimation
  //cut on RelIso
  Float_t cutXmin = 0.0;  Float_t cutX0   = 0.1;  Float_t cutX1   = 0.2;  Float_t cutXmax = 20.0;  //0., 0.1, 0.1, 1.
  //cut on d0Sign
  Float_t cutYmin = 0.0;  Float_t cutY0   = 3;  Float_t cutY1   = 4;  Float_t cutYmax = 20.0; //0., 0.2, 0.2, 1.  //0.0 3.0 5.0 15.0
  //Float_t cutYmin = 0.0;  Float_t cutY0   = 0.2;  Float_t cutY1   = 0.2;  Float_t cutYmax = 2.0; //0., 0.2, 0.2, 1.  //0.0 3.0 5.0 15.0
  Int_t region    = 1;   // regions go from 1 to 4
  
  //ABCDEstimation abcdEstim(TString("ABCD"), TString("ABCD"), TString("RelIso"), TString("IP Significance"), QCDNbinsX, QCDbinsX, QCDNbinsY, QCDbinsY, NbinsX, binsX);

 ////////////////////////////////////
 /// QCDShapeEstimation
 ////////////////////////////////////
  if(verbose>0) cout<<" - Configuration of QCDShapeEstimation ..."<<endl;
  //QCDShapeEstimation qcdShapeEstim (100, 0., 1., string("MET"));
  QCDShapeEstimation qcdShapeEstim (NbinsX, binsX, XaxisLabel);
  vector<float> relIsoVec;
  float relIsoVecWidth = 0.1;
  float relIsoVecNof = 20;
  for(int i=0;i<relIsoVecNof;i++) relIsoVec.push_back(relIsoVecWidth*i); 
  qcdShapeEstim.SetControlShape(relIsoVec);

 ////////////////////////////////////
 /// Plots used by BkgEstimationSummary
 ////////////////////////////////////
  if(verbose>0) cout<<" - Create histos used by BkgEstimationSummary ..."<<endl;
  TH1F* hData = new TH1F("hData","Data",NbinsX, binsX);
  TH1F* hQCDMC = new TH1F("hQCDMC","QCD - MC",NbinsX, binsX);
  TH1F* hWJetMC = new TH1F("hWJetMC","WJets - MC",NbinsX, binsX);
  TH1F* hTtJetMC = new TH1F("TtJetMC","TtJets - MC",NbinsX, binsX);
 
  //////////////////////////////////////////////
  //loop on datasets
  if(verbose>0) cout<<" - Loop over datasets ... "<<datasets.size()<<" datasets !"<<endl;
  for(unsigned int d=0;d<datasets.size();d++){
  //for(unsigned int d=0;d<1;d++){
  	if(verbose>1) cout<<"   Dataset "<<d<<": "<<datasets[d].Name()<<endl;
  	//open files and load
	TChain* eventTree = datasets[d].eventTree();
	TChain* runTree = datasets[d].runTree();

  	///////////////////////////////////////
  	//Branches and TCLonesArray
  	///////////////////////////////////////
  	//jets
  	jets_br = (TBranch *) eventTree->GetBranch ("Jets");
  	tcjets = new TClonesArray ("TopTree::TRootJet", 0);
  	jets_br->SetAddress (&tcjets);
  	//muons
  	muons_br = (TBranch *) eventTree->GetBranch ("Muons");
  	tcmuons = new TClonesArray ("TopTree::TRootMuon", 0);
  	muons_br->SetAddress (&tcmuons);
  	//electrons
  	electrons_br = (TBranch *) eventTree->GetBranch ("Electrons");
  	tcelectrons = new TClonesArray ("TopTree::TRootElectron", 0);
  	electrons_br->SetAddress (&tcelectrons);
  	//mets
  	mets_br = (TBranch *) eventTree->GetBranch ("MET");
  	tcmets = new TClonesArray ("TopTree::TRootMET", 0);
  	mets_br->SetAddress (&tcmets);
	//GenEvent
	genEvt_br = (TBranch *) eventTree->GetBranch ("GenEvent");
  	tcgenEvt = new TClonesArray ("TopTree::TRootGenEvent", 0);
  	if(genEvt_br) genEvt_br->SetAddress (&tcgenEvt);
	//HLT info
	run_br = (TBranch *) runTree->GetBranch("runInfos");
	run_br->SetAddress(&runInfos);
	event_br = (TBranch *) eventTree->GetBranch("Event");
	event_br->SetAddress(&event);
	////////////////////////////////////////
	
	int itrigger = -1;
	runTree->GetEvent(d);//it's "d"
	for(int t=0;t<runInfos->hltNames().size();t++){
		if(runInfos->hltNames(t)==string("HLT_Mu9"))
			itrigger = t;
	}

	//loop on events
  	if(verbose>1) cout<<"	Loop over events "<<endl;
	for(unsigned int ievt=0;ievt<eventTree->GetEntries();ievt++){
		//if(ievt == 200) return 0 ;
		TLorentzVector JESbalance(0,0,0,0);
  		
		if(verbose>2) cout<<"Event "<<ievt<<endl;
	        eventTree->GetEvent(ievt);
	        
		// clear vectors	
		init_jets.clear();
		init_muons.clear();
		init_electrons.clear();
		mets.clear();
	  	
		TRootGenEvent* genEvt;
		
		//fill vectors
		for(int i=0;i<tcjets->GetEntriesFast();i++) init_jets.push_back(* (TRootJet*) tcjets->At(i));
		for(int i=0;i<tcmuons->GetEntriesFast();i++) init_muons.push_back(* (TRootMuon*) tcmuons->At(i));
		for(int i=0;i<tcelectrons->GetEntriesFast();i++) init_electrons.push_back(* (TRootElectron*) tcelectrons->At(i));
		for(int i=0;i<tcmets->GetEntriesFast();i++) mets.push_back(* (TRootMET*) tcmets->At(i));
  	        if(tcgenEvt->GetEntriesFast()==1) genEvt = (TRootGenEvent*) tcgenEvt->At(0);
	
		if(verbose>3) cout<<"Nof muons "<<init_muons.size()<<endl;
		if(verbose>3) cout<<"Nof electrons "<<init_electrons.size()<<endl;
		if(verbose>3) cout<<"Nof jets "<<init_jets.size()<<endl;
		
		bool isSemiLep = false;
		if(genEvt_br) isSemiLep = genEvt->isSemiLeptonic();
		bool isLeptonic = false;
		if(genEvt_br) isLeptonic = genEvt->isFullLeptonic();
		bool isHadronic = false;
		if(genEvt_br) isHadronic = genEvt->isFullHadronic();

		if(doJES){
			vector<TRootJet> jets_JES;
			jets_JES = JES(init_jets, JESFactor);
			for(int i=0;i<init_jets.size();i++){
				//JESbalance+=((TLorentzVector)jets_JES[i]-(TLorentzVector)init_jets[i]);
				double dPx = jets_JES[i].Px()-init_jets[i].Px();
				double dPy = jets_JES[i].Py()-init_jets[i].Py();
				TLorentzVector a(dPx,dPy,0,sqrt(dPx*dPx+dPy*dPy));
				JESbalance+=a;
				
			}
			//after balance calculation
			init_jets = jets_JES;
		}
	
		/////////////////////////////
		//   Selection
		/////////////////////////////
		//Declare selection instance	
		Selection selection(init_jets, init_muons, init_electrons, mets);
		//selecTable.Fill(selection, d);//WARNING: to be implemented
		//selecTableCR.Fill(selection, d); //WARNING: to be defined & implemented

		//variable on which on want to make bgk estimation  & search for New Physics
		float variable = 0;
		//can be changed
		if(isMET) variable = mets[0].Et();
		if(isHT){
			variable = 0;
			if(selection.GetSelectedJets(JetsPtCut4SR,JetsEtaCutSR).size()>=4)
				for(unsigned int i=0;i<4;i++)
					variable+=selection.GetSelectedJets(JetsPtCut4SR,JetsEtaCutSR)[i].Pt();
		}
		if(isMttbar){
			variable = 0;
			if(selection.GetSelectedJets(JetsPtCut4SR,JetsEtaCutSR).size()>=4 && selection.GetSelectedMuons(MuonPtCutSR,MuonEtaCutSR,MuonRelIsoCutSR).size()>=1)	
			variable = Mttbar(selection.GetSelectedMuons(MuonPtCutSR,MuonEtaCutSR,MuonRelIsoCutSR)[0],selection.GetSelectedJets(JetsPtCut4SR,JetsEtaCutSR),mets[0]);
			
		}
	
		if(doJES) {
			double dPx = mets[0].Px()-JESbalance.Px();
			double dPy = mets[0].Py()-JESbalance.Py();
			TLorentzVector newMET(dPx,dPy,0,sqrt(dPx*dPx+dPy*dPy)) ;
			variable = newMET.Et();
		}
		//cout<<"variable: "<< mets[0].Et()<<" new "<< variable<<endl;
		//variable = selection.GetSelectedJets(30,2.4)[0].Pt();


		/////////////////////////////
		//   TtJets estimation
		/////////////////////////////
  		
		bool trigged = false;
		//cout<<"nof triggers "<<runInfos->hltNames().size()<<endl;
		if(event->trigHLT(itrigger)) trigged = true; // correspond to HLT_Mu9
		if(!trigged) continue;
		//the different regions
		bool isQCD = false;
		bool isInQCDCR = false;
		bool isInABCD = false;
		bool isInCR = false;
		bool isInSR = false;
		float RelIso = 100.;
		//////////////////////////

		
		vector<TRootJet> selectedJets = selection.GetSelectedJets(JetsPtCut4CR,JetsEtaCutSR);
		//btag jets
		vector<TRootJet> btagJets;
		for(unsigned int i=0;i<selectedJets.size();i++){
			if(btaggingAlgo == 0 && selectedJets[i].btag_trackCountingHighEffBJetTags()>WorkingPoint[0][btaggingWP]) btagJets.push_back(selectedJets[i]);
			if(btaggingAlgo == 1 && selectedJets[i].btag_trackCountingHighPurBJetTags()>WorkingPoint[1][btaggingWP]) btagJets.push_back(selectedJets[i]);
			if(btaggingAlgo == 2 && selectedJets[i].btag_jetProbabilityBJetTags()>WorkingPoint[2][btaggingWP]) btagJets.push_back(selectedJets[i]);
			if(btaggingAlgo == 3 && selectedJets[i].btag_jetBProbabilityBJetTags()>WorkingPoint[3][btaggingWP]) btagJets.push_back(selectedJets[i]);
			if(btaggingAlgo == 4 && selectedJets[i].btag_simpleSecondaryVertexBJetTags()>WorkingPoint[4][btaggingWP]) btagJets.push_back(selectedJets[i]);
		}
		////////////////////


	        ///////////////////////////////////////////
		//fill booleans for the different region
	        ///////////////////////////////////////////
	
		selecTable.Fill(d,0,1.);
		if(trigged) selecTable.Fill(d,1,1.);
		if(trigged & selection.GetSelectedMuons(MuonPtCutSR,MuonEtaCutSR,MuonRelIsoCutSR).size()==1) selecTable.Fill(d,2,1.);
		if(trigged & selection.GetSelectedMuons(MuonPtCutSR,MuonEtaCutSR,MuonRelIsoCutSR).size()==1 && selection.GetSelectedElectrons(ElectronPtCut,ElectronEtaCut,ElectronRelIsoCut).size()==0) selecTable.Fill(d,3,1.);
		if(trigged && selection.GetSelectedMuons(MuonPtCutSR,MuonEtaCutSR,MuonRelIsoCutSR).size()==1 && selection.GetSelectedElectrons(ElectronPtCut,ElectronEtaCut,ElectronRelIsoCut).size()==0 && selection.GetSelectedJets(JetsPtCut4SR,JetsEtaCutSR).size()>=4 && selection.GetSelectedJets(JetsPtCut3SR,JetsEtaCutSR).size()>=3){
			selecTable.Fill(d,4,1.);
			RelIso = 1/selection.GetSelectedMuons(MuonPtCutSR,MuonEtaCutSR,MuonRelIsoCutSR)[0].relativeIso03();
			isInSR = true;
			if(datasets[d].Name()=="TTJets") mcplotMET.Fill(variable,METCutsList[0].first);
			if(datasets[d].Name()=="TTJets") mcplotMET2.Fill(variable,METCutsList[0].first);
			if(datasets[d].Name()=="TTJets") {
				mcplotMETttbarComponents.Fill(variable,ttbarComponents[0].first);
				if(isSemiLep) mcplotMETttbarComponents.Fill(variable,ttbarComponents[2].first); 
				if(isLeptonic) mcplotMETttbarComponents.Fill(variable,ttbarComponents[4].first); 
				//if(isHadronic) mcplotMETttbarComponents.Fill(variable,ttbarComponents[6].first); 
			}
			if(mets[0].Et()>METCut) selecTable.Fill(d,5,1.);
		}
		if(trigged && selection.GetSelectedMuonsInvIso(MuonPtCutSR,MuonEtaCutSR,MuonRelIsoCutSR).size()==1  && selection.GetSelectedJets(JetsPtCut4SR,JetsEtaCutSR).size()>=4){
			isInQCDCR = true;	
			if(!isInSR) RelIso = 1/selection.GetSelectedMuonsInvIso(MuonPtCutSR,MuonEtaCutSR,MuonRelIsoCutSR)[0].relativeIso03();
			if(!isInSR && isMttbar) variable = Mttbar(selection.GetSelectedMuonsInvIso(MuonPtCutSR,MuonEtaCutSR,MuonRelIsoCutSR)[0],selection.GetSelectedJets(JetsPtCut4SR,JetsEtaCutSR),mets[0]);
		}
		if(trigged && selection.GetSelectedMuonsQCDEstim(MuonPtCutSR,MuonEtaCutSR).size()==1 && selection.GetSelectedJets(JetsPtCut4SR,JetsEtaCutSR).size()>=4){
			isInABCD = true;
		}
		//Control Region for TtJetEstimation
		selecTableCR.Fill(d,0,1.);
		if(trigged) selecTableCR.Fill(d,1,1.);
		if(trigged & selection.GetSelectedMuons(MuonPtCutCR,MuonEtaCutCR,MuonRelIsoCutCR).size()==1) selecTableCR.Fill(d,2,1.);
		if(trigged & selection.GetSelectedMuons(MuonPtCutCR,MuonEtaCutCR,MuonRelIsoCutCR).size()==1 && selection.GetSelectedElectrons(ElectronPtCut,ElectronEtaCut,ElectronRelIsoCut).size()==0) selecTableCR.Fill(d,3,1.);
		if(trigged && selection.GetSelectedMuons(MuonPtCutCR,MuonEtaCutCR,MuonRelIsoCutCR).size()==1 && selection.GetSelectedElectrons(ElectronPtCut,ElectronEtaCut,ElectronRelIsoCut).size()==0 && selection.GetSelectedJets(JetsPtCut4CR,JetsEtaCutCR).size()>=4 && selection.GetSelectedJets(JetsPtCut3CR,JetsEtaCutCR).size()>=3){
			selecTableCR.Fill(d,4,1.);
			if(btagJets.size()>1){
				selecTableCR.Fill(d,5,1.);
				if(datasets[d].Name()=="TTJets") mcplotMET.Fill(variable,METCutsList[1].first);
				if(datasets[d].Name()=="TTJets") mcplotMET2.Fill(variable,METCutsList[1].first);
				
				
				TRootMuon muon = selection.GetSelectedMuons(MuonPtCutCR,MuonEtaCutCR,0.1)[0];
				//compute Mbl
				float Mbl = -9999;
				float DRblMin = 9999;
				for(unsigned int bj=0;bj<btagJets.size();bj++){
					if(bj>1) break;
					float DeltaR = ROOT::Math::VectorUtil::DeltaR (btagJets[bj],muon);
					if (DeltaR<DRblMin){
						DRblMin = DeltaR;
						TLorentzVector BL = btagJets[bj]+selection.GetSelectedMuons(MuonPtCutCR,MuonEtaCutCR,MuonRelIsoCutCR)[0];
						Mbl = BL.M();
					}
				}

				float HTBB = btagJets[0].Pt()+btagJets[1].Pt();
				float DRBB = ROOT::Math::VectorUtil::DeltaR(btagJets[0],btagJets[1]);
				float DEtaBB = fabs(btagJets[0].Eta()-btagJets[1].Eta());
				float DPhiBB = fmod(fabs(btagJets[0].Phi()-btagJets[1].Phi()),M_PI);
				float DPhiMETlepton = fmod(fabs(mets[0].Phi()-muon.Phi()),M_PI);
				if(Mbl<MblCut && HTBB<HTBBCut && DRBB>DRBBCut){
				//if(DRBB>DRBBCut){
					isInCR = true;
					if(isMttbar) variable = Mttbar(selection.GetSelectedMuons(MuonPtCutCR,MuonEtaCutCR,MuonRelIsoCutCR)[0],selection.GetSelectedJets(JetsPtCut4CR,JetsEtaCutCR),mets[0]);
					if(datasets[d].Name()=="TTJets"){ 
						mcplotMETttbarComponents.Fill(variable,ttbarComponents[1].first); 
						if(isSemiLep) mcplotMETttbarComponents.Fill(variable,ttbarComponents[3].first); 
						if(isLeptonic) mcplotMETttbarComponents.Fill(variable,ttbarComponents[5].first); 
						//if(isHadronic) mcplotMETttbarComponents.Fill(variable,ttbarComponents[7].first); 
					}
				}
				msplotMET.Fill(variable,datasets[d]);
				msplotMbl.Fill(Mbl,datasets[d]);
				msplotDRBB.Fill(DRBB,datasets[d]);
				msplotDEtaBB.Fill(DEtaBB,datasets[d]);
				msplotDPhiBB.Fill(DPhiBB,datasets[d]);
				msplotHTBB.Fill(HTBB,datasets[d]);
				msplotDRbl.Fill(DRblMin,datasets[d]);
				msplotDPhiMETlepton.Fill(DPhiMETlepton,datasets[d]);
				
				if(datasets[d].Name()=="TTJets"){ 
					if(Mbl<MblCut) {
						mcplotMET.Fill(variable,METCutsList[2].first);
						mcplotMET2.Fill(variable,METCutsList[2].first);
					}
					if(Mbl<MblCut && DRBB>DRBBCut){
						mcplotMET.Fill(variable,METCutsList[3].first);
					}
					if(Mbl<MblCut &&  DRBB>DRBBCut && HTBB<HTBBCut){
					//if(DRBB>DRBBCut){
						mcplotMET.Fill(variable,METCutsList[4].first);
					}
					if(DRBB>DRBBCut) mcplotMET2.Fill(variable,METCutsList[3].first);
					if(HTBB<HTBBCut) mcplotMET2.Fill(variable,METCutsList[4].first);
				}
				if(Mbl<MblCut) selecTableCR.Fill(d,6,1.);
				if(Mbl<MblCut && DRBB>DRBBCut)  selecTableCR.Fill(d,7,1.);
				if(Mbl<MblCut &&  DRBB>DRBBCut && HTBB<HTBBCut){
					selecTableCR.Fill(d,8,1.);
					if(mets[0].Et()>METCut) selecTableCR.Fill(d,9,1.);
				}
				
				
				//study the effects on the cuts applied
				cimpeMbl.Fill(Mbl,variable);
				cimpeHTBB.Fill(HTBB,variable);
				cimpeDRBB.Fill(DRBB,variable);
				//for LMX
				if(datasets[d].Name().find("LM")<datasets[d].Name().size()){
					unsigned int iDat = 0;
					for(unsigned int x=0;x<LMX.size();x++) if(LMX[x].Name() == datasets[d].Name()) iDat = x;
					
					cimpeMbl.Fill(Mbl,variable, true, iDat);
					cimpeHTBB.Fill(HTBB,variable, true, iDat);
					cimpeDRBB.Fill(DRBB,variable, true, iDat);
				}
			}
		}
		///////////////////////////////////////////
		
		////////////////////////////////////
		//  WJetEstimation: norm & shape
		///////////////////////////////////
		
		//Signal Region
		if(isInSR){
			//for WJetEstimation
			if(isMttbar) variable = Mttbar(selection.GetSelectedMuons(MuonPtCutSR,MuonEtaCutSR,MuonRelIsoCutSR)[0],selection.GetSelectedJets(JetsPtCut4SR,JetsEtaCutSR),mets[0]);
			int injets = 0;
			if(selection.GetSelectedJets(JetsPtCut4SR,JetsEtaCutSR).size()==4) injets = 0;
			if(selection.GetSelectedJets(JetsPtCut4SR,JetsEtaCutSR).size()==5) injets = 1;
			if(selection.GetSelectedJets(JetsPtCut4SR,JetsEtaCutSR).size()>=6) injets = 2;
			if(btagJets.size()==0) N0_SR[injets][d]++;
			if(btagJets.size()==1) N1_SR[injets][d]++; 
			if(btagJets.size()==2) N2_SR[injets][d]++; 
			if(btagJets.size()>=3) N3_SR[injets][d]++; 
			//Iterative Shape Extrapolator
			if(btagJets.size()==0) histo_0bj_SR->Fill(variable,datasets[d].NormFactor()*Luminosity);
			if(btagJets.size()==1) histo_1bj_SR->Fill(variable,datasets[d].NormFactor()*Luminosity);
			histo_inclusif_SR->Fill(variable,datasets[d].NormFactor()*Luminosity);
			if(datasets[d].Name()=="WJets") histo_Wjets_SR->Fill(variable,datasets[d].NormFactor()*Luminosity);
			if(datasets[d].Name()=="TTJets") histo_ttjets_SR->Fill(variable,datasets[d].NormFactor()*Luminosity);
		}
		
		//Control Region
		if(isInCR){
			//Iterative Shape Extrapolator
			if(btagJets.size()==0) histo_0bj_CR->Fill(variable,datasets[d].NormFactor()*Luminosity);
			if(btagJets.size()==1) histo_1bj_CR->Fill(variable,datasets[d].NormFactor()*Luminosity);
			histo_inclusif_CR->Fill(variable,datasets[d].NormFactor()*Luminosity);
			if(datasets[d].Name()=="WJets") histo_Wjets_CR->Fill(variable,datasets[d].NormFactor()*Luminosity);
			if(datasets[d].Name()=="TTJets") histo_ttjets_CR->Fill(variable,datasets[d].NormFactor()*Luminosity);
		}
		///////////////////////
	
		///////////////////////
		//  QCDShapeEstimation
		///////////////////////
		if(datasets[d].Name()=="QCD") isQCD = true;
		//if(selection.GetSelectedMuonsInvIso(30,2.1,0.1).size()>0) 
		qcdShapeEstim.Fill( variable, datasets[d].NormFactor()*Luminosity,RelIso, isQCD, isInSR, isInQCDCR);
		//if(isQCD && trigged && isInSR) cout<<"QCD event trigged"<<endl;
		///////////////////////
		
		///////////////////////
		//  ABCDEstimation
		///////////////////////
		if(isInABCD) {
			float d0 = selection.GetSelectedMuonsQCDEstim(MuonPtCutSR,MuonEtaCutSR)[0].d0();
			float d0Signvalue = 100.;
			if(selection.GetSelectedMuonsQCDEstim(MuonPtCutSR,MuonEtaCutSR)[0].d0error()>0) d0Signvalue = selection.GetSelectedMuonsQCDEstim(MuonPtCutSR,MuonEtaCutSR)[0].d0()/selection.GetSelectedMuonsQCDEstim(MuonPtCutSR,MuonEtaCutSR)[0].d0error();
			RelIso = 100.;
			if(selection.GetSelectedMuonsQCDEstim(MuonPtCutSR,MuonEtaCutSR)[0].relativeIso03()>0) RelIso = 1/selection.GetSelectedMuonsQCDEstim(MuonPtCutSR,MuonEtaCutSR)[0].relativeIso03();
			abcdEstim.Fill(RelIso,d0Signvalue,datasets[d].NormFactor()*Luminosity, isQCD);
			//abcdEstim.Fill(RelIso,d0,datasets[d].NormFactor()*Luminosity, isQCD);
			if(datasets[d].Name()=="QCD") h_ABCD_QCD->Fill(RelIso,d0Signvalue);
			if(datasets[d].Name()=="WJets") h_ABCD_WJets->Fill(RelIso,d0Signvalue);
			if(datasets[d].Name()=="TTJets") h_ABCD_TTJets->Fill(RelIso,d0Signvalue);
		}
		//per jets & b-jet multiplicity
		float X = -999.;
		float Y = -9999.;
		if(isInABCD){
			float relIso = 100.;
			if(selection.GetSelectedMuonsQCDEstim(MuonPtCutSR,MuonEtaCutSR)[0].relativeIso03()>0) relIso = 1/selection.GetSelectedMuonsQCDEstim(MuonPtCutSR,MuonEtaCutSR)[0].relativeIso03();
			float d0Signvalue = 100.;
			if(selection.GetSelectedMuonsQCDEstim(MuonPtCutSR,MuonEtaCutSR)[0].d0error()>0) d0Signvalue = selection.GetSelectedMuonsQCDEstim(MuonPtCutSR,MuonEtaCutSR)[0].d0()/selection.GetSelectedMuonsQCDEstim(MuonPtCutSR,MuonEtaCutSR)[0].d0error();
			X = relIso;
			Y = d0Signvalue;
		}
		int nofJets = -1;
		int nofBJets = -1;
		if(selection.GetSelectedJets(JetsPtCut4SR,JetsEtaCutSR).size()>=6) nofJets = 6;
		else nofJets = selection.GetSelectedJets(JetsPtCut4SR,JetsEtaCutSR).size();
		if(btagJets.size()>=3) nofBJets = 3;
		else nofBJets = btagJets.size();

		//EstimJetsBJets.Fill( nofJets, nofBJets, X, Y, datasets[d], Luminosity, isInABCD, isInSR);
		///////////////////////

		///////////////////////
		//  TtJetEstimation
		///////////////////////
  		//Bgk
  		if(datasets[d].Name()=="QCD" || datasets[d].Name()=="WJets" || datasets[d].Name()=="SingleTop"){
			int idataset = -1;
			for(unsigned int dd=0;dd<datasetsBkg.size();dd++)
				if(datasets[d].Name() == datasetsBkg[dd].Name()) idataset = dd;
			ttjEstimation.FillByDatasetsBkg(variable,idataset);
  			if(isInSR) ttjEstimation.FillSRByDatasetsBkg(variable,idataset);
			if(isInCR) ttjEstimation.FillCRByDatasetsBkg(variable,idataset);
  		}
		//TtJets
		if(datasets[d].Name()=="TTJets"){
  			ttjEstimation.FillByDatasetsTtJets(variable,0);
  			if(isInSR) ttjEstimation.FillSRByDatasetsTtJets(variable,0);
  			if(isInCR) ttjEstimation.FillCRByDatasetsTtJets(variable,0);
		}
  		//SUSY
		//if(datasets[d].Name()=="SUSY"){
		if(datasets[d].Name().find("LM")<datasets[d].Name().size()){
			int id = 0;
			for(unsigned int ii=0;ii<datasetsNP.size();ii++) if(datasets[d].Name()==datasetsNP[ii].Name()) id = ii;
			ttjEstimation.FillByDatasetsNP(variable,id);
  			if(isInSR) ttjEstimation.FillSRByDatasetsNP(variable,id);
  			if(isInCR) ttjEstimation.FillCRByDatasetsNP(variable,id);
		}
		///////////////////////
	        if(isInSR){
		 	//Data = all datasets
			hData->Fill(variable,datasets[d].NormFactor()*Luminosity);
			if(datasets[d].Name()=="QCD") hQCDMC->Fill(variable,datasets[d].NormFactor()*Luminosity);
			if(datasets[d].Name()=="WJets") hWJetMC->Fill(variable,datasets[d].NormFactor()*Luminosity);
			if(datasets[d].Name()=="TTJets") hTtJetMC->Fill(variable,datasets[d].NormFactor()*Luminosity);
		}
	}//loop on events
  }//loop on datasets

  //Once everything is filled ...
  if(verbose>0) cout<<" We ran over all the data ;-)"<<endl;
 
  
 ////////////////////////////////////
 /// Selection Table
 ////////////////////////////////////
  if(verbose>0) cout<<" - Write SelectionTable"<<endl;
  selecTable.TableCalculator();
  selecTable.Write(selecTableFileName);
  selecTableCR.TableCalculator();
  selecTableCR.Write(selecTableCRFileName);

 ////////////////////////////////////
 /// ABCDEstimation
 ////////////////////////////////////
 if(verbose>0) cout<<" - Perform ABCD method ..."<<endl;
 abcdEstim.ComputeEstimate(cutXmin, cutX0, cutX1, cutXmax, cutYmin, cutY0, cutY1, cutYmax, region);
 if(verbose>0) abcdEstim.Print2DEstimate();
 //if(verbose>0) abcdEstim.Print3DEstimate();

 ////////////////////////////////////
 /// QCDShapeEstimation
 ////////////////////////////////////
 //call the function in that order
 if(verbose>0) cout<<" - Perform QCD shape estimation ..."<<endl;
 //qcdShapeEstim.Normalize(qcdShapeEstim.GetQCDMCShapeSR()->Integral());
 qcdShapeEstim.Normalize(abcdEstim.GetEstimate(region));
 qcdShapeEstim.AddNormalisationError(10.);//warning: add 10% from ABCD
 qcdShapeEstim.Draw();
 if(verbose>0) qcdShapeEstim.PrintResults();

  /*
  for(int i=0;i<NofJetBins;i++){
 	for(unsigned int j=0;j<4;j++){ 
		MultiJets_Estimated_Nj_SR[i][j] = EstimJetsBJets.GetQCDEstimate(i,j);
		cout<<"QCD estim "<< EstimJetsBJets.GetQCDEstimate(i,j)<<" "<<MultiJets_Estimated_Nj_SR[i][j]<<endl;
	}
   }
   */

 ////////////////////////////////////
 /// WJets Estimation = NofEvents
 ////////////////////////////////////
 
 if(verbose>0) cout<<" - Perform WJets estimation ..."<<endl;
 //SR
 //Important - Scale numbers given as input
 for(int i=0;i<NofJetBins;i++){
 	Scale(datasets, N0_SR[i], Luminosity);
 	Scale(datasets, N1_SR[i], Luminosity);
 	Scale(datasets, N2_SR[i], Luminosity);
 	Scale(datasets, N3_SR[i], Luminosity);
 }
  if(doWJetEstimation){
  	wjEstimation_SR.FillInputs(N0_SR, N1_SR, N2_SR, N3_SR, MultiJets_Estimated_Nj_SR);
  	if(verbose>0)  wjEstimation_SR.PrintInputs();
  	wjEstimation_SR.Estimation(false, true);
  	if(verbose>0)  wjEstimation_SR.PrintResults();

  	//compute the errors with Pseudo-exp
  	int NWjetsPseudoExp = 1000;
  	WJetEstPseudoExp wjPseudoExp(NWjetsPseudoExp, &(wjEstimation_SR));
  	wjPseudoExp.RollTheDice(false); 
	TFile* fout = new TFile("PseudoExp.root","RECREATE");
	wjPseudoExp.Write(fout,string(""));
	fout->Write();
	fout->Close();
  	//EstimJetsBJets.FillWJetsEstim(wjEstimation_SR, 0., 0.1, 0.1, 1., 0., 0.2, 0.2, 1., 1);
  }

 ////////////////////////////////////
 /// WJets Estimation - shape
 ////////////////////////////////////
  if(verbose>0) cout<<" - Perform WJets shape estimation ..."<<endl;
  float RatioNbEvts_Wlike_1bjet = 0;
  float RatioNbEvts_TTlike_0bjet = 0;
  float Nevts0bj = 0;
  float Nevts1bj = 0;
  //SR
  WJetShapeEstimation itshExtr_SR;
  if(doWJetEstimation){
  	RatioNbEvts_Wlike_1bjet = wjEstimation_SR.GetNWEventsRatioEstimatedNbjets(1);
  	RatioNbEvts_TTlike_0bjet = wjEstimation_SR.GetNttEventsRatioEstimatedNbjets(0);
  	Nevts0bj = wjEstimation_SR.NofEventsNbjets(0);
  	Nevts1bj = wjEstimation_SR.NofEventsNbjets(1);
  	itshExtr_SR.Initialisation(histo_inclusif_SR, histo_0bj_SR, histo_1bj_SR, histo_Wjets_SR, histo_ttjets_SR, RatioNbEvts_Wlike_1bjet, RatioNbEvts_TTlike_0bjet, Nevts0bj, Nevts1bj, Iterations);
  	itshExtr_SR.IterativeHistoSubstraction();
  	itshExtr_SR.Scale(wjEstimation_SR.GetNttEventsEstimated(), wjEstimation_SR.GetNWEventsEstimated());
  	if(verbose>0) itshExtr_SR.PrintResults();
  	itshExtr_SR.AddNormalisationError(12.);//Warning: add 12% from WJetEstimation
  	itshExtr_SR.Draw();
  }

 ////////////////////////////////////
 /// TtJet Estimation
 ////////////////////////////////////
  if(verbose>0) cout<<" - Perform TTJets Estimation ..."<<endl;
  // 1- Scale to Lumi
  ttjEstimation.Scale(Luminosity);
  // 2- Compute observables
  //make it easier ... !! and comment
  ttjEstimation.ComputeObservableSR(true, true, true, 0); //do NP
  ttjEstimation.ComputeObservableCR(true, true, false, 0); // WARNING ... w/o NP  
  //Expected Estimation w/o SUSY contamination: WARNING: had possibility to choose !!
  ttjEstimation.ComputeObservableExpected(true, true); 
  ttjEstimation.ComputeObservableObserved(true, true, true, 0); 
  // 3- Normalisation
  ttjEstimation.DefineNR(NRegion);
  ttjEstimation.ComputeNormalisation(qcdShapeEstim.GetQCDShapeEstimated(), itshExtr_SR.GetHistoWLikeEstimated());
  if(verbose>0){
  	ttjEstimation.PrintResults();
	float RelError = (ttjEstimation.NofEventsEstimated(METCut).second/ttjEstimation.NofEventsEstimated(METCut).first)*100.;	
  	cout<<"TTJetsResults nofEvtsEstimated: "<<ttjEstimation.NofEventsEstimated(METCut).first<<" +/- "<<ttjEstimation.NofEventsEstimated(METCut).second<<" RelError = "<<RelError<<" %"<<endl; 
  	cout<<"TTJetsResults nofEvtsExpected: "<<ttjEstimation.NofEventsExpected(METCut).first<<" +/- "<<ttjEstimation.NofEventsExpected(METCut).second<<endl;
  	cout<<"TTJetsResults nofEvtsObserved: "<<ttjEstimation.NofEventsObserved(METCut).first<<" +/- "<<ttjEstimation.NofEventsObserved(METCut).second<<endl;
  }
  //here
  cout<<"ComputePlotSignificance"<<endl;
  ttjEstimation.ComputePlotSignificance(qcdShapeEstim.GetQCDShapeEstimated(), itshExtr_SR.GetHistoWLikeEstimated());
  ttjEstimation.ComputeEstimationLMXContamination(qcdShapeEstim.GetQCDShapeEstimated(), itshExtr_SR.GetHistoWLikeEstimated());
  cout<<"end"<<endl;
  //
  ttjEstimation.Draw();




  //sensitivity estimated & expected depending on the cut on MET for different LMX
  /*
  TH1F** hSignEstim = new TH1F*[LMX.size()];
  TH1F** hSignExpec = new TH1F*[LMX.size()];
  for(unsigned int l=0;l<LMX.size();l++){
  	char name[100];
	sprintf(name,"hSignEstim_%d",l);
  	hSignEstim[l] = new TH1F(name,"",NbinsX, binsX);
	hSignEstim[l]->SetLineColor(LMX[l].Color());
	sprintf(name,"hSignExpec_%d",l);
  	hSignExpec[l] = new TH1F(name,"",NbinsX, binsX);
	hSignExpec[l]->SetLineColor(LMX[l].Color());
  	
  }
  for(int i=0;i<NbinsX;i++){
  	ttjEstimation.ComputeObservableExpected(true, true); 
  	ttjEstimation.DefineNR(pair<float,float>(0,binsX[i]));
  	for(unsigned int l=0;l<LMX.size();l++){
		//recompute
  		ttjEstimation.ComputeObservableSR(true, true, true, l); // add this LMX
  		ttjEstimation.ComputeObservableCR(true, true, true, l); // add this LMX
  		ttjEstimation.ComputeObservableObserved(true, true, true, l); 
		ttjEstimation.ComputeNormalisation(qcdShapeEstim.GetQCDShapeEstimated(), itshExtr_SR.GetHistoWLikeEstimated());
		//
		pair<float,float> sign;	
		//Estimated
		sign = ttjEstimation.SignificanceEstimated(binsX[i]);
		hSignEstim[l]->SetBinContent(i,sign.first);
		hSignEstim[l]->SetBinError(i,sign.second);
		//Expected
		sign = ttjEstimation.SignificanceExpected(binsX[i]);
		hSignExpec[l]->SetBinContent(i,sign.first);
		hSignExpec[l]->SetBinError(i,sign.second);
  	}
  }
  TCanvas cc;
  cc.cd();
  for(unsigned int l=0;l<LMX.size();l++){
  	if(l==0) hSignEstim[l]->Draw();
  	else hSignEstim[l]->Draw("same");
	hSignExpec[l]->Draw("same");
  }
  */

  
  //////////////////
  /// Others plots
  ////////////////////
  msplotMET.Draw(string("Met"));
  msplotDRBB.Draw(string("DRBB"));
  msplotDEtaBB.Draw(string("DEtaBB"));
  msplotDPhiBB.Draw(string("DPhiBB"));
  msplotHTBB.Draw(string("HTBB"));
  msplotMbl.Draw(string("Mbl"));
  msplotDRbl.Draw(string("DRbl"));
  msplotDPhiMETlepton.Draw(string("DPhiMETlepton"));

  mcplotMET.Draw(string("MET_Cuts"),4.);
  cout<<"Tail of MET for different cuts"<<endl;
  mcplotMET.TailEstimation(100.);
  mcplotMET2.Draw(string("MET_Cuts_OneByOne"),4.);

  mcplotMETttbarComponents.Draw(string("MET"),4.);
  
  cimpeMbl.Draw(string("CimpeMbl"));
  cimpeHTBB.Draw(string("CimpeHTBB"));
  cimpeDRBB.Draw(string("CimpeDRBB"));
  
  //Get info from CR cuts
  cout<<"------------------------------------------------------"<<endl;
  float EffGoal = 0.9;
  float Chi2Goal = 0.0001;
  cout<<" Results about CR cuts"<<endl;
  pair<float,float> GetInfoRes;
  GetInfoRes = cimpeMbl.GetInformation(MblCut);
  cout<<"Cut on Mbl: "<<MblCut<<" - X2 = "<<GetInfoRes.first<<" - Eff = "<<GetInfoRes.second<<endl;
  cout<<"For a goal of Eff = "<<EffGoal<<" Cut = "<<cimpeMbl.GetCutForAGivenEfficiency(EffGoal).first <<" Chi2 = "<<cimpeMbl.GetCutForAGivenEfficiency(EffGoal).second <<endl;
  cout<<"For a goal of Chi2 = "<<Chi2Goal<<" Cut = "<<cimpeMbl.GetCutForAGivenChi2(Chi2Goal).first<<"  eff = "<<cimpeMbl.GetCutForAGivenChi2(Chi2Goal).second<<endl;
  GetInfoRes = cimpeDRBB.GetInformation(DRBBCut);
  cout<<"Cut on DRBB: "<<DRBBCut<<" - X2 = "<<GetInfoRes.first<<" - Eff = "<<GetInfoRes.second<<endl;
  cout<<"For a goal of Eff = "<<EffGoal<<" Cut = "<<cimpeDRBB.GetCutForAGivenEfficiency(EffGoal).first<<" Chi2 = "<<cimpeDRBB.GetCutForAGivenEfficiency(EffGoal).second<<endl;
  cout<<"For a goal of Chi2 = "<<Chi2Goal<<" Cut = "<<cimpeDRBB.GetCutForAGivenChi2(Chi2Goal).first<<" eff = "<<cimpeDRBB.GetCutForAGivenChi2(Chi2Goal).second<<endl;
  GetInfoRes = cimpeHTBB.GetInformation(HTBBCut);
  cout<<"Cut on HTBB: "<<HTBBCut<<" - X2 = "<<GetInfoRes.first<<" - Eff = "<<GetInfoRes.second<<endl;
  cout<<"For a goal of Eff = "<<EffGoal<<" Cut = "<<cimpeHTBB.GetCutForAGivenEfficiency(EffGoal).first<<" eff = "<<cimpeHTBB.GetCutForAGivenEfficiency(EffGoal).second<<endl;
  cout<<"For a goal of Chi2 = "<<Chi2Goal<<" Cut = "<<cimpeHTBB.GetCutForAGivenChi2(Chi2Goal).first<<" Chi2 = "<<cimpeHTBB.GetCutForAGivenChi2(Chi2Goal).second<<endl;
  cout<<"------------------------------------------------------"<<endl;


  /////////////////////////////////
  //  Summary plot
  /////////////////////////////////
  //BkgEstimationSummary bkgSum(TString("BgkSum"), XaxisLabel, TString("# events"), hData, qcdShapeEstim.GetQCDShapeEstimated(), hQCDMC, itshExtr_SR.GetHistoWLikeEstimated(),  hWJetMC, ttjEstimation.PlotEstimated(),  hTtJetMC);
  cout<<"hData "<<hData->Integral()<<endl;
  cout<<"QCD MC "<<hQCDMC->Integral()<<endl;
  cout<<"Wj  MC "<<hWJetMC->Integral()<<endl;
  cout<<"ttj MC "<<hTtJetMC->Integral()<<endl;
  cout<<"QCD estim "<<qcdShapeEstim.GetQCDShapeEstimated()->Integral()<<endl;
  cout<<"Wj estim "<<itshExtr_SR.GetHistoWLikeEstimated()->Integral()<<endl;
  cout<<"ttj estim "<<ttjEstimation.PlotEstimated()->Integral()<<endl;
  BkgEstimationSummary bkgSum(TString("BgkSum"), XaxisLabel, TString("# events"), hData, qcdShapeEstim.GetQCDShapeEstimated(), hQCDMC, itshExtr_SR.GetHistoWLikeEstimated(),  hWJetMC, ttjEstimation.PlotEstimated(),  hTtJetMC);
  //BkgEstimationSummary bkgSum(TString("BgkSum"), XaxisLabel, TString("# events"), hData, 0, 0, itshExtr_SR.GetHistoWLikeEstimated(),  hWJetMC, ttjEstimation.PlotEstimated(),  hTtJetMC);
  bkgSum.Draw();


  ////////////////////////////////
  // Estimation per #jets & #b-jets
  ////////////////////////////////
  //EstimJetsBJets.Print(bkgEstimPerJetsBJetsFileName);

  /////////////////////////////////////////
  //loop over lumi 
  // estimate statistical errors
  /////////////////////////////////////////
  int nPointsLumi = 100;//10
  int LumiMax = 100;//100
  float PreviousTemp = Luminosity;

  //copy of bkg estimation
  WJetEstimation wjEstimCopy(wjEstimation_SR);
  ABCDEstimation abcdEstimCopy(abcdEstim);
  TtJetEstimation ttjEstimCopy(ttjEstimation);

  //TTtJets
  TH1F* h_RelErrorEstimation_TtJets = new TH1F("h_RelErrorEstimation_TtJets","",nPointsLumi,0, LumiMax);
  h_RelErrorEstimation_TtJets->GetXaxis()->SetTitle("Luminosity [pb^{-1}]");
  h_RelErrorEstimation_TtJets->GetYaxis()->SetTitle("Relative error [%]");
  //WJets
  TH1F* h_RelErrorEstimation_WJets = new TH1F("h_RelErrorEstimation_WJets","",nPointsLumi,0, LumiMax);
  h_RelErrorEstimation_WJets->GetXaxis()->SetTitle("Luminosity [pb^{-1}]");
  h_RelErrorEstimation_WJets->GetYaxis()->SetTitle("Relative error [%]");
  //ABCD
  TH1F* h_RelErrorEstimation_QCD = new TH1F("h_RelErrorEstimation_QCD","",nPointsLumi,0, LumiMax);
  h_RelErrorEstimation_QCD->GetXaxis()->SetTitle("Luminosity [pb^{-1}]");
  h_RelErrorEstimation_QCD->GetYaxis()->SetTitle("Relative error [%]");
  //
  TFile foutWJPS ("WJetsPseudoExp.root","RECREATE");
  for(int n=1;n<nPointsLumi+1;n++){
  	float Lumi = n*LumiMax/nPointsLumi;
  	float RescaleFactor = (float) Lumi/PreviousTemp;
	float RelError = 0;
	
	//TtJets
	ttjEstimCopy.ReScale(RescaleFactor);
  	ttjEstimCopy.ComputeObservableSR(true, true, false, 0); //not do NP
  	ttjEstimCopy.ComputeObservableCR(true, true, false, 0); //not do NP  
  	ttjEstimCopy.ComputeObservableExpected(true, true); 
  	ttjEstimCopy.ComputeObservableObserved(true, true, false, 0);//not do NP 
  	ttjEstimCopy.ComputeNormalisation(qcdShapeEstim.GetQCDShapeEstimated(), itshExtr_SR.GetHistoWLikeEstimated());
	RelError = (ttjEstimCopy.NofEventsEstimated(METCut).second/ttjEstimCopy.NofEventsEstimated(METCut).first)*100.;	
	h_RelErrorEstimation_TtJets->SetBinContent(n,RelError);
	
	//WJets
	/*
	wjEstimCopy.ReScale(RescaleFactor);
	int NWjetsPseudoExp = 1000;
	WJetEstPseudoExp wjPseudoExpCopy(NWjetsPseudoExp, &(wjEstimation_SR));
	wjPseudoExpCopy.RollTheDice(false,false);
	RelError = 0;
	wjEstimCopy.Estimation(Iterations, false);
	if(wjEstimCopy.GetNWEventsEstimated()) RelError = wjPseudoExpCopy.GetNWError()/wjEstimCopy.GetNWEventsEstimated();
	h_RelErrorEstimation_WJets->SetBinContent(n,RelError);
	char dirname[100];
	sprintf(dirname,"%d",n);
	foutWJPS.mkdir(dirname);
	foutWJPS.cd(dirname);
	wjPseudoExpCopy.Write(&foutWJPS,string(dirname));
	*/

	//QCD
	abcdEstimCopy.ReScale(RescaleFactor);
	abcdEstimCopy.ComputeEstimate(cutXmin, cutX0, cutX1, cutXmax, cutYmin, cutY0, cutY1, cutYmax, region);
        RelError = 0;
	cout<<"QCD errors: "<<abcdEstimCopy.GetEstimateError(region)<<" "<<abcdEstimCopy.GetEstimate(region)<<endl;
	if(abcdEstimCopy.GetEstimate(region)) RelError = abcdEstimCopy.GetEstimateError(region)/abcdEstimCopy.GetEstimate(region);
	h_RelErrorEstimation_QCD->SetBinContent(n,RelError);
	
	PreviousTemp = Lumi;
  }
  foutWJPS.Close();
  /////////////////////////////////////////

  ///////////////////
  // Writting
  //////////////////
  if(verbose>1) cout<<" - Writting outputs on files ..."<<endl;
  TFile* fout = new TFile(rootFileName.c_str(),"RECREATE");
  //SR
  if(doWJetEstimation){
  	wjEstimation_SR.Write(fout,string("SR"));
  	itshExtr_SR.Write(fout,string("SR"));
  }
  histo_Wjets_SR->Write();
  histo_ttjets_SR->Write();
  ttjEstimation.Write(fout);
  abcdEstim.Write(fout);
  h_ABCD_QCD->Write(); 
  h_ABCD_WJets->Write(); 
  h_ABCD_TTJets->Write(); 
  TCanvas cABCD("cABCD");
  cABCD.cd();
  h_ABCD_QCD->Draw();
  h_ABCD_WJets->SetLineColor(3);
  h_ABCD_TTJets->SetLineColor(4);
  h_ABCD_WJets->Draw("same");
  h_ABCD_TTJets->Draw("same");
  TLegend leg_ABCD(0.7,0.7,0.98,0.98);
  leg_ABCD.AddEntry(h_ABCD_QCD,"QCD","l");
  leg_ABCD.AddEntry(h_ABCD_WJets,"W+jets","l");
  leg_ABCD.AddEntry(h_ABCD_TTJets,"tt+jets","l");
  leg_ABCD.Draw("same");
  cABCD.Write();
  qcdShapeEstim.Write(fout);
  bkgSum.Write(fout);
  //others plots
  msplotMET.Write(fout);
  msplotMbl.Write(fout);
  msplotDRBB.Write(fout);
  msplotDEtaBB.Write(fout);
  msplotDPhiBB.Write(fout);
  msplotHTBB.Write(fout);
  msplotDRbl.Write(fout);
  msplotDPhiMETlepton.Write(fout);
  mcplotMET.Write(fout);
  mcplotMET2.Write(fout);
  mcplotMETttbarComponents.Write(fout);
  //
  cimpeMbl.Write(fout);
  cimpeHTBB.Write(fout);
  cimpeDRBB.Write(fout);
  //
  //Stat errors
  fout->cd();
  fout->mkdir("StatErrors");
  fout->cd("StatErrors");
  h_RelErrorEstimation_TtJets->Write();
  h_RelErrorEstimation_WJets->Write();
  h_RelErrorEstimation_QCD->Write();
  //
  //cc.Write();
  //

  if(verbose>1) cout<<" - Writting  the file ..."<<endl;
  fout->Write();

  cout<<"**********************************************************************"<<endl;
  cout<<"           End of the program !!" <<endl;
  cout<<" 		doesn't crashed yet ;-) "<<endl;
  cout<<"**********************************************************************"<<endl;

}
