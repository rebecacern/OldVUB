
#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

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


void initFromByHandValues(VJetEstimation *vje, UInt_t NofJetBins, UInt_t NofBtagWorkingPoint, Double_t **initForMinVJet);
void initFromMCValues(VJetEstimation *vje, UInt_t NofJetBins, UInt_t NofBtagWorkingPoint);
/**
 Results from BTV-11-003, https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG link to https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/eff_b_c-ttbar_payload.txt 
 In Range : 0<=|eta|<=2.4 , 30<=p_T<=200  (approximately)
 */
/** TCHE */
Float_t eff_b_ttbar_FtCM(Float_t /*b-disc*/ x) {
  return -3.67153247396e-07*x*x*x*x +  -2.81599797034e-05*x*x*x +  0.00293190163243*x*x +  -0.0849600849778*x +  0.928524440715 ;
};
Float_t eff_b_ttbar_FtCM_pluserr(Float_t /*b-disc*/ x) {
  return 3.03337430722e-06*x*x*x*x + -0.000171604835897*x*x*x + 0.00474711667943*x*x + -0.0929933040514*x + 0.978347619293 ;
};
Float_t eff_b_ttbar_MC(Float_t /*b-disc*/ x) {
  return 3.90732786802e-06*x*x*x*x +  -0.000239934437355*x*x*x +  0.00664986827287*x*x +  -0.112578996016*x +  1.00775721404 ;
};
Float_t eff_c_ttbar_MC(Float_t /*b-disc*/ x) {
  return 0.343760640168*exp(-0.00315525164823*x*x*x + 0.0805427315196*x*x + -0.867625139194*x + 1.44815935164 ) ;
};

int main(int argc, const char* argv[]) {
  
  clock_t start = clock();
  
  cout << "********************************" << endl;
  cout << " Backgrounds for stops searches " << endl;
  cout << "********************************" << endl;
  
    //SetStyle if needed
    //setTDRStyle(); 
  setMyStyle();
  
    /////////////////////
    // Configuration
    /////////////////////
  
    //xml file
    //  std::string xmlfile ="config/myRefSelconfig.xml";
  std::string baseDirectory = "." ;
  if (argc >= 2) {
    baseDirectory = std::string(argv[1]);
  }
  std::string outputDirectory = "." ;
  if (argc >= 3) {
    outputDirectory = std::string(argv[2]);
  }
  {
    char timestamp[51];
    struct tm *timeinfo;
    time_t currTime;
    time(&currTime);
    timeinfo = localtime(&currTime);
    strftime(timestamp, 50, "%Y%m%d_%H%M%S", timeinfo);
      //    delete timeinfo;
    outputDirectory = outputDirectory + "/" + timestamp;
    system( (std::string("mkdir -p ")+outputDirectory+" && cp "+baseDirectory+"/macros/StopBckg_fixedVarStudies.cc "+outputDirectory+"/ ;").c_str());
  }
  std::string xmlfile = baseDirectory + "/config/StopBckg.xml";
  
    //Input ROOT file
  std::string inputRootFileName = baseDirectory + "/StopBckg_Output.root";
    //Output ROOT file
  std::string outputRootFileName = outputDirectory + "/StopBckg_Output2.root";
  
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
  
  Bool_t runOnEvents_ = kFALSE;
  Bool_t reloadFromPreviousRun_ = kTRUE;
  
  TTreeLoader treeLoader;
  cout << " - Load datasets ..." << endl;
  std::vector < Dataset* > datasets;
  treeLoader.LoadDatasets (datasets, xmlfile.c_str());
  for(unsigned int i=0;i<datasets.size();i++) {
    new ((*tcdatasets)[i]) Dataset(*datasets[i]);
  }
  
  float oldLuminosity = anaEnv.Luminosity;	// in 1/pb
  float Luminosity = 5000.;
    //  Double_t eventsPrescaleFactor = 50. ;
  
  /**
   Uncertainty variables
   */
    //int nsigma = 1;
  
  TFile *fin = new TFile (inputRootFileName.c_str(), "READ");
  
  /*
   for (unsigned int d = 0; d < datasets.size (); d++)
   {
   //cout << "luminosity of dataset "<< d << " is " << datasets[d]->EquivalentLumi() << endl;
   //    if(Luminosity > datasets[d]->EquivalentLumi() ) Luminosity = datasets[d]->EquivalentLumi();
   std::string dataSetName = datasets[d]->Name();
   if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
   {
   Luminosity = datasets[d]->EquivalentLumi();
   }
   
   }
   */
  if(Luminosity != oldLuminosity) cout << "changed analysis environment luminosity to "<< Luminosity << endl;
  
  
  
    ////////////////////////////////////
    /// Normal Plots (TH1F* and TH2F*)
    ////////////////////////////////////
  
  map<std::string,TH1*> histo1D;
  map<std::string,TH2*> histo2D;
  
  
  
    ////////////////////////////////////
    /// V-jets estimation
    ////////////////////////////////////
  
  if (verbose > 0)
    cout << " - V-jets estimation instantiation ..." << endl;
  const UInt_t NofBtagWorkingPoint = 2;
  float* BtagWorkingPoint;
  BtagWorkingPoint = new float[NofBtagWorkingPoint]; BtagWorkingPoint[0]=3.2;
  const UInt_t NofJets = 3;       // Min. mult. of jets
  const UInt_t NofJetBins = 4;    // Nb of bins for mult. of jets
                                  //  double** EffXbq;   // fraction of b-tag jets, by multiplicity ([i][j] is b-mult "j" in "i+NofJets" jets) ?
  histo1D["JetMultiplicities_Vlike"]  = new TH1F("JetMultiplicities_Vlike", "JetMultiplicities_Vlike", 21, 0, 21) ;
  histo1D["JetMultiplicities_TTlike"] = new TH1F("JetMultiplicities_TTlike", "JetMultiplicities_TTlike" ,21, 0, 21) ;
  
  int NofDatasets = datasets.size();
  std::vector<int> iDTTLike; // Indices of TT-Like datasets
  std::vector<int> iDVLike;  // 
  std::vector<int> iDVbLike; // 
  for (Int_t i=0; i<NofDatasets; i++) {
    if ((datasets[i]->Name().find("TT_")==0) || (datasets[i]->Name().find("Stop_")==0)) {
      iDTTLike.push_back(i);
    } 
    if (datasets[i]->Name().find("W_Jets")==0) {
      iDVLike.push_back(i);
    }
  }
  printf("V-like datasets : ");
  for (std::vector<int>::iterator iter=iDVLike.begin(); iter!=iDVLike.end(); iter++) {
    if (iter != iDVLike.begin()) {
      printf(", ");
    }
    printf("%s", datasets[(*iter)]->Name().c_str());
  }
  printf("\n");
  printf("TT-like datasets : ");
  for (std::vector<int>::iterator iter=iDTTLike.begin(); iter!=iDTTLike.end(); iter++) {
    if (iter != iDTTLike.begin()) {
      printf(", ");
    }
    printf("%s", datasets[(*iter)]->Name().c_str());
  }
  printf("\n");
  printf("Vb-like datasets : ");
  for (std::vector<int>::iterator iter=iDVbLike.begin(); iter!=iDVbLike.end(); iter++) {
    if (iter != iDVbLike.begin()) {
      printf(", ");
    }
    printf("%s", datasets[(*iter)]->Name().c_str());
  }
  printf("\n");
  
  VJetEstimation vJetEst;//(NofBtagWorkingPoint, BtagWorkingPoint, NofJets, NofJetBins, NofDatasets, iDTTLike, iDVLike, iDVbLike);
                         //  VJetEstimation vJetEstSTasV(NofBtagWorkingPoint, BtagWorkingPoint, NofJets, NofJetBins, NofDatasets, iDTTLike, iDVLike, iDVbLike);
                         //  VJetEstimation vJetEstSTasTT(NofBtagWorkingPoint, BtagWorkingPoint, NofJets, NofJetBins, NofDatasets, iDTTLike, iDVLike, iDVbLike);
                         //  VJetEstimation vJetEst(NofBtagWorkingPoint, BtagWorkingPoint, NofJets, NofJetBins, EffXbq, NofDatasets, iDTTLike, iDVLike, iDVbLike);
  VJetEstimation* vje = NULL;
  if (reloadFromPreviousRun_) {
    char name[200];
    fin->ls();
    fin->Map();
    sprintf(name, "VJetEstimation-TCHE_LM--tt_W"/*, NofJets+NofJetBins-1*/);
    printf("%s\n", name);
    printf("Loading\n");
    if (fin->Get(name)==NULL)
      printf("TObject not found !!!\n");
      //vJetEst = (*((VJetEstimation*) fin->Get(name)));
    vje = (VJetEstimation*) (fin->Get(name));
    printf("Loaded\n");
  }
  
    //   fin->Close();
  
    //unsigned int btagAlgo = 0;
  /* btagAlgo
   0) btag_trackCountingHighEffBJetTags()
   1) btag_trackCountingHighPurBJetTags()
   2) btag_jetProbabilityBJetTags()
   3) btag_jetBProbabilityBJetTags()
   4) btag_simpleSecondaryVertexHighEffBJetTags()
   5) btag_simpleSecondaryVertexHighPurBJetTags()
   6) btag_combinedSecondaryVertexBJetTags()
   */
  Bool_t setFromPAS_11_001 = kFALSE;
  Bool_t isOnMC = kTRUE;
  
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
  
    //  vje->Dump();
    // fin->Close();
  
  
  
  std::ofstream stopReport; 
  stopReport.open((outputDirectory + "/StopBckg.tex").c_str()) ;
  stopReport << "\\documentclass{article}" << std::endl;
  stopReport << "\\usepackage{supertabular}" << std::endl;
  stopReport << "\\usepackage{pdflscape}" << std::endl;
  stopReport << "\\begin{document}" << std::endl;
  
  std::ofstream errCompReport; 
  errCompReport.open((outputDirectory + "/StopErrComp.tex").c_str()) ;
  
  
    // V-jets estimation
  printf("Extracting MC values for the efficiency\n");
  vje->ComputeEffFromMC();
  vje->ComputeEffbqFromMC(); // From ttbar sample
  
    //BTV-11-001
    // System8 method // p_T,rel
    //Muon-jet p_T 50-80 GeV, SF(= e_data/e_mc) : 20-240 GeV
  std::vector<Double_t> eff_b(NofBtagWorkingPoint,0.), eff_b_err(NofBtagWorkingPoint,0.), eff_b_SF(NofBtagWorkingPoint,0.), eff_b_SF_err(NofBtagWorkingPoint,0.),
  eff_uds(NofBtagWorkingPoint,0.), eff_uds_err(NofBtagWorkingPoint,0.), eff_uds_SF(NofBtagWorkingPoint,0.), eff_uds_SF_err(NofBtagWorkingPoint,0.),
  eff_udsc(NofBtagWorkingPoint,0.), eff_udsc_err(NofBtagWorkingPoint,0.), eff_udsc_SF(NofBtagWorkingPoint,0.), eff_udsc_SF_err(NofBtagWorkingPoint,0.);
    //TCHE-loose
  eff_b[0] = 0.77;//0.76
  eff_b_err[0] = 0.01;
  eff_b_SF[0] = 0.95;
  eff_b_SF_err[0] = 0.01+0.10;
  eff_uds[0] = 0.128;
  eff_uds_err[0] = 0.001+0.026;
  eff_udsc_SF[0] = eff_uds_SF[0] = 1.11 ;
  eff_udsc_SF_err[0] = eff_uds_SF_err[0] = 0.01+0.12;
  if (setFromPAS_11_001==kFALSE) {
    if (isOnMC) {
      eff_b_SF[0] = 1.;
      eff_uds_SF[0] = 1.;
      eff_udsc_SF[0] = 1.;
    }
    eff_b[0] = vje->GetPredEb(0);
    eff_b_err[0] = vje->GetPredEbErr(0);
    eff_uds[0] = vje->GetPredEuds(0);
    eff_uds_err[0] = vje->GetPredEudsErr(0);
    eff_udsc[0] = vje->GetPredEudsc(0);
    eff_udsc_err[0] = vje->GetPredEudscErr(0);
  }  
  
    //TCHE-medium
  eff_b[1] = 0.63;//0.63
  eff_b_err[1] = 0.02;
  eff_b_SF[1] = 0.94;
  eff_b_SF_err[1] = 0.01+0.09;
  eff_uds[1] = 0.0175;
  eff_uds_err[1] = 0.0003+0.0038 ;
  eff_udsc_SF[1] = eff_uds_SF[1] = 1.21;
  eff_udsc_SF_err[1] = eff_uds_SF_err[1] = 0.02+0.17 ; 
  if (setFromPAS_11_001==kFALSE) {
    if (isOnMC) {
      eff_b_SF[1] = 1.;
      eff_uds_SF[1] = 1.;
      eff_udsc_SF[1] = 1.;
    }
    eff_b[1] = vje->GetPredEb(1);
    eff_b_err[1] = vje->GetPredEbErr(1);
    eff_uds[1] = vje->GetPredEuds(1);
    eff_uds_err[1] = vje->GetPredEudsErr(1);
    eff_udsc[1] = vje->GetPredEudsc(1);
    eff_udsc_err[1] = vje->GetPredEudscErr(1);
  }  
  
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
  
  /*
   for (UInt_t d =0; d<vje->GetNbOfDatasets(); d++) {
   vje->ReScaleInputs(d, Ntot, 5000./Luminosity);
   }
   */
  vje->SumOverAllInputs();
  vje->PrintInputs();
  
  /*
   //OLD Method  
   vJetEst.SetInitialValues((Double_t**) initForMinVJet);
   */
  
  
  
  VJetEstimation *vJetEst__tt__W = (VJetEstimation*) fin->Get("VJetEstimation-TCHE_LM--tt_W");
  VJetEstimation *vJetEst__tt__W__Z__SingleT__QCD = (VJetEstimation*) fin->Get("VJetEstimation-TCHE_LM--tt_W_Z_SingleT_QCD");
  VJetEstimation *vJetEst__tt__W_Z = (VJetEstimation*) fin->Get("VJetEstimation-TCHE_LM--tt_W-Z");
  VJetEstimation *vJetEst__tt__W_Z__QCD = (VJetEstimation*) fin->Get("VJetEstimation-TCHE_LM--tt_W-Z_QCD");
  VJetEstimation *vJetEst__tt__W_Z__SingleT__QCD = (VJetEstimation*) fin->Get("VJetEstimation-TCHE_LM--tt_W-Z_SingleT_QCD");
  VJetEstimation *vJetEst__tt_SingleT__W_Z__QCD = (VJetEstimation*) fin->Get("VJetEstimation-TCHE_LM--tt-SingleT_W-Z_QCD");
  VJetEstimation *vJetEst__tt_SingleT__W_Z = (VJetEstimation*) fin->Get("VJetEstimation-TCHE_LM--tt-SingleT_W-Z");
  
    // V-jets estimation
  std::vector<VJetEstimation*> vJetEst_MC_methods;
  if (vJetEst__tt__W != NULL) {
    vJetEst_MC_methods.push_back(vJetEst__tt__W);
  } else { printf("Missing object\n"); }
  /**/
  if (vJetEst__tt__W__Z__SingleT__QCD != NULL) {
    vJetEst_MC_methods.push_back(vJetEst__tt__W__Z__SingleT__QCD);
  } else { printf("Missing object\n"); }
  if (vJetEst__tt__W_Z != NULL) {
    vJetEst_MC_methods.push_back(vJetEst__tt__W_Z);
  } else { printf("Missing object\n"); }
  if (vJetEst__tt__W_Z__QCD != NULL) {
    vJetEst_MC_methods.push_back(vJetEst__tt__W_Z__QCD);
  } else { printf("Missing object\n"); }
  if (vJetEst__tt__W_Z__SingleT__QCD != NULL) {
    vJetEst_MC_methods.push_back(vJetEst__tt__W_Z__SingleT__QCD);
  } else { printf("Missing object\n"); }
  /**/
  if (vJetEst__tt_SingleT__W_Z__QCD != NULL) {
    vJetEst_MC_methods.push_back(vJetEst__tt_SingleT__W_Z__QCD);
  } else { printf("Missing object\n"); }
  if (vJetEst__tt_SingleT__W_Z != NULL) {
    vJetEst_MC_methods.push_back(vJetEst__tt_SingleT__W_Z);
  } else { printf("Missing object\n"); }
  
  
  
    //fixedVar.push_back(0);
  Bool_t doMinos = false;
  TFile *fout = new TFile (outputRootFileName.c_str(), "RECREATE");
  printf("Before loop\n");
  
  Int_t index_method=0;
  for (std::vector<VJetEstimation*>::iterator iterVje=vJetEst_MC_methods.begin(); iterVje!=vJetEst_MC_methods.end(); iterVje++) {
    std::string objectName = "";
    if (*iterVje == vJetEst__tt__W) {
      objectName="vJetEst__tt__W";
    } else if (*iterVje == vJetEst__tt__W__Z__SingleT__QCD) {
      objectName="vJetEst__tt__W__Z__SingleT__QCD";
    } else if (*iterVje == vJetEst__tt__W_Z) {
      objectName="vJetEst__tt__W_Z";
    } else if (*iterVje == vJetEst__tt__W_Z__QCD) {
      objectName="vJetEst__tt__W_Z__QCD";
    } else if (*iterVje == vJetEst__tt__W_Z__SingleT__QCD) {
      objectName="vJetEst__tt__W_Z__SingleT__QCD";
    } else if (*iterVje == vJetEst__tt_SingleT__W_Z__QCD) {
      objectName="vJetEst__tt_SingleT__W_Z__QCD";
    } else if (*iterVje == vJetEst__tt_SingleT__W_Z) {
      objectName="vJetEst__tt_SingleT__W_Z";
    } else {
      printf("Houston, we have a problem !\n");
      exit(1);
    }    
    vje = *iterVje; 
    if (vje==NULL) {
      continue;
    }
    printf("\n\nChanging to VJetEstimation object %s\n\n", objectName.c_str());
    stopReport << std::endl << "\\section{\\verb|" << objectName.c_str() << "| object} " << std::endl << std::endl;
    /*    vJetEst.UnBinnedMaximumNjetsLikelihoodEst(index_wp, fixedVar, doMinos, verbose);
     vJetEst.FillSummaryHistos();
     vJetEst.PrintResults();
     vJetEst.PrintResults_LatexFormat(stopReport);
     */
    /*  vJetEst.Write(fout, "_allJetsInclusive_", (verbose>0));
     fout->cd();*/
    
      //  for (int njet=NofJets; njet<NofJets+NofJetBins; njet++) {
    std::vector<std::string> fixedVarName(3, std::string());
    fixedVarName[0] = std::string("eb");
    fixedVarName[1] = std::string("eudsc");
    fixedVarName[2] = std::string("euds");
    for (UInt_t i=0; i<pow(2.,3); i++) { // Creating the combinations of fixed variables (3, at most)
      char fixedVarCombination[100] = "__fixed" ;
      std::vector<Int_t> fixedVar;    /*  0 : b-tag eff    ;    1 : mistag eff (udsc tagged as b)    ;    2 : mistag eff light (uds tagged as b) */
      UInt_t j=i;
      Int_t n=0;
      while (j!=0) {
        if ((j&1)==1) {
          sprintf(fixedVarCombination, "%s_%s", fixedVarCombination, fixedVarName[n].c_str()) ;
          fixedVar.push_back(n);
        }
        n++;
        j = j >> 1;
      }
      if (!(fixedVar.size()==3
            || (fixedVar.size()==2 && fixedVar[0]==1 && fixedVar[1]==2)
            || (fixedVar.size()==1 && fixedVar[0]==2)
            || (fixedVar.size()==1 && fixedVar[0]==1)
            || (fixedVar.size()==0)))
        continue;
      printf("\n\n HERE !!!!!!! INTERESTING\n\n\n");
      printf("");
      for (ULong_t configErr=0; configErr<pow(3.,(Int_t) fixedVar.size()); configErr++) {
        std::vector<Double_t> k_err(3,0.); // 3 variables can be fixed
                                           //for (ULong_t fVar = 0; fVar < fixedVar.size() || (fVar==0 && fixedVar.size()==0); fVar++) { // Listing the fixed variable to change their value within their errors
                                           //  for (ULong_t systTrial = 0; systTrial<pow(2, fixedVar.size()); systTrial++) {
        k_err = std::vector<Double_t>(3,0.); // 3 variables can be fixed (at most)
        ULong_t toGetLastDigit = configErr;
        for (UInt_t errScan=0; errScan<fixedVar.size(); errScan++) { // 3 points for each variation of the value : (-1;0;1)*error
          
          k_err[fixedVar[errScan]]= ((Double_t) (toGetLastDigit % 3))-1. ;
          toGetLastDigit = toGetLastDigit / 3 ; // division of integers.
        }
        for (Int_t ii=0; ii<3; ii++) {
          printf("k_err[%d]=% 2.0lf ", ii, k_err[ii]);
        }
        printf("\tfixedVar = {");
        for (std::vector<Int_t>::const_iterator iterv=fixedVar.begin(); iterv!=fixedVar.end(); iterv++) {
          printf(" %d", *iterv);
        }
        printf(" }");
        printf("\n");
        for (UInt_t index_wp=0; index_wp< NofBtagWorkingPoint ; index_wp++) {
          for (UInt_t njet=0; njet<NofJetBins; njet++) {
            fprintf(stderr, "\rProcessing : % 6.2f %%            ", 100.*(((Float_t)index_method)+((i+((configErr+(index_wp+((Float_t)(njet))/NofJetBins)/NofBtagWorkingPoint)/(pow(3.,(Int_t) fixedVar.size()))))/pow(2.,3))/(vJetEst_MC_methods.size())));
            cout << std::endl << std::endl;
            stopReport << std::endl << std::endl;
            errCompReport << std::endl << std::endl;
            bool saveWS = (k_err[0]==0 && k_err[1]==0 && k_err[2]==0); //Only saves the nominal (or floating)
                                                                       //(Re)set the correct initial values (from "by hand" or "MC")
                                                                       //initFromByHandValues(vje, NofJetBins, NofBtagWorkingPoint, initForMinVJet);
            initFromMCValues(vje, NofJetBins, NofBtagWorkingPoint);
              //vje->PrintInputs();
            for (std::vector<Int_t>::iterator iterf=fixedVar.begin(); iterf!=fixedVar.end(); iterf++) {
              for (UInt_t njets = 0; njets<NofJetBins; njets++) {
                switch (*iterf) { // Set fixed nominal values to ?????? 
                  case 0:
                    vje->SetInitialValues_Eb(njets, index_wp, eff_b[index_wp]*eff_b_SF[index_wp]+k_err[0]*(eff_b_err[index_wp]+eff_b_SF_err[index_wp]));
                    break;
                  case 1:
                    vje->SetInitialValues_Eudsc(njets, index_wp, eff_udsc[index_wp]*eff_udsc_SF[index_wp]+k_err[1]*(eff_udsc_err[index_wp]+eff_udsc_SF_err[index_wp]));
                    break;
                  case 2:
                    vje->SetInitialValues_Euds(njets, index_wp, eff_uds[index_wp]*eff_uds_SF[index_wp]+k_err[2]*(eff_uds_err[index_wp]+eff_uds_SF_err[index_wp]));
                    break;
                }
              }
            }
            vje->UnBinnedMaximumLikelihoodEst(index_wp, njet, fixedVar, doMinos, saveWS, verbose);
            if (saveWS) {
              vje->FillSummaryHistos();
            }
            errCompReport << "Indices of fixed variables :" ;
            for (std::vector<Int_t>::iterator iterf=fixedVar.begin(); iterf!=fixedVar.end(); iterf++) {
              errCompReport << " " << (*iterf) ;
            }
            errCompReport << endl;
            errCompReport << "Values for error propagation (in number of errors) :" ;
            for (std::vector<Double_t>::iterator iterf=k_err.begin(); iterf!=k_err.end(); iterf++) {
              errCompReport << " " << (*iterf) ;
            }
            errCompReport << endl;
            cout << "Indices of fixed variables :" ;
            for (std::vector<Int_t>::iterator iterf=fixedVar.begin(); iterf!=fixedVar.end(); iterf++) {
              cout << " " << (*iterf) ;
            }
            cout << endl;
            cout << "Values for error propagation (in number of errors) :" ;
            for (std::vector<Double_t>::iterator iterf=k_err.begin(); iterf!=k_err.end(); iterf++) {
              cout << " " << (*iterf) ;
            }
            cout << endl;
            
            errCompReport << "Predicted Ntt : " << vje->GetEstNtt(index_wp) << " +/- " << vje->GetEstNttErr(index_wp) << endl;
            errCompReport << "Predicted Nv : " << vje->GetEstNv(index_wp) << " +/- " << vje->GetEstNvErr(index_wp) << endl;
          }
        }
        vje->PrintResults();
        
        stopReport << "\\begin{landscape}" << std::endl;     
        vje->PrintResults_LatexFormat(stopReport);
        stopReport << "\\end{landscape}" << std::endl;     
        stopReport << "\\pagebreak" << std::endl << std::endl;
        
        /*
         char saveName[100]; sprintf(saveName, "VJetEstimation_%djets", njet+NofJets);
         //      vJetEst.Write(fout, directorySuffix, (verbose>0));
         fout->cd();
         vJetEst.Write(saveName);
         */
      }
      std::string dirName = objectName + fixedVarCombination;
      system( (std::string("mkdir -p ")+outputDirectory+"/"+dirName+"--Reloaded && mv "+outputDirectory+"/hNbjetsSummary_*.pdf VJetEstimation_RooFit_WS_*.root "+outputDirectory+"/"+dirName+"--Reloaded/ ;").c_str());
      
    }
  }
  
  printf("After loop\n");
  stopReport << "\\end{document}" << std::endl;
  stopReport.close();
  
  errCompReport.close();
  
  printf("fout dir creationg\n");
    //  vJetEst.Write(fout, "", (verbose>0));
  fout->cd();
  
  TDirectory* th1dir = fout->mkdir("1D_histograms");
  
  fout->cd();
  th1dir->cd();
  
  printf("In the correct directory, preparing to write canvases for histo1D\n");
  if (runOnEvents_) {
    for(std::map<std::string,TH1*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++) {
      TH1 *temp = it->second;
        //int N = temp->GetNbinsX();
        //temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
        //temp->SetBinContent(N+1,0);
      temp->Write();
      TCanvas* tempCanvas = TCanvasCreator((TH1F*)temp, it->first);
      tempCanvas->SaveAs( (outputDirectory + "/"+it->first+".root").c_str() );
      
    }
  }
  char name[100];
  
  
  
  vje->FillSummaryHistos();
  printf("After filling summary histos\n");
  /*
   vector<vector<TCanvas*> > v = vje->GETtCanva_Nbjets_Summary();
   for (std::vector<std::vector<TCanvas*> >::iterator iter_1=v.begin(); iter_1!=v.end(); iter_1++) {
   for (std::vector<TCanvas*>::iterator iter_2=iter_1->begin(); iter_2!=iter_1->end(); iter_2++) {
   if(true)cout<<"Writing summary histograms for X jets, wp nr Y"<<endl;
   sprintf(name, "macros/TestReload__%s.pdf", (*iter_2)->GetName());
   (*iter_2)->Print(name);
   }
   }
   */
  
  
  
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
  }
  
  
  
  
  
  printf("Before deleting\n");
  delete vje;
  printf("Before closing input file\n");
  fin->Close();
  printf("Before closing output file\n");
  fout->Close();
  
  printf("It took us %7lgs to run the program\n", ((Double_t)clock() - start) / CLOCKS_PER_SEC);
  
}





void initFromByHandValues(VJetEstimation *vje, UInt_t NofJetBins, UInt_t NofBtagWorkingPoint, Double_t **initForMinVJet){
    // Init from the vector of values ("by hand")
  std::vector<Double_t> init;
  std::vector< std::vector<Double_t> > init2;
    //  Double_t initForMinVJetArray[NofJetBins][2+3*NofBtagWorkingPoint]
  init.clear();  for (UInt_t i=0; i<NofJetBins; i++) { init.push_back(initForMinVJet[i][0]); }
  vje->SetInitialValues_Nttlike(init);
  init.clear();  for (UInt_t i=0; i<NofJetBins; i++) { init.push_back(initForMinVJet[i][1]); }
  vje->SetInitialValues_Nvlike(init);
  init2.clear();
  for (UInt_t i=0; i<NofJetBins; i++) { 
    std::vector<Double_t> tmpInit ;
    for (UInt_t j=0; j<NofBtagWorkingPoint; j++) { tmpInit.push_back(initForMinVJet[i][2+3*j]); }
    init2.push_back(tmpInit);
  }
  vje->SetInitialValues_Eb(init2); //[njets][btagIdx]
  init2.clear();
  for (UInt_t i=0; i<NofJetBins; i++) { 
    std::vector<Double_t> tmpInit ;
    for (UInt_t j=0; j<NofBtagWorkingPoint; j++) { tmpInit.push_back(initForMinVJet[i][2+1+3*j]); }
    init2.push_back(tmpInit);
  }
  vje->SetInitialValues_Eudsc(init2); //[njets][btagIdx]
  init2.clear();
  for (UInt_t i=0; i<NofJetBins; i++) { 
    std::vector<Double_t> tmpInit ;
    for (UInt_t j=0; j<NofBtagWorkingPoint; j++) { tmpInit.push_back(initForMinVJet[i][2+2+3*j]); }
    init2.push_back(tmpInit);
  }
  vje->SetInitialValues_Euds(init2); //[njets][btagIdx]
};


void initFromMCValues(VJetEstimation *vje, UInt_t NofJetBins, UInt_t NofBtagWorkingPoint) {
    // Init from MC ("automated")
  std::vector<Double_t> init;
  std::vector< std::vector<Double_t> > init2;
    //  Double_t initForMinVJetArray[NofJetBins][2+3*NofBtagWorkingPoint]
    // For Ntt and Nv initial values, the numbers should be the same across the working points (see the Fill method)
  init.clear();  for (UInt_t i=0; i<NofJetBins; i++) { init.push_back(vje->GetPredNtt(0,i)); }
  printf("Ntt init : ");
  for (std::vector<Double_t>::iterator iter=init.begin(); iter!=init.end(); iter++) {
    printf(" %lf", *iter);
  }
  printf("\n");
  vje->SetInitialValues_Nttlike(init);
  init.clear();  for (UInt_t i=0; i<NofJetBins; i++) { init.push_back(vje->GetPredNv(0,i)); }
  printf("Nv init : ");
  for (std::vector<Double_t>::iterator iter=init.begin(); iter!=init.end(); iter++) {
    printf(" %lf", *iter);
  }
  printf("\n");
  vje->SetInitialValues_Nvlike(init);
  init2.clear();
  for (UInt_t i=0; i<NofJetBins; i++) { 
    std::vector<Double_t> tmpInit ;
    for (UInt_t j=0; j<NofBtagWorkingPoint; j++) { tmpInit.push_back(vje->GetPredEb(j,i)); }
    init2.push_back(tmpInit);
  }
  printf("e_b init : \n");
  for (std::vector< std::vector<Double_t> >::iterator iter1=init2.begin(); iter1!=init2.end(); iter1++) {
    for (std::vector<Double_t>::iterator iter=iter1->begin(); iter!=iter1->end(); iter++) {
      printf(" %lf", *iter);
    }
    printf("\n");
  }
  printf("\n");
  vje->SetInitialValues_Eb(init2); //[njets][btagIdx]
  init2.clear();
  for (UInt_t i=0; i<NofJetBins; i++) { 
    std::vector<Double_t> tmpInit ;
    for (UInt_t j=0; j<NofBtagWorkingPoint; j++) { tmpInit.push_back(vje->GetPredEudsc(j,i)); }
    init2.push_back(tmpInit);
  }
  printf("e_udsc init : \n");
  for (std::vector< std::vector<Double_t> >::iterator iter1=init2.begin(); iter1!=init2.end(); iter1++) {
    for (std::vector<Double_t>::iterator iter=iter1->begin(); iter!=iter1->end(); iter++) {
      printf(" %lf", *iter);
    }
    printf("\n");
  }
  printf("\n");
  vje->SetInitialValues_Eudsc(init2); //[njets][btagIdx]
  init2.clear();
  for (UInt_t i=0; i<NofJetBins; i++) { 
    std::vector<Double_t> tmpInit ;
    for (UInt_t j=0; j<NofBtagWorkingPoint; j++) { tmpInit.push_back(vje->GetPredEuds(j,i)); }
    init2.push_back(tmpInit);
  }
  printf("e_uds init : \n");
  for (std::vector< std::vector<Double_t> >::iterator iter1=init2.begin(); iter1!=init2.end(); iter1++) {
    for (std::vector<Double_t>::iterator iter=iter1->begin(); iter!=iter1->end(); iter++) {
      printf(" %lf", *iter);
    }
    printf("\n");
  }
  printf("\n");
  vje->SetInitialValues_Euds(init2); //[njets][btagIdx]
};
