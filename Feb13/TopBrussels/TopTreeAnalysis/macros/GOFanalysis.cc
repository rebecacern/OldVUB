#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include "TF1.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TBranch.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include <cmath>
#include <fstream>
 
//user code
#include "TopTreeProducer/interface/TRootGenEvent.h"
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "../Selection/interface/SelectionTable.h"
#include "../BkgEstimationMethods/interface/TtJetEstimation.h"
#include "../BkgEstimationMethods/interface/TtJetEstPseudoExp.h"
#include "../BkgEstimationMethods/interface/BkgEstimationSummary.h"
#include "../Tools/interface/PlottingTools.h"
#include "../Tools/interface/MultiSamplePlot.h"
#include "../Tools/interface/MultiCutPlot.h"
#include "../Tools/interface/CutImpactEvaluation.h"
#include "../Tools/interface/TemplateComparator.h"
#include "../Tools/interface/TTreeLoader.h"
#include "../Tools/interface/AnalysisEnvironmentLoader.h"
#include "../Reconstruction/interface/MEzCalculator.h"
#include "../Reconstruction/interface/Observables.h"
#include "../Reconstruction/interface/PlotObservables.h"
#include "../Reconstruction/interface/TTreeObservables.h"
#include "../Reconstruction/interface/ObservablesRanker.h"
#include "../Reconstruction/interface/MakeBinning.h"
#include "../Content/interface/AnalysisEnvironment.h"
#include "../Content/interface/Dataset.h"
#include "../Content/interface/MCExpectation.h"
#include "../Content/interface/MCObsExpectation.h"
#include "../Content/interface/Container.h"
#include "../StatProcedure/interface/EventCombinedWeightCalculator.h"
#include "../StatProcedure/interface/SampleCombinedWeightCalculator.h"
#include "../StatProcedure/interface/MCPseudoExp.h"
#include "../StatProcedure/interface/WeightProbaCalculator.h"
#include "Style.C"

using namespace std;
using namespace TopTree;

void
Scale (vector < Dataset > datasets, double *&numbers, float Luminosity)
{
  for (unsigned int i = 1; i < datasets.size (); i++) {
    numbers[i] = numbers[i] * datasets[i].NormFactor () * Luminosity;
  }
}


int main (int argc, char *argv[])
{
  Bool_t DEBUG = (argc > 2 ? (atoi (argv[2]) == 1) : 0);
  cout << "**********************************************************************" << endl;
  cout << "Begining of the program for the Goodness-Of-Fit test of the Standard Model / Hypothesis test for t' search !" << endl;
  cout << "**********************************************************************" << endl;

  cout << "This macro starts from TFiles with TTrees with selected event observables and performs the statistical procedure." << endl;
  cout << "Note that a binning is needed in any case, this can also be made with this macro (different config file parameters)." << endl;
  cout << "The Monte Carlo expectation can be made first with this macro, starting from the TTrees; but can also be made with FullChainEstimation." << endl;
  

  //SetStyle if needed
  //setTDRStyle(); 
  //setMyStyle();
 
bool Tprime = true, scanTprimespace = false; //forceSMonly = true;
if(Tprime) cout<<"----- TPRIME ANALYSIS -----"<<endl;
else{
   scanTprimespace = false;
   cout<<"----- GOF-TEST ANALYSIS -----"<<endl;
}

float xsectrange[] = {2.94,3.0}; //{1,15} //{0.5,4} // {0.1,6} // {0.1,15} // {0.1,10}
float xsectstep = 0.2; //1
unsigned int numberofmasses = 1;
//float massestprimeHardCoded[6] = {350,400,450,500,550,600};//{300,350,400,450}; //only to be used when looping over plane but the t' samples are not added in the config
///float massestprimeHardCoded[4] = {350,400,450,500};
float massestprimeHardCoded[4] = {350,400,450,500};
string SignificanceCombination("SUM_LRatioMinded"); // to combine squared bin significances to get a combined weight for the event
//string SignificanceCombination("SUM"); // to combine squared bin significances to get a combined weight for the event
bool doPoissonPseudos = false;
bool doSMPoissonPseudos = false;
if(!doPoissonPseudos) doSMPoissonPseudos = false;
int nPoissonPseudos = 1000; //10000 per point
vector<float> massestprime;
vector<float> crosssectionstprime; //in pb
vector<float> NofEventsTprimeTopTree; //= {22166,20701,17817,20601};
bool do2Dhisto = false; //should stay false...
bool likelihood = true; //if true: inverses how the p-values are calculated... (test statistic S+B less than for B on average). For V-values this should be false.
bool doWithSpecial2Dbinning = false; //false when making the first HT binning!!!, true when makein Mtop binning //true when making data Mtop histos for the HTbins

bool doJESsystematics = false;
int doJESShift = 0; // 0: off 1: minus 2: plus. Decides only the files to write!! Reading is done for minus and plus at the same time for the moment, when you put doJESsystematics = true...
bool doBkgNormsystematics = false;
float TTbarError = 0.10, SingleTopError = 0.30;//, WJetsError = 0.30, ZJetsError = 0.30; ???? //relative
bool doLumisystematics = false;
float LumiError = 0;
if(doLumisystematics) LumiError = 0.06; //relative
bool doBtagsystematics = false;
bool doDataDummy = false; //this will abort the code after a certain point, to prevent a crash which causes a file to be not closed properly... Note this has nothing to do when running the code to obtain the observed test statistic!! Only for plots... scanTprimespace = true zetten denk ik, geen echte scan, maar toch
if(doDataDummy) SignificanceCombination = "SUM"; //very bad coding

//xml file
//const char *xmlfile = "../config/myTprimeconfig.xml"; //for real GOF-test: "../config/myNPconfig.xml"
const char *xmlfile = "../config/myTprimeconfig.xml";
//const char *xmlfile = "../config/myNPconfig.xml";
string thisROUNDdir = "./ROUND19/"; //directory for (most of??) the files of this round... (binning is still in config...) //don't forget the slash at the end
  
//Configuration output format
TTree *configTree_dummy = new TTree ("configTree_dummy", "Dummy configuration Tree");
TClonesArray *tcdatasets_dummy = new TClonesArray ("Dataset", 1000);
configTree_dummy->Branch ("Datasets_dummy", "TClonesArray", &tcdatasets_dummy);
TClonesArray *tcAnaEnv_dummy = new TClonesArray ("AnalysisEnvironment", 1000);
configTree_dummy->Branch ("AnaEnv", "TClonesArray", &tcAnaEnv_dummy);
TClonesArray *tcMCObsExp_dummy = new TClonesArray ("MCObsExpectation", 1000);
configTree_dummy->Branch ("MCObsExp", "TClonesArray", &tcMCObsExp_dummy);
/*TTree *configTree = new TTree ("configTree", "configuration Tree");
TClonesArray *tcdatasets = new TClonesArray ("Dataset", 1000);
configTree->Branch ("Datasets", "TClonesArray", &tcdatasets);
TClonesArray *tcAnaEnv = new TClonesArray ("AnalysisEnvironment", 1000);
configTree->Branch ("AnaEnv", "TClonesArray", &tcAnaEnv);
TClonesArray *tcMCObsExp = new TClonesArray ("MCObsExpectation", 1000);
configTree->Branch ("MCObsExp", "TClonesArray", &tcMCObsExp);*/

////////////////////////////////////
/// AnalysisEnvironment 
////////////////////////////////////
AnalysisEnvironment anaEnv;
cout << "Loading environment ..." << endl;
AnalysisEnvironmentLoader anaLoad (anaEnv, xmlfile);
new ((*tcAnaEnv_dummy)[0]) AnalysisEnvironment (anaEnv);
int verbose = anaEnv.Verbose;
float Luminosity = anaEnv.Luminosity;	// in 1/pb
int nbins = anaEnv.nbins;
int EventsPerBin=anaEnv.eventsperbin;
//  int nPseudos = anaEnv.nPseudoSession;
cout << "The results will be obtained for a luminosity of " << Luminosity << " x 1/pb" <<endl;


/////////////////////
// Load Datasets
/////////////////////
TTreeLoader treeLoader;
if (verbose > 0)
    cout << " - Load datasets ..." << endl;
vector < Dataset* > datasets;
vector < Dataset* > datasets_allconfigMC;
treeLoader.LoadDatasets (datasets, xmlfile);
for (unsigned int i = 0; i < datasets.size (); i++){
    new ((*tcdatasets_dummy)[i]) Dataset (*datasets[i]);
}

//this is something needed to run systematics on data... I know, chaotic... should be changed in the future...
datasets_allconfigMC = datasets;
for (unsigned int d = 0; d < datasets_allconfigMC.size (); d++){
     if(datasets_allconfigMC[d]->Title () == "Data"){
        datasets.clear();
	datasets.push_back(datasets_allconfigMC[d]);
	datasets_allconfigMC.erase(datasets_allconfigMC.begin()); //not very good, only if 'data' is the very first dataset in the config file...
	break;
     }
}
//after this, in 'datasets' only the data is present (well, unless you don't run on data but only MC), while datasets_allconfigMC contains all other MC samples that had "add=1" in config file...
cout<<"   MC datasets: "<<endl;
for (unsigned int d = 0; d < datasets_allconfigMC.size (); d++){
     cout<<"        "<<datasets_allconfigMC[d]->Title()<<", xsectionerror = "<<datasets_allconfigMC[d]->XsectionError()<<endl;
}
cout<<"--------------------------------"<<endl;
cout<<"   'datasets': "<<endl;
bool NPsamplesAdded = false;
for (unsigned int d = 0; d < datasets.size (); d++) {
     string dataSetTitle = datasets[d]->Title();
     string dataSetName = datasets[d]->Name();
     cout<<"        "<<dataSetTitle<< endl;
     if(dataSetName.find("prime")<=dataSetName.size() || dataSetTitle.find("prime")<=dataSetTitle.size()){
        NPsamplesAdded = true;
	massestprime.push_back(datasets[d]->Mass());
	crosssectionstprime.push_back(datasets[d]->Xsection());
	NofEventsTprimeTopTree.push_back(datasets[d]->eventTree()->GetEntriesFast() / datasets[d]->PreSelEfficiency());
	//cout<<"     nofevents t' samples? "<<datasets[d]->eventTree()->GetEntriesFast()<<endl;
     }
}

  
//////////////////////////////
//for statistical procedure
//////////////////////////////
TString Combination("sum");  // to combine weights to get V value
vector<float> FractionHWEvts = anaEnv.FractionHWEvts;
unsigned int nofFractionsHWEvts = FractionHWEvts.size();
//vector<unsigned int> sensitive_Obs; //for historical reasons; this is in principle not needed, but this vector should be kept empty and passed as an argument to one of the StatProcedure classes (should be removed/adapted in future)
vector< pair< float , float > > Vvect; //first entry: fraction x, second entry: V value (vector without adding uncertainty to estimation)
vector< pair< float , float > > FracNPEventsInHWSubsampleVect; //first entry: fraction x, second entry: fraction NP events in highest weight subsample
vector< pair< float , float > > SignNPEventsInHWSubsampleVect; //first entry: fraction x, second entry: signal significance in highest weight subsample (number of NP events / sqrt(number of background events))

//vector<float> Vvect_poisson;

if(massestprime.size()==0 && crosssectionstprime.size()==0){
  //putting dummy values...
  for(unsigned int mi=0;mi<numberofmasses;mi++){
    massestprime.push_back(massestprimeHardCoded[mi]);
    NofEventsTprimeTopTree.push_back(-1);
    crosssectionstprime.push_back(-1);  
  }
}

if(massestprime.size()==NofEventsTprimeTopTree.size() && massestprime.size()==crosssectionstprime.size()){
  cout<<endl;
  for(unsigned int miter=0;miter<massestprime.size();miter++){
    cout<<"massestprime["<<miter<<"] = "<<massestprime[miter]<<endl;
    cout<<"  NofEventsTprimeTopTree["<<miter<<"] = "<<NofEventsTprimeTopTree[miter]<<endl;
    cout<<"  crosssectionstprime["<<miter<<"] = "<<crosssectionstprime[miter]<<endl;
  }
  cout<<endl;
}

bool firstiteration = true;
int nHTbins = 12;
//'superloop' for scan of tprime space: t' masses
//cout<<"massestprime.size() = "<<massestprime.size()<<endl;
for(unsigned int miter=0;miter<massestprime.size();miter++){

 float xsect = xsectrange[0]; 
 if(scanTprimespace){
   cout<<"... Tprime scanning modus ..."<<endl;
   crosssectionstprime.clear();
   while(xsect>=xsectrange[0] && xsect<=xsectrange[1]){
     crosssectionstprime.push_back(xsect);
     xsect = xsect + xsectstep;
   }
 }
 else{
   cout<<"... No scanning modus ..."<<endl;
 }	  
 //'superloop' for scan of tprime space: cross sections
 for(unsigned int ixsect=0;ixsect<crosssectionstprime.size();ixsect++){   
  cout<<"  -->  At point ("<<massestprime[miter]<<" GeV, "<<crosssectionstprime[ixsect]<<" pb)"<<endl;
  
  //Output files... all "dummy's" shouldn't be saved in the end (obsolete)... for the moment just ignore.
  string rootFileName ("");
  string rootStatOutputFileName ("StatOutput.root");
  string rootStatVtreeFileName ("Vtree.root");
  ostringstream masstprime_sstream, crosssections_sstream, lumi_sstream; //nPseudoExp_sstream;
  masstprime_sstream << massestprime[miter];
  crosssections_sstream << crosssectionstprime[ixsect];
  lumi_sstream << anaEnv.Luminosity;
  
  if(anaEnv.isMC){
    if(!NPsamplesAdded){
      if(anaEnv.nPseudoExp<1){
        if(anaEnv.MCRound<1)
	   rootFileName = "SM_Binning" + lumi_sstream.str() + "_MC_Dummy.root";
	else{
           if(!doJESsystematics) rootFileName = thisROUNDdir + "SM_Lumi" + lumi_sstream.str() + "_MCObsExpFile.root";
	   else if(doJESShift == 1) rootFileName = thisROUNDdir + "SM_Lumi" + lumi_sstream.str() + "_MCObsExpFile_JESminus.root";
	   else if(doJESShift == 2) rootFileName = thisROUNDdir + "SM_Lumi" + lumi_sstream.str() + "_MCObsExpFile_JESplus.root";
	}       
	rootStatOutputFileName = "SM_MC_Dummy" + rootStatOutputFileName;
      }
      else if(anaEnv.nPseudoExp>=1){
        if(!Tprime){ //so you are doing the SM pseudoexperiments for the GOF-test: V tree only file that is important
          rootFileName = "GOF_SM_Dummy.root";
	  rootStatOutputFileName = "GOF_SM_Dummy_" + rootStatOutputFileName;
          rootStatVtreeFileName = thisROUNDdir + "GOF_SM_Lumi" + lumi_sstream.str() + "_" + SignificanceCombination + "_" + rootStatVtreeFileName;
        }
	else if(Tprime){ //so you are doing the SM pseudoexperiments for the Tprime analysis: V tree only file that is important
	  rootFileName = "SMTprime_Dummy.root";
	  rootStatOutputFileName = "SMTprime_SM_Dummy_" + rootStatOutputFileName;
          if(SignificanceCombination == "SUM_LRatioMinded") rootStatVtreeFileName = thisROUNDdir + "TpAnalysis_SMonly_" + masstprime_sstream.str() + "GeV_" + crosssections_sstream.str() + "pb_Lumi" + lumi_sstream.str() + "_" + SignificanceCombination + "_" + rootStatVtreeFileName;
	  else rootStatVtreeFileName = thisROUNDdir + "TpAnalysis_SMonly_Lumi" + lumi_sstream.str() + "_" + SignificanceCombination + "_" + rootStatVtreeFileName; 
	}
      }
    }
    else if(NPsamplesAdded){
      if(Tprime){
        if(anaEnv.nPseudoExp<1){
           if(!doJESsystematics) rootFileName = thisROUNDdir + "SMTprime_" + masstprime_sstream.str() + "GeV_" + crosssections_sstream.str() + "pb_Lumi" + lumi_sstream.str() + "_" + "MCObsExpFile.root";
           else if(doJESShift == 1) rootFileName = thisROUNDdir + "SMTprime_" + masstprime_sstream.str() + "GeV_" + crosssections_sstream.str() + "pb_Lumi" + lumi_sstream.str() + "_" + "MCObsExpFile_JESminus.root";
	   else if(doJESShift == 2) rootFileName = thisROUNDdir + "SMTprime_" + masstprime_sstream.str() + "GeV_" + crosssections_sstream.str() + "pb_Lumi" + lumi_sstream.str() + "_" + "MCObsExpFile_JESplus.root";
	   rootStatOutputFileName = "SMTprime_MC_Dummy" + rootStatOutputFileName;
        }
	if(anaEnv.nPseudoExp>=1){ //so you are doing SM+NP pseudoexperiments for the Tprime analysis
	   //rootFileName = masstprime_sstream.str() + "GeV_" + crosssections_sstream.str() + "pb_" + rootFileName;
           rootFileName = "SMTprime_Dummy.root";
	   rootStatOutputFileName = "SMTprime_Dummy_" + rootStatOutputFileName;
           rootStatVtreeFileName = thisROUNDdir + "SMTprime_" + masstprime_sstream.str() + "GeV_" + crosssections_sstream.str() + "pb_Lumi" + lumi_sstream.str() + "_" + SignificanceCombination + "_" + rootStatVtreeFileName;	
	}
      }
    }
  }
  if(!anaEnv.isMC){ // this means: you run on data...
     rootFileName = "Data_Dummy.root";
     rootStatOutputFileName = thisROUNDdir + "Data_Lumi" + lumi_sstream.str() + "_" + rootStatOutputFileName; //should be changed!!
     rootStatVtreeFileName = thisROUNDdir + "Data_Lumi" + lumi_sstream.str() + "_" + masstprime_sstream.str() + "GeV_" + crosssections_sstream.str() + "pb_" + rootStatVtreeFileName;
  }


  //Configuration output format
  TTree *configTree = new TTree ("configTree", "configuration Tree");
  TClonesArray *tcdatasets = new TClonesArray ("Dataset", 1000);
  configTree->Branch ("Datasets", "TClonesArray", &tcdatasets);
  TClonesArray *tcAnaEnv = new TClonesArray ("AnalysisEnvironment", 1000);
  configTree->Branch ("AnaEnv", "TClonesArray", &tcAnaEnv);
  TClonesArray *tcMCObsExp = new TClonesArray ("MCObsExpectation", 1000);
  configTree->Branch ("MCObsExp", "TClonesArray", &tcMCObsExp);


  for (unsigned int d = 0; d < datasets.size (); d++) {
      cout << "   Dataset " << d << ": " << datasets[d]->Title () <<" with events  "<< datasets[d]->NofEvtsToRunOver ()<<endl;      
      string dataSetName = datasets[d]->Name();
      if(scanTprimespace){
        if(dataSetName.find("prime")<=dataSetName.size() && datasets[d]->Mass()==massestprime[miter]){
	  cout<<"!! Resetting cross section AND equivalent luminosity of current Tprime sample!!"<<endl;
	  datasets[d]->SetXsection(crosssectionstprime[ixsect]);
	  float tempEqLumi = NofEventsTprimeTopTree[miter] /  crosssectionstprime[ixsect];
	  cout<<"   NofEventsTprimeTopTree[miter] = "<<NofEventsTprimeTopTree[miter]<<endl;
	  datasets[d]->SetEquivalentLuminosity(tempEqLumi);
	  cout<<"   resetted xsection: "<<datasets[d]->Xsection()<<", resetted eqLumi: "<<datasets[d]->EquivalentLumi()<<endl;
          break;
	}
      } 
  }

  Vvect.clear();
  FracNPEventsInHWSubsampleVect.clear();
  SignNPEventsInHWSubsampleVect.clear();

	
  //nof selected events
  float NEvtsData = 0;
  Double_t *nEvents = new Double_t[datasets.size()];

  vector<int > nofSelEvts_vector; //(z) added


  ///////////////
  // Observables
  Observables obs;


  vector <string > lstVar;
  lstVar = obs.ListOfVariables ();
	    //Loading of the binning
  map < string, TAxis * > mapAxis;
  TFile fBinning (anaEnv.binningFile.c_str (), "READ");
  float nBinsBinningFile = 0;
  //Binning for the given observable
  for (vector < string >::iterator iter = lstVar.begin (); iter != lstVar.end (); iter++) {
    TAxis *axis = NULL;
    char tAxisName[150];

    sprintf (tAxisName, "Binning_%s_SM", iter->c_str ());
    fBinning.GetObject (tAxisName, axis);
  
    if (axis == NULL)
      {
      cout<<" While reading the  "<<tAxisName<<"  the axis is NULL ..... "<<endl;
      mapAxis.insert (make_pair (*iter, (TAxis *) NULL));}
    else{
      mapAxis.insert (make_pair (*iter, new TAxis (*axis))); 
	/*for( map<string, TAxis*>::iterator ii=mapAxis.begin(); ii!=mapAxis.end(); ++ii)
		{ 
		cout << (*ii).first << ": " <<(*ii).second->GetNbins()<< " Title  "<<(*ii).second->GetName()<<" min "<<(*ii).second->GetXmin()<<"  max   "<<(*ii).second->GetXmax()<<endl;
		}
	*/
      nBinsBinningFile = axis->GetNbins ();
      cout<<" nBinsBinningFile = "<<nBinsBinningFile<<endl;
    }
      	//delete axis;
  }  
  fBinning.Close ();

  vector<pair<string,float> > BinContentObsAll_;
  vector<pair<string,float> > BinContentObsPerSample_;  
  TFile *fout = new TFile (rootFileName.c_str (), "RECREATE");
  TFile* foutStat = new TFile(rootStatOutputFileName.c_str(),"RECREATE");

  cout<<"Now creating directories in the Stat root file..."<<endl;
  TDirectory * mytdirFraction;
  TDirectory * mytdirVariable;
  for(unsigned int k=0;k<nofFractionsHWEvts;k++){
  	ostringstream Fractionstrstream;
        Fractionstrstream << FractionHWEvts[k];
  	mytdirFraction = foutStat->mkdir("x="+TString(Fractionstrstream.str()));
	for(unsigned int v=0;v<lstVar.size (); v++){
	//	if (!anaEnv.runOnObs ((int) v))
      	//		continue;
		mytdirFraction->cd();
		char varname[50];
		sprintf(varname,lstVar[v].c_str());
		//cout<<"... In fraction directory in TFile: creating directory for variable "<<lstVar[v]<<endl;
		mytdirVariable = mytdirFraction->mkdir(varname);
	}
  }
  TDirectory * mytdirall = foutStat->mkdir("AllEvents");
  for(unsigned int v=0;v<lstVar.size (); v++){
	//	if (!anaEnv.runOnObs ((int) v))
      	//		continue;
		mytdirall->cd();
		char varname[50];
		sprintf(varname,lstVar[v].c_str());
		//cout<<"... In allevents directory in TFile: creating directory for variable "<<lstVar[v]<<endl;
		mytdirVariable = mytdirall->mkdir(varname);
   }
  cout<<"  Directories created"<<endl;
  int* ievtmin = new int[datasets.size()];
  for(unsigned int i=0;i<datasets.size();i++) ievtmin[i] = 0;

  char name[100];
  

  ////////////////////////////////////
  //  MCExpectation (if MC)
  ////////////////////////////////////
  MCExpectation *MCExp = NULL;
  vector<MCObsExpectation *> MCObsExp;
  vector<MCObsExpectation *> MCObsExp_NP;
  map<int,vector<MCObsExpectation*> > MCObsExp_bins;
  map<int,vector<MCObsExpectation*> > MCObsExp_NP_bins;
  
  //JES varied 'histograms'
  vector<MCObsExpectation *> MCObsExp_JESminus;
  vector<MCObsExpectation *> MCObsExp_JESplus;
  vector<MCObsExpectation *> MCObsExp_NP_JESminus;
  vector<MCObsExpectation *> MCObsExp_NP_JESplus;
  map<int,vector<MCObsExpectation*> > MCObsExp_JESminus_bins;  
  map<int,vector<MCObsExpectation*> > MCObsExp_JESplus_bins;
  map<int,vector<MCObsExpectation*> > MCObsExp_NP_JESminus_bins;
  map<int,vector<MCObsExpectation*> > MCObsExp_NP_JESplus_bins;  
     
  if (anaEnv.isMC && anaEnv.nPseudoExp<1) {
    //to be adapted?
    //create as first step
    MCExp = new MCExpectation ();
    MCExp->SetLuminosity (Luminosity);
    MCExp->SetLabel (string ("CrossSection"));
  }
  else {
  
   if(!doWithSpecial2Dbinning){
  //  string fin_name("../config/");
    string fin_name = thisROUNDdir + "SM_Lumi" + lumi_sstream.str() + "_MCObsExpFile.root";
    //if(SignificanceCombination=="SUM") fin_name = fin_name + "SM_Lumi" + lumi_sstream.str() + "_MCObsExpFile.root";
    //if(SignificanceCombination=="SUM_LRatioMinded") fin_name = fin_name + "SMTprime_" + masstprime_sstream.str() + "GeV_" + crosssections_sstream.str() + "pb_Lumi" + lumi_sstream.str() + "_" + "MCObsExpFile.root";
    cout<<"Reading input MCObsExp file for SM only... --> "<<fin_name<<endl;
    TFile fin(fin_name.c_str (), "READ");
    TClonesArray *tcMCObsExpTemp = new TClonesArray ("MCObsExpectation", 0);
    TTree *tree = (TTree *) fin.Get ("configTree");
    TBranch *b3 = 0;
    if (tree) {     
      b3 = tree->GetBranch ("MCObsExp");
      if (b3)
	b3->SetAddress (&tcMCObsExpTemp);
      tree->GetEntry (0);
    }
    MCObsExpectation *temp3 = 0;
    for (unsigned x = 0; x < tcMCObsExpTemp->GetEntries (); x++) {
      temp3 = (MCObsExpectation *) tcMCObsExpTemp->At (x);
      MCObsExp.push_back ((MCObsExpectation *) temp3->Clone ());
    }
//    delete temp3;
//    tcMCObsExpTemp->Delete ();
//    tree->Delete ();
    
    //regarding MCFile of SM+NP //name: TO BE CHECKED!!!!!
    if(SignificanceCombination=="SUM_LRatioMinded"){
////    string fin_NPname("../config/");
      string fin_NPname = thisROUNDdir + "SMTprime_" + masstprime_sstream.str() + "GeV_" + crosssections_sstream.str() + "pb_Lumi" + lumi_sstream.str() + "_MCObsExpFile.root";    
      cout<<"Reading input MCObsExp file for SM+NP... (SignificanceCombination = "<<SignificanceCombination<<") --> "<<fin_NPname<<endl;
      TFile fin_NP(fin_NPname.c_str(), "READ");
      //if(fin_NP.IsOpen()) cout<<"File Open"<<endl;
      //else{cout<<"file not open?"<<endl;}
      TClonesArray *tcMCObsExpTemp_NP = new TClonesArray ("MCObsExpectation", 0);//tcMCObsExpTemp = new TClonesArray ("MCObsExpectation", 0);
      TTree * tree_NP = (TTree *) fin_NP.Get ("configTree");//tree = (TTree *) fin_NP.Get ("configTree");
      TBranch *b3_NP = 0;//b3 = 0;
      if (tree_NP) {     
        b3_NP = tree_NP->GetBranch ("MCObsExp");
        if (b3_NP)
	  b3_NP->SetAddress (&tcMCObsExpTemp_NP);
        tree_NP->GetEntry (0);
      }
      MCObsExpectation *temp3_NP = 0;//temp3 = 0;
      for (unsigned x = 0; x < tcMCObsExpTemp_NP->GetEntries (); x++) {
        temp3_NP = (MCObsExpectation *) tcMCObsExpTemp_NP->At (x);
        MCObsExp_NP.push_back ((MCObsExpectation *) temp3_NP->Clone ());
      }
      delete temp3_NP;
      tcMCObsExpTemp_NP->Delete ();
      tree_NP->Delete ();
     }
   } //end if(!doWithSpecial2Dbinning)
   else if(doWithSpecial2Dbinning){

    int nHTbins = 12;
    for(int HTbini=1;HTbini<nHTbins+1;HTbini++){
       string fin_name("");
       stringstream sstr_HTbini;
       sstr_HTbini << HTbini;
       fin_name = "./ROUND19/HTbin" + sstr_HTbini.str() + "/SM_Lumi" + lumi_sstream.str() + "_MCObsExpFile.root";
       cout<<"Reading input MCObsExp file for SM only, HT bin "<<HTbini<<"... --> "<<fin_name<<endl;
       TFile fin(fin_name.c_str (), "READ");
       TClonesArray *tcMCObsExpTemp = new TClonesArray ("MCObsExpectation", 0);
       TTree *tree = (TTree *) fin.Get ("configTree");
       TBranch *b3 = 0;
       if (tree) {     
          b3 = tree->GetBranch ("MCObsExp");
          if (b3)
	     b3->SetAddress (&tcMCObsExpTemp);
          tree->GetEntry (0);
       }
       MCObsExpectation *temp3 = 0;
       vector<MCObsExpectation*> MCObsExp_bini;
       for (unsigned x = 0; x < tcMCObsExpTemp->GetEntries (); x++) {
           temp3 = (MCObsExpectation *) tcMCObsExpTemp->At (x);
           MCObsExp_bini.push_back ((MCObsExpectation *) temp3->Clone ()); 
       }
       MCObsExp_bins[HTbini] = MCObsExp_bini;
       //delete temp3;//added
       //tcMCObsExpTemp_NP->Delete ();//added
       //tree_NP->Delete ();//added
              
       string fin_NPname = "./ROUND19/HTbin" + sstr_HTbini.str() + "/SMTprime_" + masstprime_sstream.str() + "GeV_" + crosssections_sstream.str() + "pb_Lumi" + lumi_sstream.str() + "_MCObsExpFile.root";    
       cout<<"Reading input MCObsExp file for SM+NP, HT bin "<<HTbini<<"... --> "<<fin_NPname<<endl;
       TFile fin_NP(fin_NPname.c_str(), "READ");
       TClonesArray *tcMCObsExpTemp_NP = new TClonesArray ("MCObsExpectation", 0);//tcMCObsExpTemp = new TClonesArray ("MCObsExpectation", 0);
       TTree * tree_NP = (TTree *) fin_NP.Get ("configTree");//tree = (TTree *) fin_NP.Get ("configTree");
       TBranch *b3_NP = 0;//b3 = 0;
       if (tree_NP) {     
          b3_NP = tree_NP->GetBranch ("MCObsExp");
          if (b3_NP)
	     b3_NP->SetAddress (&tcMCObsExpTemp_NP);
          tree_NP->GetEntry (0);
       }
       MCObsExpectation *temp3_NP = 0;//temp3 = 0;
       vector<MCObsExpectation*> MCObsExp_NP_bini;
       for (unsigned x = 0; x < tcMCObsExpTemp_NP->GetEntries (); x++) {
          temp3_NP = (MCObsExpectation *) tcMCObsExpTemp_NP->At (x);
          MCObsExp_NP_bini.push_back ((MCObsExpectation *) temp3_NP->Clone ());
       }
       MCObsExp_NP_bins[HTbini] = MCObsExp_NP_bini;
       delete temp3_NP;
       tcMCObsExpTemp_NP->Delete ();
       tree_NP->Delete ();
       
       
       if(doJESsystematics){
        //JES minus
       string fin_JESminusname = "./ROUND19/HTbin" + sstr_HTbini.str() + "/SM_Lumi" + lumi_sstream.str() + "_MCObsExpFile_JESminus.root";
       cout<<"Reading input MCObsExp file for SM only (JES minus), HT bin "<<HTbini<<"... --> "<<fin_JESminusname<<endl;
       TFile fin_JESminus(fin_JESminusname.c_str (), "READ");
       TClonesArray *tcMCObsExpTemp_JESminus = new TClonesArray ("MCObsExpectation", 0);
       TTree *tree_JESminus = (TTree *) fin_JESminus.Get ("configTree");
       TBranch *b3_JESminus = 0;
       if (tree_JESminus) {     
          b3_JESminus = tree_JESminus->GetBranch ("MCObsExp");
          if (b3_JESminus)
	     b3_JESminus->SetAddress (&tcMCObsExpTemp_JESminus);
          tree_JESminus->GetEntry (0);
       }
       MCObsExpectation *temp3_JESminus = 0;
       vector<MCObsExpectation*> MCObsExp_JESminus_bini;
       for (unsigned x = 0; x < tcMCObsExpTemp_JESminus->GetEntries (); x++) {
           temp3_JESminus = (MCObsExpectation *) tcMCObsExpTemp_JESminus->At (x);
           MCObsExp_JESminus_bini.push_back ((MCObsExpectation *) temp3_JESminus->Clone ()); 
       }
       MCObsExp_JESminus_bins[HTbini] = MCObsExp_JESminus_bini;

       
         string fin_JESminusNPname = "./ROUND19/HTbin" + sstr_HTbini.str() + "/SMTprime_" + masstprime_sstream.str() + "GeV_" + crosssections_sstream.str() + "pb_Lumi" + lumi_sstream.str() + "_MCObsExpFile_JESminus.root";    
         cout<<"Reading input MCObsExp file for SM+NP (JES minus), HT bin "<<HTbini<<"... --> "<<fin_JESminusNPname<<endl;
         TFile fin_JESminusNP(fin_JESminusNPname.c_str(), "READ");
         TClonesArray *tcMCObsExpTemp_NP_JESminus = new TClonesArray ("MCObsExpectation", 0);//tcMCObsExpTemp = new TClonesArray ("MCObsExpectation", 0);
         TTree * tree_NP_JESminus = (TTree *) fin_JESminusNP.Get ("configTree");//tree = (TTree *) fin_NP.Get ("configTree");
         TBranch *b3_NP_JESminus = 0;//b3 = 0;
         if (tree_NP_JESminus) {     
            b3_NP_JESminus = tree_NP_JESminus->GetBranch ("MCObsExp");
            if (b3_NP_JESminus)
	       b3_NP_JESminus->SetAddress (&tcMCObsExpTemp_NP_JESminus);
            tree_NP_JESminus->GetEntry (0);
         }
         MCObsExpectation *temp3_NP_JESminus = 0;//temp3 = 0;
         vector<MCObsExpectation*> MCObsExp_NP_JESminus_bini;
         for (unsigned x = 0; x < tcMCObsExpTemp_NP_JESminus->GetEntries (); x++) {
            temp3_NP_JESminus = (MCObsExpectation *) tcMCObsExpTemp_NP_JESminus->At (x);
            MCObsExp_NP_JESminus_bini.push_back ((MCObsExpectation *) temp3_NP_JESminus->Clone ());
         } 
         MCObsExp_NP_JESminus_bins[HTbini] = MCObsExp_NP_JESminus_bini;
         delete temp3_NP_JESminus;
         tcMCObsExpTemp_NP_JESminus->Delete ();
         tree_NP_JESminus->Delete ();    

//for debugging
      /*TH1F* estim_SM=0;
      TH1F* estim_NP=0;
      TH1F* estim_SM_JESminus=0;
      TH1F* estim_SM_JESplus=0;
      TH1F* estim_NP_JESminus=0;
      TH1F* estim_NP_JESplus=0;
      cout<<"OUTSIDE 1"<<endl;
      //copy from MCObsExpectation to the vector
      for (unsigned int l=0;l<MCObsExp_JESminus_bins[HTbini].size();l++){
		if (MCObsExp_JESminus_bins[HTbini][l]->GetHistoSMProcesses()->GetName()=="SMProcess_MassHadTop"){
		   cout<<endl<<"  found MCSMObsExp plots....... "<<MCObsExp_JESminus_bins[HTbini][l]->GetHistoSMProcesses()->GetName()<<"   "<<l<<endl;           
		   estim_SM = MCObsExp_bins[HTbini][l]->GetHistoSMProcesses();
		   estim_NP = MCObsExp_NP_bins[HTbini][l]->GetHistoAll();
		   cout<<"OUTSIDE 2"<<endl;
		     cout<<"INSIDE 3"<<endl;
		     estim_SM_JESminus = MCObsExp_JESminus_bins[HTbini][l]->GetHistoSMProcesses();
		     cout<<"OUTSIDE 3.1"<<endl;
		     estim_SM_JESplus = MCObsExp_JESplus_bins[HTbini][l]->GetHistoSMProcesses();
		     cout<<"OUTSIDE 3.2"<<endl;
		     estim_NP_JESminus = MCObsExp_NP_JESminus_bins[HTbini][l]->GetHistoAll();
		     cout<<"OUTSIDE 3.3"<<endl;
		     estim_NP_JESplus = MCObsExp_NP_JESplus_bins[HTbini][l]->GetHistoAll();
		   cout<<"OUTSIDE 4"<<endl;
		   break;
	        }
      }*/
    
       //JES plus
         string fin_JESplusname = "./ROUND19/HTbin" + sstr_HTbini.str() + "/SM_Lumi" + lumi_sstream.str() + "_MCObsExpFile_JESplus.root";
         cout<<"Reading input MCObsExp file for SM only (JES plus), HT bin "<<HTbini<<"... --> "<<fin_JESplusname<<endl;
     	 TFile fin_JESplus(fin_JESplusname.c_str (), "READ");
    	 TClonesArray *tcMCObsExpTemp_JESplus = new TClonesArray ("MCObsExpectation", 0);
    	 TTree *tree_JESplus = (TTree *) fin_JESplus.Get ("configTree");
    	 TBranch *b3_JESplus = 0;
    	 if (tree_JESplus) {     
           b3_JESplus = tree_JESplus->GetBranch ("MCObsExp");
           if (b3_JESplus)
	     b3_JESplus->SetAddress (&tcMCObsExpTemp_JESplus);
           tree_JESplus->GetEntry (0);
         }
         MCObsExpectation *temp3_JESplus = 0;
	 vector<MCObsExpectation*> MCObsExp_JESplus_bini;
         for (unsigned x = 0; x < tcMCObsExpTemp_JESplus->GetEntries (); x++) {
           temp3_JESplus = (MCObsExpectation *) tcMCObsExpTemp_JESplus->At (x);
           MCObsExp_JESplus_bini.push_back ((MCObsExpectation *) temp3_JESplus->Clone ());
         }
	 MCObsExp_JESplus_bins[HTbini] = MCObsExp_JESplus_bini;
       
         string fin_JESplusNPname = "./ROUND19/HTbin" + sstr_HTbini.str() + "/SMTprime_" + masstprime_sstream.str() + "GeV_" + crosssections_sstream.str() + "pb_Lumi" + lumi_sstream.str() + "_MCObsExpFile_JESplus.root";    
         cout<<"Reading input MCObsExp file for SM+NP (JES plus), HT bin "<<HTbini<<"... --> "<<fin_JESplusNPname<<endl;
         TFile fin_JESplusNP(fin_JESplusNPname.c_str(), "READ");
         TClonesArray *tcMCObsExpTemp_NP_JESplus = new TClonesArray ("MCObsExpectation", 0);//tcMCObsExpTemp = new TClonesArray ("MCObsExpectation", 0);
         TTree * tree_NP_JESplus = (TTree *) fin_JESplusNP.Get ("configTree");//tree = (TTree *) fin_NP.Get ("configTree");
         TBranch *b3_NP_JESplus = 0;//b3 = 0;
         if (tree_NP_JESplus) {     
            b3_NP_JESplus = tree_NP_JESplus->GetBranch ("MCObsExp");
            if (b3_NP_JESplus)
	       b3_NP_JESplus->SetAddress (&tcMCObsExpTemp_NP_JESplus);
            tree_NP_JESplus->GetEntry (0);
         }
         MCObsExpectation *temp3_NP_JESplus = 0;//temp3 = 0;
         vector<MCObsExpectation*> MCObsExp_NP_JESplus_bini;
         for (unsigned x = 0; x < tcMCObsExpTemp_NP_JESplus->GetEntries (); x++) {
            temp3_NP_JESplus = (MCObsExpectation *) tcMCObsExpTemp_NP_JESplus->At (x);
            MCObsExp_NP_JESplus_bini.push_back ((MCObsExpectation *) temp3_NP_JESplus->Clone ());
         } 
         MCObsExp_NP_JESplus_bins[HTbini] = MCObsExp_NP_JESplus_bini;
         delete temp3_NP_JESplus;
         tcMCObsExpTemp_NP_JESplus->Delete ();
         tree_NP_JESplus->Delete ();
       } 
            
    }
 
  }
  
 }  

  vector<pair<float,float > > historanges;
  ///////////////////// choosing the Variables to run on...
  if(!Tprime){  ////GOF-test
     lstVar.clear();
     lstVar.push_back("ET1oET4");
     lstVar.push_back("AllJetsPtMET");
     lstVar.push_back("MET");
     lstVar.push_back("PtMuon");
  }
  if(Tprime){
     lstVar.clear();
   /*  lstVar.push_back("MassHadTop");
     lstVar.push_back("MassTtbar");
     lstVar.push_back("HT4jets");*/
     lstVar.push_back("HT4jetsMuonMET");
     historanges.push_back(make_pair(0,3000));
     lstVar.push_back("MassHadTop");
     historanges.push_back(make_pair(0,3000));
     //just for kinematic plots after b-tagging
     lstVar.push_back("MET");
     historanges.push_back(make_pair(40,1000));
     lstVar.push_back("PtMuon");
     historanges.push_back(make_pair(35,500));
     lstVar.push_back("PtQuark1");
     historanges.push_back(make_pair(35,1000));
     lstVar.push_back("PtQuark2");
     historanges.push_back(make_pair(35,1000));
     lstVar.push_back("PtLepbquark");
     historanges.push_back(make_pair(35,1000));
     lstVar.push_back("PtHadbquark");
     historanges.push_back(make_pair(35,1000));     
     
     ///////////////lstVar.push_back("EtaHadTop"); //not used anymore
  }
  cout<<"We will perform the statistical method for "<<lstVar.size()<<" variables:"<<endl;
  for(unsigned int vari=0;vari<lstVar.size();vari++){
     cout<<lstVar[vari]<<endl;
  }
  
      //if you do this there will be problems if the TAxes of the other variables are NULL etc?

 map<int,TH1F* > DataHistos_Special2Dbinning;
 
 if(!doPoissonPseudos){
  ////////////////////////////////////
  //    Loop on  Pseudo-experiments
  ////////////////////////////////////
  int nRuns = 1;
  if (anaEnv.nPseudoExp != 0)
    nRuns = anaEnv.nPseudoExp;
  for (int iPseudoExp = 0; iPseudoExp < nRuns; iPseudoExp++) {

    cout << " ** PseudoExp no: " << iPseudoExp << "/" << nRuns-1 << " **" << endl;
    //////////////////////////////////////// 
    // For the StatProcedure
    //////////////////////////////////////// 
    vector < vector < float > > VarValues; //first dim: events ; second dim: variables
    vector < unsigned int > NPbooleans;  //to indicate if an event is NP (1) or not (0)
    vector < vector < pair < string, float > > > VarValues_WithString; //first dim: events ; second dim: variables with string...
    
    unsigned int isNP = 0;
    vector < TH1F * > ObsEstimated;
    vector < TH1F * > ObsData;

    Double_t *nEvents = new Double_t[datasets.size ()];         
    

    //For Background summary
    TH1F **hData = new TH1F *[lstVar.size ()];

    //for data... with the special 2D binning in HT bins
    map < int, TAxis * > AxesDataHistos_Special2Dbinning;
    map<int,TH1F* > DataHistos_Special2Dbinning;
    int nBinsBinningFile = 0;
//    TFile currentfBinning(TString("./ROUND11/Binning_12_Bins_TTJetsSemiMuFlat_MassHadTop_HTbin")+TString(sstr_HTbini.str())+".root", "READ");
    TFile* currentfBinning = 0;
    if(doWithSpecial2Dbinning && !anaEnv.isMC){
      for(int HTbini=1;HTbini<nHTbins+1;HTbini++){
       stringstream sstr_HTbini;
       sstr_HTbini << HTbini;
       
       currentfBinning = TFile::Open(TString("./ROUND19/Binning_14_Bins_TTJetsSemiMuFlat_MassHadTop_HTbin")+TString(sstr_HTbini.str())+".root", "READ");
       TAxis *axisi = NULL;
       char tAxisiName[150];
       sprintf (tAxisiName, "Binning_%s_SM", "MassHadTop"); //hardcoded
       cout<<" tAxisiName = "<<tAxisiName<<endl;
       currentfBinning->GetObject (tAxisiName, axisi);  
       AxesDataHistos_Special2Dbinning.insert (make_pair (HTbini, new TAxis (*axisi))); 
	/*for( map<string, TAxis*>::iterator ii=AxesDataHistos_Special2Dbinningmap.begin(); ii!=AxesDataHistos_Special2Dbinningmap.end(); ++ii)
		{ 
		cout << (*ii).first << ": " <<(*ii).second->GetNbins()<< " Title  "<<(*ii).second->GetName()<<" min "<<(*ii).second->GetXmin()<<"  max   "<<(*ii).second->GetXmax()<<endl;
		}
	*/
        nBinsBinningFile = axisi->GetNbins ();
        cout<<" nBinsCurrentBinningFile = "<<nBinsBinningFile<<endl;
      	//delete axis;             
        DataHistos_Special2Dbinning[HTbini] = new TH1F (TString ("hData_") + TString(sstr_HTbini.str()), "Data", nBinsBinningFile, axisi->GetXbins ()->fArray);   
        cout<<"DataHistos_Special2Dbinning["<<HTbini<<"]->GetName() = "<<DataHistos_Special2Dbinning[HTbini]->GetName()<<endl;
//	currentfBinning->Close(); //if closed, crashes
      }
    }
//    for(int HTbini=1;HTbini<nHTbins+1;HTbini++){
//       cout<<"DataHistos_Special2Dbinning["<<HTbini<<"]->GetName() = "<<DataHistos_Special2Dbinning[HTbini]->GetName()<<endl;
//    }
    
    //LOOP OVER VARIABLES ... only for shape estimation
    if (verbose > 0)
      cout << "Run on " << lstVar.size () << " variables " << endl;
   if(!do2Dhisto){
    for (unsigned int i = 0; i < lstVar.size (); i++) {
      /*if (!anaEnv.runOnObs ((int) i)){
	if(anaEnv.isMC && anaEnv.nPseudoExp<1) MCObsExp.push_back(new MCObsExpectation (nbins, obs.RangeVariables ()[i].first, obs.RangeVariables ()[i].second, lstVar[i], Luminosity));
	continue;
      }*/    //why was this??
      if (verbose > 0)
	cout << "For the " << i << " st variable " << lstVar[i] << " with nbins read out from config file... " <<nbins<<" and mapaxis.size....." <<mapAxis.size()<< endl;

      TAxis *axis = mapAxis[lstVar[i]];
      // cout<<"  AND FOR EVERY VARIABLE THE AXIS RANGE  "<<lstVar[i]<<"   "<<axis->GetNbins()<<"  "<<axis->GetXbins ()->fArray<<"  "<<axis->GetXbins ()<<"  "<<axis->GetName()<<endl;
  
      string decision="Bins";
     

      if (verbose > 0)
	cout << " - Create histos used by BkgEstimationSummary ..." << endl;
      if (axis == NULL) {
        cout<<"i = "<<i<<endl;
        cout<<"Warning: axis == NULL"<<endl;
	if(anaEnv.isMC && anaEnv.nPseudoExp<1){
	  cout<<"nbins = "<<nbins<<", range obs: "<<historanges[i].first<<" -> "<<historanges[i].second<<", variable "<<lstVar[i]<<endl;
	  MCObsExp.push_back(new MCObsExpectation (datasets, nbins, historanges[i].first, historanges[i].second, lstVar[i], Luminosity));	  	  
	}
	hData[i] = new TH1F (TString ("hData_") + lstVar[i], "Data", nbins, historanges[i].first, historanges[i].second);
      }
      else {
	if (axis->GetNbins () == 0) {
	  cout<<"Warning: axis->GetNbins () == 0"<<endl;
	  if(anaEnv.isMC && anaEnv.nPseudoExp<1) MCObsExp.push_back(new MCObsExpectation (datasets, axis, lstVar[i], Luminosity));
	  hData[i] = new TH1F (TString ("hData_") + lstVar[i], "Data", nbins, historanges[i].first, historanges[i].second);
	}
	else {
	  if(anaEnv.isMC && anaEnv.nPseudoExp<1)
            MCObsExp.push_back(new MCObsExpectation (datasets, axis,lstVar[i], Luminosity));	  
	  hData[i] = new TH1F (TString ("hData_") + lstVar[i], "Data", axis->GetNbins (), axis->GetXbins ()->fArray);
        }
      }

    }				//end of loop over variables
   }
   else{
       TAxis *axis_x = mapAxis[lstVar[0]];
       TAxis *axis_y = mapAxis[lstVar[1]];
       if (axis_x == NULL || axis_y == NULL) {
         cout<<"Warning: axis == NULL"<<endl;
       }
       else {
          if(anaEnv.isMC && anaEnv.nPseudoExp<1)
            MCObsExp.push_back(new MCObsExpectation (datasets, axis_x, axis_y,lstVar[0],lstVar[1], Luminosity));	  
	  //hData[i] = new TH1F (TString ("hData_") + lstVar[i], "Data", axis->GetNbins ()
       }
   }


   
    float oldfillweight = 0, fillweight = 0., fillweight_TTJets = 0; //(z) initialized to 0

    
//   if(anaEnv.nPseudoExp>0 && doPoissonPseudos){
//         cout<<"... We will NOT loop over events!"<<endl;
//   }
//   else{
    //////////////////////////////////////////////
    //loop on datasets
    if (verbose > 0)
      cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
    for (unsigned int d = 0; d < datasets.size (); d++) {
      cout << "   Dataset " << d << ": " << datasets[d]->Title () <<" with events  "<< datasets[d]->NofEvtsToRunOver ()<<endl;      
      string dataSetName = datasets[d]->Name();
      if(dataSetName.find("prime")<=dataSetName.size() && datasets[d]->Mass()!=massestprime[miter]){
         continue; //i.e. do not consider this tprime sample if the mass in the loop does not correspond to the mass of the t' sample...
      }
      cout << "datasets NormFactor: " << datasets[d]->NormFactor () << endl;
     
      ///////////////////////////////////////
      BinContentObsPerSample_.clear();
    
      int nofSelEvts = 0;

      if(anaEnv.nPseudoExp == 0){ //or at least, also when it is not data...
      	//fillweight = datasets[d]->NormFactor () * Luminosity;
        oldfillweight = datasets[d]->NormFactor () * Luminosity; //* F; //replaced by a scaling of the histo's in a later stage, when you read them... BTW, for plotmaker.cc of so, make sure you scale there too if you want to plot the rescaled histo's
      }
      //else if(datasets[d]->Name() == "Data"){oldfillweight = 1;}
      else
      	oldfillweight = 1; //when data or pseudoexperiment... wait, now with F, fillweight != 1 for pseudoexps... //update: this F-thing is abandonded?     
    //  }
    //  else
    //  	fillweight = F;
      cout<<"for dataset  "<<datasets[d]->Title()<<"  oldfillweight = "<<oldfillweight<<" and normfactor = "<<datasets[d]->NormFactor()<<endl;

      if (datasets[d]->Name()=="TTJets") fillweight_TTJets = datasets[d]->NormFactor()*Luminosity; //for binning...
      
      MCPseudoExp pseudoExp;
      int nPseudos=anaEnv.nPseudoExp;
      /*int PseudoSessions=anaEnv.nPseudoSession;
      bool manyPseudoExp;
      if (nPseudos <1 && PseudoSessions==0) manyPseudoExp=0;
      if (nPseudos <1 && PseudoSessions!=0) manyPseudoExp=1;
      if (nPseudos >0 && PseudoSessions==0) manyPseudoExp=0;
      if (nPseudos >0 && PseudoSessions!=0) manyPseudoExp=1;
      if (nPseudos >0 && PseudoSessions!=0 && datasets[d]->Name()=="Data") manyPseudoExp=0;*/

     //to be used when runonTTrees is true
      string filename;
      string treesdirname("./");
      TFile* f1 = 0;
      TTree* seleventtree = 0;
      if(anaEnv.runonTTrees){
        filename = treesdirname + datasets[d]->Name() + "_tree.root"; //previously: datasets[d]->Title() !!
        f1 = new TFile (filename.c_str (), "READ");
      	seleventtree = (TTree*) f1->Get("OBS");
      }
      //float met = 0;
      //seleventtree->SetBranchAddress("MET",&met);
      //
     
      if(!anaEnv.runonTTrees){
        cout<<"Running on toptrees! WARNING: nothing is done."<<endl;
      	//pseudoExp.CalculatePseudoExp(datasets[d], Luminosity, ievtmin[d], manyPseudoExp);
      }
      else {
          cout<<"Running on TTrees!"<<endl;		
	  float eqlumi = datasets[d]->EquivalentLumi();
	  if(anaEnv.nPseudoExp>0 && datasets[d]->Name() != "Data"){
	   if(anaEnv.Systematics == 2){
	      TRandom3 randomgeneratorGaus(0);		
	      float LuminositySmeared = randomgeneratorGaus.Gaus(Luminosity,anaEnv.LuminosityError * Luminosity); //smearing lumi... warning: maybe lumi syst not needed when using data-driven norm factors??  
	      pseudoExp.CalculatePseudoExp(seleventtree, LuminositySmeared, datasets[d], anaEnv); //smearing of cross section inside function...
	   }
	   else pseudoExp.CalculatePseudoExp(seleventtree, Luminosity, datasets[d], anaEnv);
          }
	  else{
	   //pseudoExp.CalculatePseudoExp(seleventtree, Luminosity, eqlumi); //maybe not even needed, because you clear the indices when not (in principle)
          }
      }



   //   PlotObservables Plots (obs.ListOfVariables (), obs.RangeVariables ()); //keep this for the moment 
      //TTreeObservables TtreeObs(obs.ListOfVariables ());
    
      vector<int> indices;
      if ((anaEnv.isMC==1 && anaEnv.nPseudoExp<1 && anaEnv.runonTTrees) || (datasets[d]->Name()=="Data" && anaEnv.runonTTrees) || (datasets[d]->Name()=="Data1" && anaEnv.runonTTrees) || (datasets[d]->Name()=="Data2" && anaEnv.runonTTrees)){ //anaEnv.runonTTrees just to be sure
        indices.clear();
     	nEvents[d] = 0;
   
    	if (verbose > 1) cout << "	Loop over events " << endl;

    	int  nevt = seleventtree->GetEntries();

  	//cout<<" --------------------------------------------->>>>>>>>>>>>>>>>Will loop for "<<datasets[d]->Name()<<" for   "<<nevt<<"  events...... "<<endl;
    
    	for (int ii=0; ii<nevt; ii++){
        	indices.push_back(ii);
        }
      }
      else
	indices = pseudoExp.GetRandomNumbers();
         
  int n = 0;
     /////////////////////////////////////////////loop on events
      nEvents[d] = 0;
      //nEvtsControlSample[d] = 0;
      cout << "	Loop over events, indices.size() =  "<< indices.size() << endl;
   //cout<<" --------------------------------------------->>>>>>>>>>>>>>>>Will loop for "<<datasets[d]->Name()<<" for   "<<ievtmin[d]<<"  events...... "<<endl;
      float sumscaleFactors = 0;
      int dataeventcount = 0;
      int Nevents_supposedtobeinHTbin12 = 0;
      for (vector<int>::iterator iter=indices.begin(); iter!=indices.end(); iter++) {
        int ievt = *iter;
	nEvents[d]++;

        if (anaEnv.nPseudoExp!=0 && ievt % 5000 == 0)
          cout << "Processing the " << ievt << "th event" << " from "<<indices.size()<<" for Pseudo "<<iPseudoExp<<flush << "\r";
       //  if (ievt %1 ==0)
       //  cout<< "Starting from the "<<ievt<<"th event"<<" from "<<indices.size()<<flush<<"\r";
        if (anaEnv.nPseudoExp==0 && ievt % 5000 == 0)
	  cout << "Processing the " << ievt << "th event" << " from "<<indices.size()<<flush << "\r";

	//vectors needed to store variables
        vector < float > ValuesInSR;
	vector < pair  < string , float > > ValuesInSR_WithString;
	ValuesInSR_WithString.clear();
	ValuesInSR.clear();	
	

	//load event... NOTE: different when running on the toptrees, and on TTrees of variables of selected events
        if(!anaEnv.runonTTrees){	      
	   cout<<"Running on toptrees! WARNING: nothing is done."<<endl;
	}
	
	//when running on TTrees rather than toptrees...

        /* //to be place in runonTTrees...?
	        
	*/
        else if(anaEnv.runonTTrees){		  
	  ///////////////////////
	  //LOOP OVER VARIABLES
	  int nvar = lstVar.size ();
	  float variable[nvar]; //TEMPORARY?
	  string variable_string;
	  bool HTcut = true; //'temporary' on true, as if this was not present, REMEMBER THAT YOU HAVE PUT THIS HERE!!!!
	  for (unsigned int v = 0; v < lstVar.size (); v++) {	     
	          variable_string = lstVar[v].c_str();
		  seleventtree->SetBranchAddress(variable_string.c_str(),&variable[v]); //is this the way to do it?	     
	  }
	  //piece for the scaleFactor of the event...
	  float scaleFactor = 1.;
	  seleventtree->SetBranchAddress("scaleFactor",&scaleFactor);	  
	  
	  seleventtree->GetEntry(ievt);
	  fillweight = oldfillweight*scaleFactor;
	  sumscaleFactors = sumscaleFactors + scaleFactor;
	
		float xbins[14]; //this one obtained from HT4jetsMuonMET, 12 bins, TTJets semi mu flat	  
	       //12 HT bins
 
 		/*xbins[0] = 0;
		xbins[1] = 150;
 		xbins[2] = 312.441;
 		xbins[3] = 339.367;
 		xbins[4] = 361.02;
 		xbins[5] = 383.141;
 		xbins[6] = 407.255;
 		xbins[7] = 431.133;
 		xbins[8] = 458.443;
 		xbins[9] = 492.924;
 		xbins[10] = 530.696;
 		xbins[11] = 588.808;
 		xbins[12] = 677.971;
 		xbins[13] = 767.134;*/
 
		xbins[0] = 0;  //left edge of 1st bin (but is dummy bin, right?)
 		xbins[1] = 150;
 		xbins[2] = 332.682;
 		xbins[3] = 359.467;
		xbins[4] = 383.009;
 		xbins[5] = 405.223;
 		xbins[6] = 428.57;
 		xbins[7] = 455.44;
 		xbins[8] = 484.705;
 		xbins[9] = 519.487;
 		xbins[10] = 562.545;
 		xbins[11] = 618.977;
 		xbins[12] = 716.949;
 		xbins[13] = 814.921;  //right edge of 12th bin
		
 //		cout<<" xinbs: ";
//		for(int jj=0;jj<axis->GetNbins ()+1;jj++) cout<<xbins[jj]<<",  ";
//		cout<<endl;
	  
	  int currentHTbin = 0;
	  //this is for the '2D' binning 
 	  if((doWithSpecial2Dbinning && anaEnv.isMC==1 && anaEnv.MCRound==0) || (doWithSpecial2Dbinning && anaEnv.isMC==1 && anaEnv.MCRound==1)){	
	  	HTcut = false;
//		TAxis *axis = mapAxis["HT4jetsMuonMET"]; //hardcoded...
//		double* xbins = axis->GetXbins ()->fArray; //this one obtained from HT4jetsMuonMET, 8 bins, TTJets semi mu flat
	        
		bool inthisbin = false; 
		
		for (unsigned int v = 0; v < lstVar.size (); v++) { //really have to loop over all variables for this first... REMEMBER THAT YOU HAVE PUT THIS HERE!!!!	  	   
		   inthisbin = false;
		   if(thisROUNDdir == "./ROUND19/HTbin1/") inthisbin = (variable[v]>=xbins[1] && variable[v]<xbins[2]);
		   if(thisROUNDdir == "./ROUND19/HTbin2/") inthisbin = (variable[v]>=xbins[2] && variable[v]<xbins[3]);
		   if(thisROUNDdir == "./ROUND19/HTbin3/") inthisbin = (variable[v]>=xbins[3] && variable[v]<xbins[4]);
		   if(thisROUNDdir == "./ROUND19/HTbin4/") inthisbin = (variable[v]>=xbins[4] && variable[v]<xbins[5]);
		   if(thisROUNDdir == "./ROUND19/HTbin5/") inthisbin = (variable[v]>=xbins[5] && variable[v]<xbins[6]);
		   if(thisROUNDdir == "./ROUND19/HTbin6/") inthisbin = (variable[v]>=xbins[6] && variable[v]<xbins[7]);
		   if(thisROUNDdir == "./ROUND19/HTbin7/") inthisbin = (variable[v]>=xbins[7] && variable[v]<xbins[8]);
		   if(thisROUNDdir == "./ROUND19/HTbin8/") inthisbin = (variable[v]>=xbins[8] && variable[v]<xbins[9]);
		   if(thisROUNDdir == "./ROUND19/HTbin9/") inthisbin = (variable[v]>=xbins[9] && variable[v]<xbins[10]);
		   if(thisROUNDdir == "./ROUND19/HTbin10/") inthisbin = (variable[v]>=xbins[10] && variable[v]<xbins[11]);
		   if(thisROUNDdir == "./ROUND19/HTbin11/") inthisbin = (variable[v]>=xbins[11] && variable[v]<xbins[12]);
		   if(thisROUNDdir == "./ROUND19/HTbin12/") inthisbin = (variable[v]>=xbins[12] && variable[v]<xbins[13]+1000000);
		   //if(lstVar[v]=="HT4jetsMuonMET" && variable[v]>=xbins[8] && variable[v]<xbins[9]+10000000) HTcut = true;
	//  	   if(lstVar[v]=="HT4jetsMuonMET" && variable[v]>=xbins[1] && variable[v]<xbins[2]) HTcut = true;
//	if(lstVar[v]=="HT4jetsMuonMET") cout<<" variable = "<<variable[v]<<endl;
	           if(lstVar[v]=="HT4jetsMuonMET" && inthisbin){//inthisbin
		     HTcut = true;
		     n++;
		     //cout<<"  HTcut = true"<<endl;
		   }

		}
		//if(datasets[d]->Name()=="Tprime350") cout<<"  n = "<<n<<endl;
		//if(datasets[d]->Name()=="TTbarJets_SemiMuon_Analysis") cout<<"  n = "<<n<<endl;
	  }
	  else if (doWithSpecial2Dbinning && anaEnv.isMC==0){
	    //data
	    for (unsigned int v = 0; v < lstVar.size (); v++) {
	      if(lstVar[v]=="HT4jetsMuonMET"){
	      //SOME HARDCODING!!!
//	        cout<<"      variable[v] = "<<variable[v]<<endl;
		for(int xbinsi=1;xbinsi<14;xbinsi++){
//		  cout<<"		xbins["<<xbinsi<<"] = "<<xbins[xbinsi]<<",  xbins["<<xbinsi+1<<"] = "<<xbins[xbinsi+1]<<endl;
		  if(variable[v]>=xbins[xbinsi] && variable[v]<xbins[xbinsi+1]){
		  	currentHTbin = xbinsi;
			//HTcut = true; //well, you know, it should already be true
			break;
		  }
		  if(variable[v]>xbins[13]){
		  	currentHTbin = 12;
			//HTcut = true; //well, you know, it should already be true
			break;
		  }
		}
//		cout<<"         -> currentHTbin = "<<currentHTbin<<endl;	      
	      }
	    }	  
	  }
	  
	  //for the test of the plots with a HTcut > 400 GeV, tbc
	  /*for (unsigned int v = 0; v < lstVar.size (); v++) {
	      if(lstVar[v]=="HT4jetsMuonMET" && variable[v]>400){		  
			HTcut = true; //well, you know, it should already be true
			break;
//		cout<<"         -> currentHTbin = "<<currentHTbin<<endl;	      
	      }
	  }*/
	    
	    
	    
	  //if(METcut){ // I switched off the MET cut (temporary?)	  
	  if(HTcut){
	   if(!do2Dhisto){
	    for (unsigned int v = 0; v < lstVar.size (); v++) {
	      //cout<<"lstVar["<<v<<"] = "<<lstVar[v]<<endl;
	      variable_string = lstVar[v];
              BinContentObsPerSample_.push_back( pair< string,float > (lstVar[v],variable[v]));	  	    
	      ValuesInSR_WithString.push_back(pair < string, float > (variable_string,variable[v]));	  
	      hData[v]->Fill (variable[v], fillweight);
	      int N = hData[v]->GetNbinsX();
	      hData[v]->SetBinContent(N,hData[v]->GetBinContent(N)+hData[v]->GetBinContent(N+1));// same as in MultiSamplePlot
	      hData[v]->SetBinContent(N+1,0);// same as in MultiSamplePlot
	     
	      //if(lstVar[v]=="MET" && variable[v]<40){
	      //	nEvtsControlSample_MET[d]++;
	      //}
	      if(lstVar[v]=="ET1oET4" && variable[v]<1){
	     	cout<<"ET1oET4 = "<<variable[v]<<", not possible???"<<endl;
	      }
	      if(lstVar[v]=="HT4jetsMuonMET" && variable[v]>716.949){
	         Nevents_supposedtobeinHTbin12++;
	      }
	     
	      //when you want to (re)do the MC with the selected events in the TTree
	      if(anaEnv.isMC && anaEnv.nPseudoExp<1){
		  //cout<<"     datasets["<<d<<"]->Name() = "<<datasets[d]->Name()<<endl;
		  //cout<<"     variable["<<v<<"] ="<<variable[v]<<endl;
		  //cout<<"     fillweight ="<<fillweight<<endl;
		     MCObsExp[v]->Fill (datasets[d], variable[v], fillweight);
	      }
	      if(doWithSpecial2Dbinning && anaEnv.isMC==0){
	      	  //for data with special 2D binning
		  if(lstVar[v]=="MassHadTop"){
		    dataeventcount++;
//		    cout<<"      currentHTbin = "<<currentHTbin<<", variable[v] = "<<variable[v]<<", fillweight = "<<fillweight<<endl;
///		    cout<<"       "<<DataHistos_Special2Dbinning[currentHTbin]->GetName()<<endl;
		    DataHistos_Special2Dbinning[currentHTbin]->Fill(variable[v], fillweight); //fillweight should be 1 for data...	          
		    int Nbins = DataHistos_Special2Dbinning[currentHTbin]->GetNbinsX();
		    //cout<<"	Nbins = "<<Nbins<<endl;
		    //cout<<"	DataHistos_Special2Dbinning[currentHTbin]->GetBinContent(Nbins) = "<<DataHistos_Special2Dbinning[currentHTbin]->GetBinContent(Nbins)<<", DataHistos_Special2Dbinning[currentHTbin]->GetBinContent(Nbins+1) = "<<DataHistos_Special2Dbinning[currentHTbin]->GetBinContent(Nbins+1)<<endl;
		    DataHistos_Special2Dbinning[currentHTbin]->SetBinContent(Nbins,DataHistos_Special2Dbinning[currentHTbin]->GetBinContent(Nbins)+DataHistos_Special2Dbinning[currentHTbin]->GetBinContent(Nbins+1));
		    DataHistos_Special2Dbinning[currentHTbin]->SetBinContent(Nbins+1,0);
		  }
	      }
	    }/////////////end of looping for variables
	   }
	   else{
	      //when you want to (re)do the MC with the selected events in the TTree
	      if(anaEnv.isMC && anaEnv.nPseudoExp<1){
		     MCObsExp[0]->Fill (datasets[d], variable[0], variable[1], fillweight);
	      }
	   }
	   VarValues_WithString.push_back(ValuesInSR_WithString);
	   isNP = 0;
	   if(dataSetName.find("prime")<=dataSetName.size()) isNP = 1; 
	   NPbooleans.push_back(isNP);
	   nofSelEvts++;	   
	  }  //////end of HTcut	  
	}/////////////////////end of runonTTrees
	
      }	////////////loop on 'indices' ('events')
      cout<<" "<<endl;
      cout<<" dataeventcount = "<<dataeventcount<<endl;
      cout<<" sumscaleFactors = "<<sumscaleFactors<<endl;
      cout<<" XS/PreSelEff/nEvents : " << datasets[d]->Xsection () << "/" << datasets[d]->PreSelEfficiency () << "/" << nEvents[d] << endl;
      cout<<"   VarValues.size() = "<<VarValues.size()<<endl;
      cout<<"   NPbooleans.size() = "<<NPbooleans.size()<<" and total BinContentObs size for this dataset  "<<BinContentObsPerSample_.size()<<endl;
      cout<<" <---> Nevents_supposedtobeinHTbin12 = "<<Nevents_supposedtobeinHTbin12<<endl;
      string dsTitle = datasets[d]->Title ();
    		
      sort( BinContentObsPerSample_.begin(), BinContentObsPerSample_.end());
	
      string bin_decision = "Bins";
     /* if(d>0 && anaEnv.nPseudoExp==0){   
        cout<<"  Will try to create and fill .... "<<dsTitle<<"   "<<BinContentObsPerSample_.size()<<endl;
        TString dss= datasets[d]->Title();
        if (anaEnv.MCRound==1 && bin_decision == "Bins")   Plots.WriteNBins(dss,BinContentObsPerSample_,anaEnv.nbins,fillweight,true);
        if (anaEnv.MCRound==1 && bin_decision == "Events")     Plots.WriteEvtsPerBin(dss,BinContentObsPerSample_,anaEnv.eventsperbin,fillweight,true);
      }
     */
      BinContentObsAll_ = BinContentObsPerSample_; //first copy to another vector
      BinContentObsPerSample_.clear(); //and then clear, because will be needed empty for the next dataset

      indices.clear();
      cout<<"   DONE FOR DATASET........................ "<<dsTitle<<endl; 
      nofSelEvts_vector.push_back(nofSelEvts);
      cout<<"   real number of selected events of this dataset: "<<nofSelEvts_vector[d]<<endl<<endl;
      cout<<"   nEvents["<<d<<"] = "<<nEvents[d]<<endl;
      //if(f1->IsOpen()) cout<<"!!! f1 is open !!!"<<endl;
      f1->Close();     
    } //////loop on datasets

    BinContentObsPerSample_.clear();  
    
    
    //regarding binning   
    if (anaEnv.isMC && anaEnv.MCRound==0){
      string bin_decision = "Bins";	
      cout<<" ... Will try to create the Binning according to the decision \""<<bin_decision<<"\" ..."<<endl;	
      MakeBinning NewBins;
      cout<<"BinContentObsAll_.size() = "<<BinContentObsAll_.size()<<endl;
      NewBins.Binning(bin_decision,BinContentObsAll_,anaEnv.nbins,EventsPerBin,fillweight_TTJets);
      cout<<"HERE 1"<<endl;
    }    
    
   //Once everything is filled ...
    if (verbose > 0)
      cout << " We ran over all the data ;-)" << endl;     
     
   string lstChosenVar[lstVar.size ()]; //needed for the variable names in a later loop over ObsData; only related to the subdirectory structure of the Stat output file
   int kk =0;
  //LOOP OVER VARIABLES
   cout<<"lstVar.size () = "<<lstVar.size ()<<endl;
   for (unsigned int v = 0; v < lstVar.size (); v++) {
      /*if (!anaEnv.runOnObs ((int) v))
	continue;
      */
      
      /////////////////////////////////
      //  Summary plot
      /////////////////////////////////
      if (verbose > 0)
	cout << " - Summary plots ..." << endl;
/*
      cout<<"MCObsExp.size() = "<<MCObsExp.size()<<endl;
      for (int l=0;l<MCObsExp.size();l++){
        //TH1F* dummyhisto = 0;
	//dummyhisto = MCObsExp[l]->GetHistoSMProcesses();
	//cout<<"HERE1"<<endl;
	//cout<<"MCObsExp["<<l<<"]->GetHistoSMProcesses()->GetName() = "<<MCObsExp[l]->GetHistoSMProcesses()->GetName()<<endl;
	if (MCObsExp[l]->GetHistoSMProcesses()->GetName()=="SMProcess_"+lstVar[v]){
	   cout<<"  found MCObsExp plots....... "<<MCObsExp[l]->GetHistoSMProcesses()->GetName()<<"   "<<l<<endl;
	   TH1F* htemp = MCObsExp[l]->GetHistoSMProcesses(); //added this
	   ObsEstimated.push_back(htemp); //actually not really used anymore... +- everything goes directly from MCObsExp
	   break;
	}
      }
      //      ObsEstimated.push_back(MCObsExp[v]->GetHistoSMProcesses());
*/
      
      ObsData.push_back (hData[v]);
      lstChosenVar[kk] = lstVar[v];
    
      kk++;
      /////////////////////////////////////////
   }	//end of loop on variables
   
   //writing observable plots in Stat Procedure output file
   if (verbose > 0)
      cout << " - Start writing observable plots in Stat Procedure output file ..." << endl;
   if(anaEnv.nPseudoExp > 0){
      foutStat->cd();
      cout<<"ObsData.size() = "<<ObsData.size()<<endl;
      for(unsigned int i=0;i<ObsData.size();i++){
	 mytdirall->cd(lstChosenVar[i].c_str());
	 TH1F hdatatemp;
	 TH1F hestimtemp;
	 hdatatemp = *ObsData[i];
	 
	 cout<<" classname: "<<hdatatemp.ClassName()<<", nbins "<<hdatatemp.GetNbinsX()<<", integral "<<hdatatemp.Integral()<<endl;	 
	 hdatatemp.Write();  
	 
/*	 hestimtemp = *ObsEstimated[i];
	 hestimtemp.Write(); // always the same when MC expectation as estimation, but always different when data-driven estimation...
*/       
         /*TH1F hestimtemp_StackSM; // SM processes stacked
         hestimtemp_StackSM = MCObsExp[i]->GetTHStackSMProcesses();
	 hestimtemp_StackSM.Write();
	 TH1F hestimtemp_hW_LF; // W_LF
	 hestimtemp_hW_LF = *(MCObsExp[i]->Get_hW_LF());
	 hestimtemp_hW_LF.Write();*/
       }
   }

   if(doDataDummy){
     foutStat->Close(); //only when you want to obtain the data observable distributions, while scanning some dummy small plane...
     exit(1);
   }

  ///////////////////
  // Writing
  //////////////////
   if (verbose > 0)
      cout << " - Start writing the module outputs in the output file ..." << endl;
   char nom[200] = "";
   sprintf (nom, "_exp%d_", iPseudoExp);

   //loop on variables
   if(!do2Dhisto){
     int it = 0;
     for (unsigned int v = 0; v < lstVar.size (); v++) {
      //if (anaEnv.runOnObs ((int) v)) {
	 sprintf (nom, "%s_exp%d_", lstVar[v].c_str (), iPseudoExp);
	 cout<<" v = "<<v<<endl;
	 if (anaEnv.isMC && anaEnv.nPseudoExp<1) { //first: only anaEnv.isMC 
	   MCObsExp[v]->SetColors (datasets);
	   MCObsExp[v]->Compute ();
	   new ((*tcMCObsExp)[it]) MCObsExpectation (*MCObsExp[v]);
	 }
	 it++;
      //}
     }
   }
   else{
     new ((*tcMCObsExp)[0]) MCObsExpectation (*MCObsExp[0]); //in the 2D case, there is only 1 MCObsExp...
   }
// } //end 'if not Poission pseudo's, but real pseudo's'
  ///////////////////////////////
  ///  STATISTICAL PROCEDURE
  //////////////////////////////
  
   if(anaEnv.nPseudoExp > 0 && !doWithSpecial2Dbinning){
     if (verbose > 0)
        cout << " - Statistical procedure ..." << endl;
     vector<float> Svect;
     SampleCombinedWeightCalculator mySampleCombinedWeightCalculator;////// For the V-value
     if(verbose>0) 
        cout << endl << "Starting calculation of weights..." << endl;
     // LOOP OVER SELECTED EVENTS

     int nVar = VarValues_WithString.size();
     //cout<<"lstVar size "<<lstVar.size()<<" VarValues_WithString.size() = "<<VarValues_WithString.size()<<" VarValues "<<VarValues.size()<<"   obsData  "<<ObsData.size()<<" obsEstimated  "<<ObsEstimated.size()<<endl;
//     for (int k=0;k<ObsData.size();++k){
//  	cout<<ObsData[k]->GetName()<<"  "<<ObsEstimated[k]->GetName()<<"   "<<k<<endl;
//     }
     VarValues.clear();

     vector <float> dummy; //working with this is a bit cumbersome, but it is a temporary fix    

     sort(VarValues_WithString.begin(),VarValues_WithString.end()); //wait, why is this??

     for(unsigned int i=0;i<VarValues_WithString.size();i++){
//	cout<<"VarValues_WithString["<<i<<"].size() = "<<VarValues_WithString[i].size()<<endl;
	for ( int ii=0;ii<VarValues_WithString[i].size();ii++){	
////////cout<<"  Varvalues......"<<i<<"   "<<VarValues_WithString[i][ii].second<<endl;	   
	    dummy.push_back(VarValues_WithString[i][ii].second);
        }
        VarValues.push_back(dummy);
	dummy.clear();
     }
     for(unsigned int i=0;i<VarValues_WithString.size();i++){
	/*for(unsigned int j=0; j<nVar;j++){ 
		cout << " VarValues["<<i<<"]["<<j<<"]: "<< VarValues[i][j];
	}*/
  //cout<<"lstVar size "<<lstVar.size()<<" VarValues_WithString.size() = "<<VarValues_WithString.size()<<" VarValues "<<VarValues.size()<<"   obsData  "<<ObsData.size()<<" obsEstimated  "<<ObsEstimated.size()<<endl;
        EventCombinedWeightCalculator myEventCombinedWeightCalculator(SignificanceCombination);////// For the squared - significance and S-product
   	myEventCombinedWeightCalculator.SetObsExp(MCObsExp,MCObsExp_NP);
	//if(anaEnv.Systematics == 0 || anaEnv.Systematics == 2) myEventCombinedWeightCalculator.CalculateWeight(VarValues_WithString[i],ObsData,ObsEstimated);
	float LumiError = anaEnv.LuminosityError;	
	bool doLumiSyst = false, doXsectionSyst = false, doJESSyst = false; //this is configuration if anaEnv.Systematics == 0
	if(anaEnv.Systematics == 1){
	       doLumiSyst = true; doXsectionSyst = true; doJESSyst = true;
	}
	if(anaEnv.Systematics == 2){
	       doLumiSyst = false; doXsectionSyst = false; doJESSyst = true;
	}
//	myEventCombinedWeightCalculator.SetSystematics(LumiError, MCObsExp_JESSyst, doLumiSyst, doXsectionSyst, doJESSyst); //OUTDATED anyway...
	myEventCombinedWeightCalculator.CalculateWeight(VarValues_WithString[i],ObsData,datasets_allconfigMC); //not 'datasets'...
	if(!myEventCombinedWeightCalculator.isUnderflowWarning()){
	    Svect.push_back(myEventCombinedWeightCalculator.GetCombinedWeight());
	    if(verbose>3) cout << "Weight event " << i << ": S = " << Svect[i] << endl;
	}
	else{
	    cout<<"WARNING: myEventCombinedWeightCalculator.isUnderflowWarning() = true. This problem has to be solved!"<<endl;
	    Svect.push_back(-9999);  // or -999? By setting it to 0 or a negative number, these events will not end up in the highest-weight subsample 
	}
     }


     if(Svect.size()!=NPbooleans.size()) cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                  ---------------->>>>>>>>>>>> WARNING: vectors Svect and NPbooleans not equal size!"<<endl;
   
     vector< eventinfo > Svect_eventsinfos; //the eventinfo struct is defined in SampleCombinedWeightCalculator.h
     eventinfo w;

     cout<<"Svect.size() = "<<Svect.size()<<endl;
     for(unsigned int i=0;i<Svect.size();i++){
////   	cout<<" Svect["<<i<<"] = "<<Svect[i]<<endl;
	w.weight = Svect[i];
	w.isnp = NPbooleans[i];
	//cout<<" VarValues["<<i<<"].size() = "<<VarValues[i].size()<<endl;
	w.variableValues = &VarValues[i];
	Svect_eventsinfos.push_back(w);
     }
     
     foutStat->cd();
     for(unsigned int k=0;k<nofFractionsHWEvts;k++){
        mySampleCombinedWeightCalculator.CalculateV(Svect_eventsinfos, Combination, FractionHWEvts[k]);   // after this step, Svect_eventsinfos is sorted (according to the weigths), and the V value is calculated
        ostringstream Fractionstrstream;
        Fractionstrstream << FractionHWEvts[k];
        mytdirFraction = foutStat->GetDirectory("x="+TString(Fractionstrstream.str()));
        mytdirFraction->cd();
	
        //writing highest-weight events histograms			
        for(unsigned int i=0;i<ObsData.size();i++){  // why does this in an if(ips==0) loop give no such histograms written in the file	
           TH1F *hdatatemp=new TH1F(ObsData[i]->GetName()+TString("_highestSweightsEvts_fr")+TString(Fractionstrstream.str()), ObsData[i]->GetTitle()+TString(" (highest-S-weight subsample [x = ")+TString(Fractionstrstream.str())+TString("])"), ObsData[i]->GetNbinsX(), ObsData[i]->GetXaxis()->GetXbins ()->fArray);
	   hdatatemp->GetYaxis()->SetTitle("#evts");
	  //cout<<"mySampleCombinedWeightCalculator.GetnEventsInHighestWeightSubsample() = "<<mySampleCombinedWeightCalculator.GetnEventsInHighestWeightSubsample()<<endl;
	   for(unsigned int j=0; j<mySampleCombinedWeightCalculator.GetnEventsInHighestWeightSubsample();j++){
		hdatatemp->Fill((*mySampleCombinedWeightCalculator.GetEventsInfos()[j].variableValues)[i]);
	   }
	   //cout<< "NAME OF CURRENT DIRECTORY: " << gDirectory->GetName() <<endl;           
	   int N = hdatatemp->GetNbinsX();
	   hdatatemp->SetBinContent(N,hdatatemp->GetBinContent(N)+hdatatemp->GetBinContent(N+1));// same as in MultiSamplePlot
	   hdatatemp->SetBinContent(N+1,0);// same as in MultiSamplePlot
	   hdatatemp->SetMarkerStyle(20);
	   hdatatemp->Sumw2(); //I don't do calculations with these histo's (doesn't matter), only for plotting (and by calling this and the the line above, you should get data points with error bars)
	   mytdirFraction->cd(lstChosenVar[i].c_str());
	   hdatatemp->Write();
	   delete hdatatemp;
        }
	TH1F *hSHWevts = new TH1F(TString("highestSweights_fr")+TString(Fractionstrstream.str()), TString("highest-S-weights [x = ")+TString(Fractionstrstream.str())+TString("]"),5000,0,50); //you can rebin/zoom in later
	hSHWevts->GetXaxis()->SetTitle("S");
	hSHWevts->GetYaxis()->SetTitle("#evts");
	for(unsigned int j=0; j<mySampleCombinedWeightCalculator.GetnEventsInHighestWeightSubsample();j++){
	   hSHWevts->Fill(mySampleCombinedWeightCalculator.GetEventsInfos()[j].weight);	   	   
        }
	mytdirFraction->cd();
	hSHWevts->Write();
	delete hSHWevts;//still have to make sure this does not make it crash!
        if(verbose>0) cout<<"FractionHWEvts["<<k<<"]: "<<FractionHWEvts[k]<<endl;
        if(verbose>0) cout << "  V = " << mySampleCombinedWeightCalculator.GetV() << endl << endl;
        Vvect.push_back(pair< float , float >(FractionHWEvts[k],mySampleCombinedWeightCalculator.GetV()));
        float fracNP = -9999, signNP = -9999;
        fracNP = mySampleCombinedWeightCalculator.GetFracNPEventsInHWSubsample();
	signNP = mySampleCombinedWeightCalculator.GetSignalSignificanceInHWSubSample();
        FracNPEventsInHWSubsampleVect.push_back(pair< float , float >(FractionHWEvts[k],fracNP));
	SignNPEventsInHWSubsampleVect.push_back(pair< float , float >(FractionHWEvts[k],signNP));  
     }  //end loop on fraction
 
     if(verbose>0) cout<<"   Number of events in final sample = "<<mySampleCombinedWeightCalculator.GetnEventsInFinalSample()<<endl; //final sample = sample in SR
     if(verbose>0) cout<<"   Number of NP events in final sample = "<<mySampleCombinedWeightCalculator.GetnNPEventsInFinalSample()<<endl;
     
    // loop to free memory allocated to hData //worked (and was necessary) in 22X version, should be checked in 35X version
     if(verbose>0) cout<<"Releasing memory of hData..."<<endl; 
     for (unsigned int v = 0; v < lstVar.size();v++){//obs.ListOfVariables ().size (); v++) {
        //if(!anaEnv.runOnObs((int)v)) continue; 
        if(hData[v]) delete hData[v];  // if(hData[v]!=NULL)?
     }
   } //end if(anaEnv.nPseudoExp > 0) statement
   
     
   delete MCExp;
   MCExp = NULL;
    

   ObsEstimated.clear();
   ObsData.clear();
   cout << "** End pseudoexperiment  **" << endl << endl;
   
   //Q value for data
   if(!anaEnv.isMC && likelihood && doWithSpecial2Dbinning){
     /*for(int HTbin=1;HTbin<13;HTbin++){
      for(int MassHadTopbin=2;MassHadTopbin<14;MassHadTopbin++){
        cout<<"	 ";
	cout<<"       DataHistos_Special2Dbinning["<<HTbin<<"]->GetBinContent("<<MassHadTopbin<<") = "<<DataHistos_Special2Dbinning[HTbin]->GetBinContent(MassHadTopbin)<<endl;
      }     
     }*/
      EventCombinedWeightCalculator myEventCombinedWeightCalculator_Data(SignificanceCombination);
      myEventCombinedWeightCalculator_Data.SetObsExp_towards2D(MCObsExp_bins,MCObsExp_NP_bins);
      //if(doJESsystematics) myEventCombinedWeightCalculator_Data.SetObsExp_Systematics_towards2D(MCObsExp_JESminus_bins, MCObsExp_JESplus_bins, MCObsExp_NP_JESminus_bins, MCObsExp_NP_JESplus_bins); //no, not for data...
      myEventCombinedWeightCalculator_Data.SetDataHistos_Special2Dbinning(DataHistos_Special2Dbinning);
      cout<<" Q = "<<myEventCombinedWeightCalculator_Data.GetQValue_PoissonPseudos_towards2D("MassHadTop",doSMPoissonPseudos)<<endl;
      Vvect.push_back(pair<float,float>(1.,myEventCombinedWeightCalculator_Data.GetQValue_PoissonPseudos_towards2D("MassHadTop",doSMPoissonPseudos))); //is in this case not pseudos...
      
      //if data, write the histograms per HTbin in a file...
      if(firstiteration){
        cout<<" First iteration for data (in the loop over the plane) --> writing file with histos"<<endl;
        TFile * fData = new TFile("./ROUND19/plotsData_905ipb_HTbins.root","RECREATE");
        for(int HTbini=1;HTbini<nHTbins+1;HTbini++){
           stringstream sstr_HTbini;
           sstr_HTbini << HTbini;
           string label = "Data_MassHadTop_HTbin" + sstr_HTbini.str();
           DataHistos_Special2Dbinning[HTbini]->GetName();
           cout<<"     ...writing histo for HTbin "<<HTbini<<endl  ; 
           DataHistos_Special2Dbinning[HTbini]->Write(label.c_str());       
        }
        firstiteration = false;
        cout<<" Closing file..."<<endl;
        //fData->Close();
      }
   }
   //if data, write the histograms per HTbin in a file...
   
 }// loop over pseudo-experiments
}  //end if(!doPoissonPseudos)
else if(doPoissonPseudos){
   if(!do2Dhisto && !doWithSpecial2Dbinning){ //Vvect_poisson.clear();
     for(int ppseudo=0;ppseudo<nPoissonPseudos;ppseudo++){
      EventCombinedWeightCalculator myEventCombinedWeightCalculator_poisson(SignificanceCombination);
      myEventCombinedWeightCalculator_poisson.SetObsExp(MCObsExp,MCObsExp_NP);    
      if(!likelihood) Vvect.push_back(pair<float,float>(1.,myEventCombinedWeightCalculator_poisson.GetVValue_PoissonPseudos(lstVar[0],doSMPoissonPseudos)));
      else if(likelihood) Vvect.push_back(pair<float,float>(1.,myEventCombinedWeightCalculator_poisson.GetQValue_PoissonPseudos(lstVar[0],doSMPoissonPseudos)));
     }     
   }
   else if(!do2Dhisto && doWithSpecial2Dbinning){
     for(int ppseudo=0;ppseudo<nPoissonPseudos;ppseudo++){
      EventCombinedWeightCalculator myEventCombinedWeightCalculator_poisson(SignificanceCombination);
  //    myEventCombinedWeightCalculator_poisson.SetObsExp_towards2D(MCObsExp_bin1,MCObsExp_bin2,MCObsExp_bin3,MCObsExp_bin4,MCObsExp_bin5,MCObsExp_bin6,MCObsExp_bin7,MCObsExp_bin8,MCObsExp_NP_bin1,MCObsExp_NP_bin2,MCObsExp_NP_bin3,MCObsExp_NP_bin4,MCObsExp_NP_bin5,MCObsExp_NP_bin6,MCObsExp_NP_bin7,MCObsExp_NP_bin8);    
      myEventCombinedWeightCalculator_poisson.SetObsExp_towards2D(MCObsExp_bins,MCObsExp_NP_bins);
      if(doJESsystematics) myEventCombinedWeightCalculator_poisson.SetObsExp_Systematics_towards2D(MCObsExp_JESminus_bins, MCObsExp_JESplus_bins, MCObsExp_NP_JESminus_bins, MCObsExp_NP_JESplus_bins);
      if(doLumisystematics) myEventCombinedWeightCalculator_poisson.SetLumi_Systematics(LumiError);
      //if(doBkgNormsystematics) myEventCombinedWeightCalculator_poisson.SetBkgNorm_Systematics(TTbarError,SingleTopError);
      
      //Mass had top hardcoded!!!! whatever, this code is screwed up anyway...
//      cout<<" * pseudo "<< ppseudo <<" *"<<endl;
      if(likelihood){
         /*if(myEventCombinedWeightCalculator_poisson.GetQValue_PoissonPseudos_towards2D("MassHadTop",doSMPoissonPseudos) == -999999){ 
	     ppseudo--; //redo pseudoexperiment when one of the expected bin contents is negative or 0
             cout<<"ENCOUNTERED NEGATIVE OR ZERO EXPECTED BIN CONTENT ... redoing pseudoexperiment"<<endl;
	 }
	 else*/
	 Vvect.push_back(pair<float,float>(1.,myEventCombinedWeightCalculator_poisson.GetQValue_PoissonPseudos_towards2D("MassHadTop",doSMPoissonPseudos)));
      }
//      cout<<endl;
     } 
   }
   else if(do2Dhisto){
     for(int ppseudo=0;ppseudo<nPoissonPseudos;ppseudo++){
      EventCombinedWeightCalculator myEventCombinedWeightCalculator_poisson(SignificanceCombination);
      myEventCombinedWeightCalculator_poisson.SetObsExp(MCObsExp,MCObsExp_NP);    
      Vvect.push_back(pair<float,float>(1.,myEventCombinedWeightCalculator_poisson.GetVValue_PoissonPseudos_2D(lstVar[0],lstVar[1],doSMPoissonPseudos)));
     }
   }
}



//OUTDATED OR WHAT???
/*//for data
if(!doPoissonPseudos && !anaEnv.isMC && likelihood && doWithSpecial2Dbinning){
   EventCombinedWeightCalculator myEventCombinedWeightCalculator_Data(SignificanceCombination);
   cout<<"   HEY 1 "<<endl;
   myEventCombinedWeightCalculator_Data.SetObsExp_towards2D(MCObsExp_bins,MCObsExp_NP_bins);
   cout<<"   HEY 2 "<<endl;
   //if(doJESsystematics) myEventCombinedWeightCalculator_Data.SetObsExp_Systematics_towards2D(MCObsExp_JESminus_bins, MCObsExp_JESplus_bins, MCObsExp_NP_JESminus_bins, MCObsExp_NP_JESplus_bins); //no, not for data...
   cout<<"                          <---> "<<DataHistos_Special2Dbinning[1]->GetName();
   myEventCombinedWeightCalculator_Data.SetDataHistos_Special2Dbinning(DataHistos_Special2Dbinning);
   cout<<"   HEY 3"<<endl;
   cout<<" Q = "<<myEventCombinedWeightCalculator_Data.GetQValue_PoissonPseudos_towards2D("MassHadTop",doSMPoissonPseudos)<<endl;
   Vvect.push_back(pair<float,float>(1.,myEventCombinedWeightCalculator_Data.GetQValue_PoissonPseudos_towards2D("MassHadTop",doSMPoissonPseudos))); //is in this case not pseudos...
}*/
 
 
 ///////////////////////////////
 ///  STATISTICAL PROCEDURE (to write V values)
 //////////////////////////////
 if(anaEnv.nPseudoExp > 0){
  //in principle calculations related to V_alfa can or even should be (re)done in a seperate macro, because all the V values will be stored
   WeightProbaCalculator myWeightProbaCalculator;
   float alfa= 0.1;
  //writing this vector to a TFile... in TTree format for the moment; can be changed in the future
   cout<<"Writing V vector in TTree format to TFile..."<<endl;  // experimenting with trees to write and read a collection of numbers (the V values)!
   TFile treefile(rootStatVtreeFileName.c_str(),"RECREATE");
   float Vtemp, Valfatemp;
   for(int k=0;k<nofFractionsHWEvts;k++){
     ostringstream Fractionstrstream, alfastrstream;
     Fractionstrstream << FractionHWEvts[k];
     TTree *Vtree = new TTree("Vtree_fr"+TString(Fractionstrstream.str()),"Tree of V values [fr="+TString(Fractionstrstream.str())+"]");
     Vtree->Branch("V",&Vtemp,"Vleaf");
     vector<float> VvectCurrentFraction;
     for(int i=0;i<Vvect.size();i++){	
	if(Vvect[i].first == FractionHWEvts[k]){
/////	   cout<<"Vvect  "<<Vvect[i].first<<"  "<<Vvect[i].second<<"   "<<Vvect.size()<<endl;
	   VvectCurrentFraction.push_back(Vvect[i].second);
	}
     }

     sort(VvectCurrentFraction.begin(),VvectCurrentFraction.end());
	
     for(int i=0;i<VvectCurrentFraction.size();i++){
	//cout<<"   V values output  "<<VvectCurrentFraction[i]<<"  and fraction  "<<(k+1)*0.1<<endl;
	Vtemp = VvectCurrentFraction[i];
	Vtree->Fill();
     }	
     Vtree->Write();
     delete Vtree;
   }
   
   treefile.Close();
 
  //related to filling of histograms of fractions of NP events in highest weight subsample etc...
   cout<<"Updating stat output file"<<endl;
   foutStat->cd();
   for(int k=0;k<nofFractionsHWEvts;k++){
     ostringstream Fractionstrstream;
     Fractionstrstream << FractionHWEvts[k];
     TH1F* FracNPHisto = new TH1F("Fraction of events in highest weight subsample that are NP "+TString(Fractionstrstream.str()),"Fraction of events in highest weight subsample that are NP "+TString(Fractionstrstream.str()),100,0,1);
     FracNPHisto->GetXaxis()->SetTitle("fraction NP events");
     TH1F* SignNPHisto = new TH1F("Signal significance of highest weight subsample "+TString(Fractionstrstream.str()),"Signal significance of highest weight subsample "+TString(Fractionstrstream.str()),100,0,1);
     SignNPHisto->GetXaxis()->SetTitle("s/sqrt(b)");
     for(int i=0;i<FracNPEventsInHWSubsampleVect.size();i++){
	if(FracNPEventsInHWSubsampleVect[i].first == FractionHWEvts[k]){
	   FracNPHisto->Fill(FracNPEventsInHWSubsampleVect[i].second);
	   SignNPHisto->Fill(SignNPEventsInHWSubsampleVect[i].second);
	}
     }
     mytdirFraction = foutStat->GetDirectory("x="+TString(Fractionstrstream.str()));
     mytdirFraction->cd();
     FracNPHisto->Write();
     SignNPHisto->Write();
   }
   foutStat->Close();
 
 } //end if(anaEnv.nPseudoExp > 0) statement
 
 

  ///////////////////
  // Writing
  //////////////////
  if(anaEnv.nPseudoExp > 0){
/*    if (verbose > 0)
      cout << " - Writing configTree in Stat output file ..." << endl;
    //add configuration
    foutStat->cd ();
    cout<<"HEREEEEE 1"<<endl;
    configTree->Fill ();
    cout<<"HEREEEEE 2"<<endl;
    configTree->Write ();
    cout<<"HEREEEEE 3"<<endl;
    foutStat->Close();*/
  }

//only to test... but...??? this (dummy) writing is needed, otherwise there it gives "Error in <TClass::New>: cannot create object of class TH1" when filling configtree... what??? why??? 
  if(do2Dhisto && (anaEnv.nPseudoExp < 1)){
     string Var_x = lstVar[0];
     string Var_y = lstVar[1];

     TH2F* estim_SM_2D=0;    	
//     cout<<"MCObsExp.size() = "<<MCObsExp.size()<<endl;
     for (unsigned int l=0;l<MCObsExp.size();l++){
   	        cout<<"MCObsExp["<<l<<"]->GetHistoSMProcesses_2D()->GetName() = "<<MCObsExp[l]->GetHistoSMProcesses_2D()->GetName()<<endl;
		if (MCObsExp[l]->GetHistoSMProcesses_2D()->GetName()=="SMProcess_"+Var_x+"_vs_"+Var_y){ 
	//	   cout<<endl<<"  found MCObsExp plots....... "<<MCObsExp[l]->GetHistoSMProcesses_2D()->GetName()<<"   "<<l<<endl;           
		   estim_SM_2D = MCObsExp[l]->GetHistoSMProcesses_2D();
		 //  estim_NP_2D = MCNPObsExp_[l]->GetHistoAll_2D();
		   break;
	        }
     }
//     cout<<"estim_SM_2D->GetNbinsX() = "<<estim_SM_2D->GetNbinsX()<<", estim_SM_2D->GetNbinsY() = "<<estim_SM_2D->GetNbinsY()<<endl;     
   cout<<" Bin content bin "<<estim_SM_2D->GetNbinsX()<<","<<2<<": "<<estim_SM_2D->GetBinContent(estim_SM_2D->GetNbinsX(),2)<<endl;
   TFile * foutpseudo = new TFile(TString(thisROUNDdir)+"Pseudoplotmaker.root","RECREATE");
   foutpseudo->cd();
   estim_SM_2D->Write();
   cout<<"Output written in "<< foutpseudo->GetName()<<endl;
   foutpseudo->Close();
  }
  
  
  if (verbose > 0)
    cout << " - Writing outputs on files ..." << endl;
  //add configuration
  fout->cd ();
  cout<<"HERE 1..."<<endl;
  configTree->Fill ();
  cout<<"HERE 2..."<<endl;
  configTree->Write ();
  cout<<"HERE 3..."<<endl;
  //
  if (verbose > 0)
    cout << " - Writing  the file ..." << endl;
  fout->Write ();

  //delete
  delete configTree;
  delete tcdatasets;
  delete tcAnaEnv;
  delete tcMCObsExp;
  fout->Close();
  
  if(!scanTprimespace) break;
  
 }//end superloop over cross sections
 
 if(!scanTprimespace) break;
 
} //end superloop over t' masses

  

  cout << "**********************************************************************" << endl;
  cout << "           End of the program !!" << endl;
  cout << "**********************************************************************" << endl;
  return 0;
}
