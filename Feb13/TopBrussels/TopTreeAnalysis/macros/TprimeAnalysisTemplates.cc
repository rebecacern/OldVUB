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
  cout << "Begining of the new macro for the t' analysis !" << endl;
  cout << "**********************************************************************" << endl;  

  //SetStyle if needed
  setTDRStyle(); 
  //setMyStyle();

  float massestprimeHardCoded[6] = {350,400,450,500,550,600};
  vector<float> massestprime;
  vector<float> crosssectionstprime; //in pb
  vector<float> NofEventsTprimeTopTree; //= {22166,20701,17817,20601};
  
  //xml file
  const char *xmlfile = "../config/myTprimeconfig.xml";
  string thisROUNDdir = "./ROUND24/"; //directory for (most of??) the files of this round... (binning is still in config...) //don't forget the slash at the end

  int doJESShift = 0; // 0: off 1: minus 2: plus.
  int doJERShift = 0; // 0: off (except nominal scalefactor for jer) 1: minus 2: plus
  int dobTagEffShift = 0; //0: off (except nominal scalefactor for btag eff) 1: minus 2:plus
  int domisTagEffShift = 0; //0: off (except nominal scalefactor for mistag eff) 1: minus 2: plus
  string treesdirname,postfix;
  if(doJESShift==0 && doJERShift==0 && dobTagEffShift==0 && domisTagEffShift==0)
  {
    treesdirname = thisROUNDdir + "Trees_Nominal/";
    postfix = "";
    //postfix = "_correctedfixedbinsplots";
  }
  else if(doJESShift==1)
  {
    treesdirname = thisROUNDdir + "./Trees_JESminus/";
    postfix = "_JESminus";
  }
  else if(doJESShift==2)
  {
    treesdirname = thisROUNDdir + "./Trees_JESplus/";
    postfix = "_JESplus";
  }
  else if(doJERShift==1)
  {
    treesdirname = thisROUNDdir + "./Trees_JERminus/";
    postfix = "_JERminus";
  }
  else if(doJERShift==2)
  {
    treesdirname = thisROUNDdir + "./Trees_JERplus/";
    postfix = "_JERplus";
  }
  else if(dobTagEffShift==1)
  {
    treesdirname = thisROUNDdir + "./Trees_bTagEffSFminus/";
    postfix = "_bTagEffSFminus";
  }
  else if(dobTagEffShift==2)
  {
    treesdirname = thisROUNDdir + "./Trees_bTagEffSFplus/";
    postfix = "_bTagEffSFplus";
  }
  else if(domisTagEffShift==1)
  {
    treesdirname = thisROUNDdir + "./Trees_misTagEffSFminus/";
    postfix = "_misTagEffSFminus";
  }
  else if(domisTagEffShift==2)
  {
    treesdirname = thisROUNDdir + "./Trees_misTagEffSFplus/";
    postfix = "_misTagEffSFplus";
  }
  
  
  string xvariable = "HT4jetsMuonMET", yvariable = "MassHadTop"; //these are the two variables for which the 2D plane is made
  string XaxisLabel = "Reconstructed top mass (GeV/c^{2})", YaxisLabel = "#Events";
  int nbinsxvariable = 12, nbinsyvariable = 14;
  int nfixedbins = 80;
  
  AnalysisEnvironment anaEnv;
  cout << "Loading environment ..." << endl;
  AnalysisEnvironmentLoader anaLoad (anaEnv, xmlfile);
  int verbose = anaEnv.Verbose;
  float Luminosity = anaEnv.Luminosity;	// in 1/pb
  //int nbins = anaEnv.nbins;
 
  cout << "The results will be obtained for a luminosity of " << Luminosity << " x 1/pb" <<endl;

  /////////////////////
  // Load Datasets
  /////////////////////
  TTreeLoader treeLoader;
  if (verbose > 0)
    cout << " - Load datasets ..." << endl;
  vector < Dataset* > datasets;
  treeLoader.LoadDatasets (datasets, xmlfile);
  
  vector < Dataset* > datasets_forMSPlot_variables_fixedbins;
  for(unsigned int d=0;d<datasets.size();d++){
    if(!((datasets[d]->Name().find("Tprime")<=datasets[d]->Name().size()) && !(datasets[d]->Name().find("NP_overlay")==0)))
      datasets_forMSPlot_variables_fixedbins.push_back(datasets[d]);
  }
	  
  
  //Output files... all "dummy's" shouldn't be saved in the end (obsolete)... for the moment just ignore.
  string rootFileName = thisROUNDdir + "TprimeAnalysisTemplates" + postfix + ".root";
  TFile* fout = new TFile(rootFileName.c_str(),"RECREATE");
  ostringstream lumi_sstream;
  lumi_sstream << anaEnv.Luminosity;
  map<int,MultiSamplePlot*> MSPlot_xvariableBins;
  map<string,MultiSamplePlot*> MSPlot_variables_fixedbins;
  
  
  vector <string > lstVar;
  vector<pair<float,float > > historanges;
  ///////////////////// choosing the Variables to run on...
  lstVar.clear();
  lstVar.push_back("HT4jetsMuonMET");
  historanges.push_back(make_pair(0,3000));
  lstVar.push_back("MassHadTop");
  historanges.push_back(make_pair(0,3000));
  //just for kinematic plots after b-tagging
  lstVar.push_back("MET");
  historanges.push_back(make_pair(0,1000));
  lstVar.push_back("PtMuon");
  historanges.push_back(make_pair(0,500));
  lstVar.push_back("PtQuark1");
  historanges.push_back(make_pair(0,1000));
  lstVar.push_back("PtQuark2");
  historanges.push_back(make_pair(0,1000));
  lstVar.push_back("PtLepbquark");
  historanges.push_back(make_pair(0,1000));
  lstVar.push_back("PtHadbquark");
  historanges.push_back(make_pair(0,1000));  
     
  int nvar = lstVar.size ();
  
  MSPlot_variables_fixedbins["HT4jetsMuonMET"] = new MultiSamplePlot(datasets_forMSPlot_variables_fixedbins, "HT",nfixedbins, historanges[0].first, historanges[0].second, "HT (GeV)");
  MSPlot_variables_fixedbins["MassHadTop"] = new MultiSamplePlot(datasets_forMSPlot_variables_fixedbins, "MassHadTop",nfixedbins, historanges[1].first, historanges[1].second, "Reconstructed top mass (GeV)");    
  MSPlot_variables_fixedbins["MET"] = new MultiSamplePlot(datasets_forMSPlot_variables_fixedbins, "MET",nfixedbins, historanges[2].first, historanges[2].second, "Missing transverse energy (GeV)");    
  MSPlot_variables_fixedbins["PtQuark1"] = new MultiSamplePlot(datasets_forMSPlot_variables_fixedbins, "",nfixedbins, historanges[3].first, historanges[3].second, "Pt jet 1 (MVA) (GeV/c)");
  MSPlot_variables_fixedbins["PtQuark2"] = new MultiSamplePlot(datasets_forMSPlot_variables_fixedbins, "",nfixedbins, historanges[4].first, historanges[4].second, "Pt jet 2 (MVA) (GeV/c)");    
  MSPlot_variables_fixedbins["PtLepbquark"] = new MultiSamplePlot(datasets_forMSPlot_variables_fixedbins, "",nfixedbins, historanges[5].first, historanges[5].second, "Pt leptonic jet(MVA) (GeV/c)");
  MSPlot_variables_fixedbins["PtHadbquark"] = new MultiSamplePlot(datasets_forMSPlot_variables_fixedbins, "",nfixedbins, historanges[6].first, historanges[6].second, "Pt hadronic jet (MVA) (GeV/c)");    
  
  map<string,vector<float> > VariableValuesMap;
  map<string,int> nbinsMap;
  nbinsMap[xvariable] = nbinsxvariable;
  nbinsMap[yvariable] = nbinsyvariable;  

  ////////// BINNING
  //Loading of the binning
  map < int, TAxis * > mapAxis; // the int designates the HT bin
  map < int, vector < pair < TH1F*,Dataset* > > > Histos_xvariableBins; // the int designates the HT bin
  vector < pair < TH1F*,Dataset* > > Histos_datasets;
  TFile fBinning (anaEnv.binningFile.c_str (), "READ");
  float nBinsBinningFile = 0;
  unsigned int xarraysize = nbinsMap[xvariable] + 2;
  //double xbins[xarraysize];
  double *xbins = NULL;
  
  //Binning for HT
  TAxis *xaxis = NULL;
  char txAxisName[150];
  sprintf (txAxisName, "Binning_%s_SM", xvariable.c_str ());
  if (!(anaEnv.isMC && anaEnv.MCRound==0))
  {
    fBinning.GetObject (txAxisName, xaxis);
    xbins = xaxis->GetXbins ()->fArray; //this is the HT binning
    for (unsigned int b = 0; b < xarraysize; b++)
    {        
	cout<<" xbins["<<b<<"] = "<<xbins[b]<<endl;
    }
  }

  //Binning for MassHadTop
  TAxis *axis = NULL;
  char tAxisName[150];
  for(int k = 1;k < nbinsMap[xvariable]+1; k++) //k < nbinsMap[xvariable]+1
  {    
    axis = NULL;
    //char tAxisName[150];
    ostringstream k_sstream;
    k_sstream << k;
    
    sprintf (tAxisName, "Binning_%s_SM_HTbin%i", yvariable.c_str (),k);
    cout<<" - Getting "<<tAxisName<<"..."<<endl;
    fBinning.GetObject (tAxisName, axis);
  
    if (axis == NULL)
    {
      cout<<"  While reading the  "<<tAxisName<<"  the axis is NULL ..... "<<endl;
      mapAxis.insert (make_pair (k, (TAxis *) NULL));
      //cout<<"range obs: "<<historanges[i].first<<" -> "<<historanges[i].second<<endl;
    }
    else
    {
      mapAxis.insert (make_pair (k, new TAxis (*axis))); 
	//for( map<string, TAxis*>::iterator ii=mapAxis.begin(); ii!=mapAxis.end(); ++ii)
	//	{ 
	//	cout << (*ii).first << ": " <<(*ii).second->GetNbins()<< " Title  "<<(*ii).second->GetName()<<" min "<<(*ii).second->GetXmin()<<"  max   "<<(*ii).second->GetXmax()<<endl;
	//	}
	
      nBinsBinningFile = axis->GetNbins ();
      //cout<<" nBinsBinningFile = "<<nBinsBinningFile<<endl;
       cout<<"  Axis inserted in map"<<endl;       
    }
    /*for (unsigned int b = 0; b < xarraysize; b++)
    {        
	cout<<"    ybins["<<b<<"] = "<<(axis->GetXbins ()->fArray)[b]<<endl;
    }*/
    
    if (verbose > 0)
	cout << " - Creating histos..." << endl;
    if (axis == NULL) {
        cout<<"Warning: axis == NULL"<<endl;
	//cout<<"nbins = "<<nbins<<", range obs: "<<historanges[i].first<<" -> "<<historanges[i].second<<", variable "<<lstVar[i]<<endl;	
	vector < pair < TH1F*,Dataset* > > Histos_datasets;
	for (unsigned int d = 0; d < datasets.size (); d++)
	{  
	   Histos_datasets.push_back(make_pair (new TH1F (TString ("h")+yvariable+"_HTbin"+k_sstream.str()+"_"+datasets[d]->Name(), "data", nfixedbins, historanges[1].first, historanges[1].second),datasets[d])); //to be fixed!!
           (Histos_datasets.back()).first->GetXaxis()->SetTitle(XaxisLabel.c_str());
	   (Histos_datasets.back()).first->GetYaxis()->SetTitle(YaxisLabel.c_str());
	}
	Histos_xvariableBins[k] = Histos_datasets;
    }
    else {	  
        Histos_datasets.clear();
	for (unsigned int d = 0; d < datasets.size (); d++)
	{  
	   Histos_datasets.push_back(make_pair (new TH1F (TString ("h")+yvariable+"_HTbin"+k_sstream.str()+"_"+datasets[d]->Name(), "data", axis->GetNbins (), axis->GetXbins ()->fArray),datasets[d]));
           (Histos_datasets.back()).first->GetXaxis()->SetTitle(XaxisLabel.c_str());
	   (Histos_datasets.back()).first->GetYaxis()->SetTitle(YaxisLabel.c_str());
	}
	Histos_xvariableBins[k] = Histos_datasets;
	//cout<<" name = "<<(Histos_xvariableBins[k][0].first)->GetName()<<endl;	
    } 
      	//delete axis;
  }
  //fBinning.Close (); //crash later on... actually this is why you first had to put the axes in the ma, then close the file, and make the histos with the aid of the map... I am trying to do it too quickly here... 
  
      
  //////////////// LOOPING OVER DATASETS
  float oldfillweight = 0, fillweight = 0;
  vector<float> eventweightvector;
  int counter[nbinsMap[xvariable]];
  for(int k = 1;k < nbinsMap[xvariable]+1; k++){
    counter[k] = 0;
  }
  
  for (unsigned int d = 0; d < datasets.size (); d++)
  {    
      cout << "   Dataset " << d << ": " << datasets[d]->Title () << endl;      
      string dataSetName = datasets[d]->Name();

      if(!(datasets[d]->Title () == "Data" || datasets[d]->Title () == "data" || datasets[d]->Title () == "DATA"))
        oldfillweight = datasets[d]->NormFactor () * Luminosity;
      else
      	oldfillweight = 1;
      cout<<"for dataset  "<<datasets[d]->Title()<<"  oldfillweight = "<<oldfillweight<<" and normfactor = "<<datasets[d]->NormFactor()<<endl;
      
      string filename;
      TFile* f1 = 0;
      TTree* seleventtree = 0;
      filename = treesdirname + datasets[d]->Name() + "_tree.root"; //previously: datasets[d]->Title() !!
      f1 = new TFile (filename.c_str (), "READ");
      seleventtree = (TTree*) f1->Get("OBS");  

      int  nevt = seleventtree->GetEntries();
      cout << "	Loop over " << nevt << " events " << endl;

      for (unsigned int ievt = 0; ievt < nevt; ievt++)
      {
      
	cout << "Processing the " << ievt << "th event" << " from "<<nevt<<flush << "\r";	
	
		  
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
	fillweight = oldfillweight*scaleFactor; //this is the weight for the events that are filled in the histos
	
	////cout<<"Looping over variables"<<endl;
	for (unsigned int v = 0; v < lstVar.size (); v++) {
	      //cout<<"lstVar["<<v<<"] = "<<lstVar[v]<<endl;
	      variable_string = lstVar[v];
	      //vector<float> ValuesThisSample;
              VariableValuesMap[variable_string].push_back(variable[v]);	  	    	  
	}/////////////end of loop on variables
	
	/*//b-tagging efficiencies...; be careful, can already be put in the trees now (not now)
	if(datasets[d]->Name()=="TTbarJets_SemiMuon" || (datasets[d]->Name()).find("Tprime")<=(datasets[d]->Name()).size()){//check again!
	    //cout<<"Applying b-efficiency scale factor to TTbarJets_SemiMuon"<<endl;
	    float btageff_SF = 0.95;
	    fillweight = fillweight*btageff_SF;		 
	}
	if(datasets[d]->Name()=="TTbarJets_Other"){
	    //cout<<"Applying b-efficiency scale factor to TTbarJets_Other"<<endl;
	    float btageff_SF = 0.95;
	    fillweight = fillweight*btageff_SF;	 
	}
	if((datasets[d]->Name()).find("ST_")==0){
	    //cout<<"Applying b-efficiency scale factor to Single top"<<endl;
	    float btageff_SF = 0.95;
	    fillweight = fillweight*btageff_SF;	 
	}*/
	/*		
	if(datasets[d]->Name()=="WJets_Nobtagging"){ //be careful, can already be put in the trees now (since Sep 5th 2011, but maybe I'll change it again)
	    //cout<<"RESCALING WJETS (without btagging)"<<endl;
	    float btageff = 1.;
	    if(doJESShift==0) btageff= 0.14859;
	    else if(doJESShift==1) btageff= 0.14898; //JES minus
	    else if(doJESShift==2) btageff= 0.13874; //JES plus
	    
	    fillweight = fillweight*btageff; 
	}*/
	
	eventweightvector.push_back(fillweight);
	
	if(anaEnv.MCRound!=0)
	{
	   float xvariable_value = -9999;
	   for (unsigned int v = 0; v < lstVar.size (); v++) {
	     if(!((datasets[d]->Name().find("Tprime")<=datasets[d]->Name().size()) && !(datasets[d]->Name().find("NP_overlay")==0)))
	     {	     
	        if(lstVar[v]=="HT4jetsMuonMET")
		   MSPlot_variables_fixedbins["HT4jetsMuonMET"]->Fill(variable[v], datasets[d], true, Luminosity*fillweight/oldfillweight); //not fillweight!, but scaleFactor(*extrascalefactor)		    
    	        if(lstVar[v]=="MassHadTop")
		   MSPlot_variables_fixedbins["MassHadTop"]->Fill(variable[v], datasets[d], true, Luminosity*fillweight/oldfillweight);
	        if(lstVar[v]=="MET")
		   MSPlot_variables_fixedbins["MET"]->Fill(variable[v], datasets[d], true, Luminosity*fillweight/oldfillweight);
		if(lstVar[v]=="PtQuark1")
		   MSPlot_variables_fixedbins["PtQuark1"]->Fill(variable[v], datasets[d], true, Luminosity*fillweight/oldfillweight);
                if(lstVar[v]=="PtQuark2")
		   MSPlot_variables_fixedbins["PtQuark2"]->Fill(variable[v], datasets[d], true, Luminosity*fillweight/oldfillweight);
	        if(lstVar[v]=="PtLepbquark")
		   MSPlot_variables_fixedbins["PtLepbquark"]->Fill(variable[v], datasets[d], true, Luminosity*fillweight/oldfillweight);
		if(lstVar[v]=="PtHadbquark")
		   MSPlot_variables_fixedbins["PtHadbquark"]->Fill(variable[v], datasets[d], true, Luminosity*fillweight/oldfillweight); 
	     }
	     
	     //cout<<"HT = "<<variable[0]<<endl;
	     if(lstVar[v]==xvariable){
	        xvariable_value = variable[v]; //hmm, then the xvariable has to come first in the lstVar list, because to be used in loop later on... bad coding!!
	     }
	     if(lstVar[v]==yvariable){
	       for(int k = 1;k < nbinsMap[xvariable]+1; k++)//k < nbinsMap[xvariable]+1
  	       {		    
		  //cout<<"---> HT bin "<<k<<endl;
       		  double rightcriterium;
       		  if(k==nbinsMap[xvariable]){
         		rightcriterium = xbins[k+1] + 1000000;
		  }
       		  else{
         		rightcriterium = xbins[k+1];
		  }
       
                  //cout<<" range of current HT bin = ["<<xbins[k]<<" , "<<rightcriterium<<"["<<endl;
       		  //cout<<variable[v]<<endl;
		  if(xvariable_value>=xbins[k] && xvariable_value<rightcriterium) //then the event is in the 1st bin of the base variable (HT)
       		  {
		        counter[k]++;
          		(Histos_xvariableBins[k][d].first)->Fill(variable[v],fillweight);
			int N = (Histos_xvariableBins[k][d].first)->GetNbinsX();
			(Histos_xvariableBins[k][d].first)->SetBinContent(N,(Histos_xvariableBins[k][d].first)->GetBinContent(N)+(Histos_xvariableBins[k][d].first)->GetBinContent(N+1)); // same as in MultiSamplePlot
			(Histos_xvariableBins[k][d].first)->SetBinContent(N+1,0); // same as in MultiSamplePlot
			//cout<<" name = "<<(Histos_xvariableBins[k][d].first)->GetName()<<endl;
       		  }
	       }
	     }
	   }
	
	}

      } //end loop on events
  } //end loop on datasets

  //regarding binning   
  if (anaEnv.isMC && anaEnv.MCRound==0){	
      cout<<endl<<" ... Creating the binning"<<endl;	
      MakeBinning NewBins;
      NewBins.Binning_forTprimeAnalysis(VariableValuesMap,nbinsMap,eventweightvector);
  }
  else{  
    fout->cd();
    cout<<" - Writing histograms"<<endl;
    for(int k = 1;k < nbinsMap[xvariable]+1; k++) //k < nbinsMap[xvariable]+1
    {  
       cout<<"counter["<<k<<"] = "<<counter[k]<<endl;
       MSPlot_xvariableBins[k] = new MultiSamplePlot(Histos_xvariableBins[k]); 

       //for (unsigned int d = 0; d < datasets.size (); d++)
       //{
	//     (Histos_xvariableBins[k][d].first)->Write(); 
       //}	       
    }
    for(map<int,MultiSamplePlot*>::const_iterator it = MSPlot_xvariableBins.begin(); it != MSPlot_xvariableBins.end(); it++)
    {
         MultiSamplePlot *temp = it->second;
         int name_int = it->first;
	 ostringstream name_sstream;
	 name_sstream << name_int;
         temp->Draw(false, "HTbin"+name_sstream.str(), true, true, true, true, true, 1); //ttbar stond op true, maar door iets in MS werd het toch niet samen genomen
         temp->Write(fout, "HTbin"+name_sstream.str(), false, "");
     }
  }

  //fout->Close();
  cout<<"Output written in "<<fout->GetName()<<endl;
  
if (!(anaEnv.isMC && anaEnv.MCRound==0)){ 
  //////////////////////////////////////////////////////////////////////  
  //real 2D to 1D converted templates; following GOFSuperplotmaker_1DHistoConverter.cc
  cout<<" - Producing 2D to 1D converted templates"<<endl;
  vector<TH1F* > hHT_Mtop_1Dconverted; //vector of dataset histograms, for each entry in the vector there is the histogram of a dataset which should be 1D converted from 2D
  vector<TH1F* > hHT_Mtop_1Dconverted_JESminus; //not for data, but maybe I put a dummy in for the moment...
  vector<TH1F* > hHT_Mtop_1Dconverted_JESplus; //not for data  
  for(unsigned int d=0;d<datasets.size();d++){
     string datasetname = datasets[d]->Name();
     string histoName = datasetname + "_Mtop";
     TH1F* htemp = new TH1F(histoName.c_str(),histoName.c_str(),nbinsMap[xvariable]*nbinsMap[yvariable],0,nbinsMap[xvariable]*nbinsMap[yvariable]);
     hHT_Mtop_1Dconverted.push_back(htemp);
    
     string histoName_JESminus = datasetname + "_Mtop_JESminus";
     TH1F* htemp_JESminus = new TH1F(histoName_JESminus.c_str(),histoName_JESminus.c_str(),nbinsMap[xvariable]*nbinsMap[yvariable],0,nbinsMap[xvariable]*nbinsMap[yvariable]);
     hHT_Mtop_1Dconverted_JESminus.push_back(htemp_JESminus);
    
     string histoName_JESplus = datasetname + "_Mtop_JESplus";
     TH1F* htemp_JESplus = new TH1F(histoName_JESplus.c_str(),histoName_JESplus.c_str(),nbinsMap[xvariable]*nbinsMap[yvariable],0,nbinsMap[xvariable]*nbinsMap[yvariable]);
     hHT_Mtop_1Dconverted_JESplus.push_back(htemp_JESplus);
  }
  
  TH1F* hHT_Mtop_1Dconverted_TotalSMExpected = new TH1F("TotalSMExpected_Mtop","TotalSMExpected",nbinsMap[xvariable]*nbinsMap[yvariable],0,nbinsMap[xvariable]*nbinsMap[yvariable]);
  //hHT_Mtop_1Dconverted_TotalSMExpected->Sumw2(); //hmm?
  TH1F* hHT_Mtop_1Dconverted_TotalSMExpected_copy = new TH1F("TotalSMExpected_Mtop","TotalSMExpected",nbinsMap[xvariable]*nbinsMap[yvariable],0,nbinsMap[xvariable]*nbinsMap[yvariable]);

  
  int b=1,b_JESminus=1,b_JESplus=1; //0 is underflowbin; k is the bin number of the 1D converted histogram
  int b_remember=1, b_JESminus_remember=1,b_JESplus_remember=1;
  int nbins_MassHadTop = 0;
  for(int k = 1;k < nbinsMap[xvariable]+1; k++)
  {
     for(unsigned int d=0;d<datasets.size();d++)
     {
        cout<<"d = "<<d<<", "<<datasets[d]->Name()<<endl;
	//nominal
	TH1F* htemp =0;
	htemp = (TH1F*) Histos_xvariableBins[k][d].first;
	b=b_remember;
	cout<<"b before loop over Mtop bins = "<<b<<endl;
        nbins_MassHadTop = htemp->GetNbinsX(); //should be 15
	cout<<"nbins_MassHadTop = "<<nbins_MassHadTop<<endl;	
	
	for(int j=2;j<nbins_MassHadTop+1;j++){
	    float bincontentj = 0;
	    bincontentj = htemp->GetBinContent(j);
	    cout<<" b = "<<b<<endl;
	    cout<<"      bincontent "<<j<<" = "<<bincontentj<<endl;
	    hHT_Mtop_1Dconverted[d]->SetBinContent(b,bincontentj);
	    b++;
	}	
     }
     
     cout<<" b_remember before update = "<<b_remember<<endl;
     b_remember = b_remember + nbins_MassHadTop - 1;
     cout<<" b_remember after update = "<<b_remember<<endl;
     b_JESminus_remember = b_JESminus_remember + nbins_MassHadTop - 1;
     b_JESplus_remember = b_JESplus_remember + nbins_MassHadTop - 1;
  
  }
  
  fout->cd();
 
  int d_ttbarsemimu = -1,d_ttbarother = -1,d_WJets = -1, d_ZJets = -1, d_ST_t_t = -1, d_ST_tW_t = -1, d_ST_s_t = -1, d_ST_t_tbar = -1, d_ST_tW_tbar = -1, d_ST_s_tbar = -1;
  for(unsigned int d=0;d<datasets.size();d++){
     if(datasets[d]->Name()=="TTbarJets_SemiMuon"){
        d_ttbarsemimu = d;
     }
     if(datasets[d]->Name()=="TTbarJets_Other"){
        d_ttbarother = d;
     }
     if(datasets[d]->Name().find("WJets")==0){
        d_WJets = d;
     }
     if(datasets[d]->Name().find("ZJets")==0){
        d_ZJets = d;
     }
     if(datasets[d]->Name().find("ST_tChannel_t")==0){
        d_ST_t_t = d;
     }
     if(datasets[d]->Name().find("ST_tChannel_tbar")==0){
        d_ST_t_tbar = d;
     }
     if(datasets[d]->Name().find("ST_tWChannel_t")==0){
        d_ST_tW_t = d;
     }
     if(datasets[d]->Name().find("ST_tWChannel_tbar")==0){
        d_ST_tW_tbar = d;
     }
     if(datasets[d]->Name().find("ST_sChannel_t")==0){
        d_ST_s_t = d;
     }
     if(datasets[d]->Name().find("ST_sChannel_tbar")==0){
        d_ST_s_tbar = d;
     }
     
     if(!(datasets[d]->Name()=="Data" || datasets[d]->Name()=="data" || datasets[d]->Name()=="DATA") && !(datasets[d]->Name().find("Tprime")<=datasets[d]->Name().size()))
        hHT_Mtop_1Dconverted_TotalSMExpected->Add(hHT_Mtop_1Dconverted[d]);
  }
  
  for(int bini=0;bini<hHT_Mtop_1Dconverted_TotalSMExpected->GetNbinsX()+1;bini++)
  {
     hHT_Mtop_1Dconverted_TotalSMExpected_copy->SetBinContent(bini,hHT_Mtop_1Dconverted_TotalSMExpected->GetBinContent(bini));  
  }

  if(d_ttbarsemimu != -1 && d_ttbarother != -1){
    hHT_Mtop_1Dconverted[d_ttbarsemimu]->Add(hHT_Mtop_1Dconverted[d_ttbarother]); //add the TTJets histos
    hHT_Mtop_1Dconverted_JESminus[d_ttbarsemimu]->Add(hHT_Mtop_1Dconverted_JESminus[d_ttbarother]); //add the TTJets histos
    hHT_Mtop_1Dconverted_JESplus[d_ttbarsemimu]->Add(hHT_Mtop_1Dconverted_JESplus[d_ttbarother]); //add the TTJets histos
  }
  if(d_WJets != -1 && d_ZJets != -1){
    hHT_Mtop_1Dconverted[d_WJets]->Add(hHT_Mtop_1Dconverted[d_ZJets]); //add the ZJets histos
    hHT_Mtop_1Dconverted_JESminus[d_WJets]->Add(hHT_Mtop_1Dconverted_JESminus[d_ZJets]); //add the ZJets histos
    hHT_Mtop_1Dconverted_JESplus[d_WJets]->Add(hHT_Mtop_1Dconverted_JESplus[d_ZJets]); //add the ZJets histos
  }
  if(d_ST_t_t != -1 && d_ST_t_tbar != -1){
    hHT_Mtop_1Dconverted[d_ST_t_t]->Add(hHT_Mtop_1Dconverted[d_ST_t_tbar]); //add the  histos
    hHT_Mtop_1Dconverted_JESminus[d_ST_t_t]->Add(hHT_Mtop_1Dconverted_JESminus[d_ST_t_tbar]); //add the  histos
    hHT_Mtop_1Dconverted_JESplus[d_ST_t_t]->Add(hHT_Mtop_1Dconverted_JESplus[d_ST_t_tbar]); //add the  histos
  }
  if(d_ST_tW_t != -1 && d_ST_tW_tbar != -1){
    hHT_Mtop_1Dconverted[d_ST_tW_t]->Add(hHT_Mtop_1Dconverted[d_ST_tW_tbar]); //add the  histos
    hHT_Mtop_1Dconverted_JESminus[d_ST_tW_t]->Add(hHT_Mtop_1Dconverted_JESminus[d_ST_tW_tbar]); //add the  histos
    hHT_Mtop_1Dconverted_JESplus[d_ST_tW_t]->Add(hHT_Mtop_1Dconverted_JESplus[d_ST_tW_tbar]); //add the  histos
  }
  if(d_ST_s_t != -1 && d_ST_s_tbar != -1){
    hHT_Mtop_1Dconverted[d_ST_s_t]->Add(hHT_Mtop_1Dconverted[d_ST_s_tbar]); //add the  histos
    hHT_Mtop_1Dconverted_JESminus[d_ST_s_t]->Add(hHT_Mtop_1Dconverted_JESminus[d_ST_s_tbar]); //add the  histos
    hHT_Mtop_1Dconverted_JESplus[d_ST_s_t]->Add(hHT_Mtop_1Dconverted_JESplus[d_ST_s_tbar]); //add the  histos
  }
  
  vector<pair<TH1F*,Dataset*> > templates1D_forMS;
  cout<<"datasets.size() = "<<datasets.size()<<endl;
  for(unsigned int d=0;d<datasets.size();d++){
    string histoname;
    cout<<" d = "<<d<<endl;
    /*histoname = hHT_Mtop_1Dconverted[d]->GetName();  
    hHT_Mtop_1Dconverted[d]->GetYaxis()->SetTitle("#Events");
    hHT_Mtop_1Dconverted[d]->Write(histoname.c_str());*/
   if(!(datasets[d]->Name()=="TTbarJets_Other") && !(datasets[d]->Name().find("ZJets")==0) && !(datasets[d]->Name().find("ST_tChannel_tbar")==0) && !(datasets[d]->Name().find("ST_tWChannel_tbar")==0) && !(datasets[d]->Name().find("ST_sChannel_tbar")==0))
   {
    histoname = hHT_Mtop_1Dconverted[d]->GetName();
    //cout<<" histoname = "<<histoname<<endl;
    if(datasets[d]->Name()=="TTbarJets_SemiMuon")
    {
      histoname = "TTJets_Mtop" + postfix; //because 'other' was added
    }
    else if(datasets[d]->Name().find("WJets")==0)
    {
  ///    histoname = "WJetsZJets_Mtop"; //because 'ZJets' was added
       histoname ="WJets_Mtop" + postfix;
    }
    else if(datasets[d]->Name().find("ST_tChannel_t")==0)
    {
       histoname ="ST_tChannel_Mtop" + postfix;
    }
    else if(datasets[d]->Name().find("ST_tWChannel_t")==0)
    {
       histoname ="ST_tWChannel_Mtop" + postfix;
    }
    else if(datasets[d]->Name().find("ST_sChannel_t")==0)
    {
       histoname ="ST_sChannel_Mtop" + postfix;
    }
    else if(datasets[d]->Name().find("Tprime")<=datasets[d]->Name().size()){ //the name has NP or NP_overlay in it, should be gone in name of template histo that is written in file...
       ostringstream mass_sstream;
       mass_sstream << datasets[d]->Mass();
       histoname = "Tprime" + mass_sstream.str() + "_Mtop" + postfix;
    }
    else if(datasets[d]->Name()=="Data" || datasets[d]->Name()=="data" || datasets[d]->Name()=="DATA")
    {
       histoname = histoname;
    }
    else{
       histoname = histoname + postfix;
    }
    hHT_Mtop_1Dconverted[d]->GetYaxis()->SetTitle("#Events");
    if(!(doJESShift!=0 && (datasets[d]->Name()=="Data" || datasets[d]->Name()=="data" || datasets[d]->Name()=="DATA"))){
       cout<<"   -> Writing "<<histoname<<"..."<<endl; 
       hHT_Mtop_1Dconverted[d]->Write(histoname.c_str());
    }
    templates1D_forMS.push_back(make_pair(hHT_Mtop_1Dconverted[d],datasets[d])); //tbc?
  
   }
  }
  
  if(doJESShift==0) hHT_Mtop_1Dconverted_TotalSMExpected_copy->Write();
  
  MultiSamplePlot* templates1D_MS = new MultiSamplePlot(templates1D_forMS);
  templates1D_MS->Draw(false, "Templates", true, true, true, true, true, 1);// templates1D_MS->Draw(false, "Templates", true, true, true, true, true, 1);
  templates1D_MS->Write(fout, "Templates", false, "");
  
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot_variables_fixedbins.begin(); it != MSPlot_variables_fixedbins.end(); it++)
  {
    MultiSamplePlot *temp = it->second;
    string name = it->first;
    temp->Draw(false, name, true, true, true, true, true,1);//temp->Draw(false, name, true, true, true, true, true,5);
    temp->Write(fout, name, false, "");
  }

  cout<<"Output written in "<< fout->GetName()<<endl;
  fout->Close();
 
} //end if(!(anaEnv.isMC && anaEnv.MCRound==0))
 
  cout << "**********************************************************************" << endl;
  cout << "           End of the program !!" << endl;
  cout << "**********************************************************************" << endl;
  return 0;

}
