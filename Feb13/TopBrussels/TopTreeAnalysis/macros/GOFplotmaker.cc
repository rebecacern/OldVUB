//This macro makes use of the MultiSamplePlot class
//you have to put all wanted MC files to add=1 in the config!! (and data =0) + Make sure the correct MCFile is put in the config

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


int
main (int argc, char *argv[])
{

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptFit(1111);

  bool Tprime = true;
  const char *xmlfile = 0;
  if(!Tprime) xmlfile = "../config/myNPconfig.xml";
  if(Tprime) xmlfile = "../config/myTprimeconfig.xml"; 
  AnalysisEnvironment anaEnv;
  cout << "Loading environment ..." << endl;
  AnalysisEnvironmentLoader anaLoad (anaEnv, xmlfile);
  float WHFkFactor = 1.5;
  
  /////////////////////
  // Load Datasets	//in plotreader this was not needed
  /////////////////////
  TTreeLoader treeLoader;
  cout << " - Load datasets ..." << endl;
  vector < Dataset* > datasets;
  treeLoader.LoadDatasets (datasets, xmlfile);

  vector<MCObsExpectation *> MCObsExp; 
  //read one
  TFile fin (anaEnv.MCExpFilename.c_str (), "READ"); //the histos will be read from the MC file that was produced with GOFanalysis.cc (or FullChainEstimation.cc) in the format of MCObsExpectation class
  //TFile fin ("FullChainEstimation.root", "READ");
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
    delete temp3;
    tcMCObsExpTemp->Delete ();
    tree->Delete ();
    
   vector <string > lstVar; 
   Observables obs;
   lstVar = obs.ListOfVariables ();
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
     /*lstVar.push_back("MassHadTop");
     lstVar.push_back("MassTtbar");
     lstVar.push_back("HT4jets");*/
     lstVar.push_back("MassHadTop");
     lstVar.push_back("HT4jetsMuonMET"); //was commented...
     lstVar.push_back("MET");
     lstVar.push_back("PtMuon");
     lstVar.push_back("PtQuark1");
     lstVar.push_back("PtQuark2");
     lstVar.push_back("PtLepbquark");
     lstVar.push_back("PtHadbquark");
     //lstVar.push_back("EtaHadTop");
   }
   

   vector < vector < TH1F* > > hAllObservables; //for each observable, you have a vector of the different contributions of the processes (datasets)
   cout<<"lstVar.size () = "<<lstVar.size ()<<endl;  
    
   //LOOP OVER VARIABLES
   for (unsigned int v = 0; v < lstVar.size (); v++) { //so this is a loop over the observables you chose...
      /*cout<<"!anaEnv.runOnObs ((int) v) = "<<!anaEnv.runOnObs ((int) v)<<endl; //WHY always true?????
      if (!anaEnv.runOnObs ((int) v))
	continue; */
	
      cout<<"MCObsExp.size() = "<<MCObsExp.size()<<endl;
      for (int l=0;l<MCObsExp.size();l++){ //normally this is a loop over all '42' observables...
	if (MCObsExp[l]->GetHistoSMProcesses()->GetName()=="SMProcess_"+lstVar[v]){ //and here you require the observable to be the one you chose to plot
	   cout<<"  found MCObsExp plots....... "<<MCObsExp[l]->GetHistoSMProcesses()->GetName()<<"   "<<l<<endl;
	   vector < TH1F* > hDatasets; // all variable histos of one dataset (process) for a certain observable
	   for(unsigned int d=0;d<datasets.size();d++){ //loop over datasets (processes)	   	
		cout<<"d = "<<d<<", "<<datasets[d]->Name()<<endl;
		hDatasets.push_back(MCObsExp[l]->GetHistogram(datasets[d]->Name()));
		/*cout<<"'number of bins' = "<<hDatasets[d]->GetNbinsX()<<endl;
		for(int bini=0;bini<hDatasets[d]->GetNbinsX()+1;bini++){
		   cout<<" hDatasets["<<d<<"], content bin "<<bini<<" = "<<hDatasets[d]->GetBinContent(bini)<<endl;
		}*/
		 int nbins = hDatasets[d]->GetNbinsX();
		 cout<<"  nbins = "<<nbins<<endl;
		 /*//this is for the purpose of plotting...
		  const TArrayD* thisarray = hDatasets[d]->GetXaxis()->GetXbins();
		  int thisarraysize = sizeof(thisarray)/sizeof(thisarray[0]);
		  thisarray[thisarraysize-1] = 0//2*thisarray[thisarraysize-1] - thisarray[0]; //this is just my choice to obtain a reasonable last bin when plotted
		  hDatasets[d]->GetXaxis()->Set(nbins,thisarray);
		  hDatasets[d]->SetBinContent(nbins,hDatasets[d]->GetBinContent(nbins)); //well...*/
	/*	 float beginrange = hDatasets[d]->GetXaxis()->GetBinLowEdge(1);;// first: = 0;		 
		 float endrange = hDatasets[d]->GetXaxis()->GetBinUpEdge(nbins);//first I had GetBinLowEdge(nbins), with the old binning...
		 cout<<"beginrange = "<<beginrange<<", endrange = "<<endrange<<endl;
		 hDatasets[d]->GetXaxis()->SetLimits(beginrange,endrange);		 
	*/	 
		 
		 if(datasets[d]->Name()=="WJets_Nobtagging"){
		 	cout<<"RESCALING WJETS (of not btagging)"<<endl;
			float btageff = 0.15876;
			hDatasets[d]->Scale(btageff*WHFkFactor);		 
		 }
		 if(datasets[d]->Name()=="TTbarJets_SemiMuon" || datasets[d]->Name()=="TTbarJets_Other"){
		 	cout<<"Applying b-efficiency scale factor to TTbarJets"<<endl;
			float btageff_SF = 0.95;
			hDatasets[d]->Scale(btageff_SF);
	         }
		 if(datasets[d]->Name()=="ST_tChannel" || datasets[d]->Name()=="ST_tWChannel"){
		 	cout<<"Applying b-efficiency scale factor to Single top"<<endl;
			float btageff_SF = 0.95;
			hDatasets[d]->Scale(btageff_SF);
		 }
		 
		 if(lstVar[v]=="ET1oET4") hDatasets[d]->GetXaxis()->SetTitle("PT1oPT4");
	         if(lstVar[v]=="AllJetsPtMET") hDatasets[d]->GetXaxis()->SetTitle("HTMET (GeV)");
	         if(lstVar[v]=="MET") hDatasets[d]->GetXaxis()->SetTitle("Missing Transverse Energy (GeV)");
	         if(lstVar[v]=="PtMuon") hDatasets[d]->GetXaxis()->SetTitle("Pt Muon (GeV/c)");
		 
	         if(lstVar[v]=="MassHadTop") hDatasets[d]->GetXaxis()->SetTitle("Mass hadronic top (GeV/c^2)");
	         if(lstVar[v]=="MassTtbar") hDatasets[d]->GetXaxis()->SetTitle("Mass ttbar system (GeV/c^2)");
	         if(lstVar[v]=="HT4jets") hDatasets[d]->GetXaxis()->SetTitle("HT4jets (GeV/c)");
		 if(lstVar[v]=="HT4jetsMuonMET") hDatasets[d]->GetXaxis()->SetTitle("HT (GeV/c)");
		 if(lstVar[v]=="PtQuark1") hDatasets[d]->GetXaxis()->SetTitle("Pt jet 1 (GeV/c)");
		 if(lstVar[v]=="PtQuark2") hDatasets[d]->GetXaxis()->SetTitle("Pt jet 2 (GeV/c)");
	         if(lstVar[v]=="PtLepbquark") hDatasets[d]->GetXaxis()->SetTitle("Pt leptonic b jet (MVA) (GeV/c)");
		 if(lstVar[v]=="PtHadbquark") hDatasets[d]->GetXaxis()->SetTitle("Pt hadronic b jet (MVA) (GeV/c)");		 
		 //if(lstVar[v]=="EtaHadTop") hDatasets[d]->GetXaxis()->SetTitle("Pseudorapidity hadronic top");
	  
		 if(d==datasets.size()-1) hDatasets[d]->GetYaxis()->SetTitle("#Events"); //ok, the if is maybe not needed
	   }
	   hAllObservables.push_back(hDatasets);	   
	   break;
	}
      }
   }
     
//   TFile * fout = new TFile("./ROUND19/HTcut400_75bins/plotmaker_854ipb_fixed75bins_dataMC_Tprime400Scaled_HTcut400.root","RECREATE");
   TFile * fout = new TFile("./ROUND19/plotmaker_854ipb_fixed75bins_dataMC_Tprime400Scaled_WHFkFactor.root","RECREATE");

//   fout->cd();
   //TFile *fData = TFile::Open("~/CMSSW/CMSSW_3_5_8_patch1/src/38X/TopTreeAnalysis/macros/GOF_StatOutput_Data_20Bins.root");
   //TFile *fData = TFile::Open("GOF_StatOutput_Data_Tprime.root"); //this contains the four initial observables for the S2-method for Tprime search
//   TFile *fData = TFile::Open("./ROUND10/Data_Lumi308.5_StatOutput.root"); //this contains all initial observables for the S2-method for Tprime search

   TFile *fData = TFile::Open("./ROUND19/Data_Lumi853.8_StatOutput_fixed75bins.root"); //this contains all initial observables for the S2-method for Tprime search
//   TFile *fData = TFile::Open("./ROUND19/HTcut400_75bins/Data_Lumi853.8_StatOutput_fixed75bins_HTcut400.root");
   
   TDirectory *dir;
   cout<<" hAllObservables.size() = "<<hAllObservables.size()<<endl;
   for(unsigned int v=0;v<hAllObservables.size();v++){ //in principle, loop over the variables you chose
	 cout<<"  iter "<<v<<endl;
	 vector<pair<TH1F*,Dataset*> > vec; //vector of dataset with corresponding histos, for one specific variable   
   	 for(unsigned int d=0;d<datasets.size();d++){	 
		pair<TH1F*,Dataset*> thispair;
		thispair.first = hAllObservables[v][d]; //histogram of the v-th variable, d-th dataset
		thispair.second = datasets[d];
		vec.push_back(thispair); //there is a certain risk in the matching between histo and datasets, be careful...!!
	 }
	 MultiSamplePlot myMultiSamplePlot(vec);
	 string label = lstVar[v], directory = "AllEvents/" + lstVar[v], histokey = "hData_" + lstVar[v] + ";1";
//	 string label = lstVar[v], directory = "x=1/" + lstVar[v], histokey = "hData_" + lstVar[v] + ";1";
	 dir = fData->GetDirectory(directory.c_str());
	 dir->cd();
  	 TH1F* DataHisto = (TH1F*) dir->Get(histokey.c_str());
	  int nbins = DataHisto->GetNbinsX();	  
/*	  float beginrange = DataHisto->GetXaxis()->GetBinLowEdge(1);;// first: = 0;
	  float endrange = DataHisto->GetXaxis()->GetBinUpEdge(nbins); //first I had GetBinLowEdge(nbins), with the old binning...
	  DataHisto->GetXaxis()->SetLimits(beginrange,endrange);
*/
	  DataHisto->GetYaxis()->SetTitle("#Events"); //doesn't work on it's own...
	 myMultiSamplePlot.AddDataHisto(DataHisto);
	 myMultiSamplePlot.Draw(false,label,true,false,false,false,true); //merge TTJets //no merge () // no merge () // no merge // merge ST
	 myMultiSamplePlot.Write(fout,label,false,"");	 
   }
   cout<<"Output written in "<< fout->GetName()<<endl;;

}
