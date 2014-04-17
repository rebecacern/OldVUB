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
#include "TMultiGraph.h"
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

// TH1F* RefVDistribution_ = new TH1F("V distribution MC","Monte Carlo V distribution ",50,0,3000);
// RefVDistribution_->GetXaxis()->SetTitle("V");
 
 //xml file
  const char *xmlfile = "../config/myNPconfig.xml";
 //const char *xmlfile = "../config/myTprimeconfig.xml";
 
 /////////////////////////////////////
 /// AnalysisEnvironment 
 ////////////////////////////////////
 AnalysisEnvironment anaEnv;
 cout << "Loading environment ..." << endl;
 AnalysisEnvironmentLoader anaLoad (anaEnv, xmlfile);  
  
 vector<float> FractionHWEvts = anaEnv.FractionHWEvts; 
 unsigned int nofFractionsHWEvts = FractionHWEvts.size();
 float FractionHWEvts_array[nofFractionsHWEvts];
 float FractionHWEvts_array_syst[nofFractionsHWEvts]; //just as a trick to superimpose graphs...

 float pdataArray[nofFractionsHWEvts];
 float pdataArray_syst[nofFractionsHWEvts];

 TFile *_treeNullHypoFile = TFile::Open("ROUND16_nosyst/GOF_Vtree_StatOutput_SM_1000PsExps.root");  
 TFile *_treeNullHypoFile_syst = 0;//TFile::Open("ROUND13/GOF_Vtree_StatOutput_ROUND13_10_Bins_Corr_PF_nosyst_1000PsExps.root"); 
 TFile *dataVtreeFile = new TFile("ROUND16_nosyst/GOF_Vtree_StatOutput_Data.root");
 TFile *dataVtreeFile_syst = 0;//new TFile("ROUND13/GOF_Vtree_StatOutput_Data_10_Bins_Corr_PF_JESSyst.root");
 
 _treeNullHypoFile->cd();
 for(unsigned int k=0;k<nofFractionsHWEvts;k++){
   ostringstream Fractionstrstream;
   Fractionstrstream << FractionHWEvts[k];
   FractionHWEvts_array[k] = FractionHWEvts[k];
   FractionHWEvts_array_syst[k] = FractionHWEvts[k]+0.001;//just as a trick to superimpose graphs...
   cout<<"For x = "<<Fractionstrstream.str()<<endl;
 //  TFile *_treeNullHypoFile = TFile::Open("mergedVtrees_GOF_Round2_100PsExps.root");
   //TFile *_treeNullHypoFile = TFile::Open("GOF_38X_Round5_METcut40_8Bins_nosyst/mergedVtrees_GOF_Round5_METcut40_8Bins_nosyst_100PsExps.root");

   TTree *T = (TTree*) _treeNullHypoFile->Get("Vtree_fr"+TString(Fractionstrstream.str()));
   vector<float> MCVvect;
   float V;
   T->SetBranchAddress("V",&V);
   for(int i=0;i<T->GetEntriesFast();i++){
 		T->GetEntry(i);
 		MCVvect.push_back(V);
//		RefVDistribution_->Fill(V);
   }
   
   TTree *T_syst = 0;
   vector<float> MCVvect_syst;
   float V_syst;
   if(anaEnv.Systematics != 0){
     T_syst = (TTree*) _treeNullHypoFile_syst->Get("Vtree_fr"+TString(Fractionstrstream.str()));
     T_syst->SetBranchAddress("V",&V_syst);
     for(int i=0;i<T_syst->GetEntriesFast();i++){
 		T_syst->GetEntry(i);
 		MCVvect_syst.push_back(V_syst);
//		RefVDistribution_->Fill(V);
     }
   }
  
  float Vdata = -999;
  TTree *T_data = (TTree*) dataVtreeFile->Get("Vtree_fr"+TString(Fractionstrstream.str()));
  T_data->SetBranchAddress("V",&Vdata);
  for(int i=0;i<T_data->GetEntriesFast();i++){ //normally 1 value
 		T_data->GetEntry(i);
  }
  cout<<"  V data (no syst) = "<<Vdata<<endl;
  
  float Vdata_syst= -999;
  TTree *T_data_syst = 0;
  if(anaEnv.Systematics != 0){     
     T_data_syst = (TTree*) dataVtreeFile_syst->Get("Vtree_fr"+TString(Fractionstrstream.str()));
     T_data_syst->SetBranchAddress("V",&Vdata_syst);
     for(int i=0;i<T_data_syst->GetEntriesFast();i++){ //normally 1 value
 		T_data_syst->GetEntry(i);
     }
     cout<<"  V data (syst) = "<<Vdata_syst<<endl;
  }
  
  pdataArray[k] = -999;
  pdataArray_syst[k] = -999;
  
  WeightProbaCalculator myWeightProbaCalculator;
  myWeightProbaCalculator.CalculateProba(MCVvect,Vdata);
  pdataArray[k] = myWeightProbaCalculator.GetProba();
  cout<<"  pvalue data (no syst): "<<pdataArray[k]<<endl;
  
  if(anaEnv.Systematics != 0){
    WeightProbaCalculator myWeightProbaCalculator_syst;
    myWeightProbaCalculator_syst.CalculateProba(MCVvect_syst,Vdata_syst);
    pdataArray_syst[k] = myWeightProbaCalculator_syst.GetProba();
    cout<<"  pvalue data (with syst): "<<pdataArray_syst[k]<<endl;
  }
  
  
//  TFile * VdistrFile = new TFile("VdistrFile_GOFRound5_METcut40_8Bins_nosyst_SignCombSUM_x0p1.root","RECREATE");
//  RefVDistribution_->Write();
 }
 TFile *_output = TFile::Open(dataVtreeFile->GetName(),"UPDATE"); //or no _syst...
 _output->cd();
 TCanvas* c1 = new TCanvas("c1","",500,500);
// TMultiGraph* mg = new TMultiGraph(); 
 TGraph* graph_pvaluevsFraction = new TGraph(nofFractionsHWEvts,FractionHWEvts_array,pdataArray);
 graph_pvaluevsFraction->SetMarkerColor(kBlue); //kBlue no syst, kGreen syst
 graph_pvaluevsFraction->SetMarkerStyle(21);
 graph_pvaluevsFraction->GetXaxis()->SetTitle("x");
 graph_pvaluevsFraction->GetYaxis()->SetTitle("p-value data");
 graph_pvaluevsFraction->GetYaxis()->SetRangeUser(0,1);
 graph_pvaluevsFraction->SetTitle("");
 graph_pvaluevsFraction->Draw("AP");
// mg->Add(graph_pvaluevsFraction,"AP");
// mg->Draw("AP");

 if(anaEnv.Systematics != 0){
   TGraph* graph_pvaluevsFraction_syst = new TGraph(nofFractionsHWEvts,FractionHWEvts_array_syst,pdataArray_syst);
   graph_pvaluevsFraction_syst->SetMarkerColor(kGreen); //kBlue no syst, kGreen syst
   graph_pvaluevsFraction_syst->SetMarkerStyle(21);
   graph_pvaluevsFraction_syst->GetXaxis()->SetTitle("x");
   graph_pvaluevsFraction_syst->GetYaxis()->SetTitle("p-value data");
   graph_pvaluevsFraction_syst->GetYaxis()->SetRangeUser(0,1);
   graph_pvaluevsFraction_syst->SetTitle("");
   graph_pvaluevsFraction_syst->Draw("P"); //SAMEAP?
 }
/* 
 TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
// TLegend* leg = new TLegend(0.30,0.80,0.99,0.99);
 leg->AddEntry("graph_pvaluevsFraction","systematics","P");
 leg->Draw("SAME");
 */ //doesn't work
 
 c1->Write("pvalue_vs_x_superimposed",TObject::kOverwrite);
 
 cout<<endl<<"Output written to "<<_output->GetName()<<endl;
 _output->Close() ;


}
