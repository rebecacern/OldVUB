#include "TStyle.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

// Root headers
#include "TArrow.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include  "TGraphErrors.h"
#include "TMath.h"
#include "TMarker.h"
#include "TPaveStats.h"
#include "TRandom3.h"

// User headers
#include "TopTreeAnalysis/BkgEstimationMethods/interface/VJetEstimation.h"

// RooFit headers
#include "RooAddition.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooCatType.h"
#include "RooCBShape.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFormula.h"
#include "RooFormulaVar.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooPolynomial.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"

#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooProdPdf.h"

#include "RooMCStudy.h"
#include "RooMinuit.h"
#include "RooMinimizer.h"
#include "RooNLLVar.h"

#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooTable.h"
#include "Roo1DTable.h"

#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/ConfInterval.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ConfidenceBelt.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/ModelConfig.h"

#include "RooStats/HypoTestInverterOriginal.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"
#include "RooStats/HybridCalculatorOriginal.h"

#include "Style.C"

using namespace RooFit;
using namespace RooStats ;
using namespace std;
//using namespace TopTree;

int main (int argc, char *argv[])
{

  clock_t start = clock();

  cout << "*************************************************************" << endl;
  cout << " Beginning of the program for the VJetEstimation analysis ! " << endl;
  cout << "*************************************************************" << endl;

  //SetStyle if needed
  //setTDRStyle();
  setGregStyle();
  //setMyStyle();
  bool verbose = true;
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////Configuration ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //Input ROOT file
  char inputname[100];
  /** TO BE ADAPTED **/
  TFile *fin = TFile::Open("StopBckg_Output.root");

  //Retrieve VJetEStimation object from file
  /** TO BE ADAPTED **/
  sprintf(inputname,"VJetEstimation-TCHE_LM--tt_W");
  VJetEstimation* vje = (VJetEstimation*) fin->Get(inputname);

  //Output ROOT file
  string postfix = "_Analysis";
  string channelpostfix = "_SemiMuon";
  string comment = "_Constants_euds";
  string rootFileName ("VJetEstimation"+postfix+channelpostfix+comment+".root");
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");
  TDirectory *myDir = 0;

//  UInt_t btag_wp_idx = 0; // B-tagging working point index
//  UInt_t njets = 4;       // Nb of selected jets
  const int NbOfCPU = 1;
  const int NbOfPE = 10000;
  const int pullNbOfBins = 100;
//  const UInt_t NbOfJetBins = vje->GetNbOfJetsBins();
  const UInt_t  NbOfJets = 3;
  const UInt_t  JetIdx = NbOfJets-vje->GetNjets(0);
  const UInt_t  NbOfBtagWorkingPoint = vje->GetNbOfBtagWorkingPoint();
  const UInt_t  BtagWorkingPoint     = 0;
  const Float_t Nttlike = vje->GetPredNtt(0,0);
  const Float_t Nvlike  = vje->GetPredNv(0,0);

  const double IntLumi = vje->GetIntLumi();

  vector<int>::iterator Idx;
  vector<int> FixedVarIdx;
  //FixedVarIdx.push_back(0);
  //FixedVarIdx.push_back(1);
  FixedVarIdx.push_back(2);

  TString wp[3]  = {"loose","medium","tight"};
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////// Analysis ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // C r e a t e   o b s e r v a b l e
  // ---------------------------------

  RooRealVar Ntt("Ntt","N_{t#bar{t}-like}",Nttlike,0.0,Nttlike*3);
  RooRealVar Nv("Nv","N_{V-like}",Nvlike,0.0,Nvlike*4);

  RooRealVar* eb[NbOfBtagWorkingPoint];
  Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 0 );
  for(int i=0; i<NbOfBtagWorkingPoint; i++){
  	eb[i] = new RooRealVar( "eb_"+wp[i],"#epsilon_{b-tag} ("+wp[i]+")",vje->GetPredEb(i,JetIdx),0.0,1.0);
	if(Idx != FixedVarIdx.end()) eb[i]->setConstant(kTRUE);
	cout<<"eb ("<<wp[i]<<")"<<eb[i]->getVal()<<endl;
  }

  RooRealVar* eudsc[NbOfBtagWorkingPoint];
  Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 1 );
  for(int i=0; i<NbOfBtagWorkingPoint; i++){
	eudsc[i] = new RooRealVar( "eudsc_"+wp[i],"#epsilon_{mis-tag} ("+wp[i]+")",vje->GetPredEudsc(i,JetIdx),0.0,1.0);
	if(Idx != FixedVarIdx.end()) eudsc[i]->setConstant(kTRUE);
	cout<<"eudsc ("<<wp[i]<<")"<<eudsc[i]->getVal()<<endl;
  }

  RooRealVar* euds[NbOfBtagWorkingPoint];
  Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 2 );
  for(int i=0; i<NbOfBtagWorkingPoint; i++){
	euds[i]  = new RooRealVar( "euds_"+wp[i],"#epsilon^{,}_{mis-tag} ("+wp[i]+")",vje->GetPredEuds(i,JetIdx),0.0,1.0);
	if(Idx != FixedVarIdx.end()) euds[i]->setConstant(kTRUE);
	cout<<"euds ("<<wp[i]<<")"<<euds[i]->getVal()<<endl;
  }

  // C r e a t e   c o n s t a n t s
  // -------------------------------

  RooConstVar n("n","number of selected jets",NbOfJets) ;

  RooConstVar e0bq("e0bq","e0bq",vje->GetTTEff0bq(JetIdx));
  RooConstVar e1bq("e1bq","e1bq",vje->GetTTEff1bq(JetIdx));
  RooConstVar e2bq("e2bq","e2bq",vje->GetTTEff2bq(JetIdx));
//  RooConstVar e3bq("e3bq","e3bq",vje->GetTTEff3bq(JetIdx));
  
  // C o n s t r u c t   a   c a t e g o r y   w i t h   l a b e l s    a n d   i n d e c e s
  // -----------------------------------------------------------------------------------------

  RooCategory nbjets("nbjets","Number of b-jets");
  nbjets.defineType("N0bjet", 0);
  nbjets.defineType("N1bjet", 1);
  nbjets.defineType("N2bjets",2);
  nbjets.defineType("N3bjets",3);

  RooCategory WPCat("WPCat","B-tagging working point") ;
  for(int i=0; i<NbOfBtagWorkingPoint; i++) WPCat.defineType(wp[i]) ;

  // C o n s t r u c t   f o r m u l a s
  // -------------------------------------------------------------------------------
  
  RooFormulaVar *p0bjets_tt[NbOfBtagWorkingPoint];
  RooFormulaVar *p1bjets_tt[NbOfBtagWorkingPoint];
  RooFormulaVar *p2bjets_tt[NbOfBtagWorkingPoint];
  RooFormulaVar *p3bjets_tt[NbOfBtagWorkingPoint];

  for(int i=0; i<NbOfBtagWorkingPoint; i++){
	p0bjets_tt[i] = new RooFormulaVar("p0bjets_tt_"+wp[i],"p0bjets_tt_"+wp[i],"(1-@0)*(1-@0)*pow((1-@1),@2-2)*e2bq+(1-@0)*pow((1-@1),@2-1)*e1bq+pow((1-@1),@2)*e0bq",RooArgList((*eb[i]),(*eudsc[i]),n,e0bq,e1bq,e2bq));
	p1bjets_tt[i] = new RooFormulaVar("p1bjets_tt_"+wp[i],"p1bjets_tt_"+wp[i],"(2*@0*(1-@0)*pow(1-@1,@2-2)+(1-@0)*(1-@0)*(@2-2)*@1*pow(1-@1,@2-3))*e2bq+(@0*pow(1-@1,@2-1)+(1-@0)*(@2-1)*@1*pow(1-@1,@2-2))*e1bq+(@2*@1*pow(1-@1,@2-1))*e0bq",RooArgList((*eb[i]),(*eudsc[i]),n,e0bq,e1bq,e2bq));
	p2bjets_tt[i] = new RooFormulaVar("p2bjets_tt_"+wp[i],"p2bjets_tt_"+wp[i],"(@0*@0*pow(1-@1,@2-2)+2*@0*(1-@0)*(@2-2)*@1*pow(1-@1,@2-3)+(1-@0)*(1-@0)*((@2-2)*(@2-3)/2)*@1*@1*pow(1-@1,@2-4))*e2bq+(@0*(@2-1)*@1*pow(1-@1,@2-2)+(1-@0)*((@2-1)*(@2-2)/2)*@1*@1*pow(1-@1,@2-3))*e1bq+((@2*(@2-1)/2)*@1*@1*pow(1-@1,@2-2))*e0bq",RooArgList((*eb[i]),(*eudsc[i]),n,e0bq,e1bq,e2bq));
	p3bjets_tt[i] = new RooFormulaVar("p3bjets_tt_"+wp[i],"p3bjets_tt_"+wp[i],"(@0*@0*(@2-2)*@1*pow(1-@1,@2-3)+2*@0*(1-@0)*((@2-2)*(@2-3)/2)*@1*@1*pow(1-@1,@2-4)+(@2>4 ? pow((1-@0),2)*((@2-2)*(@2-3)*(@2-4)/6)*pow(@1,3)*pow((1-@1),@2-5) : 0 ))*e2bq+(@0*((@2-1)*(@2-2)/2)*@1*@1*pow(1-@1,@2-3)+(1-@0)*((@2-1)*(@2-2)*(@2-3)/6)*pow(@1,3)*pow(1-@1,@2-4))*e1bq+((@2*(@2-1)*(@2-2)/6)*pow(@1,3)*pow(1-@1,@2-3))*e0bq",RooArgList((*eb[i]),(*eudsc[i]),n,e0bq,e1bq,e2bq));
  }

  RooFormulaVar *p0bjets_v[NbOfBtagWorkingPoint]; 
  RooFormulaVar *p1bjets_v[NbOfBtagWorkingPoint]; 
  RooFormulaVar *p2bjets_v[NbOfBtagWorkingPoint]; 
  RooFormulaVar *p3bjets_v[NbOfBtagWorkingPoint]; 

  for(int i=0; i<NbOfBtagWorkingPoint; i++){
	p0bjets_v[i] = new RooFormulaVar("p0bjets_v_"+wp[i], "p0bjets_v_"+wp[i] ,"pow(1-@0,@1)",RooArgList((*euds[i]),n));
	p1bjets_v[i] = new RooFormulaVar("p1bjets_v_"+wp[i], "p1bjets_v_"+wp[i] ,"@1*@0*pow(1-@0,@1-1)",RooArgList((*euds[i]),n));
	p2bjets_v[i] = new RooFormulaVar("p2bjets_v_"+wp[i], "p2bjets_v_"+wp[i] ,"(@1*(@1-1)/2)*@0*@0*pow(1-@0,@1-2)",RooArgList((*euds[i]),n));
	p3bjets_v[i] = new RooFormulaVar("p3bjets_v_"+wp[i], "p3bjets_v_"+wp[i] ,"((@1)*(@1-1)*(@1-2)/6)*pow(@0,3)*pow(1-@0,@1-3)",RooArgList((*euds[i]),n));
  }

  // C o n s t r u c t   p . d . f 's
  // -------------------------------------------------------------------------------

  RooGenericPdf *pbjets_tt[NbOfBtagWorkingPoint];
  RooGenericPdf *pbjets_v[NbOfBtagWorkingPoint];
  RooExtendPdf  *pbjets_v_ext[NbOfBtagWorkingPoint];
  RooExtendPdf  *pbjets_tt_ext[NbOfBtagWorkingPoint];

  for(int i=0; i<NbOfBtagWorkingPoint; i++){
	pbjets_tt[i] = new RooGenericPdf("pbjets_tt_"+wp[i],"pbjets_tt_"+wp[i],"(@0==0)*@1+(@0==1)*@2+(@0==2)*@3+(@0==3)*@4",RooArgList(nbjets,*p0bjets_tt[i],*p1bjets_tt[i],*p2bjets_tt[i],*p3bjets_tt[i]));
	pbjets_v[i]  = new RooGenericPdf("pbjets_v_"+wp[i],"pbjets_v_"+wp[i],"((@0)==0)*@1+(@0==1)*@2+(@0==2)*@3+(@0==3)*@4",RooArgList(nbjets,*p0bjets_v[i],*p1bjets_v[i],*p2bjets_v[i],*p3bjets_v[i]));

	pbjets_tt_ext[i] = new RooExtendPdf("pbjets_tt_ext_"+wp[i],"pbjets_tt_ext_"+wp[i],*pbjets_tt[i],Ntt);
	pbjets_v_ext[i]  = new RooExtendPdf("pbjets_v_ext_"+wp[i] ,"pbjets_v_ext_"+wp[i] ,*pbjets_v[i],Nv);
  }

  RooAddPdf* model[NbOfBtagWorkingPoint];
  for(int i=0; i<NbOfBtagWorkingPoint; i++){
	model[i] = new RooAddPdf("model"+wp[i],"model"+wp[i],RooArgList(*pbjets_tt_ext[i],*pbjets_v_ext[i]));//,RooArgList(Ntt,Nv));
  }

  RooSimultaneous simPdf("simPdf","simultaneous pdf",WPCat) ;
  for(int i=0; i<NbOfBtagWorkingPoint; i++) simPdf.addPdf(*model[i],wp[i]) ;

  // C r e a t e   d a t a s e t 
  // -------------------------------------------------------------------------------

  // Generate dataset for use in fitting below
  //RooDataSet* data = simPdf.generate(RooArgList(nbjets,WPCat),NbOfBtagWorkingPoint*(NbOfTTlike[NbOfJets-3]+NbOfVlike[NbOfJets-3])) ;
  //data->Print("v");
  //Roo1DTable* table = data->table(RooArgSet(nbjets,WPCat));
  //table->Print("v");
  

  // F i t   t h e   d a t a   a n d   c o n s t r u c t   t h e   l i k e l i h o o d   f u n c t i o n
  // ----------------------------------------------------------------------------------------------

  //RooFitResult* fit_result = simPdf.fitTo(*data,Save(1),Extended(1),Minos(1),Strategy(2));
	//fit_result->Print("v");

  RooMCStudy* mcstudy = new RooMCStudy(simPdf,RooArgList(nbjets,WPCat),Binned(kTRUE),Silence(),Extended(kTRUE),FitOptions(Save(0),Minos(1),Extended(1),Strategy(2),PrintEvalErrors(-1)));

  // Generate and fit 1000 samples of Poisson(nExpected) events
  mcstudy->generateAndFit(NbOfPE,(int)(NbOfBtagWorkingPoint*(Nttlike+Nvlike))) ;


  // P l o t   f i t   r e s u l t s 
  // ---------------------------------------------------------------------
  fout->cd();
  myDir = fout->mkdir("Pulls");

  TCanvas* canvas = 0;
  TH1* hPull = 0;
  TF1* fit = 0;
  TPaveStats* stat = 0;

  // Canva for Ntt pull distribution
  canvas = new TCanvas("myCanva_Nttpull","",600,600);
  canvas->cd();
  hPull = mcstudy->fitParDataSet().createHistogram("Nttpull",pullNbOfBins) ;
  hPull->Draw("E1");
  hPull->SetMarkerStyle(20);
  hPull->Fit("gaus","IEM");
  fit = (TF1*)hPull->FindObject("gaus");
  fit->SetLineWidth(3);
  fit->SetLineColor(kBlue);
  hPull->GetXaxis()->SetTitle("Pull");
  hPull->GetYaxis()->SetTitle("Number of events");
  hPull->GetYaxis()->SetTitleOffset(1.35);
  gPad->Update();
  stat = (TPaveStats*)hPull->FindObject("stats");
  stat->SetX1NDC(0.60);
  stat->SetX2NDC(0.98);
  stat->SetY1NDC(0.81);
  stat->SetY2NDC(0.98);
  stat->Draw("same");
  canvas->Write();
  canvas->Delete();

  // Canva for Nv pull distribution
  canvas = new TCanvas("myCanva_Nvpull","",600,600);
  canvas->cd();
  hPull = mcstudy->fitParDataSet().createHistogram("Nvpull",pullNbOfBins) ;
  hPull->Draw("E1");
  hPull->SetMarkerStyle(20);
  hPull->Fit("gaus","IEM");
  fit = (TF1*)hPull->FindObject("gaus");
  fit->SetLineWidth(3);
  fit->SetLineColor(kBlue);
  hPull->GetXaxis()->SetTitle("Pull");
  hPull->GetYaxis()->SetTitle("Number of events");
  hPull->GetYaxis()->SetTitleOffset(1.35);
  gPad->Update();
  stat = (TPaveStats*)hPull->FindObject("stats");
  stat->SetX1NDC(0.60);
  stat->SetX2NDC(0.98);
  stat->SetY1NDC(0.81);
  stat->SetY2NDC(0.98);
  stat->Draw("same");
  canvas->Write();
  canvas->Delete();

/*  
  RooPlot* frame_Nttpull = mcstudy->plotPull(Ntt,-4,4,100) ;
  frame_Nttpull->SetName("myRooPlot_Nttpull");
  frame_Nttpull->Write();
  RooPlot* frame_Ntt = mcstudy->plotParam(Ntt);//,Binning(100,NbOfTTlike[NbOfJets-3]*(1-0.2),NbOfTTlike[NbOfJets-3]*(1+0.2))) ;
  frame_Ntt->SetName("myRooPlot_Ntt");
  frame_Ntt->Write();

  RooPlot* frame_Nvpull  = mcstudy->plotPull(Nv,-4,4,100) ;
  frame_Nvpull->SetName("myRooPlot_Nvpull");
  frame_Nvpull->Write();
  RooPlot* frame_Nv = mcstudy->plotParam(Nv);//,Binning(100,NbOfVlike[NbOfJets-3]*(1-0.2),NbOfVlike[NbOfJets-3]*(1+0.2))) ;
  frame_Nv->SetName("myRooPlot_Nv");
  frame_Nv->Write();

  RooPlot* frame = 0;

  for(int i=0; i<NbOfBtagWorkingPoint; i++){
  	if(!eb[i]->isConstant()){
		  if(SavePull){
			  frame = mcstudy->plotPull(*eb[i],-4,4,100) ;
			  frame->SetName("myRooPlot_Ebpull_"+wp[i]);
			  frame->Write();
		  }
		  frame = mcstudy->plotParam(*eb[i]);//,Binning(100)) ;
		  frame->SetName("myRooPlot_Eb_"+wp[i]);
		  frame->Write();
  	}
  	if(!eudsc[i]->isConstant()){
		  if(SavePull){
			  frame = mcstudy->plotPull(*eudsc[i],-4,4,100) ;
			  frame->SetName("myRooPlot_Eudscpull_"+wp[i]);
			  frame->Write();
		  }
		  frame = mcstudy->plotParam(*eudsc[i]);//,Binning(100)) ;
		  frame->SetName("myRooPlot_Eudsc_"+wp[i]);
		  frame->Write();
	  }
  	if(!euds[i]->isConstant()){
		  if(SavePull){
			  frame = mcstudy->plotPull(*euds[i],-4,4,100) ;
			  frame->SetName("myRooPlot_Eudspull_"+wp[i]);
			  frame->Write();
		  }
		  frame = mcstudy->plotParam(*euds[i]);//,Binning(100)) ;
		  frame->SetName("myRooPlot_Euds_"+wp[i]);
		  frame->Write();
	  }
  }
*/
  
  // Closing files
  fin->Close();
  fout->Close();

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;

  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;

  return 0;
}

