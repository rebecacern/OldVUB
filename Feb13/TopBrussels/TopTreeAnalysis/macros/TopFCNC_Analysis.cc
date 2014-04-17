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
  cout << " Beginning of the program for the FCNC analysis ! " << endl;
  cout << "*************************************************************" << endl;

  //SetStyle if needed
  //setTDRStyle();
  setGregStyle();
  //setMyStyle();


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////// Configuration ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //Output ROOT file
  string channelpostfix = "_diMuon";
  string postfix = "_Analysis";
  string rootFileName ("TopFCNC"+postfix+channelpostfix+".root");

  TFile      *fin   =  TFile::Open("TopFCNC_EventSelection_diMu.root");
  TFile      *fout  = new TFile (rootFileName.c_str(), "RECREATE");
  TDirectory *myDir = 0;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////// Histograms /////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////
  ////////////////// 1D histograms  //////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  myDir = fin->GetDirectory("MultiSamplePlot_Mzq_mm_ch");
  TH1F* Mzq_mm_ch_ttbar_fcnc = (TH1F*)myDir->Get("Mzq_mm_ch_ttbar_fcnc");
  TH1F* Mzq_mm_ch_ttjets     = (TH1F*)myDir->Get("Mzq_mm_ch_ttjets");
  TH1F* Mzq_mm_ch_zjets      = (TH1F*)myDir->Get("Mzq_mm_ch_zjets");

  myDir = fin->GetDirectory("MultiSamplePlot_MET_mm_ch");
  TH1F* MET_mm_ch_ttbar_fcnc = (TH1F*)myDir->Get("MET_mm_ch_ttbar_fcnc");
  TH1F* MET_mm_ch_ttjets     = (TH1F*)myDir->Get("MET_mm_ch_ttjets");
  TH1F* MET_mm_ch_zjets      = (TH1F*)myDir->Get("MET_mm_ch_zjets");

  cout << " - Declared histograms ..." <<  endl;
	
  ////////////////////////////////////////////////////////////////////
  ////////////////// Plots  //////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  string pathPNG = "TopFCNC_Analysis"+postfix+channelpostfix;
  pathPNG += "_MSPlots/"; 	
//  pathPNG = pathPNG +"/"; 	
  mkdir(pathPNG.c_str(),0777);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////// Analysis //////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  float N_fcnc   = Mzq_mm_ch_ttbar_fcnc->Integral() ;
  float N_ttjets = Mzq_mm_ch_ttjets->Integral() ;
  float N_zjets  = Mzq_mm_ch_zjets->Integral() ;

  int RebinFact = 5;
  int InterpolFact = 5;

  Mzq_mm_ch_ttbar_fcnc->Rebin(RebinFact);
  Mzq_mm_ch_ttjets->Rebin(RebinFact);
  Mzq_mm_ch_zjets->Rebin(RebinFact);

  MET_mm_ch_ttbar_fcnc->Rebin(RebinFact);
  MET_mm_ch_ttjets->Rebin(RebinFact);
  MET_mm_ch_zjets->Rebin(RebinFact);

  // C r e a t e   p d f ' s   a n d   d a t a s e t s
  // -------------------------------------------------------

  // Declare variables
  RooRealVar mzq("mzq","m_{Zq}",50,300,"[GeV/c^2]") ;
  RooRealVar MET("MET","\\slashE_{T}",50,300,"[GeV]") ;
  RooRealVar n_fcnc("n_fcnc","number of tt+jets events (FCNC)",N_fcnc,0.,2*N_fcnc) ;
  RooRealVar n_tt("n_tt","number of tt+jets events",N_ttjets,0.,2*N_ttjets) ;
  RooRealVar n_z("n_z","number of Z+jets events",N_zjets,0,2*N_zjets) ;

  RooRealVar m0("m0","",170,0,500);
  RooRealVar sigmaCB("sigmaCB","",30,0,100);
  RooRealVar alpha("alpha","",-0.4,-1.,1.);
  RooRealVar n("n","",5,0,20);
  RooRealVar mean("mean","",170,120,200);
  RooRealVar sigma("sigma","",10,0,50);

  RooCBShape CB("CB","CB", mzq, m0, sigmaCB, alpha, n);
  RooGaussian Gaus("Gaus", "Gaus", mzq, mean, sigma);
  RooRealVar frac("frac","frac",0.5,0.,1.) ;
  RooAddPdf sig("sig","Signal",RooArgList(CB,Gaus),frac) ;

  cout << "/****************************************************/" << endl;
  cout << "/* Initial numbers :                                */" << endl;
  cout << "/* - n(FCNC)    :  " << N_fcnc <<   "                         */" << endl;
  cout << "/* - n(tt+jets) :  " << N_ttjets << "                         */" << endl;
  cout << "/* - n(z+jets)  :  " << N_zjets <<  "                         */" << endl;
  cout << "/****************************************************/" << endl;
  
  RooDataHist histmc_mzq_mm_ch_ttbar_fcnc("histmc_mzq_mm_ch_ttbar_fcnc", "", RooArgList(mzq), Mzq_mm_ch_ttbar_fcnc) ;
  RooDataHist histmc_mzq_mm_ch_ttjets    ("histmc_mzq_mm_ch_ttjets",     "", RooArgList(mzq), Mzq_mm_ch_ttjets) ;
  RooDataHist histmc_mzq_mm_ch_zjets     ("histmc_mzq_mm_ch_zjets",      "", RooArgList(mzq), Mzq_mm_ch_zjets) ;

  RooDataHist histmc_MET_mm_ch_ttbar_fcnc("histmc_MET_mm_ch_ttbar_fcnc", "", RooArgList(MET), MET_mm_ch_ttbar_fcnc) ;
  RooDataHist histmc_MET_mm_ch_ttjets    ("histmc_MET_mm_ch_ttjets",     "", RooArgList(MET), MET_mm_ch_ttjets) ;
  RooDataHist histmc_MET_mm_ch_zjets     ("histmc_MET_mm_ch_zjets",      "", RooArgList(MET), MET_mm_ch_zjets) ;

  // Represent data in dh as pdf in mzq
  RooHistPdf histpdf_mzq_mm_ch_ttbar_fcnc("histpdf_mzq_mm_ch_ttbar_fcnc","",mzq,histmc_mzq_mm_ch_ttbar_fcnc,InterpolFact) ;
  RooHistPdf histpdf_mzq_mm_ch_ttjets    ("histpdf_mzq_mm_ch_ttjets",    "",mzq,histmc_mzq_mm_ch_ttjets,InterpolFact) ;
  RooHistPdf histpdf_mzq_mm_ch_zjets     ("histpdf_mzq_mm_ch_zjets",     "",mzq,histmc_mzq_mm_ch_zjets,InterpolFact) ;

  // Represent data in dh as pdf in MET
  RooHistPdf histpdf_MET_mm_ch_ttbar_fcnc("histpdf_MET_mm_ch_ttbar_fcnc","",MET,histmc_MET_mm_ch_ttbar_fcnc,InterpolFact) ;
  RooHistPdf histpdf_MET_mm_ch_ttjets    ("histpdf_MET_mm_ch_ttjets",    "",MET,histmc_MET_mm_ch_ttjets,InterpolFact) ;
  RooHistPdf histpdf_MET_mm_ch_zjets     ("histpdf_MET_mm_ch_zjets",     "",MET,histmc_MET_mm_ch_zjets,InterpolFact) ;

  RooAddPdf model_mzq("model_mzq","model_mzq",RooArgList(histpdf_mzq_mm_ch_ttbar_fcnc,histpdf_mzq_mm_ch_ttjets,histpdf_mzq_mm_ch_zjets),RooArgList(n_fcnc,n_tt,n_z));
  RooAddPdf model_mzq_backgd("model_mzq_backgd","model_mzq_backgd",RooArgList(histpdf_mzq_mm_ch_ttjets,histpdf_mzq_mm_ch_zjets),RooArgList(n_tt,n_z));
  model_mzq.Print();
  RooFormulaVar nb("nb","n_tt+n_z",RooArgList(n_tt,n_z));
  RooExtendPdf model_mzq_backgd_ext("model_mzq_backgd_ext","B-only model",model_mzq_backgd,nb);

  RooAddPdf model_MET("model_MET","model_MET",RooArgList(histpdf_MET_mm_ch_ttbar_fcnc,histpdf_MET_mm_ch_ttjets,histpdf_MET_mm_ch_zjets),RooArgList(n_fcnc,n_tt,n_z));
  model_MET.Print();
/*
  // Define category to distinguish physics and control samples events
  RooCategory vars("vars","vars") ;
  vars.defineType("MZQ") ;
  vars.defineType("MissET") ;

  // C o n s t r u c t   a   s i m u l t a n e o u s   p d f
  // -----------------------------------------------------------------------------------

  // Construct a simultaneous pdf using category sample as index
  RooSimultaneous simPdf("simPdf","simultaneous pdf",vars) ;

  // Associate model with the physics state and model_ctl with the control state
  simPdf.addPdf(model_mzq,"MZQ") ;
  simPdf.addPdf(model_MET,"MissET") ;

*/
  // S a m p l e ,   f i t   e x t e n d e d   m o d e l 
  // ---------------------------------------------------------------------

  // Generate a data sample of expected number events in X from model
  // = model.expectedEvents() = nsig+nbkg
  RooDataSet *data_mzq = model_mzq.generate(mzq) ; data_mzq->Print() ; 
  RooDataSet *data_MET = model_MET.generate(MET) ;
/*
  RooDataSet *data = (RooDataSet*)data_mzq->Clone();

  vars.setLabel("MZQ");
  data->addColumn(vars);
  vars.setLabel("MissET");
  data_MET->addColumn(vars);

  data->append(*data_MET);
  // Construct combined dataset
  //RooDataSet data_comb("data_comb","combined data",x,Index(sample),Import("physics",*data),Import("control",*data_ctl)) ;

  // P e r f o r m   a   s i m u l t a n e o u s   f i t
  // ---------------------------------------------------

  // Perform simultaneous fit of model to data and model_ctl to data_ctl
  simPdf.fitTo(*data) ;
*/
  sig.fitTo(histmc_mzq_mm_ch_ttbar_fcnc);
/*
  // Fit model to data, extended ML term automatically included
  RooAbsReal* nll1 = model_mzq.createNLL(*data_mzq) ;
  RooAbsReal* nll2 = model_MET.createNLL(*data_MET) ;
  RooAddition nllsum("nllsum","nllsum",RooArgSet(*nll1,*nll2)) ;
  RooMinuit m(nllsum) ; // etc   
  m.migrad() ;
  m.hesse() ; 
  //model_mzq.fitTo(*data_mzq) ; 
  //model_MET.fitTo(*data_MET) ; 
*/
//  RooMCStudy* mcstudy = new RooMCStudy(model,mzq,Binned(kTRUE),Silence(),Extended(kTRUE),FitOptions(Save(0),Minos(1),Extended(1),Strategy(2),PrintEvalErrors(-1)));
//  mcstudy->generateAndFit(1000) ;

  // prepare the calculator
  HybridCalculatorOriginal myhc(*data_mzq, model_mzq, model_mzq_backgd_ext,0,0);
  myhc.SetTestStatistic(2);
  myhc.SetNumberOfToys(10);
  myhc.UseNuisance(false);                            

  // run the hypothesis-test invertion
  HypoTestInverterOriginal myInverter(myhc,n_fcnc);
  myInverter.SetTestSize(0.10);
  myInverter.UseCLs(true);
  // myInverter.RunFixedScan(5,1,6);
  // scan for a 95% UL
  myInverter.RunAutoScan(3.,5,myInverter.Size()/2,0.005);  
  // run an alternative autoscan algorithm 
  // myInverter.RunAutoScan(1,6,myInverter.Size()/2,0.005,1);  
  //myInverter.RunOnePoint(3.9);


  HypoTestInverterResult* results = myInverter.GetInterval();

  HypoTestInverterPlot myInverterPlot("myInverterPlot","",results);
/*
  RooWorkspace* w = new RooWorkspace();
  ModelConfig modelConfig("FCNC",w);
  modelConfig.SetPdf(simPdf);
  modelConfig.SetParametersOfInterest(RooArgSet(n_fcnc));
  RooDataSet ValToScan("ValToScan","",n_fcnc);
  n_fcnc = 0;
  ValToScan.add(n_fcnc);
  n_fcnc = 1;
  ValToScan.add(n_fcnc);
  n_fcnc = 2;
  ValToScan.add(n_fcnc);
  //modelConfig.SetObservables(RooArgSet(mzq,MET));
  //modelConfig.SetNuisanceParameters(RooArgSet(n_tt,n_z));
  w->Print();

  //////// show use of Feldman-Cousins
  RooStats::FeldmanCousins fc(*data,modelConfig);
  fc.SetTestSize(.05); // set size of test
  fc.UseAdaptiveSampling(true);
  fc.FluctuateNumDataEntries(false); // number counting analysis: dataset always has 1 entry with N events observed
  //fc.SetNBins(20); // number of points to test per parameter
  fc.SetParameterPointsToTest(ValToScan);
  // use the Feldman-Cousins tool
  PointSetInterval* interval = (PointSetInterval*)fc.GetInterval();
  
  std::cout << "is this point in the interval? " << 
    interval->IsInInterval(RooArgSet(n_fcnc)) << std::endl;

  std::cout << "interval is ["<<  
    interval->LowerLimit(n_fcnc)  << ", "  <<
    interval->UpperLimit(n_fcnc) << "]" << endl;

*/
  // P l o t   p d f ' s   a n d   f i t   r e s u l t s 
  // ---------------------------------------------------------------------

  // Plot data and histogram pdf overlaid
  RooPlot* frame_mzq_pdf_1 = mzq.frame(Title("Histogram pdf (FCNC)"),Bins(50)) ;
  histmc_mzq_mm_ch_ttbar_fcnc.plotOn(frame_mzq_pdf_1) ;
  histpdf_mzq_mm_ch_ttbar_fcnc.plotOn(frame_mzq_pdf_1) ;
  sig.plotOn(frame_mzq_pdf_1,LineStyle(kDotted)) ;

  RooPlot* frame_mzq_pdf_2 = mzq.frame(Title("Histogram pdf (tt+jets)"),Bins(50)) ;
  histmc_mzq_mm_ch_ttjets.plotOn(frame_mzq_pdf_2) ;
  histpdf_mzq_mm_ch_ttjets.plotOn(frame_mzq_pdf_2) ;

  RooPlot* frame_mzq_pdf_3 = mzq.frame(Title("Histogram pdf (Z+jets)"),Bins(50)) ;
  histmc_mzq_mm_ch_zjets.plotOn(frame_mzq_pdf_3) ;
  histpdf_mzq_mm_ch_zjets.plotOn(frame_mzq_pdf_3) ;  

  RooPlot* frame_MET_pdf_1 = MET.frame(Title("Histogram pdf (FCNC)"),Bins(100)) ;
  histmc_MET_mm_ch_ttbar_fcnc.plotOn(frame_MET_pdf_1) ;
  histpdf_MET_mm_ch_ttbar_fcnc.plotOn(frame_MET_pdf_1) ;  

  RooPlot* frame_MET_pdf_2 = MET.frame(Title("Histogram pdf (tt+jets)"),Bins(100)) ;
  histmc_MET_mm_ch_ttjets.plotOn(frame_MET_pdf_2) ;
  histpdf_MET_mm_ch_ttjets.plotOn(frame_MET_pdf_2) ;

  RooPlot* frame_MET_pdf_3 = MET.frame(Title("Histogram pdf (Z+jets)"),Bins(100)) ;
  histmc_MET_mm_ch_zjets.plotOn(frame_MET_pdf_3) ;
  histpdf_MET_mm_ch_zjets.plotOn(frame_MET_pdf_3) ;  

  // Plot data and PDF overlaid, use expected number of events for p.d.f projection normalization
  // rather than observed number of events (==data->numEntries())
  RooPlot* frame_mzq_fit = mzq.frame(Title("Extended ML fit")) ;
  data_mzq->plotOn(frame_mzq_fit) ;
  model_mzq.plotOn(frame_mzq_fit,Normalization(1.0,RooAbsReal::RelativeExpected)) ;
  //simPdf.plotOn(frame_mzq_fit,Slice(vars,"mzq"),ProjWData(vars,*data)) ;
  // Overlay the background component of model with a dashed line
  model_mzq.plotOn(frame_mzq_fit,Components(RooArgSet(histpdf_mzq_mm_ch_ttjets,histpdf_mzq_mm_ch_zjets)),LineStyle(kDashed),Normalization(1.0,RooAbsReal::RelativeExpected)) ;
  //simPdf.plotOn(frame_mzq_fit,Slice(vars,"mzq"),Components(RooArgSet(histpdf_mzq_mm_ch_ttjets,histpdf_mzq_mm_ch_zjets)),ProjWData(vars,*data),LineStyle(kDashed)) ; 
  // Overlay the background+sig2 components of model with a dotted line
  //model.plotOn(xframe,Components(RooArgSet(bkg,sig2)),LineStyle(kDotted),Normalization(1.0,RooAbsReal::RelativeExpected)) ;

  RooPlot* frame_MET_fit = MET.frame(Title("Extended ML fit")) ;
  data_MET->plotOn(frame_MET_fit) ;
  model_MET.plotOn(frame_MET_fit,Normalization(1.0,RooAbsReal::RelativeExpected)) ;
  // Overlay the background component of model with a dashed line
  model_MET.plotOn(frame_MET_fit,Components(RooArgSet(histpdf_MET_mm_ch_ttjets,histpdf_MET_mm_ch_zjets)),LineStyle(kDashed),Normalization(1.0,RooAbsReal::RelativeExpected)) ;

  // Save histograms, canvas, etc.
  TCanvas* canva_mzq = new TCanvas("canva_mzq","",600,600) ;
  canva_mzq->Divide(2,2);
  canva_mzq->cd(1) ; gPad->SetLeftMargin(0.15) ; frame_mzq_pdf_1->GetYaxis()->SetTitleOffset(1.4) ; frame_mzq_pdf_1->Draw() ;
  canva_mzq->cd(2) ; gPad->SetLeftMargin(0.15) ; frame_mzq_pdf_2->GetYaxis()->SetTitleOffset(1.4) ; frame_mzq_pdf_2->Draw() ;
  canva_mzq->cd(3) ; gPad->SetLeftMargin(0.15) ; frame_mzq_pdf_3->GetYaxis()->SetTitleOffset(1.4) ; frame_mzq_pdf_3->Draw() ;
  canva_mzq->cd(4) ; gPad->SetLeftMargin(0.15) ; frame_mzq_fit->GetYaxis()->SetTitleOffset(1.4)   ; frame_mzq_fit->Draw() ;

  TCanvas* canva_MET = new TCanvas("canva_MET","",600,600) ;
  canva_MET->Divide(2,2);
  canva_MET->cd(1) ; gPad->SetLeftMargin(0.15) ; frame_MET_pdf_1->GetYaxis()->SetTitleOffset(1.4) ; frame_MET_pdf_1->Draw() ;
  canva_MET->cd(2) ; gPad->SetLeftMargin(0.15) ; frame_MET_pdf_2->GetYaxis()->SetTitleOffset(1.4) ; frame_MET_pdf_2->Draw() ;
  canva_MET->cd(3) ; gPad->SetLeftMargin(0.15) ; frame_MET_pdf_3->GetYaxis()->SetTitleOffset(1.4) ; frame_MET_pdf_3->Draw() ;
  canva_MET->cd(4) ; gPad->SetLeftMargin(0.15) ; frame_MET_fit->GetYaxis()->SetTitleOffset(1.4)   ; frame_MET_fit->Draw() ;

  canva_mzq->SaveAs("canva_mzq.eps");
  canva_MET->SaveAs("canva_MET.eps");

  TCanvas* intervalCanvas = new TCanvas("intervalCanvas","",600,600) ;
  TGraphErrors* gr1 = myInverterPlot.MakePlot();
  gr1->Draw("ALP");

  double ulError = results->UpperLimitEstimatedError();
  double upperLimit = results->UpperLimit();
  std::cout << "The computed upper limit is: " << upperLimit << std::endl;
  std::cout << "an estimated error on this upper limit is: " << ulError << std::endl;

  intervalCanvas->SaveAs("intervalCanvas.eps");
/*
  // make a canvas for plots
  TCanvas* intervalCanvas =  new TCanvas("intervalCanvas");

  RooDataHist* parameterScan = (RooDataHist*) fc.GetPointsToScan();
  TH1F* hist = (TH1F*) parameterScan->createHistogram("n_fcnc",30);
  hist->Draw();

  RooArgSet* tmpPoint;
  // loop over points to test
  for(Int_t i=0; i<parameterScan->numEntries(); ++i){
    //    cout << "on parameter point " << i << " out of " << parameterScan->numEntries() << endl;
     // get a parameter point from the list of points to test.
    tmpPoint = (RooArgSet*) parameterScan->get(i)->clone("temp");

    TMarker* mark = new TMarker(tmpPoint->getRealValue("n_fcnc"), 1, 25);
    if (interval->IsInInterval( *tmpPoint ) ) 
      mark->SetMarkerColor(kBlue);
    else
      mark->SetMarkerColor(kRed);

    mark->Draw("s");
    //delete tmpPoint;
    //    delete mark;
  }
  intervalCanvas->SaveAs("intervalCanvas.eps");
*/
  fout->cd();
  canva_mzq->Write();
  canva_MET->Write();
  intervalCanvas->Write();

  // Closing files
  fin->Close();
  fout->Close();

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;

  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;

  return 0;
}

