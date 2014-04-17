//Rebeca Gonzalez Suarez
//rebeca@cern.ch

#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TFile.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "setTDRStyle.C"
using namespace std;

void plotmaker(int mode = 0){
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0);
  setTDRStyle();
  gROOT->SetBatch(1);
  
  labelcms  = new TPaveText(0.12,0.88,0.5,0.94,"NDCBR");
  labelcms->SetTextAlign(12);
  labelcms->SetTextSize(0.045);
  labelcms->SetFillColor(kWhite);
  labelcms->AddText("CMS Preliminary, #sqrt{s} = 8 TeV");
  labelcms->SetBorderSize(0);
    
  labelcms2  = new TPaveText(0.12,0.85,0.5,0.88,"NDCBR");
  labelcms2->SetTextAlign(12);
  labelcms2->SetTextSize(0.045);
  labelcms2->SetFillColor(kWhite);
  
  if (mode == 0) labelcms2->AddText("5.1 fb^{-1}, e#mu channel  ");
  if (mode == 1) labelcms2->AddText("5.1 fb^{-1}, #mu#mu channel  ");
  if (mode == 2) labelcms2->AddText("5.1 fb^{-1}, ee channel  ");
  
  labelcms2->SetBorderSize(0);
  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(600);
  gStyle->SetCanvasDefW(600);
  gStyle->SetLabelFont(18,"");
  
  gStyle->SetTitleXOffset(1.2);//1.5
  gStyle->SetTitleYOffset(1.2);//1.7
  
  labelKStest = new TPaveText(0.12,0.80,0.6,0.84,"NDCBR");
  labelKStest->SetTextAlign(12);
  labelKStest->SetTextSize(0.035);
  labelKStest->SetFillColor(kWhite);
  labelKStest->SetBorderSize(0);

  
  char myRootFile[300];
  double lumi = 1000;
  
  if (mode == 0 )        lumi = 5085.246; 
  else if ( mode == 1)   lumi = 1000;
  else if ( mode == 2)   lumi = 1000;
  
  
  sprintf(myRootFile,"results/an_%dpb_%d.root", lumi, mode);
  
  TFile *_file0 = TFile::Open(myRootFile);
  cout << myRootFile << endl;
  
  const int nProcess = 8;
  const int nPlots = 29;
  TString processName[nProcess] =  { "twdr", "st", "tt","di", "zjets", "wjets",  "qcd_mu", "data"};
  TString processTitle[nProcess] = { "tW", "t/s-channel", "t#bar{t}", "WW/WZ/ZZ", "Z/#gamma*+jets", "W+jets",  "QCD", "data"};
  Color_t color[nProcess] =        { kWhite, kMagenta-10, kRed+1, kYellow-10,  kAzure-2, kGreen-3, 40, kBlack};
  

  TString cutLabel[nPlots] =     { "cuts", "met", "mll", "njets",  "njetsbt", "ptsys", "ht", "pt_leading", "nvertex", "pt_max","pt_min","nvertex_2lep", "nvertex_final_purw", "nvertex_final_3D",  "eta_jet", "mll_zgamma", "met_zgamma","eta_jet_30","eta_jet_50","eta_jet_70", "eta_jet_90","eta_jet_110","eta_jet1","eta_jet1_30","eta_jet1_50","eta_jet1_70","eta_jet1_90","eta_jet1_110" };
  int rebinHisto[nPlots] =       { 1, 4, 4, 1, 1, 4, 12, 4, 1, 2, 2,  1, 1, 1, 2, 2,2,2,2,2,2,2,2,2,2,2,2,2,2};
  TString cutTitle[nPlots] =     { "Analysis Cut", "E_{T}^{miss}","Inv. Mass", "# of jets", "# of jets(bt)" , "P_{T} system[GeV]","H_{T}[GeV]","P_{T} of the leading jet", "# of vertex", "p_T ofthe first lepton [GeV]", "p_T  of the second lepton [GeV]", "# ofvertex", "# of vertex after cuts PU reweighted", "# of vertex aftercuts 3D PU reweighted", "eta jets", "mll in zgamma control region","MET in zgamma  control region", "eta with pt30", "eta with pt50","eta with pt70", "eta with pt90", "eta with pt110","eta secondleading jet", "eta secondleading jet with pt30", "eta second leading jet with pt50","eta second leading jet with pt70", "eta second leading jet with pt90", "eta second leading jet with pt110"};


  TString modeString[3] = {"0", "1", "2"};
  
  TString plotExtension = "plot_";
  
  TH1F*  h [nPlots][nProcess];
  TH1F*  hextra [nPlots];
  THStack* hStack[nPlots];
  TGraphAsymmErrors *GE[nPlots];
  
  for (const int iVariable = 0; iVariable < nPlots; iVariable++){
    leg = new TLegend(0.7,0.7,0.94,0.94);
    leg ->SetFillStyle(1001);
    leg ->SetFillColor(kWhite);
    leg ->SetBorderSize(1);
    hStack[iVariable] = new THStack(cutLabel[iVariable],cutLabel[iVariable]);
    for (int iProcess = 0; iProcess < nProcess; iProcess++){
      h[iVariable][iProcess] = (TH1F*) _file0->Get(cutLabel[iVariable]+ "_" + processName[iProcess]);
      h[iVariable][iProcess]->Rebin(rebinHisto[iVariable]);
      h[iVariable][iProcess]->SetFillColor(color[iProcess]);
      h[iVariable][iProcess]->SetLineColor(kBlack);
      h[iVariable][iProcess]->SetLineWidth(1);
   
      
    }
    
    h[iVariable][5]->Add(h[iVariable][1]);
    h[iVariable][5]->Add(h[iVariable][3]);
    h[iVariable][5]->Add(h[iVariable][6]);
    
    hStack[iVariable]->Add(h[iVariable][5]);
    hStack[iVariable]->Add(h[iVariable][4]);
    hStack[iVariable]->Add(h[iVariable][2]);
    hStack[iVariable]->Add(h[iVariable][0]);
    
    if (mode == 0) leg->AddEntry(h[iVariable][7],  processTitle[7], "p");
    if (mode == 1) leg->AddEntry(h[iVariable][7],  processTitle[7], "p");
    if (mode == 2) leg->AddEntry(h[iVariable][7], processTitle[7], "p");
    
    leg->AddEntry(h[iVariable][0], processTitle[0], "f");
    leg->AddEntry(h[iVariable][2], processTitle[2], "f");
    leg->AddEntry(h[iVariable][4], processTitle[4], "f");
    leg->AddEntry(h[iVariable][5], "Other", "f");
    
    h[iVariable][7]->SetMarkerStyle(20);
    h[iVariable][7]->SetMarkerSize(1.2);
    h[iVariable][7]->SetLineWidth(1);
    h[iVariable][7]->SetMarkerColor(kBlack);
    h[iVariable][7]->SetLineColor(kBlack);
    
    
    
    double max = TMath::Max(hStack[iVariable]->GetMaximum(), h[iVariable][7]->GetMaximum());
    TCanvas *c1 = new TCanvas();
    hStack[iVariable]->Draw("histo");
    hStack[iVariable]->SetMaximum(max * 1.2);
    hStack[iVariable]->SetMinimum(1);
    hStack[iVariable]->GetXaxis()->SetTitle(cutTitle[iVariable]);
    hStack[iVariable]->GetYaxis()->SetTitle("events / 4.9 fb^{-1}");
    
    hStack[iVariable]->GetYaxis()->CenterTitle(); 
    
    if (iVariable == 0){
      hStack[iVariable]->GetXaxis()->SetBinLabel(2,"Lepton Selection");
      hStack[iVariable]->GetXaxis()->SetBinLabel(3,"m_{ll}");
      hStack[iVariable]->GetXaxis()->SetBinLabel(4,"E_{T}^{miss}");
      hStack[iVariable]->GetXaxis()->SetBinLabel(5,"1 jet");
      hStack[iVariable]->GetXaxis()->SetBinLabel(6,"b-tag");
      hStack[iVariable]->GetXaxis()->SetBinLabel(7,"H_{T}");
      hStack[iVariable]->GetXaxis()->SetRangeUser(1,6);
    }
    
    if (iVariable == 5) hStack[iVariable]->GetYaxis()->SetRangeUser(1,100);
    h[iVariable][7]->Draw("e, sames");

    
    leg->Draw();
    labelcms->Draw();
    labelcms2->Draw();
    
    c1->SaveAs("plots/" + plotExtension + modeString[mode] + "_" + cutLabel[iVariable] + ".png");
    c1->SaveAs("plots/pdf/" + plotExtension + modeString[mode] + "_" + cutLabel[iVariable] + ".pdf");
    
    c1->SetLogy();
    hStack[iVariable]->SetMaximum(max * 10);
    c1->SaveAs("plots/" + plotExtension + modeString[mode] + "_" + cutLabel[iVariable] + "_log.png");
    c1->SaveAs("plots/pdf/" + plotExtension + modeString[mode] + "_" + cutLabel[iVariable] + "_log.pdf");
    
    
    c1->SetLogy(0);
   
    TH1F* hextra[iVariable] = (TH1F*) h[iVariable][5]->Clone();
    hextra[iVariable]->Add(h[iVariable][4]);
    hextra[iVariable]->Add(h[iVariable][2]);
    hextra[iVariable]->Add(h[iVariable][0]);
    TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.5)");
    setex2->Draw();
    TGraphAsymmErrors *GE[iVariable] = (TGraphAsymmErrors*) hextra[iVariable];
   
    GE[iVariable]->SetFillColor(28);
    GE[iVariable]->SetFillStyle(3018);
    GE[iVariable]->SetMarkerSize(0);
    GE[iVariable]->SetLineWidth(0);
    GE[iVariable]->SetLineColor(kWhite);
    leg->AddEntry( GE[iVariable], "uncertainty", "f");
    hStack[iVariable]->SetMaximum(max*1.2);
    GE[iVariable]->Draw("sames, e2");
    TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0)");
    setex1->Draw();
    h[iVariable][7]->Draw("e, sames");
    
    c1->SaveAs("plots/error_" + plotExtension + modeString[mode] + "_" + cutLabel[iVariable] + ".png");
    c1->SaveAs("plots/pdf/error_" + plotExtension + modeString[mode] + "_" + cutLabel[iVariable] + ".pdf");
     
    c1->SetLogy();
    hStack[iVariable]->SetMaximum(max * 10);
    c1->SaveAs("plots/error_" + plotExtension+ modeString[mode] + "_" + cutLabel[iVariable] + "_log.png");
    c1->SaveAs("plots/pdf/error_" + plotExtension + modeString[mode] + "_" + cutLabel[iVariable] + "_log.pdf");
    
    
  }
  
}
