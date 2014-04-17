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

void controlplots(int mode = 0){
 
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0);
  setTDRStyle();
  gROOT->SetBatch(1);

  labelcms  = new TPaveText(0.12,0.88,0.5,0.94,"NDCBR");
  labelcms->SetTextAlign(12);
  labelcms->SetTextSize(0.045);
  labelcms->SetFillColor(kWhite);
  labelcms->AddText("CMS Preliminary, #sqrt{s} = 7 TeV");
  labelcms->SetBorderSize(0);
  
  labelcms2  = new TPaveText(0.12,0.85,0.5,0.88,"NDCBR");
  labelcms2->SetTextAlign(12);
  labelcms2->SetTextSize(0.045);
  labelcms2->SetFillColor(kWhite);
  
  if (mode == 0) labelcms2->AddText("4.9 fb^{-1}, e#mu channel  ");
  if (mode == 1) labelcms2->AddText("4.9 fb^{-1}, #mu#mu channel  ");
  if (mode == 2) labelcms2->AddText("4.9 fb^{-1}, ee channel  ");
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
  
  char myRootFile[300];
  double lumi = 0;
  
  if (mode == 0 )        lumi = 4904.338;
  else if ( mode == 1)   lumi = 4919.924;
  else if ( mode == 2)   lumi = 4895.249;

  sprintf(myRootFile,"results/sf_an_%dpb_%d.root", lumi, mode);
  
  
  TFile *_file0 = TFile::Open(myRootFile);
  cout << myRootFile << endl;

  const int nProcess = 8;
  TString processName[nProcess] =  { "twdr", "st", "tt","di", "zjets", "wjets",  "qcd_mu", "data"};
  TString processTitle[nProcess] = { "tW", "t/s-channel", "t#bar{t}", "WW/WZ/ZZ", "Z/#gamma*+jets", "W+jets",  "QCD", "data"};
  Color_t color[nProcess] =        {kWhite, kMagenta-10, kRed+1, kYellow-10,  kAzure-2, kGreen-3, 40, kBlack};
  TString modeString[3] = {"0", "1", "2"};

  TString cutLabel = "R";
 
  
  TH1F*  h [nProcess];
  TH1F*  histo[nProcess];
  THStack* hStack;
  
  leg = new TLegend(0.7,0.7,0.94,0.94);
  leg ->SetFillStyle(1);
  leg ->SetFillColor(kWhite);
  leg ->SetBorderSize(1);
  leg ->SetTextSize(0.03);
  hStack = new THStack();
  hStack2 = new THStack();
  for (int iProcess = 0; iProcess < 8; iProcess++){
    h[iProcess] = (TH1F*) _file0->Get("R_" + processName[iProcess]);
    
    h[iProcess]->SetFillColor(color[iProcess]);
    h[iProcess]->SetLineColor(kBlack);
    
    histo[iProcess] = new TH1F("histo"+processName[iProcess], "", 3, 0, 3);
    histo[iProcess]->SetLineColor(kBlack);
    histo[iProcess]->SetLineWidth(1);
    histo[iProcess]->SetFillColor(color[iProcess]);
    histo[iProcess]->SetBinContent(1, h[iProcess]->GetBinContent(28));
    histo[iProcess]->SetBinContent(2, h[iProcess]->GetBinContent(30));
    histo[iProcess]->SetBinContent(3, h[iProcess]->GetBinContent(31));
    histo[iProcess]->SetBinError(1, h[iProcess]->GetBinError(28));
    histo[iProcess]->SetBinError(2, h[iProcess]->GetBinError(30));
    histo[iProcess]->SetBinError(3, h[iProcess]->GetBinError(31));
    
    if (iProcess == 0){
      h[iProcess]->SetLineColor(kBlack);
      histo[iProcess]->SetLineColor(kBlack);
    }
   
  }  
  
  histo[5]->Add(histo[1]);
  histo[5]->Add(histo[3]);
  histo[5]->Add(histo[6]);
  
  hStack->Add(histo[5]);
  hStack->Add(histo[4]);
  hStack->Add(histo[2]);
  hStack->Add(histo[0]);
  
  if (mode == 0) leg->AddEntry(histo[7], processTitle[7], "p");
  if (mode == 1) leg->AddEntry(histo[7], processTitle[7], "p");
  if (mode == 2) leg->AddEntry(histo[7], processTitle[7], "p");
  
  leg->AddEntry(histo[0], processTitle[0], "f");
  leg->AddEntry(histo[2], processTitle[2], "f");
  leg->AddEntry(histo[4], processTitle[4], "f");
  leg->AddEntry(histo[5], "Other", "f");
 
  
  histo[7]->SetMarkerStyle(20);
  histo[7]->SetMarkerSize(1.2);
  histo[7]->SetMarkerColor(kBlack);
  histo[7]->SetLineColor(kBlack);
  histo[7]->SetLineWidth(1);
  
  double max = TMath::Max(hStack->GetMaximum(), histo[7]->GetMaximum());
  TCanvas *c1 = new TCanvas();
  hStack->Draw("histo");
  hStack->SetMaximum(max * 1.25);
  hStack->SetMinimum(1);
  hStack->GetYaxis()->SetTitle("events / 4.9 fb^{-1}");
  hStack->GetYaxis()->SetTitleOffset(1.4);
  hStack->GetYaxis()->CenterTitle(); 
  
  hStack->GetYaxis()->SetTitleOffset(1.4);
  hStack->GetYaxis()->CenterTitle(); 
  hStack->GetXaxis()->SetBinLabel(1,"1 jet 1 tag");
  hStack->GetXaxis()->SetBinLabel(2,"2 jet 1 tag");
  hStack->GetXaxis()->SetBinLabel(3,"2 jet 2 tag");
  
  histo[7]->Draw("e, sames");
  leg->Draw();
  labelcms->Draw();
  labelcms2->Draw();
  
  c1->SaveAs("plots/sf_control_summer_tt_" + modeString[mode] + "_" + cutLabel + ".png");
  c1->SaveAs("plots/pdf/sf_control_summer_tt_" + modeString[mode] + "_" + cutLabel + ".pdf");
  c1->SetLogy();
  hStack->SetMaximum(max * 10);
  c1->SaveAs("plots/sf_control_summer_tt_" + modeString[mode] + "_" + cutLabel + "_log.png");
  c1->SaveAs("plots/pdf/sf_control_summer_tt_" + modeString[mode] + "_" + cutLabel + "_log.pdf");
  
}
