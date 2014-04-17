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

void controlplotsall(bool dyonly = false, bool ttdy = false){

  //gROOT->SetStyle("Plain");
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
  
  labelcms2->AddText("4.9 fb^{-1}, ee/e#mu/#mu#mu");
  //labelcms2->AddText("2.1 fb^{-1}, ee/e#mu/#mu#mu");
  
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
 
  TString processName[8] =  { "twdr", "st", "tt","di", "zjets", "wjets",  "qcd_mu", "data"};
  TString processTitle[8] = { "tW", "t/s-channel", "t#bar{t}", "WW/WZ/ZZ", "Z/#gamma*+jets", "W+jets",  "QCD", "data"};
  Color_t color[8] =        {kWhite, kMagenta-10, kRed+1, kYellow-10,  kAzure-2, kGreen-3, 40, kBlack};
 
  TString cutLabel = "R";
  TString cutTitle = " ";
  TString nregion = "tt"; 
  
  TH1F*  h [8];
  TH1F*  histo[8];
  TH1F*  h0[8];
  TH1F*  h1[8];
  TH1F*  h2[8];
  TH1F*  histo2[8];
  THStack* hStack;
  THStack* hStack2;
  
  leg = new TLegend(0.7,0.7,0.94,0.94);
  leg ->SetFillStyle(1);
  leg ->SetFillColor(kWhite);
  leg ->SetBorderSize(1);
  leg ->SetTextSize(0.03);
  hStack = new THStack();
  hStack2 = new THStack();
  for (int iProcess = 0; iProcess < 8; iProcess++){
    
    sprintf(myRootFile,"results/sf_an_4904pb_0.root");
    
    TFile *_file0 = TFile::Open(myRootFile);
    h0[iProcess] = (TH1F*) _file0->Get("R_" + processName[iProcess]);
   
    sprintf(myRootFile,"results/sf_an_4919pb_1.root");
    
    TFile *_file1 = TFile::Open(myRootFile);
    h1[iProcess] = (TH1F*) _file1->Get("R_" + processName[iProcess]);
 
    sprintf(myRootFile,"results/sf_an_4895pb_2.root");
    
    TFile *_file2 = TFile::Open(myRootFile);
    h2[iProcess] = (TH1F*) _file2->Get("R_" + processName[iProcess]);
    /*
    if (ttdy){
      
      if (iProcess == 2){
	h0[iProcess]->SetBinContent(2, h0[iProcess]->GetBinContent(2)*1.02);
	h0[iProcess]->SetBinContent(7, h0[iProcess]->GetBinContent(7)*1.1);
	h0[iProcess]->SetBinContent(8, h0[iProcess]->GetBinContent(8)*1.12);
	
	h1[iProcess]->SetBinContent(2, h1[iProcess]->GetBinContent(2)*1.02);
	h1[iProcess]->SetBinContent(7, h1[iProcess]->GetBinContent(7)*1.06);
	h1[iProcess]->SetBinContent(8, h1[iProcess]->GetBinContent(8)*1.1);
	
	h2[iProcess]->SetBinContent(2, h2[iProcess]->GetBinContent(2)*0.98);
	h2[iProcess]->SetBinContent(7, h2[iProcess]->GetBinContent(7)*1.09);
	h2[iProcess]->SetBinContent(8, h2[iProcess]->GetBinContent(8)*1.03);
	
	
      }
      
    }*/

    h[iProcess] =  h0[iProcess];
    h[iProcess]->Add(h1[iProcess]);
    h[iProcess]->Add(h2[iProcess]);
   
    h[iProcess]->SetFillColor(color[iProcess]);
    h[iProcess]->SetLineColor(kBlack);
    h[iProcess]->SetLineWidth(1);
    
    histo[iProcess] = new TH1F("histo"+processName[iProcess], "", 3, 0, 3);
    histo[iProcess]->SetLineColor(kBlack);
    histo[iProcess]->SetFillColor(color[iProcess]);
    histo[iProcess]->SetBinContent(1, h[iProcess]->GetBinContent(28));
    histo[iProcess]->SetBinContent(2, h[iProcess]->GetBinContent(30));
    histo[iProcess]->SetBinContent(3, h[iProcess]->GetBinContent(31));
    histo[iProcess]->SetBinError(1, h[iProcess]->GetBinError(28));
    histo[iProcess]->SetBinError(2, h[iProcess]->GetBinError(30));
    histo[iProcess]->SetBinError(3, h[iProcess]->GetBinError(31));
    /*
    if (dyonly == true){
      if (iProcess == 4){
	histo[iProcess]->SetBinContent(1,76.4);
	histo[iProcess]->SetBinError(1,6.1);
      }
    }
    
    if (ttdy){
      if (iProcess == 4){
	histo[iProcess]->SetBinContent(1,76.4);
	histo[iProcess]->SetBinError(1,6.1);
      }
    }*/
    
    histo2[iProcess] = new TH1F("histo2"+processName[iProcess], "", 5, 0, 5);
    histo2[iProcess]->SetLineColor(kBlack);
    histo2[iProcess]->SetLineWidth(1);
    histo2[iProcess]->SetFillColor(color[iProcess]);
    histo2[iProcess]->SetBinContent(1, h[iProcess]->GetBinContent(2));
    histo2[iProcess]->SetBinContent(2, h[iProcess]->GetBinContent(7));
    histo2[iProcess]->SetBinContent(3, h[iProcess]->GetBinContent(8));
    histo2[iProcess]->SetBinContent(4, h[iProcess]->GetBinContent(12));
    histo2[iProcess]->SetBinContent(5, h[iProcess]->GetBinContent(15));
    histo2[iProcess]->SetBinError(1, h[iProcess]->GetBinError(2));
    histo2[iProcess]->SetBinError(2, h[iProcess]->GetBinError(7));
    histo2[iProcess]->SetBinError(3, h[iProcess]->GetBinError(8));
    histo2[iProcess]->SetBinError(4, h[iProcess]->GetBinError(12));
    histo2[iProcess]->SetBinError(5, h[iProcess]->GetBinError(15));
    
    if (iProcess == 0){
      h[iProcess]->SetLineColor(kBlack);
      histo[iProcess]->SetLineColor(kBlack);
      histo2[iProcess]->SetLineColor(kBlack);
    }
    
  }
  for (int iP = 0; iP < 7; iP ++){
    int  index = 6 - iP;
    // hStack->Add(histo[index]);
    hStack2->Add(histo2[index]);
  }
 
  histo[5]->Add(histo[1]);
  histo[5]->Add(histo[3]);
  histo[5]->Add(histo[6]);  
  
  hStack->Add(histo[5]);
  hStack->Add(histo[4]);
  hStack->Add(histo[2]);
  hStack->Add(histo[0]);
  
  leg->AddEntry(histo[7], processTitle[7], "p");
  
  leg->AddEntry(histo[0], processTitle[0], "f");
  leg->AddEntry(histo[2], processTitle[2], "f");
  leg->AddEntry(histo[4], processTitle[4], "f");
  leg->AddEntry(histo[5], "Other", "f");
     
  
  histo[7]->SetMarkerStyle(20);
  histo[7]->SetMarkerSize(1.2);
  histo[7]->SetLineWidth(1);
  histo[7]->SetMarkerColor(kBlack);
  histo[7]->SetLineColor(kBlack);
  
  double max = TMath::Max(hStack->GetMaximum(), histo[7]->GetMaximum());
  TCanvas *c1 = new TCanvas();
  hStack->Draw("histo");
  hStack->SetMaximum(max * 1.25);
  hStack->SetMinimum(1);
  hStack->GetXaxis()->SetTitle(cutTitle);
  hStack->GetYaxis()->SetTitle("events / 4.9 fb^{-1}");
  //hStack->GetYaxis()->SetTitle("events / 2.1 fb^{-1}");
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
  
  c1->SaveAs("plots/sf_control_summer_" + nregion + "_3_" + cutLabel + ".png");
  c1->SaveAs("plots/pdf/sf_control_summer_" + nregion + "_3_" + cutLabel + ".pdf");
  c1->SetLogy();
  hStack->SetMaximum(max * 20);
  c1->SaveAs("plots/sf_control_summer_" + nregion + "_3_" + cutLabel + "_log.png");  
  c1->SaveAs("plots/pdf/sf_control_summer_" + nregion + "_3_" + cutLabel + "_log.pdf"); 
  
  
}
