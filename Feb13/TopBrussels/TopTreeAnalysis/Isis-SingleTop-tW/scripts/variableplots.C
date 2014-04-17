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
#include "inputs.h"
using namespace std;

void variableplots(int mode = 0){
  
  char myTexFile[300];
  sprintf(myTexFile,"tables/variables_%d.tex", mode);
  ofstream salida(myTexFile); 
  
  
  salida << "\\documentclass{cmspaper}" << endl;
  salida << "\\usepackage{graphicx}" << endl;
  salida << "\\begin{document}" << endl;
  salida << endl;
  salida << endl;
  
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0);
  setTDRStyle();
  gROOT->SetBatch(1);
  
  labelcms  = new TPaveText(0.17,0.87,0.17,0.97,"NDCBR");
  labelcms->SetTextAlign(12);
  labelcms->SetTextSize(0.03);
  labelcms->SetFillColor(kWhite);
  labelcms->AddText("CMS Preliminary, #sqrt{s} = 7 TeV, L_{Int} = npi pb^{-1}");
  labelcms->SetBorderSize(0);
  
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
  sprintf(myRootFile,"outputs/var_%d_twdr.root", mode);
  TFile *_file0 = TFile::Open(myRootFile);
  
  sprintf(myRootFile,"outputs/var_%d_tt.root", mode);
  TFile *_file1 = TFile::Open(myRootFile);
  
  TString cutLabel[40] = { "met", "promet", "mll", "njets", "njetsbt", "ptsys", "ht", "ht_nomet", "pt_max", "pt_min", "pt_leading",
			   "oblateness", "sphericity", "aplanarity", "njetW", "sqrts", "deltaphi", "deltaR", "deltaeta", "deltaphiclosemet", 
			   "deltaphiclosejet", "deltaphijetmet", "eta_ptmax","eta_ptmin", "eta_all", "eta_jet", "phi_jet", "phi_met",
			   "total_pt", "total_eta", "njetsextra", "ptjetsextra", "etajetsextra", "mll_lepclosejet",
			   "mll_lepfarjet", "metminusht", "metminuspt" };
  int rebinHisto[40] = {2, 2, 2, 1, 1, 2, 4, 4, 2, 2, 2, 2, 2, 2, 4, 4, 2, 2, 2, 2, 4, 4, 4, 4, 8, 8, 8, 4, 4, 4, 1, 4, 4, 2, 2, 2, 2};
  TString cutTitle[40] = { "Missing E_{T}", "Projected Missing E_{T}", "Inv. Mass", "# of jets", "# of jets (bt)" , "P_{T} system", "H_{T}",
			   "H_{T} without MET", "P_{T}^{max}", "P_{T}^{min}", "P_{T} leading jet", "Oblateness", "Sphericity", "Aplanarity", "njetW", "sqrts",
			   "#Delta#Phi_{ll}", "#DeltaR_{ll}","#Delta#eta_{ll}", "#Delta#Phi_{closest lep, MET}", "#Delta#Phi_{closest lep, jet}", "#Delta#Phi_{MET, jet}",
			   "#eta first lepton", "#eta second lepton", "#eta lepton", "#eta jet", "#phi jet", "#phi MET", "p_{T} (l+l+jet)", "#eta (l+l+jet)",
			   "# of extra jets", "p_{T} of the extra jets", "#eta of the extra jets", "Inv. Mass lepton+jet (close lepton)",
			   "Inv. Mass lepton+jet (far lepton)", "MET - H_{T}^{nomet}", "MET - P_{T}^{system, nomet}"};
  
  TH1F*  h [40];
  TH1F*  htt [40];
  TString modes[3] = { "emu", "mumu", "ee"};
  
  for (int iVariable = 0; iVariable < 37; iVariable++){
    
    cout << "histo_" + cutLabel[iVariable] << endl;
    h[iVariable] = (TH1F*) _file0->Get("histo_" + cutLabel[iVariable]);
    h[iVariable]->Rebin(rebinHisto[iVariable]);
    h[iVariable]->SetLineColor(kBlack);
    h[iVariable]->SetFillColor(kBlue);
    h[iVariable]->SetFillStyle(3003);
    h[iVariable]->SetLineWidth(1);
    
    
    cout << "histo_" + cutLabel[iVariable] << endl;
    htt[iVariable] = (TH1F*) _file1->Get("histo_" + cutLabel[iVariable]);
    htt[iVariable]->Rebin(rebinHisto[iVariable]);
    htt[iVariable]->SetLineColor(kRed);
    htt[iVariable]->SetMarkerColor(kRed);
    htt[iVariable]->SetLineWidth(3);
    
    leg = new TLegend(0.65,0.75,0.90,0.90);
    leg ->SetFillStyle(0);
    leg ->SetFillColor(kWhite);
    leg ->SetBorderSize(0);
    leg ->SetTextSize(0.030);
    leg->AddEntry(h[iVariable], "tW (DR)", "f");
    leg->AddEntry(htt[iVariable], "t#bar{t}", "l");
    
    h[iVariable]->SetNormFactor(1);
    htt[iVariable]->SetNormFactor(1);
    
    double max = TMath::Max(h[iVariable]->GetMaximum(),htt[iVariable]->GetMaximum());
    
    TCanvas *c1 = new TCanvas();
    h[iVariable]->Draw("histo");
    h[iVariable]->SetMaximum(max * 1.07);
    // h[iVariable]->SetMinimum(0);
    
    htt[iVariable]->Draw("histo, sames");
    h[iVariable]->GetXaxis()->SetTitle(cutTitle[iVariable]);
    h[iVariable]->GetYaxis()->SetTitle("Normalized to 1");
    h[iVariable]->GetYaxis()->SetTitleOffset(1.4);
    h[iVariable]->GetYaxis()->CenterTitle(); 
    leg->Draw();
    cout << iVariable << endl;
    c1->SaveAs("variables/plot_" + cutLabel[iVariable]  + "_" + modes[mode] + ".png");
    // h[iVariable]->SetMinimum(1);
    c1->SetLogy();
    c1->SaveAs("variables/plot_" + cutLabel[iVariable] + "_" + modes[mode] + "_log.png");
    
    
       
    
    salida <<" \\includegraphics[width=0.5\\linewidth]{../variables/plot_" << cutLabel[iVariable] <<
      "_" << modes[mode] << ".png} \\includegraphics[width=0.5\\linewidth]{../variables/plot_" <<cutLabel[iVariable] <<
      "_" << modes[mode] << "_log.png} " << endl;
    salida <<  " \\\\  " << endl;
    
    
    
    
  }
  
  salida << endl;
  salida << endl;
  
  
  salida << "\\end{document}" << endl;

}
