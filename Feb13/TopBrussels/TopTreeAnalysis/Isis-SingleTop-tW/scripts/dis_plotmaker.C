//Isis Van Parijs


#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TFile.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "setTDRStyle.C"
//#include "../Tools/interface/MultiSamplePlot.h"
//#include "../Tools/interface/PlottingTools.h"
using namespace std;

void dis_plotmaker(int mode = 0){
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0);
  setTDRStyle();
  gROOT->SetBatch(1);   // to show histos on canvas
  
  //Making it pretty
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
  labelcms2->AddText("4.4 fb^{-1}, e#mu channel  ");
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
  double lumi = 4399; // emu channel

    
  sprintf(myRootFile,"results/an_%dpb_%d.root", lumi, mode); // take output from looper
  
  TFile *_file0 = TFile::Open(myRootFile);
  cout << myRootFile << endl;
  
  const int nProcess = 2;
  const int nPlots = 6;
  TString processName[nProcess] =  { "twdr", "tt"};
  TString processTitle[nProcess] = { "tW","t#bar{t}"};
  Color_t color[nProcess] =        { kBlue, kRed};
  

  TString cutLabel[nPlots] =     {  "met_2j2t", "mll_2j2t", "ptsys_2j2t", "ht_2j2t", "pt_leading_2j2t",  "nvertex_2j2t" };
  int rebinHisto[nPlots] =       {  4, 4,  4, 12, 4, 1};
  TString cutTitle[nPlots] =     { "E_{T}^{miss} 2j2t", "Inv. Mass 2j2t", "P_{T} system [GeV] 2j2t", "H_{T} [GeV] 2j2t","P_{T} of the leading jet 2j2t", "# ofvertex 2j2t"};


  TString modeString[1] = {"0"};
  
  TString plotExtension = "discriminatingVar_plot_"; // name of the plots
  TString plotAnalysis = "DiscriminatingVar"; // directory in plots where the plots are saved
  

  
   
  
  
  //Make plots   
  TH1F* histo_tt;
  TH1F* histo_twdr;
   

   
  
   for(int iPlots = 0; iPlots< nPlots; iPlots++)
   {
      
     leg = new TLegend(0.7,0.7,0.94,0.94);
     leg ->SetFillStyle(1001);
     leg ->SetFillColor(kWhite);
     leg ->SetBorderSize(1);
     
     
 
	
	
     histo_tt= (TH1F*) _file0->Get(cutLabel[iPlots]+ "_" + processName[1]);
	       cout << cutLabel[iPlots]+ "_" + processName[1] << endl; 
	       
	       
       histo_tt->Rebin(rebinHisto[iPlots]);
               //histo_tt->SetFillColor(color[1]);
        histo_tt->SetLineColor(color[1]);
        histo_tt->SetLineWidth(2);


       histo_twdr= (TH1F*) _file0->Get(cutLabel[iPlots]+ "_" + processName[0]);
        cout << cutLabel[iPlots]+ "_" + processName[0] << endl; 

       histo_twdr->Rebin(rebinHisto[iPlots]);
              // histo_twdr->SetFillColor(color[0]);
       histo_twdr->SetLineColor(color[0]);
       histo_twdr->SetLineWidth(2);
	      
       leg ->AddEntry(histo_tt, "tt bckgr." , "l");  
       leg ->AddEntry(histo_twdr, "twdr signal", "l");  
	
       double max = TMath::Max(histo_tt->GetMaximum(), histo_twdr->GetMaximum());
    
       TCanvas *c1 = new TCanvas();
       histo_tt->Draw("h");
       histo_tt->SetMaximum(max * 1.2);
       histo_tt->SetMinimum(1);
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
     
       
       histo_twdr->Draw("histo , sames");
       histo_twdr->SetMaximum(max * 1.2);
       histo_twdr->SetMinimum(1);

       
      leg->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + modeString[mode] + "_" + cutLabel[iPlots] + ".png");
     // c1->SaveAs("plots/pdf/" + plotExtension + modeString[mode] + "_" + cutLabel[iPlots] + ".pdf");
      
      int y = 1/max; 
      
      
      TCanvas *c2 = new TCanvas();
       histo_tt->DrawNormalized("h",1);
      // histo_tt->SetMaximum(max * 1.2);
       //histo_tt->SetMinimum(1);
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
    //   histo_tt->GetYaxis()->SetTitle("events / 4.4 fb^{-1}");
    
       
       
       histo_twdr->DrawNormalized("h,sames",1);
    //   histo_twdr->SetMaximum(max * 1.2);
     //  histo_twdr->SetMinimum(1);

       
      leg->Draw();
     // labelcms->Draw();
     // labelcms2->Draw();
    
      c2->SaveAs("plots/"  + plotAnalysis +"/"  + plotExtension + modeString[mode] + "_" + cutLabel[iPlots] + "_normalized" +  ".png");
    //  c2->SaveAs("plots/pdf/" + plotExtension + modeString[mode] + "_" + cutLabel[iPlots] + "_normalized" + ".pdf");
      

      
    } // end plots loop 
    
 
} // end constructor loop
