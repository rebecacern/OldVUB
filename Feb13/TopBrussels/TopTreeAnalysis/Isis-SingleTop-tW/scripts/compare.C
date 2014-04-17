#include "setTDRStyle.C"

void compare(){

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0);
  setTDRStyle();
  gROOT->SetBatch(1);    // om de histo's te tonen op canvas
  
  
  //______________________________ all ___________________________________________________________________

  eta_jet_twdr->SetLineColor(kBlue);
  eta_jet_tt->SetLineColor(kRed);
  
  eta_jet_twdr->SetLineWidth(2);
  eta_jet_tt->SetLineWidth(2);

 // Rebin  // om aantal bins te kiezen 
  eta_jet_twdr->Rebin(2);
  eta_jet_tt->Rebin(2);



 
  
  leg = new TLegend(0.7,0.7,0.94,0.94);  // maak legende
  leg ->SetFillStyle(1001);
  leg ->SetFillColor(kWhite);
  leg ->SetBorderSize(1);
  leg ->AddEntry(eta_jet_twdr, "tW signal", "l");  
  leg ->AddEntry(eta_jet_tt, "tt bckgr.", "l");

  TCanvas *c41 = new TCanvas();


  
  eta_jet_twdr->Draw("h");
  eta_jet_tt->Draw("h, sames");
  
  leg->Draw();
  
 // c41->SaveAs("plots/compairison/compare_tw_tt_eta_jet.png");


  TCanvas *c1 = new TCanvas();
  // This is only for visualization
  eta_jet_twdr->SetNormFactor(1);   // Zet alles genormaliseerd tot 1
  eta_jet_tt->SetNormFactor(1);

  
  eta_jet_twdr->Draw("h");
  eta_jet_tt->Draw("h, sames");
  
  leg->Draw();
  
  c1->SaveAs("plots/compairison/compare_tw_tt_eta_jet_norm.png");
  
  



//______________________________ 30 ___________________________________________________________________

  eta_jet_30_twdr->SetLineColor(kBlue);
  eta_jet_30_tt->SetLineColor(kRed);
  
  eta_jet_30_twdr->SetLineWidth(2);
  eta_jet_30_tt->SetLineWidth(2);
  
 
 // Rebin  // om aantal bins te kiezen 
  eta_jet_30_twdr->Rebin(2);
  eta_jet_30_tt->Rebin(2);
 
  
 // TCanvas *c1 = new TCanvas();
 // pt_leading_twdr->Draw("h");  // teken in histo mode
//  pt_leading_tt->Draw("h, sames");   // teken in histogram mode op hetzelfde canvas als de vorige
  
  leg = new TLegend(0.7,0.7,0.94,0.94);  // maak legende
  leg ->SetFillStyle(1001);
  leg ->SetFillColor(kWhite);
  leg ->SetBorderSize(1);
  leg ->AddEntry(eta_jet_30_twdr, "tW signal", "l");  
  leg ->AddEntry(eta_jet_30_tt, "tt bckgr.", "l");
//  leg->Draw();
  
//  c1->SaveAs("plots/compairison/compare_tw_tt_pt_leading.png");  // sla histo op als png
  
  TCanvas *c2 = new TCanvas();
  // This is only for visualization
  eta_jet_30_twdr->SetNormFactor(1);   // Zet alles genormaliseerd tot 1
  eta_jet_30_tt->SetNormFactor(1);
  
  /*
  //Scale 
  pt_leading_twdr->Scale(1); //-> Means that all the entries are multiplied by 1
  pt_leading_tt->Scale(1);
  */
  
  eta_jet_30_twdr->Draw("h");
  eta_jet_30_tt->Draw("h, sames");
  
  leg->Draw();
  
  c2->SaveAs("plots/compairison/compare_tw_tt_eta_jet_30.png");

//______________________________ 50 ___________________________________________________________________

  eta_jet_50_twdr->SetLineColor(kBlue);
  eta_jet_50_tt->SetLineColor(kRed);
  
  eta_jet_50_twdr->SetLineWidth(2);
  eta_jet_50_tt->SetLineWidth(2);
 
  
   // Rebin  // om aantal bins te kiezen 
  eta_jet_50_twdr->Rebin(2);
  eta_jet_50_tt->Rebin(2);
  
  
  leg = new TLegend(0.7,0.7,0.94,0.94);  // maak legende
  leg ->SetFillStyle(1001);
  leg ->SetFillColor(kWhite);
  leg ->SetBorderSize(1);
  leg ->AddEntry(eta_jet_50_twdr, "tW signal", "l");  
  leg ->AddEntry(eta_jet_50_tt, "tt bckgr.", "l");

  TCanvas *c5 = new TCanvas();
  // This is only for visualization
  eta_jet_50_twdr->SetNormFactor(1);   // Zet alles genormaliseerd tot 1
  eta_jet_50_tt->SetNormFactor(1);

  
  eta_jet_50_twdr->Draw("h");
  eta_jet_50_tt->Draw("h, sames");
  
  leg->Draw();
  
  c5->SaveAs("plots/compairison/compare_tw_tt_eta_jet_50.png");
 // ______________________________ 70 ___________________________________________________________________

  eta_jet_70_twdr->SetLineColor(kBlue);
  eta_jet_70_tt->SetLineColor(kRed);
  
  eta_jet_70_twdr->SetLineWidth(2);
  eta_jet_70_tt->SetLineWidth(2);
 
   // Rebin  // om aantal bins te kiezen 
  eta_jet_70_twdr->Rebin(2);
  eta_jet_70_tt->Rebin(2);
  leg = new TLegend(0.7,0.7,0.94,0.94);  // maak legende
  leg ->SetFillStyle(1001);
  leg ->SetFillColor(kWhite);
  leg ->SetBorderSize(1);
  leg ->AddEntry(eta_jet_70_twdr, "tW signal", "l");  
  leg ->AddEntry(eta_jet_70_tt, "tt bckgr.", "l");

  TCanvas *c7 = new TCanvas();
  // This is only for visualization
  eta_jet_70_twdr->SetNormFactor(1);   // Zet alles genormaliseerd tot 1
  eta_jet_70_tt->SetNormFactor(1);

  
  eta_jet_70_twdr->Draw("h");
  eta_jet_70_tt->Draw("h, sames");
  
  leg->Draw();
  
  c7->SaveAs("plots/compairison/compare_tw_tt_eta_jet_70.png");
 //   ______________________________ 90 ___________________________________________________________________

  eta_jet_90_twdr->SetLineColor(kBlue);
  eta_jet_90_tt->SetLineColor(kRed);
  
  eta_jet_90_twdr->SetLineWidth(2);
  eta_jet_90_tt->SetLineWidth(2);
 
  // Rebin  // om aantal bins te kiezen 
  eta_jet_90_twdr->Rebin(2);
  eta_jet_90_tt->Rebin(2);
  
  leg = new TLegend(0.7,0.7,0.94,0.94);  // maak legende
  leg ->SetFillStyle(1001);
  leg ->SetFillColor(kWhite);
  leg ->SetBorderSize(1);
  leg ->AddEntry(eta_jet_90_twdr, "tW signal", "l");  
  leg ->AddEntry(eta_jet_90_tt, "tt bckgr.", "l");

  TCanvas *c9 = new TCanvas();
  // This is only for visualization
  eta_jet_90_twdr->SetNormFactor(1);   // Zet alles genormaliseerd tot 1
  eta_jet_90_tt->SetNormFactor(1);

  
  eta_jet_90_twdr->Draw("h");
  eta_jet_90_tt->Draw("h, sames");
  
  leg->Draw();
  
  c9->SaveAs("plots/compairison/compare_tw_tt_eta_jet_90.png");
 //   ______________________________ 110 ___________________________________________________________________

  eta_jet_110_twdr->SetLineColor(kBlue);
  eta_jet_110_tt->SetLineColor(kRed);
  
  eta_jet_110_twdr->SetLineWidth(2);
  eta_jet_110_tt->SetLineWidth(2);
 
  // Rebin  // om aantal bins te kiezen 
  eta_jet_110_twdr->Rebin(2);
  eta_jet_110_tt->Rebin(2);
  
  leg = new TLegend(0.7,0.7,0.94,0.94);  // maak legende
  leg ->SetFillStyle(1001);
  leg ->SetFillColor(kWhite);
  leg ->SetBorderSize(1);
  leg ->AddEntry(eta_jet_110_twdr, "tW signal", "l");  
  leg ->AddEntry(eta_jet_110_tt, "tt bckgr.", "l");

  TCanvas *c11 = new TCanvas();
  // This is only for visualization
  eta_jet_110_twdr->SetNormFactor(1);   // Zet alles genormaliseerd tot 1
  eta_jet_110_tt->SetNormFactor(1);

  
  eta_jet_110_twdr->Draw("h");
  eta_jet_110_tt->Draw("h, sames");
  
  leg->Draw();
  
  c11->SaveAs("plots/compairison/compare_tw_tt_eta_jet_110.png");

}
