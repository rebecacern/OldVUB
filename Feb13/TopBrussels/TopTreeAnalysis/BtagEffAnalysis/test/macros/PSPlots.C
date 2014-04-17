{

  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.04);
  latex->SetTextAlign(31); // align right

  TFile* f = new TFile("../TTrees/DATA_6/NTupleAnalyzed_pseudoExp.root");

  ////////////////////////////////////////////////////////////
  // BEFF

  TH1D* histo = (TH1D*) f->Get("Efficiency_PseudoExps_WP_3.3");

  histo->Draw();

  //gPad->SetGrid();

  histo->SetTitle("");
  histo->GetXaxis()->SetTitle("#hat{#epsilon}_{b}");
  histo->GetYaxis()->SetTitle("pseudo experiments");

  c1->Update();
                
  TPaveStats* sb = (TPaveStats*)histo->GetListOfFunctions()->FindObject("stats");
  sb->SetX1NDC(1);
  sb->SetX2NDC(1);
  sb->SetY1NDC(1);
  sb->SetY2NDC(1);

  latex->DrawLatex(0.39, 0.865, "CMS Simulation");

  c1->SaveAs("PEXP_beff.pdf");

  histo = (TH1D*) f->Get("Pull_WP_3.3");

  histo->Draw();

  //gPad->SetGrid();

  histo->Fit("gaus");

  c1->Update();
                
  sb = (TPaveStats*)histo->GetListOfFunctions()->FindObject("stats");
  sb->SetOptStat(0);

  histo->SetTitle("");
  histo->GetXaxis()->SetTitle("Pull");
  histo->GetYaxis()->SetTitle("pseudo experiments");

  latex->DrawLatex(0.39, 0.865, "CMS Simulation");

  c1->SaveAs("PEXP_beff_pull.pdf");

  ////////////////////////////////////////////////////////////
  // XS

  histo = (TH1D*) f->Get("XS_FDATA_nDisc_0_WP_3.3");

  histo->Draw();

//gPad->SetGrid();

  histo->SetTitle("");
  histo->GetXaxis()->SetTitle("#hat{#sigma}_{t#bar{t}}");
  histo->GetYaxis()->SetTitle("pseudo experiments");

  c1->Update();
                
  TPaveStats* sb = (TPaveStats*)histo->GetListOfFunctions()->FindObject("stats");
  sb->SetX1NDC(1);
  sb->SetX2NDC(1);
  sb->SetY1NDC(1);
  sb->SetY2NDC(1);

  latex->DrawLatex(0.39, 0.865, "CMS Simulation");

  c1->SaveAs("PEXP_XS.pdf");

  histo = (TH1D*) f->Get("Pull_XS_FDATA_nDisc_0_WP_3.3");

  histo->Draw();

//gPad->SetGrid();

  histo->Fit("gaus");

  c1->Update();
                
  sb = (TPaveStats*)histo->GetListOfFunctions()->FindObject("stats");
  sb->SetOptStat(0);

  histo->SetTitle("");
  histo->GetXaxis()->SetTitle("Pull");
  histo->GetYaxis()->SetTitle("pseudo experiments");

  latex->DrawLatex(0.39, 0.865, "CMS Simulation");

  c1->SaveAs("PEXP_XS_pull.pdf");

  ////////////////////////////////////////////////////////////
  // 2D

  TH2F* nocut = (TH2F*) f->Get("MC_XS_VS_BtagEfficiency_PseudoExps");
  TH2F* cut = (TH2F*) f->Get("MC_XS_VS_BtagEfficiency_PseudoExps_FDATA_nDisc_0_WP_3.3");

  cut->SetTitle("");
  nocut->SetTitle("");

  nocut->GetYaxis()->SetLabelSize(0.04);
  cut->GetYaxis()->SetLabelSize(0.04);
  nocut->GetXaxis()->SetLabelSize(0.04);
  cut->GetXaxis()->SetLabelSize(0.04);
  //nocut->Draw();

  TCanvas* c2 = new TCanvas("c2","c2",1024,768);

  c2->Divide(2,1);

  c2->cd(1);

  nocut->Draw("colz");

  c2->Update();

  TPaveStats* sb2 = (TPaveStats*)nocut->GetListOfFunctions()->FindObject("stats");
  sb2->SetX1NDC(1);
  sb2->SetX2NDC(1);
  sb2->SetY1NDC(1);
  sb2->SetY2NDC(1);

//gPad->SetGrid();

  stringstream s; s << nocut->GetCorrelationFactor();
  TLatex* corr = new TLatex(0.7,240,("#rho = "+s.str()).c_str());

  corr->SetTextColor(kRed);
  corr->SetTextSize(0.07);

  corr->Draw("same");

  c2->cd(2);

  cut->Draw("colz");

  c2->Update();

  TPaveStats* sb3 = (TPaveStats*)cut->GetListOfFunctions()->FindObject("stats");
  sb3->SetX1NDC(1);
  sb3->SetX2NDC(1);
  sb3->SetY1NDC(1);
  sb3->SetY2NDC(1);

//gPad->SetGrid();

  stringstream t; t<< cut->GetCorrelationFactor();
  corr = new TLatex(0.7,240,("#rho = "+t.str()).c_str());
corr->SetTextColor(kRed);
  corr->SetTextSize(0.07);

  corr->Draw("same");  

  c2->Draw();

  c2->SaveAs("PEXP_2D.pdf");



  ////////////////////////////////////////////////////////////
  // 2D

  TH2F* correl = (TH2F*) f->Get("Meas_MisTagEfficiency_VS_BtagEffMeas_WP_3.3");

  correl->SetTitle("");

  correl->GetYaxis()->SetLabelSize(0.04);
  correl->GetXaxis()->SetLabelSize(0.04);
  correl->GetZaxis()->SetLabelSize(0.04);
  correl->GetZaxis()->SetTitleOffset(0.7);
  correl->GetZaxis()->SetTitle("Number of experiments");
  //nocorrel->Draw();

  TCanvas* c3 = new TCanvas("c3","c3",1024,768);

  correl->Draw("colz");

  c3->Update();

  TPaveStats* sb2 = (TPaveStats*)correl->GetListOfFunctions()->FindObject("stats");
  sb2->SetX1NDC(1);
  sb2->SetX2NDC(1);
  sb2->SetY1NDC(1);
  sb2->SetY2NDC(1);

//gPad->SetGrid();

  stringstream s; s << correl->GetCorrelationFactor();
  TLatex* corr = new TLatex(0.645,0.0515,("#rho = "+s.str()).c_str());

  corr->SetTextColor(kRed);
  corr->SetTextSize(0.04);

  //corr->Draw("same");

  c3->Draw();

  c3->SaveAs("PEXP_correl_eb_eq.pdf");

  //  gApplication->Terminate();
}
