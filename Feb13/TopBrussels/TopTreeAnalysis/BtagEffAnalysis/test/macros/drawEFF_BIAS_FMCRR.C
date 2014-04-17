{

  TFile* f = new TFile("../TTrees/DATA_6/NTupleAnalyzed.root");

  TH1F* effMC = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Sng_BtagEffMC");
  TH1F* eff = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_BtagEffMCMeasuredRR"); 
  TH1F* bias = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_BtagEffMCMeasuredRRDiff");

  TCanvas* c1 = new TCanvas("c1","c1",800,600);
  c1->Divide(2,1);
 
  c1->cd(1);
  //gPad->SetGrid();

  eff->SetLineColor(kRed);
  eff->SetMarkerColor(kRed);
  effMC->SetTitle("");
  effMC->GetXaxis()->SetTitle("TCHE discriminator");
  effMC->GetYaxis()->SetTitle("Efficiency");
  effMC->Draw();
  eff->Draw("sames");

  c1->Update();
  TPaveStats* sb = (TPaveStats*)eff->GetListOfFunctions()->FindObject("stats");
  sb->SetTextColor(kWhite);
  sb->SetLineColor(kWhite);
  sb->SetX1NDC(1.);
  sb->SetX2NDC(1.);
  sb->SetY1NDC(1.);
  sb->SetY2NDC(1.);

  c1->Update();
  sb = (TPaveStats*)effMC->GetListOfFunctions()->FindObject("stats");
  sb->SetTextColor(kWhite);
  sb->SetLineColor(kWhite);
  sb->SetX1NDC(1.);
  sb->SetX2NDC(1.);
  sb->SetY1NDC(1.);
  sb->SetY2NDC(1.);

  c1->Draw();

TLatex* text = new TLatex(0.13,0.87,"CMS Simulation");
	
	text->SetTextSize(0.05);
	
	text->SetNDC();
	
	text->Draw();

  c1->cd(2);

  //gPad->SetGrid();
  bias->SetTitle("");
  bias->GetXaxis()->SetTitle("TCHE discriminator");
  bias->GetYaxis()->SetTitle("rel. bias (%)");
  bias->GetYaxis()->SetLabelSize(0.045);
  bias->GetYaxis()->SetTitleOffset(1.1);
  bias->Draw();

  c1->Update();
  TPaveStats* sb2 = (TPaveStats*)bias->GetListOfFunctions()->FindObject("stats");
  sb2->SetTextColor(kWhite);
  sb2->SetLineColor(kWhite);
  sb2->SetX1NDC(1.);
  sb2->SetX2NDC(1.);
  sb2->SetY1NDC(1.);
  sb2->SetY2NDC(1.);

  c1->Draw();
TLatex* text = new TLatex(0.13,0.87,"CMS Simulation");
	
	text->SetTextSize(0.05);
	
	text->SetNDC();
	
	text->Draw();
  
  c1->SaveAs("effRR.pdf");
 
}
