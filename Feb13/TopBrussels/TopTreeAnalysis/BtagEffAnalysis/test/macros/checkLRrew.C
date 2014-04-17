{

  TFile* f = new TFile("../TTrees/DATA_16/NTupleAnalyzed.root");
	
	TH1F* a = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_Left1DX");
	TH1F* b = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_Right1DX");
	TH1F* c = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_Right1DXReweigh");
	
  TCanvas* c1 = new TCanvas("c1","c1");

	c1->Divide(3,1);
	
	c1->cd(1);
	
	gPad->SetGrid();
	
	a->Scale(1./a->Integral());
	a->GetYaxis()->SetRangeUser(0,1);

	a->Draw("hist");
	
	b->DrawNormalized("sames");
	c->SetMarkerColor(kRed);
	c->SetLineColor(kRed);
	c->DrawNormalized("sames");
	
	a = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Sng_Left1DX");
	b = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Sng_Right1DX");
	c = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Sng_Right1DXReweigh");
	
	c1->cd(2);
	
	gPad->SetGrid();
	
	a->Scale(1./a->Integral());
	a->GetYaxis()->SetRangeUser(0,1);
	
	a->Draw("hist");
	
	b->DrawNormalized("sames");
	c->SetMarkerColor(kRed);
	c->SetLineColor(kRed);
	c->DrawNormalized("sames");

	a = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Bkg_Left1DX");
	b = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Bkg_Right1DX");
	c = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Bkg_Right1DXReweigh");
	
	c1->cd(3);
	
	gPad->SetGrid();
	
	a->Scale(1./a->Integral());
	a->GetYaxis()->SetRangeUser(0,1);
	
	a->Draw("hist");
	
	b->DrawNormalized("sames");
	c->SetMarkerColor(kRed);
	c->SetLineColor(kRed);
	c->DrawNormalized("sames");
	

}
