{

  TFile* f = new TFile("../store/nom.root");
  TFile* g = new TFile("../store/jes.root");

	TH1F* x1 = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_Left1DX");
	TH1F* y1 = (TH1F*) g->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_Left1DX");

	TH1F* x2 = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_Left1DY");
	TH1F* y2 = (TH1F*) g->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_Left1DY");

	TH1F* x3 = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_Right1DX");
	TH1F* y3 = (TH1F*) g->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_Right1DX");
	
	TH1F* x4 = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_Right1DY");
	TH1F* y4 = (TH1F*) g->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_Right1DY");
	
	TH1F* x6 = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_Right1DXReweigh");
	TH1F* y6 = (TH1F*) g->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_Right1DXReweigh");
	
	TH1F* x5 = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_Right1DYReweigh");
	TH1F* y5 = (TH1F*) g->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_Right1DYReweigh");

  TCanvas* c1 = new TCanvas("c1","c1");
  c1->Divide(2,3);

	c1->cd(1);
	
	x1->SetTitle("Left pt bcand");

	//x1->GetYaxis()->SetRangeUser(0,1000);
	
	x1->Draw("hist");
	
	y1->SetLineColor(kRed);
	y1->SetMarkerColor(kRed);
	y1->Draw("sames");
	
	c1->cd(2);
	
	x2->SetTitle("Left tagger");
	
	//x2->GetYaxis()->SetRangeUser(0,1000);
	
	x2->Draw("hist");
	
	y2->SetLineColor(kRed);
	y2->SetMarkerColor(kRed);
	y2->Draw("sames");

	c1->cd(3);
	
	x3->SetTitle("Rihgt pt bcand");
	x3->Draw("hist");

	y3->SetLineColor(kRed);
	y3->SetMarkerColor(kRed);
	y3->Draw("sames");
	
	c1->cd(4);
	
	x4->SetTitle("Right tagger");
	x4->Draw("hist");
	
	y4->SetLineColor(kRed);
	y4->SetMarkerColor(kRed);
	y4->Draw("sames");
	
	c1->cd(5);
	
	x5->SetTitle("RightRew tagger");
	x5->Draw("hist");

	y5->SetLineColor(kRed);
	y5->SetMarkerColor(kRed);
	y5->Draw("sames");
	
}
