{

  TFile* f = new TFile("../TTrees/DATA_16/NTupleAnalyzed.root");
  TFile* g = new TFile("../TTrees/DATA_22/NTupleAnalyzed.root");
	
	TH1F* a = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_Var0");
	TH1F* b = (TH1F*) g->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_Var0");

	TH1F* c = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_ControlVar");
	TH1F* d = (TH1F*) g->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_ControlVar");
	
  TCanvas* c1 = new TCanvas("c1","c1");

	c1->Divide(2,1);
	
	c1->cd(1);
	
	gPad->SetGrid();
	
	a->Divide(b);
	
	a->Draw();
	
	a->GetXaxis()->SetRangeUser(0,300);
	a->GetYaxis()->SetRangeUser(0.5,1.5);
	
	a->Fit("pol1");
	
	c1->cd(2);
	
	gPad->SetGrid();
	
	c->Divide(d);
	
	c->Draw();
	
	c->GetXaxis()->SetRangeUser(0,300);
	c->GetYaxis()->SetRangeUser(0.5,1.5);
	
	c->Fit("pol1");

}
