{

  TFile* f = new TFile("../store/nom.root");
  TFile* g = new TFile("../store/jes.root");

  TH1F* a = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_LeftRight1DX");
  TH1F* b = (TH1F*) g->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_LeftRight1DX");

  TH2F* c = (TH2F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_SignalControlVar12_");
  TH2F* d = (TH2F*) g->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_SignalControlVar12_");
 
  TH1F* h = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Bkg_Var0");
  TH1F* i = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_ControlVarReweigh");
  
  TH1F* j = (TH1F*) g->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Bkg_Var0");
  TH1F* k = (TH1F*) g->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_ControlVarReweigh");
  
  TCanvas* c1 = new TCanvas("c1","c1");
  c1->Divide(2,3);

  c1->cd(1);
  
  a->Draw();

  c1->cd(2);
  
  b->Draw();

  c1->cd(3);
  
  c->Draw("colz");

  c->GetXaxis()->SetRangeUser(30,200);
  c->GetZaxis()->SetRangeUser(0,50);

  c1->cd(4);
  
  d->Draw("colz");

  d->GetXaxis()->SetRangeUser(30,200);
  d->GetZaxis()->SetRangeUser(0,50);

  c1->cd(5);

  h->SetLineColor(kRed);
  h->SetMarkerColor(kRed);
  h->DrawNormalized(); 
  i->SetLineColor(kBlack);
  i->DrawNormalized("sames");

  c1->cd(6);

  j->SetLineColor(kRed);
  j->SetMarkerColor(kRed);
  j->DrawNormalized(); 
  k->SetLineColor(kBlack);
  k->DrawNormalized("sames");
}
