{

  TFile* f = new TFile("../store/dev_f1.root");
  TFile* g = new TFile("../store/dev_f2.root");

  TH1F* a = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_XS");
  TH1F* b = (TH1F*) g->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_XS");

 
  
  TCanvas* c1 = new TCanvas("c1","c1");
  c1->Divide(1,1);

  c1->cd(1);
  
  a->DrawNormalized("hist");

  b->DrawNormalized("sames E");

}
