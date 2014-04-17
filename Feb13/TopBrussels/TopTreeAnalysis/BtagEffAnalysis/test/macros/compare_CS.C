{

  TFile* f = new TFile("../store/nom.root");
  TFile* g = new TFile("../store/jes.root");

  TH1F* mc = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_ControlVar");
  TH1F* mcR = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_ControlVarReweigh");
  TH1F* mcPTL = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_Pt_Left");
   TH1F* mcPTR = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_Pt_Right");
  
  TH1F* data = (TH1F*) g->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_ControlVar");
  TH1F* dataR = (TH1F*) g->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_ControlVarReweigh");
  TH1F* dataPTL = (TH1F*) g->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_Pt_Left");  
  TH1F* dataPTR = (TH1F*) g->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_Pt_Right");
 

  TCanvas* c1 = new TCanvas("c1","c1");
  c1->Divide(2,2);
  TCanvas* c2 = new TCanvas("c2","c2");
  c2->Divide(2,2);

  c1->cd(1);
  mc->SetLineColor(kBlack);
  mc->DrawNormalized("hist");
 
  data->SetLineColor(kRed);
  data->DrawNormalized("sames");

  c1->cd(2);
  TH1F* div = (TH1F*) data->Clone();

  div->Divide(mc);

  div->Draw();

  div->Fit("pol1");

  c1->cd(3);
  mcR->SetLineColor(kBlack);
  mcR->DrawNormalized("hist");
  
  dataR->SetLineColor(kRed);
  dataR->DrawNormalized("sames");
  
  c1->cd(4);
  TH1F* divR = (TH1F*) dataR->Clone();

  divR->Divide(mcR);

  divR->Draw();

  divR->Fit("pol1");


   c2->cd(1);
  mcPTL->SetLineColor(kBlack);
  mcPTL->DrawNormalized("hist");
  
  dataPTL->SetLineColor(kRed);
  dataPTL->DrawNormalized("sames");
  
  c2->cd(2);
  TH1F* divPTL = (TH1F*) dataPTL->Clone();

  divPTL->Divide(mcPTL);

  divPTL->Draw();

  divPTL->Fit("pol1");   

  c2->cd(3);
  mcPTR->SetLineColor(kBlack);
  mcPTR->DrawNormalized("hist");
  
  dataPTR->SetLineColor(kRed);
  dataPTR->DrawNormalized("sames");
  
  c2->cd(4);
  TH1F* divPTR = (TH1F*) dataPTR->Clone();

  divPTR->Divide(mcPTL);

  divPTR->Draw();

  divPTR->Fit("pol1");

  /*c1->cd(2);

  data->SetLineColor(kBlack);
  data->DrawNormalized("hist");
  
  dataR->SetLineColor(kRed);
  dataR->DrawNormalized("sames");*/
  
}
