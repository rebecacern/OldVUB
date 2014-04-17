{

TCanvas *c = new TCanvas("c","c");
c->Divide(1,2);

TFile *_file0 = TFile::Open("OutPut/Data.root");
TFile *_file1 = TFile::Open("OutPut/MC.root");

c->cd(1);
TH1D* h1 = (TH1D*) _file0->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_Pt_Control");
TH1D* h2 = (TH1D*) _file1->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_Pt_Control");

h1->SetLineColor(kRed);
h2->SetLineColor(kBlue);

h1->DrawNormalized("hist");
h2->DrawNormalized("sameshist");

gPad->SetGrid();

c->cd(2);
TH1D* h3 = (TH1D*) _file0->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_VVData");
TH1D* h4 = (TH1D*) _file1->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_VVData");

h3->SetLineColor(kRed);
h4->SetLineColor(kBlue);

h3->DrawNormalized("hist");
h4->DrawNormalized("sameshist");

gPad->SetGrid();

}
