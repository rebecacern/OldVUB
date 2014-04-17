{

TFile* _file0 = new TFile("TTrees/DATA_14/NTupleAnalyzed.root");

TH1D* histo1 =_file0->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Sng_ControlVar");
TH1D* histo2 =_file0->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Bkg_ControlVar");
TH1D* histo3 =_file0->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_ControlVar");

histo1->SetMarkerColor(kRed);
histo2->SetMarkerColor(kBlue);
histo3->SetMarkerColor(kGreen);

//histo1->DrawNormalized(); histo1->SetMarkerColor(kRed);
//histo2->DrawNormalized("sames");histo2->SetMarkerColor(kBlue);
//histo3->DrawNormalized("sames");histo3->SetMarkerColor(kGreen);

cout << "non-b pur " <<  histo2->GetEntries()/histo3->GetEntries() << endl;

cout << histo1->GetEntries()+histo2->GetEntries() << " ?= " << histo3->GetEntries() << endl;

 cout << endl << endl;

TH1D* btag =_file0->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_BtagMeasuredRR");
btag->Rebin(5);
 btag->Draw();

 float plus=0;
 float min=0;

 for(int i=0; i<=btag->GetNbinsX()+1; i++){
   plus+=btag->GetBinContent(i);
   
 }
   
 for(int i=0; i<=btag->GetNbinsX(); i++){	
   plus-=btag->GetBinContent(i);
   min+=btag->GetBinContent(i);

   cout << "Bin " << i << " Btag Cut " << btag->GetBinLowEdge(i) << " min " << min << " plus " << plus << " plus/(plus+min) " << plus/(plus+min) <<endl;
 }   




    
 }
}
