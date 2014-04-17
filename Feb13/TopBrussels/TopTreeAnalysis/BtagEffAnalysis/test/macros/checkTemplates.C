{

  TFile* f = new TFile("../FitTemplates/nDisc_0_ptbinlow_0_etabinlow_-9990__Chi2Cut_90_MLJTemplates.root ");

	TH1F* x1 = (TH1F*) f->Get("nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Bkg_Var0_TTbar");
	TH1F* y1 = (TH1F*) f->Get("nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Bkg_Var0_TTbar_bTagM");
	
	TH1F* x2 = (TH1F*) f->Get("nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Bkg_Var0_VVMC");
	TH1F* y2 = (TH1F*) f->Get("nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Bkg_Var0_VVMC_bTagM");

	cout << "ttbar non-b before btag " << x1->Integral() << endl;
	cout << "ttbar non-b after btag " << y1->Integral() << endl;
	cout << "ttbar mistag rate " << y1->Integral()/x1->Integral()  << endl;
	
	cout << "non-ttbar non-b before btag " << x2->Integral() << endl;
	cout << "non-ttbar non-b after btag " << y2->Integral() << endl;
	cout << "non-ttbar mistag rate " << y2->Integral()/x2->Integral()  << endl;

}
