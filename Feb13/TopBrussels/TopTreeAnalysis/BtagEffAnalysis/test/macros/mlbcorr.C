{

  TFile* f = new TFile("../corrHistos.root","READ");

	TH2D* bcorr = (TH2D*) f->Get("b_mlj_correlation_tche");
	TH2D* nbcorr = (TH2D*) f->Get("nonb_mlj_correlation_tche");

	TH2D* corr = (TH2D*) f->Get("mlj_correlation_tche");
	TH2D* datacorr = (TH2D*) f->Get("data_mlj_correlation_tche");

	TCanvas* c1 = new TCanvas("c1","c1",800,600);
	
	c1->cd();
	
	bcorr->Draw("colz");
	
	c1->Update();
	
	TPaveStats* sb = (TPaveStats*)bcorr->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	c1->Draw();
	
TLatex* text = new TLatex(0.13,0.87,"CMS Simulation");
	
	text->SetTextSize(0.05);
	
	text->SetNDC();
	
	text->Draw();

	TCanvas* c2 = new TCanvas("c4","c2",800,600);
	
	c2->cd();
	
	nbcorr->Draw("colz");
	
	c2->Update();
	
	sb = (TPaveStats*)nbcorr->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	c2->Draw();
	
TLatex* text = new TLatex(0.13,0.87,"CMS Simulation");
	
	text->SetTextSize(0.05);
	
	text->SetNDC();
	
	text->Draw();

	c1->SaveAs("mlbcorr_b.pdf");
	c2->SaveAs("mlbcorr_nonb.pdf");
	
	/////////////////////////////++
	
	TCanvas* c3 = new TCanvas("c3","c3",800,600);
	
	c3->Divide(2,2);
	
	c3->cd(1);
	
	corr->SetTitle("MC");
	
	corr->Draw("colz");
	
	c3->Update();
	
	TPaveStats* sb = (TPaveStats*)corr->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	c3->cd(2);
	
	TH1D* profcorr = (TH1D*) corr->ProfileX();
	
	//profcorr->GetYaxis()->SetTitle("m_{#mu j} (GeV/c^{2})");
	profcorr->SetTitle("ProfileX red=data");
	profcorr->GetYaxis()->SetTitle("TCHE discriminant");
	profcorr->GetYaxis()->SetTitleOffset(1.0);
	profcorr->Draw("hist");
	
	TH1D* profdatacorr = (TH1D*) datacorr->ProfileX();
	
	//profdatacorr->GetYaxis()->SetTitle("m_{#mu j} (GeV/c^{2})");
	profdatacorr->GetYaxis()->SetTitle("TCHE discriminant");
	profdatacorr->GetYaxis()->SetTitleOffset(1.0);
	profdatacorr->SetLineColor(kRed);
	
	profdatacorr->Draw("sames hist");
	
	
	c3->Draw();
		
	c3->cd(3);
	
	datacorr->SetTitle("data");
	
	datacorr->Draw("colz");
	
	c3->Update();
	
	sb = (TPaveStats*)datacorr->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	c3->Draw();
	
	c3->cd(4);
	
	TH1D* profcorr = (TH1D*) corr->ProfileY();
	
	profcorr->GetYaxis()->SetTitle("m_{#mu j} (GeV/c^{2})");
	profcorr->SetTitle("ProfileY red=data");
	//profcorr->GetYaxis()->SetTitle("TCHE discriminant");
	profcorr->GetYaxis()->SetTitleOffset(1.0);
	profcorr->Draw("hist");
	
	TH1D* profdatacorr = (TH1D*) datacorr->ProfileY();
	
	profdatacorr->GetYaxis()->SetTitle("m_{#mu j} (GeV/c^{2})");
	//profdatacorr->GetYaxis()->SetTitle("TCHE discriminant");
	profdatacorr->GetYaxis()->SetTitleOffset(1.0);
	profdatacorr->SetLineColor(kRed);
	
	profdatacorr->Draw("sames hist");
	
	
	c3->SaveAs("correlation_mlb_TCHE_data.pdf");
	
	///////////////////////////////
	
	cout << "Correlation b (MC)" << bcorr->GetCorrelationFactor() << endl;
	cout << "Correlation non-b (MC)" << nbcorr->GetCorrelationFactor() << endl;
	cout << "Correlation all (MC) " << corr->GetCorrelationFactor() << endl;
	cout << "Correlation all (DATA) " << datacorr->GetCorrelationFactor() << endl;

}
