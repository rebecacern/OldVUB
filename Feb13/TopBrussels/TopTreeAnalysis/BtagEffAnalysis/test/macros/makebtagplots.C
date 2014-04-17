{
	
	TFile* f = new TFile("../TTrees/DATA_6/NTupleAnalyzed.root");
	
	// S OVER S+B
	
	TH1F* sob = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1SoverSB_Var0");
	
	TCanvas* c1 = new TCanvas("c1","c1",800,600);
	
	//gPad->SetGrid();
	
	sob->SetTitle("");
	sob->GetYaxis()->SetRangeUser(0,0.8);
	sob->GetXaxis()->SetTitle("m_{lj} (GeV/c^{2})");
	sob->GetYaxis()->SetTitle("b purity");
	
	sob->Draw();
	
	c1->Update();
	
	TPaveStats* sb = (TPaveStats*)sob->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	TLatex* text = new TLatex(0.13,0.87,"CMS Simulation");
	
	text->SetTextSize(0.05);
	
	text->SetNDC();
	
	text->Draw();
	
	c1->Draw();
	
	c1->SaveAs("sovern.root");
	c1->SaveAs("sovern.pdf");
	c1->SaveAs("sovern.png");
	
	// btag eff all vs SS
	
	TH1F* all = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Sng_BtagEffAll");
	TH1F* ss = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Sng_BtagEffMC");
	
	TCanvas* c2 = new TCanvas("c2","c2",800,600);
	
	c2->cd();
	
	//gPad->SetGrid();
	
	all->SetTitle("");
	all->GetXaxis()->SetTitle("TCHE discriminant");
	all->GetYaxis()->SetTitle("Efficiency");
	
	all->Draw("hist");
	
	c2->Update();
	
	sb = (TPaveStats*)all->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	c2->Draw();
	
	ss->SetLineColor(kRed);
	ss->SetMarkerColor(kRed);
	ss->Draw("same");
	
	TLegend* leg = new TLegend(0.45,0.65,0.8,0.8);
	
	leg->AddEntry(all,"b candidate sample","l");
	leg->AddEntry(ss,"b enriched subsample","p");
	
	leg->Draw();
    
    text->Draw();
	
	c2->SaveAs("effcomp.pdf");
	
	// METHOD VS MC 
	
	all = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Sng_BtagAll");
	//all = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Sng_Left1DY");
	ss = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_BtagMCMeasured");
	
	TCanvas* c3 = new TCanvas("c3","c3",800,600);
	
	c3->cd();
	
	//gPad->SetGrid();
	
	all->SetTitle("");
	all->GetXaxis()->SetTitle("TCHE discriminant");
	all->GetYaxis()->SetTitle("a.u.");
	
	all->Scale(1./all->Integral());
	all->Draw("hist");
	
	c3->Update();
	
	sb = (TPaveStats*)all->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	c3->Draw();
	
	ss->SetLineColor(kRed);
	ss->SetMarkerColor(kRed);
	ss->DrawNormalized("same");
	
	leg = new TLegend(0.45,0.65,0.8,0.8);
	
	leg->AddEntry(all,"MC Truth","l");
	leg->AddEntry(ss,"Method","lp");
	
	leg->Draw();
    
    text->Draw();
	
	c3->SaveAs("measMC-truth.pdf");
	
	// LEFT btag vs right btag
	
	all = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_Left1DY");
	ss = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_Right1DY");
	
	TCanvas* c0 = new TCanvas("c0","c0",800,600);
	
	c0->cd();
	
	//gPad->SetGrid();
	gPad->SetLogy();
	
	all->SetTitle("");
	all->GetXaxis()->SetTitle("TCHE discriminant");
	all->GetYaxis()->SetTitle("a.u.");
	
	all->Draw("");
	
	c0->Update();
	
	sb = (TPaveStats*)all->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	c0->Draw();
	
	ss->SetLineColor(kRed);
	ss->SetMarkerColor(kRed);
	ss->Draw("same");
	
	leg = new TLegend(0.45,0.25,0.8,0.4);
	
	leg->AddEntry(all,"b enriched","lp");
	leg->AddEntry(ss,"b depleted","lp");
	
	leg->Draw();
    
    text->Draw();
	
	c0->SaveAs("left-right-btag.pdf");
	
	// LEFT-RIGHT
	
	TH1F* leftrightrew = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_LeftRight1DX");
	
	TCanvas* c4 = new TCanvas("c4","c4",800,600);
	
	c4->cd();
	
	//gPad->SetGrid();
	
	leftrightrew->SetTitle("");
	leftrightrew->GetXaxis()->SetTitle("p_{T} (GeV)");
	leftrightrew->GetYaxis()->SetTitle("weight");
	leftrightrew->Draw();
	
	c4->Update();
	
	sb = (TPaveStats*)leftrightrew->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	//sb->SetOptStat(1111111);
	
	c4->Draw();
	
    text->Draw();

	c4->SaveAs("pt-rew-ss.pdf");

	// PT COMP LEFT-RIGHT
	
	TH1F* lpt = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_Pt_Left");
	TH1F* rpt = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_Right1DX");
	TH1F* rptR = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_Pt_RightReweigh");
	
	rptR->Sumw2();

	TCanvas* c10 = new TCanvas("c10","c10",800,600);
	
	c10->cd();

	lpt->Scale(1./lpt->Integral());
	rpt->Scale(1./rpt->Integral());
	rptR->Scale(1./rptR->Integral());
	
	//gPad->SetGrid();
	
	rptR->SetTitle("");
	rptR->GetXaxis()->SetTitle("p_{T}^{b-jet} (GeV)");
	rptR->GetYaxis()->SetTitle("a.u.");
	rptR->Draw("E");

	lpt->Draw("sames hist");
	rpt->Draw("sames hist");
	rpt->SetLineStyle(kDashed);
	
	c10->Update();
	
	sb = (TPaveStats*)lpt->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);

	c10->Update();
	
	sb = (TPaveStats*)rpt->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);

	c10->Update();
	
	sb = (TPaveStats*)rptR->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);

	c10->Draw();
	
	//sb->SetOptStat(1111111);
		
	text->Draw();

	TLegend* leg0 = new TLegend(0.45,0.65,0.8,0.8);
	
	leg0->AddEntry(lpt,"b-enriched sample","l");
	leg0->AddEntry(rpt,"b-depleted sample","l");
	leg0->AddEntry(rptR,"p_{T}-reweighed b-depleted sample","lp");

	leg0->SetTextSize(0.025);

	leg0->Draw();

	c10->SaveAs("pt-ss.pdf");
	
	// CS VS SS no rew
	
	TH1F* all = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Bkg_Var0");
	TH1F* ss = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_ControlVar");
	
	TCanvas* c5 = new TCanvas("c5","c5",800,600);
	
	c5->cd();
	
	//gPad->SetGrid();
	
	all->SetTitle("");
	all->GetXaxis()->SetTitle("m_{#mu b} (GeV)");
	all->GetYaxis()->SetTitle("a.u.");
	
	all->Scale(1./all->Integral());
	ss->Scale(1./ss->Integral());
	
	all->GetYaxis()->SetRangeUser(0,0.1);
	
	all->Draw("hist");
	
	c5->Update();
	
	sb = (TPaveStats*)all->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	c5->Draw();
	
	ss->SetLineStyle(kDashed);
	ss->Draw("same hist");
	
	TLegend* leg = new TLegend(0.45,0.65,0.8,0.8);
	
	leg->AddEntry(all,"signal sample non b-jets","l");
	leg->AddEntry(ss,"control sample","l");
	leg->SetTextSize(0.025);
	leg->Draw();
    
	text->Draw();
	
	c5->SaveAs("cs-before.pdf");
	
	// CS VS SS no rew
	
	TH1F* ssr = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_ControlVarReweigh");
	
	ssr->Scale(1./ssr->Integral());

	TCanvas* c6 = new TCanvas("c6","c6",800,600);
	
	c6->cd();
	
	//gPad->SetGrid();
	
	all->Draw("hist");
	ss->Draw("same hist");
	ssr->Draw("same E");

	/*c6->Update();
	
	sb = (TPaveStats*)all->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	c6->Draw();*/
	
	TLegend* leg = new TLegend(0.45,0.65,0.8,0.8);
	
	leg->AddEntry(all,"signal sample non b-jets","l");
	leg->AddEntry(ss,"control sample","l");
	leg->AddEntry(ssr,"(p_{T},#eta) reweighed control sample","lp");
	leg->SetTextSize(0.025);

	leg->Draw();
	
	text->Draw();

	c6->SaveAs("cs-after.pdf");

	// CS VS SS pt
	
	TH2F* mapsig = (TH2F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2BkgVar12_Left");
	TH2F* mapsig2 = (TH2F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2BkgVar12_Right");
	mapsig->Add(mapsig2,1);

	cout << mapsig->GetNbinsY() << endl;
	TH1F* pt_ss = (TH1F*) mapsig->ProjectionX("projpt",0,mapsig->GetNbinsX(),"e");
	TH1F* pt_cs = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_Pt_Control");
	TH1F* pt_cs_rew = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_Pt_ControlReweigh");
	
	TCanvas* c11 = new TCanvas("c11","c11",800,600);
	
	c11->cd();
	
	//gPad->SetGrid();
	
	pt_cs->SetTitle("");
	pt_cs->GetXaxis()->SetTitle("p_{T}^{jet} (GeV)");
	pt_cs->GetYaxis()->SetTitle("a.u.");
	
	pt_ss->Scale(1./pt_ss->Integral());
	pt_cs->Scale(1./pt_cs->Integral());
	pt_cs_rew->Scale(1./pt_cs_rew->Integral());
		
	pt_cs->SetLineStyle(kDashed);
	pt_cs->Draw("hist");
	pt_ss->Draw("sames hist");
	
	c11->Update();
	
	sb = (TPaveStats*)pt_ss->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);

	c11->Update();
	
	sb = (TPaveStats*)pt_cs->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	c11->Draw();

	TLegend* leg11 = new TLegend(0.45,0.65,0.8,0.8);
	
	leg11->AddEntry(pt_ss,"signal sample non b-jets","l");
	leg11->AddEntry(pt_cs,"control sample","l");
	leg11->SetTextSize(0.025);

    	leg11->Draw();
	text->Draw();

	c11->SaveAs("cs-pt.pdf");

	pt_cs_rew->Draw("sames E");
	leg11->AddEntry(pt_cs_rew,"(p_{T},#eta) reweighed control sample","lp");
	leg11->Draw();
	
	c11->Update();
	
	sb = (TPaveStats*)pt_cs_rew->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);

	c11->Draw();

	c11->SaveAs("cs-pt-rew.pdf");

	// CS VS SS eta
	
	TH1F* eta_ss = (TH1F*) mapsig->ProjectionY("projeta",0,mapsig->GetNbinsY(),"e");
	TH1F* eta_cs = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_Eta_Control");
	TH1F* eta_cs_rew = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_Eta_ControlReweigh");
	
	TCanvas* c12 = new TCanvas("c12","c12",800,600);
	
	c12->cd();
	
	//gPad->SetGrid();
	
	eta_cs->SetTitle("");
	eta_cs->GetXaxis()->SetTitle("|#eta^{jet}|");
	eta_cs->GetYaxis()->SetTitle("a.u.");
	
	eta_ss->Scale(1./eta_ss->Integral());
	eta_cs->Scale(1./eta_cs->Integral());
	eta_cs_rew->Scale(1./eta_cs_rew->Integral());
		
	eta_cs->SetLineStyle(kDashed);
	eta_cs->Draw("hist");
	eta_ss->Draw("sames hist");
	
	c12->Update();
	
	sb = (TPaveStats*)eta_ss->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);

	c12->Update();
	
	sb = (TPaveStats*)eta_cs->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	c12->Draw();

	TLegend* leg12 = new TLegend(0.45,0.65,0.8,0.8);
	
	leg12->AddEntry(eta_ss,"signal sample non b-jets","l");
	leg12->AddEntry(eta_cs,"control sample","l");
	leg12->SetTextSize(0.025);

	leg12->Draw();
    
	text->Draw();

	c12->SaveAs("cs-eta.pdf");

	eta_cs_rew->Draw("sames E");

	leg12->AddEntry(eta_cs_rew,"(p_{T},#eta) reweighed control sample","lp");

	leg12->Draw();

	c12->Update();
	
	sb = (TPaveStats*)eta_cs_rew->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	c12->Draw();

	c12->SaveAs("cs-eta-rew.pdf");


	// 2D weights map
	
	TH2F* map2D = (TH2F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_SignalControlVar12_");
	
	TCanvas* c13 = new TCanvas("c13","c13",800,600);
	
	c13->cd();
	
	//gPad->SetGrid();
	
	map2D->SetTitle("");
	map2D->GetYaxis()->SetTitle("|#eta^{jet}|");
	map2D->GetXaxis()->SetTitle("p_{T}^{jet} (GeV)");
	map2D->GetZaxis()->SetTitle("weight");
	map2D->GetZaxis()->SetTitleOffset(.8);
	
	map2D->Draw("colz");

	map2D->GetXaxis()->SetRange(3,17);
	gPad->SetLogz();

	c13->Update();
	
	sb = (TPaveStats*)map2D->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);

	
	
	c13->Draw();

	text->Draw();
	c13->SaveAs("cs-rewmap.pdf");

	  

	
}


	// LEFT btag vs right btag
	
	all = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_Left1DY");
	
	TCanvas* c10 = new TCanvas("c10","c10",800,600);
	
	c10->cd();
	
	//gPad->SetGrid();
	gPad->SetLogy();
	
	all->SetTitle("");
	all->GetXaxis()->SetTitle("b-tagging discriminant value ");
	all->GetYaxis()->SetTitle("a.u.");
	
	all->Draw("hist");
	
	c10->Update();
	
	sb = (TPaveStats*)all->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	text->Draw();

	c10->Draw();
	
	c10->SaveAs("btag-left.root");
	c10->SaveAs("btag-left.pdf");
	c10->SaveAs("btag-left.png");
	
	// LEFT btag vs right btag
	
	all = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_Right1DY")->Clone();
	
	TCanvas* c11 = new TCanvas("c11","c11",800,600);
	
	c11->cd();
	
	//gPad->SetGrid();
	gPad->SetLogy();
	
	all->SetTitle("");
	all->GetXaxis()->SetTitle("b-tagging discriminant value ");
	all->GetYaxis()->SetTitle("a.u.");
	
	all->SetLineColor(kBlack);
	
	all->Draw("hist");
	
	c11->Update();
	
	sb = (TPaveStats*)all->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	text->Draw();

	c11->Draw();
	
	c11->SaveAs("btag-right.root");
	c11->SaveAs("btag-right.pdf");
	c11->SaveAs("btag-right.png");
	
	// LEFT btag vs right btag
	
	all = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_BtagMeasuredRR")->Clone();
	
	TCanvas* c12 = new TCanvas("c12","c12",800,600);
	
	c12->cd();
	
	//gPad->SetGrid();
	gPad->SetLogy();
	
	all->SetTitle("");
	all->GetXaxis()->SetTitle("b-tagging discriminant value ");
	all->GetYaxis()->SetTitle("a.u.");
	
	all->SetLineColor(kBlack);
	
	all->Draw("hist");
	
	c12->Update();
	
	sb = (TPaveStats*)all->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	text->Draw();

	c12->Draw();
	
	c12->SaveAs("btag-true.root");
	c12->SaveAs("btag-true.pdf");
	c12->SaveAs("btag-true.png");
	
	
	
}
