{
    	    
    TCanvas* c1 = new TCanvas("c1","c1",1000,600);
	TCanvas* c2 = new TCanvas("c2","c2",1400,600);
    
    c1->cd();

    TFile* data = new TFile("../corrHistos.root","READ");
    TFile* data2 = new TFile("../datamc-ptrew.root","READ");
    
    TH2D* data_pt_vs_mlb = (TH2D*) data->Get("MLB_VS_pT_DATA");
    TH2D* mc_pt_vs_mlb = (TH2D*) data->Get("MLB_VS_pT_MC");
    
    TProfile* data_profile = (TProfile*) data_pt_vs_mlb->ProfileX("data_profile");
    TProfile* mc_profile = (TProfile*) mc_pt_vs_mlb->ProfileX("mc_profile");
    
    data_profile->SetTitle("");
    mc_profile->SetTitle("");
    
    data_profile->GetXaxis()->SetTitle("m_{#mu j} (GeV/c^{2})");
    mc_profile->GetXaxis()->SetTitle("m_{#mu j} (GeV/c^{2})");
    
    data_profile->GetYaxis()->SetTitle("<p_{T}> (GeV/c)");
    mc_profile->GetYaxis()->SetTitle("<p_{T}> (GeV/c)");
    data_profile->GetYaxis()->SetTitleOffset(1.);
    mc_profile->GetYaxis()->SetTitleOffset(1.);

    data_profile->SetMarkerColor(kRed);
    data_profile->SetLineColor(kRed);

    data_profile->GetYaxis()->SetRangeUser(60,110);
    
    data_profile->Draw("");
    mc_profile->Draw("same hist");
    
    c1->Update();
	
	TPaveStats* sb = (TPaveStats*)data_profile->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
    
    c1->Draw();
    
    TLatex text;
    text.SetNDC(true);
    text.SetTextAlign(12);
	text.SetTextFont(42);
	text.SetTextSize(0.05);
    
    text.DrawLatex(0.13,0.9,"CMS Preliminary, 2.14 fb^{-1} at #sqrt{s} = 7 TeV");
    
    leg = new TLegend(0.25,0.7,0.6,0.8);
	
	leg->AddEntry(data_profile,"Data","p");
	leg->AddEntry(mc_profile,"Simulation","l");
	
	leg->Draw();

    
    c1->SaveAs("pt-profile.pdf");
    
    c2->Divide(2,1);
        
    TH1F* ptdata = (TH1F*) data2->Get("ControlPlots/pTEnr_DATA");
    TH1F* ptmc = (TH1F*) data2->Get("ControlPlots/pTEnr_MC");
    TH1F* ptmcrew = (TH1F*) data2->Get("ControlPlots/pTEnr_MCRew_enrdata");
    
    c2->cd(1);
    
    gPad->SetGrid();
    
    ptdata->SetMarkerColor(kRed);
    ptdata->SetLineColor(kRed);
    ptdata->GetYaxis()->SetTitle("a.u.");
    ptdata->SetTitle("");
    
    ptdata->Draw("E");
    ptmc->Draw("sames hist");
    
    text.DrawLatex(0.13,0.9,"CMS Preliminary, 2.14 fb^{-1} at #sqrt{s} = 7 TeV");

    c2->Update();
    TPaveStats* sb = (TPaveStats*)ptdata->GetListOfFunctions()->FindObject("stats");
    sb->SetTextColor(kWhite);
    sb->SetLineColor(kWhite);
    sb->SetX1NDC(1.);
    sb->SetX2NDC(1.);
    sb->SetY1NDC(1.);
    sb->SetY2NDC(1.);
    sb = (TPaveStats*)ptmc->GetListOfFunctions()->FindObject("stats");
    sb->SetTextColor(kWhite);
    sb->SetLineColor(kWhite);
    sb->SetX1NDC(1.);
    sb->SetX2NDC(1.);
    sb->SetY1NDC(1.);
    sb->SetY2NDC(1.);
    
    leg = new TLegend(0.22,0.7,0.85,0.8);
    
    double datamean = ptdata->GetMean(); datamean=datamean*10; int idatamean = datamean+0.5; datamean = (double)idatamean/10;
    double mcmean = ptmc->GetMean(); mcmean = mcmean*10; int imcmean = mcmean+0.5; mcmean = (double)imcmean/10;
    
    cout << datamean << " " << datamean << " " << endl;
    stringstream s, t;
    
    s<<datamean;
    t<<mcmean;
	
	leg->AddEntry(ptdata,("Data (<p_{T}>_{jets}="+s.str()+"GeV/c)").c_str(),"p");
	leg->AddEntry(ptmc,("Simulation (<p_{T}>_{jets}="+t.str()+"GeV/c)").c_str(),"l");
	
	leg->Draw();
    
    c2->Draw();

    
    c2->cd(2);
    
    gPad->SetGrid();
    
    ptdata->SetMarkerColor(kRed);
    ptdata->SetLineColor(kRed);
    ptdata->GetYaxis()->SetTitle("a.u.");
    ptdata->SetTitle("");
    
    ptdata->Draw("E");
    ptmcrew->Draw("sames hist");
    
    text.DrawLatex(0.13,0.9,"CMS Preliminary, 2.14 fb^{-1} at #sqrt{s} = 7 TeV");

    c2->Update();
    TPaveStats* sb = (TPaveStats*)ptdata->GetListOfFunctions()->FindObject("stats");
    sb->SetTextColor(kWhite);
    sb->SetLineColor(kWhite);
    sb->SetX1NDC(1.);
    sb->SetX2NDC(1.);
    sb->SetY1NDC(1.);
    sb->SetY2NDC(1.);
    sb = (TPaveStats*)ptmcrew->GetListOfFunctions()->FindObject("stats");
    sb->SetTextColor(kWhite);
    sb->SetLineColor(kWhite);
    sb->SetX1NDC(1.);
    sb->SetX2NDC(1.);
    sb->SetY1NDC(1.);
    sb->SetY2NDC(1.);
    
    leg = new TLegend(0.22,0.7,0.85,0.8);
    
    double mcmean2 = ptmcrew->GetMean(); mcmean2 = mcmean2*10; int imcmean2 = mcmean2+0.5; mcmean2 = (double)imcmean2/10;
    
    stringstream u;
    
    u<<mcmean2;
	
	leg->AddEntry(ptdata,("Data (<p_{T}>_{jets}="+s.str()+"GeV/c)").c_str(),"p");
	leg->AddEntry(ptmc,("Simulation (<p_{T}>_{jets}="+u.str()+"GeV/c)").c_str(),"l");
	
	leg->Draw();

    
    c2->Draw();

    c2->SaveAs("pt-rew-datamc.pdf");

}
