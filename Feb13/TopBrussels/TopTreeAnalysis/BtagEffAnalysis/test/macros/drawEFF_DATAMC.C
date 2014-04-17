{

    TFile* f = new TFile("../store/data.root");
    TFile* g = new TFile("../store/mc.root");

    TH1F* effMC = (TH1F*) g->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_BtagEffMeasuredRR");
    TH1F* eff = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_BtagEffMeasuredRR"); 
    
    TCanvas* c1 = new TCanvas("c1","c1",800,600);
    c1->Divide(2,1);
    
    c1->cd(1);
    gPad->SetGrid();
    
    eff->SetLineColor(kRed);
    eff->SetMarkerColor(kRed);
    effMC->SetTitle("");
    effMC->GetXaxis()->SetTitle("TCHE discriminant");
    effMC->GetYaxis()->SetTitle("Efficiency");
    effMC->Draw("hist");
    eff->Draw("sames");
    
    c1->Update();
    TPaveStats* sb = (TPaveStats*)eff->GetListOfFunctions()->FindObject("stats");
    sb->SetTextColor(kWhite);
    sb->SetLineColor(kWhite);
    sb->SetX1NDC(1.);
    sb->SetX2NDC(1.);
    sb->SetY1NDC(1.);
    sb->SetY2NDC(1.);
    
    c1->Update();
    sb = (TPaveStats*)effMC->GetListOfFunctions()->FindObject("stats");
    sb->SetTextColor(kWhite);
    sb->SetLineColor(kWhite);
    sb->SetX1NDC(1.);
    sb->SetX2NDC(1.);
    sb->SetY1NDC(1.);
    sb->SetY2NDC(1.);
    
    c1->Draw();
    TLatex* text = new TLatex(0.13,0.87,"CMS Preliminary, 2.14 fb^{-1} at #sqrt{s} = 7 TeV");
	
	text->SetTextSize(0.05);
	
	text->SetNDC();
	
	text->Draw();
    
    TLegend* leg = new TLegend(0.45,0.65,0.8,0.8);
    
    leg->AddEntry(effMC,"Simulation","l");
    leg->AddEntry(eff,"Data","p");
    
    leg->Draw();
    
    
    c1->cd(2);
    gPad->SetGrid();
    
    TH1F* sf = (TH1F*) eff->Clone();
    
    sf->Divide(effMC);
    
    sf->SetLineColor(kBlack);
    sf->SetMarkerColor(kBlack);
    sf->SetTitle("");
    sf->GetXaxis()->SetTitle("TCHE discriminant");
    sf->GetYaxis()->SetTitle("Scale-Factor");
    sf->GetYaxis()->SetTitleOffset(1.05);
    sf->GetYaxis()->SetLabelSize(0.035);
    
    sf->Draw();
    
    c1->SaveAs("eff-datamc.pdf");
    
    /////////
    delete eff;
    delete effMC;
    
    effMC = (TH1F*) g->Get("Variable_1_0/Discriminator_4/nDisc_4_ptbinlow_0_etabinlow_-9990_TH1Data_BtagEffMeasuredRR");
    eff = (TH1F*) f->Get("Variable_1_0/Discriminator_4/nDisc_4_ptbinlow_0_etabinlow_-9990_TH1Data_BtagEffMeasuredRR"); 
    
    TCanvas* c2 = new TCanvas("c2","c2",800,600);
    c2->Divide(2,1);
    
    c2->cd(1);
    gPad->SetGrid();
    
    eff->SetLineColor(kRed);
    eff->SetMarkerColor(kRed);
    effMC->SetTitle("");
    effMC->GetXaxis()->SetTitle("Discriminant");
    effMC->GetYaxis()->SetTitle("Efficiency");
    
    effMC->GetYaxis()->SetRangeUser(0.4,1.2);

    effMC->Draw("hist");
    eff->Draw("sames");
        
    c2->Update();
    TPaveStats* sb = (TPaveStats*)eff->GetListOfFunctions()->FindObject("stats");
    sb->SetTextColor(kWhite);
    sb->SetLineColor(kWhite);
    sb->SetX1NDC(1.);
    sb->SetX2NDC(1.);
    sb->SetY1NDC(1.);
    sb->SetY2NDC(1.);
    
    c2->Update();
    sb = (TPaveStats*)effMC->GetListOfFunctions()->FindObject("stats");
    sb->SetTextColor(kWhite);
    sb->SetLineColor(kWhite);
    sb->SetX1NDC(1.);
    sb->SetX2NDC(1.);
    sb->SetY1NDC(1.);
    sb->SetY2NDC(1.);
    
    c2->Draw();
    TLatex* text = new TLatex(0.13,0.87,"CMS Preliminary, 2.14 fb^{-1} at #sqrt{s} = 7 TeV");
	
	text->SetTextSize(0.05);
	
	text->SetNDC();
	
	text->Draw();
    
    TLegend* leg = new TLegend(0.45,0.65,0.8,0.8);
    
    leg->AddEntry(effMC,"Simulation","l");
    leg->AddEntry(eff,"Data","p");
    
    leg->Draw();
    
    
    c2->cd(2);
    gPad->SetGrid();
    
    TH1F* sf = (TH1F*) eff->Clone();
    
    sf->Divide(effMC);
    
    sf->SetLineColor(kBlack);
    sf->SetMarkerColor(kBlack);
    sf->SetTitle("");
    sf->GetXaxis()->SetTitle("Discriminant");
    sf->GetYaxis()->SetTitle("Scale-Factor");
    sf->GetYaxis()->SetTitleOffset(1.05);
    sf->GetYaxis()->SetLabelSize(0.035);
    
    sf->GetYaxis()->SetRangeUser(0.7,1.2);
    
    double sys=0.081;
    
    for (unsigned int b=1;b<sf->GetNbinsX()+1;b++) {
        double stat = sf->GetBinError(b);
        double syssytat = sqrt(pow(sys,2)+pow(stat,2));
        cout << stat << " " << sys << " " << syssytat << endl;
        sf->SetBinError(b,syssytat);
    }
    
    sf->Draw();
    
    c2->SaveAs("eff-datamc-CSV.pdf");
    c2->SaveAs("eff-datamc-CSV.root");

 
}
