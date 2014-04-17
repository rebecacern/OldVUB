{
	
  string files[4];
  string lables[8];
  int colors[8];
  
  files[0]="pseudoexp/allMC_500_pb.root";
  files[1]="pseudoexp/allMC_1000_pb.root";
//  files[2]="pseudoexp/allMC_300_pb.root";
  
  lables[0]="PseudoExp MC L=500 /pb";
  lables[1]="PseudoExp MC L=1000 /pb";
//  lables[2]="PseudoExp MC L=300 /pb";
  
  colors[0] = 1;
  colors[1] = 2;
  colors[2] = 4;
  colors[3] = 6;
  
  float WPvals[3];
  
  WPvals[0] = 0.936121;
  WPvals[1] = 0.788262;
  WPvals[2] = 0.302541;
  
  WPvals[0] = 0.844378;
  WPvals[1] = 0.724369;
  WPvals[2] = 0.242906;
  
  Int_t w = 1344;
  Int_t h = 840;

  w = 800;
  h = 600;
  
  TCanvas* c = new TCanvas("c","Btag eff",w,h);
  TCanvas* c2 = new TCanvas("c2","Btag eff unc",w,h);
  TCanvas* c3 = new TCanvas("c3","Btag pulls",w,h);
  TCanvas* c4 = new TCanvas("c4","XS",w,h);
  TCanvas* c5 = new TCanvas("c5","unc XS",w,h);
  TCanvas* c6 = new TCanvas("c6","pulls XS",w,h);
  
  
  c->Divide(2,2);
  c2->Divide(2,2);
  c3->Divide(2,2);
  c4->Divide(2,2);
  c5->Divide(2,2);
  c6->Divide(2,2);
  
  //  TCanvas* XS = new TCanvas("XS","Cross Section Results 1 (Control Sample)",w,h);
  //  XS->Divide(3,(int)(sizeof(files)/sizeof(files[0])));
  
  TCanvas* XSMC = new TCanvas("XSMC","Cross Section Results 1 (MC templates)",w,h);
  XSMC->Divide(4,(int)(sizeof(files)/sizeof(files[0])));
  
  
  
  TLegend* leg = new TLegend(0.2,0.2,0.8.,0.8);
  
  float y1=0.79;
  float y2=0.99;
  
  for (int i=0; i<(int)(sizeof(files)/sizeof(files[0])); i++) {
    
    if (files[i] == "")
      continue;
    
    TFile* f = new TFile(files[i].c_str(),"READ");
    
    f->cd();
    
    for (int j=0; j<3; j++) {
      
      //cout << i << ' ' << j << endl;
      
      TH1D* WP;
      TH1D* uWP;
      TH1D* WPPull;

      TH1D* XS;
      TH1D* uXS;
      TH1D* XSPull;
      
      if (j==0) WP = (TH1D*) f->Get("Efficiency_PseudoExps_WP_1.7");
      else if (j==1) WP = (TH1D*) f->Get("Efficiency_PseudoExps_WP_3.3");
      else if (j==2) WP = (TH1D*) f->Get("Efficiency_PseudoExps_WP_10.2");

      if (j==0) XS = (TH1D*) f->Get("XS_FDATA_nDisc_0_WP_1.7");
      else if (j==1) XS = (TH1D*) f->Get("XS_FDATA_nDisc_0_WP_3.3");
      else if (j==2) XS = (TH1D*) f->Get("XS_FDATA_nDisc_0_WP_10.2");
            
      if (j == 0) WP->SetTitle("TrackCountingHighEff Loose");
      else if (j == 1) WP->SetTitle("TrackCountingHighEff Medium");
      else if (j == 2) WP->SetTitle("TrackCountingHighEff Tight");
      
      WP->SetLineColor(colors[i]);
      XS->SetLineColor(colors[i]);

      c->cd(j+1);

      if (i == 0) WP->Draw();
      else WP->Draw("sames");

      if (i==((int)(sizeof(files)/sizeof(files[0])))-1) {
	
	TLine* line = new TLine(WPvals[j],0.75,WPvals[j],150);
	line->SetLineColor(kRed);
	line->SetLineWidth(4);
	line->Draw();
	
	TLine* line2 = new TLine(WPvals[j],0.75,WPvals[j]+0.05,1);
	line2->SetLineColor(kRed);
	line2->SetLineWidth(4);
	line2->Draw();
	
	TLine* line3 = new TLine(WPvals[j],0.75,WPvals[j]-0.05,1);
	line3->SetLineColor(kRed);
	line3->SetLineWidth(4);
	line3->Draw();
	
	TLatex* text = new TLatex(WPvals[j]+0.075,155,"Data (L=36.1 /pb)");
	text->SetTextColor(kRed);
	text->SetTextSize(0.05);
	text->Draw();
      }
      
      WP->GetYaxis()->SetRangeUser(0.75,400);
      
      if (j == 0) leg->AddEntry(WP,lables[i].c_str(),"lp");
      
      gPad->SetLogy();
      
      gPad->SetGrid();
      
      c->Update();
      
      TPaveStats* sb = (TPaveStats*)WP->GetListOfFunctions()->FindObject("stats");
      sb->SetTextColor(colors[i]);
      sb->SetLineColor(colors[i]);
      sb->SetX1NDC(.83);
      sb->SetX2NDC(.98);
      sb->SetY1NDC(y1);
      sb->SetY2NDC(y2);

      c4->cd(j+1);
      
      if (i == 0) XS->Draw();
      else XS->Draw("sames");

      gPad->SetLogy();
       
      gPad->SetGrid();
      
      c4->Update();

      c->Update();
      
      TPaveStats* sc = (TPaveStats*)XS->GetListOfFunctions()->FindObject("stats");
      sc->SetTextColor(colors[i]);
      sc->SetLineColor(colors[i]);
      sc->SetX1NDC(.83);
      sc->SetX2NDC(.98);
      sc->SetY1NDC(y1);
      sc->SetY2NDC(y2);
      
      /////////////////
	
	if (j==0)  uWP = (TH1D*) f->Get("UncEfficiency_PseudoExps_WP_1.7");
	else if (j==1) uWP = (TH1D*) f->Get("UncEfficiency_PseudoExps_WP_3.3");
	else if (j==2) uWP = (TH1D*) f->Get("UncEfficiency_PseudoExps_WP_10.2");

	if (j==0) uXS = (TH1D*) f->Get("UNC_XS_FDATA_nDisc_0_WP_1.7");
	else if (j==1) uXS = (TH1D*) f->Get("UNC_XS_FDATA_nDisc_0_WP_3.3");
	else if (j==2) uXS = (TH1D*) f->Get("UNC_XS_FDATA_nDisc_0_WP_10.2");
	
	uWP->SetTitle(WP->GetTitle());
	
	uWP->SetLineColor(colors[i]);
	uXS->SetLineColor(colors[i]);

	c2->cd(j+1);
	
	if (i == 0) uWP->Draw();
	else uWP->Draw("sames");
	
	uWP->GetYaxis()->SetRangeUser(0.75,400);
	
	gPad->SetLogy();
	
	gPad->SetGrid();
	
	c2->Update();
	
	TPaveStats* sb2 = (TPaveStats*)uWP->GetListOfFunctions()->FindObject("stats");
	sb2->SetTextColor(colors[i]);
	sb2->SetLineColor(colors[i]);
	sb2->SetX1NDC(.83);
	sb2->SetX2NDC(.98);
	sb2->SetY1NDC(y1);
	sb2->SetY2NDC(y2);

	c5->cd(j+1);
	
	if (i == 0) uXS->Draw();
	else uXS->Draw("sames");
	
	uXS->GetYaxis()->SetRangeUser(0.75,400);
	
	gPad->SetLogy();
	
	gPad->SetGrid();
	
	c5->Update();
	
	TPaveStats* sc2 = (TPaveStats*)uXS->GetListOfFunctions()->FindObject("stats");
	sc2->SetTextColor(colors[i]);
	sc2->SetLineColor(colors[i]);
	sc2->SetX1NDC(.83);
	sc2->SetX2NDC(.98);
	sc2->SetY1NDC(y1);
	sc2->SetY2NDC(y2);

		
	////////
	
	if (j==0) WPPull = (TH1D*) f->Get("Pull_WP_1.7");
	else if (j==1) WPPull = (TH1D*) f->Get("Pull_WP_3.3");
	else if (j==2) WPPull = (TH1D*) f->Get("Pull_WP_10.2");
	
	if (j==0) XSPull = (TH1D*) f->Get("Pull_XS_FDATA_nDisc_0_WP_1.7");
	else if (j==1) XSPull = (TH1D*) f->Get("Pull_XS_FDATA_nDisc_0_WP_3.3");
	else if (j==2) XSPull = (TH1D*) f->Get("Pull_XS_FDATA_nDisc_0_WP_10.2");
	
	WPPull->SetTitle(WP->GetTitle());
	
	WPPull->SetLineColor(colors[i]);

	XSPull->SetLineColor(colors[i]);
	
	c3->cd(j+1);

	if (i == 0) WPPull->Draw();
	else WPPull->Draw("sames");
	
	WPPull->GetYaxis()->SetRangeUser(0.75,400);
	
	gPad->SetLogy();
	
	gPad->SetGrid();
	
	c3->Update();
	
	TPaveStats* sb3 = (TPaveStats*)WPPull->GetListOfFunctions()->FindObject("stats");
	sb3->SetTextColor(colors[i]);
	sb3->SetLineColor(colors[i]);
	sb3->SetX1NDC(.83);
	sb3->SetX2NDC(.98);
	sb3->SetY1NDC(y1);
	sb3->SetY2NDC(y2);

	c6->cd(j+1);

	if (i == 0) XSPull->Draw();
	else XSPull->Draw("sames");
	
	XSPull->GetYaxis()->SetRangeUser(0.75,400);
	
	gPad->SetLogy();
	
	gPad->SetGrid();
	
	c6->Update();
	
	TPaveStats* sc3 = (TPaveStats*)XSPull->GetListOfFunctions()->FindObject("stats");
	sc3->SetTextColor(colors[i]);
	sc3->SetLineColor(colors[i]);
	sc3->SetX1NDC(.83);
	sc3->SetX2NDC(.98);
	sc3->SetY1NDC(y1);
	sc3->SetY2NDC(y2);
	
    }
    
    // XS RESULTS
    
    for (int j=0; j<4; j++) {
      // MC templates
      
      TH2D* nTT;
      
      if (j==0) nTT = (TH2D*) f->Get("MC_XS_VS_BtagEfficiency_PseudoExps");
      else if (j==1) nTT = (TH2D*) f->Get("MC_XS_VS_BtagEfficiency_PseudoExps_FDATA_nDisc_0_WP_1.7;1");
      else if (j==2) nTT = (TH2D*) f->Get("MC_XS_VS_BtagEfficiency_PseudoExps_FDATA_nDisc_0_WP_3.3;1");
      else if (j==3) nTT = (TH2D*) f->Get("MC_XS_VS_BtagEfficiency_PseudoExps_FDATA_nDisc_0_WP_10.2;1");
      
      XSMC->cd((i*4)+j+1);
      
      gPad->SetGrid();
      //cout << (i*3)+j+1 << endl;
      
      nTT->GetYaxis()->SetTitleOffset(0.5);
      nTT->GetYaxis()->SetTitleSize(0.1);
      nTT->GetXaxis()->SetTitleOffset(0.55);
      nTT->GetXaxis()->SetTitleSize(0.09);
      nTT->Draw("colz");
      
      XSMC->Update();
      
      sb3 = (TPaveStats*)nTT->GetListOfFunctions()->FindObject("stats");
      sb3->SetTextColor(0);
      sb3->SetLineColor(0);
      sb3->SetX1NDC(1);
      sb3->SetX2NDC(1);
      sb3->SetY1NDC(1);
      sb3->SetY2NDC(1);
      
      int pos1 = lables[i].find("L=");
      int pos2 = lables[i].find(" /pb");
      
      cout << pos1 << " " << pos2 << " " <<  lables[i].substr(pos1+2,pos2-(pos1+1)) << endl;
      
      string txt = "";
      
      if (j==1)
	txt += "WP Loose";
      else if (j==2)
	txt += "WP Medium";
      else if (j==3)
	txt += "WP Tight";
      
      txt += " L="+lables[i].substr(pos1+2,pos2-(pos1+1))+" pb^{-1}";
      label = new TLatex(0.8,260,txt.c_str());
      
      label->SetTextSize(0.07);
      label->SetTextColor(kRed);
      label->Draw("sames");
      
      stringstream corr; corr<<fabs(nTT->GetCorrelationFactor());
      TLatex* corrLabel = new TLatex(0.8,240,("#rho = "+corr.str()).c_str());
      corrLabel->SetTextSize(0.07);
      corrLabel->SetTextColor(kRed);
      corrLabel->Draw("sames");
      
    }
    
    y1-=0.2;
    y2-=0.2;	
    
    c->cd(4);
    
    leg->Draw();
    
    c2->cd(4);
    
    leg->Draw();
    
    c3->cd(4);
    
    leg->Draw();

    c4->cd(4);
    
    leg->Draw();
    
    c5->cd(4);
    
    leg->Draw();

    c6->cd(4);
    
    leg->Draw();

  }
  
  c->Draw();
  c2->Draw();
  c3->Draw();
  c4->Draw();
  c5->Draw();
  c6->Draw();
  
  //XS->Draw();
  XSMC->Draw();
  
  c->SaveAs("Pseudo-Exp-EffResults.pdf");
  c2->SaveAs("Pseudo-Exp-UNCEffResults.pdf");
  c3->SaveAs("Pseudo-Exp-PullResults.pdf");
							
  //XS->SaveAs("Pseudo-Exp-XSResults_CSTemplate.pdf");
  XSMC->SaveAs("Pseudo-Exp-XSResults_MCTemplate.pdf");
							
}
