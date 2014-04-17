{
    
    string lumiPlot="11.9";
    
    string plotName = "MLB_BTV";
    
    string filetmp="../StackPlots_Mu.root";
    
    string toSave = "NONE";
    
    double lumi_error = 0.044;
    double ttbar_error = 0.15;
    double other_error = 0.30;
    
    ttbar_error = 0;
    other_error = 0;    
    
    double btagSF=1;
    double ubtagSF=0;
    
    ifstream t(".plotopts",ios::in);
    
    while (!t.eof()) {
        
        // << t.readline() << endl;
        
        t >> filetmp >> lumiPlot >> btagSF >> ubtagSF >> plotName >> toSave;
        
    }
    
    cout << filetmp << " " << lumiPlot << " " << btagSF << " " << ubtagSF << " " << plotName << endl;
    

    TString fileName = filetmp.c_str();
    
    //ttbar_error=sqrt((ttbar_error*ttbar_error)+(ubtagSF*ubtagSF));
    //other_error=ttbar_error;
    
    //plotName="nPV";
    //fileName = "BtaggingOutput.root";
    //plotName = "Selected_Events_pT_alljets";
    //plotName = "Selected_Events_pT_jet1";
    //lotName = "Selected_Events_pT_4leadingjets";
    //plotName = "BestJetCombMLB";
    //plotName="bTagger_TCHE";
    
    //lumiPlot="2.14";
    
    // RUN ME WITH OPTIONS -n and -l!!
    
    // SET THE TDRStyle
    
    // For the canvas:
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasColor(0); // must be kWhite but I dunno how to do that in PyROOT
    gStyle->SetCanvasDefH(600); //Height of canvas
    gStyle->SetCanvasDefW(600); //Width of canvas
    gStyle->SetCanvasDefX(0);   //POsition on screen
    gStyle->SetCanvasDefY(0);
    
    
    // For the Pad:
    gStyle->SetPadBorderMode(0);
    // gStyle->SetPadBorderSize(Width_t size = 1);
    gStyle->SetPadColor(0); // kWhite
    gStyle->SetPadGridX(0); //false
    gStyle->SetPadGridY(0); //false
    gStyle->SetGridColor(0);
    gStyle->SetGridStyle(3);
    gStyle->SetGridWidth(1);
    
    // For the frame:
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameBorderSize(1);
    gStyle->SetFrameFillColor(0);
    gStyle->SetFrameFillStyle(0);
    gStyle->SetFrameLineColor(1);
    gStyle->SetFrameLineStyle(1);
    gStyle->SetFrameLineWidth(1);
    
    // For the histo:
    // gStyle->SetHistFillColor(1);
    // gStyle->SetHistFillStyle(0);
    gStyle->SetHistLineColor(1);
    gStyle->SetHistLineStyle(0);
    gStyle->SetHistLineWidth(1);
    // gStyle->SetLegoInnerR(Float_t rad = 0.5);
    // gStyle->SetNumberContours(Int_t number = 20);
    
    gStyle->SetEndErrorSize(2);
    //gStyle->SetErrorMarker(20);   /// I COMMENTED THIS OUT
    //gStyle->SetErrorX(0.);
    
    //gStyle->SetMarkerStyle(20);
    
    
    //For the fit/function:
    gStyle->SetOptFit(1011);
    gStyle->SetFitFormat("5.4g");
    gStyle->SetFuncColor(2);
    gStyle->SetFuncStyle(1);
    gStyle->SetFuncWidth(1);
    
    //For the date:
    gStyle->SetOptDate(0);
    // gStyle->SetDateX(Float_t x = 0.01);
    // gStyle->SetDateY(Float_t y = 0.01);
    
    // For the statistics box:
    gStyle->SetOptFile(0);
    gStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
    gStyle->SetStatColor(0); // kWhite
    gStyle->SetStatFont(42);
    //gStyle->SetStatFontSize(0.025);
    gStyle->SetStatFontSize(0.04);
    gStyle->SetStatTextColor(1);
    gStyle->SetStatFormat("6.4g");
    gStyle->SetStatBorderSize(1);
    gStyle->SetStatH(0.1);
    gStyle->SetStatW(0.15);
    // gStyle->SetStatStyle(Style_t style = 1001);
    // gStyle->SetStatX(Float_t x = 0);
    // gStyle->SetStatY(Float_t y = 0);
    
    // Margins:
    gStyle->SetPadTopMargin(0.07);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadLeftMargin(0.16);
    //gStyle->SetPadRightMargin(0.12);
    gStyle->SetPadRightMargin(0.03);
    
    // For the Global title:
    
    gStyle->SetOptTitle(0);
    gStyle->SetTitleFont(42);
    gStyle->SetTitleColor(1);
    gStyle->SetTitleTextColor(1);
    gStyle->SetTitleFillColor(10);
    gStyle->SetTitleFontSize(0.05);
    // gStyle->SetTitleH(0); // Set the height of the title box
    // gStyle->SetTitleW(0); // Set the width of the title box
    // gStyle->SetTitleX(0); // Set the position of the title box
    // gStyle->SetTitleY(0.985); // Set the position of the title box
    // gStyle->SetTitleStyle(Style_t style = 1001);
    // gStyle->SetTitleBorderSize(2);
    
    // For the axis titles:
    
    gStyle->SetTitleColor(1, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetTitleSize(0.06, "XYZ");
    // gStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
    // gStyle->SetTitleYSize(Float_t size = 0.02);
    gStyle->SetTitleXOffset(0.9);
    gStyle->SetTitleYOffset(1.25);
    // gStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset
    
    // For the axis labels:
    
    gStyle->SetLabelColor(1, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetLabelOffset(0.007, "XYZ");
    gStyle->SetLabelSize(0.05, "XYZ");
    
    // For the axis:
    
    gStyle->SetAxisColor(1, "XYZ");
    gStyle->SetStripDecimals(1); // kTRUE
    gStyle->SetTickLength(0.03, "XYZ");
    gStyle->SetNdivisions(510, "XYZ");
    gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    gStyle->SetPadTickY(1);
    
    // Change for log plots:
    gStyle->SetOptLogx(0);
    gStyle->SetOptLogy(0);
    gStyle->SetOptLogz(0);
    
    // Postscript options:
    gStyle->SetPaperSize(20.,20.);
    
    // CREATE THE CANVAS
    
    TCanvas* c = new TCanvas("c","c", 1000, 800);
    
    // OPEN TFILE
    TFile* f = new TFile(fileName);
    
    //TH1F* effMC = (TH1F*) g->Get("Variable_1_0/Discriminator_4/nDisc_4_ptbinlow_0_etabinlow_-9990_TH1Sng_BtagEffAll");
    
    // GET THE MC STACK
    
    THStack* stack = (THStack*) f->Get(("MultiSamplePlot_"+plotName+"/Stack_"+plotName).c_str());
    stack->SetTitle("");
    
    // GET AND ADD SEPERATE MC PLOTS
    
    cout << "TH1D* h1 = (TH1D*) f->Get((\"MultiSamplePlot_"+plotName+"/"+plotName+"_TTbarJets_Signal));" << endl;
    
    TH1D* h1 = (TH1D*) f->Get(("MultiSamplePlot_"+plotName+"/"+plotName+"_TTbarJets_Signal").c_str());
    TH1D* h2 = (TH1D*) f->Get(("MultiSamplePlot_"+plotName+"/"+plotName+"_TTbarJets_Other").c_str());
    TH1D* h3 = (TH1D*) f->Get(("MultiSamplePlot_"+plotName+"/"+plotName+"_ST_t_t").c_str());
    TH1D* h4 = (TH1D*) f->Get(("MultiSamplePlot_"+plotName+"/"+plotName+"_ST_t_tbar").c_str());   
    TH1D* h5 = (TH1D*) f->Get(("MultiSamplePlot_"+plotName+"/"+plotName+"_ST_tW_t").c_str());
    TH1D* h6 = (TH1D*) f->Get(("MultiSamplePlot_"+plotName+"/"+plotName+"_ST_tW_tbar").c_str());
    //TH1D* h7 = (TH1D*) f->Get(("MultiSamplePlot_"+plotName+"/"+plotName+"_WJets").c_str());
    TH1D* h7a = (TH1D*) f->Get(("MultiSamplePlot_"+plotName+"/"+plotName+"_W_WJets_1jets").c_str());
    TH1D* h7b = (TH1D*) f->Get(("MultiSamplePlot_"+plotName+"/"+plotName+"_W_WJets_2jets").c_str());
    TH1D* h7c = (TH1D*) f->Get(("MultiSamplePlot_"+plotName+"/"+plotName+"_W_WJets_3jets").c_str());
    TH1D* h7d = (TH1D*) f->Get(("MultiSamplePlot_"+plotName+"/"+plotName+"_W_WJets_4jets").c_str()); 
    TH1D* h7 = (TH1D*) h7a->Clone();
    h7->Add(h7b);
    h7->Add(h7c);
    h7->Add(h7d);

    //TH1D* h8 = (TH1D*) f->Get(("MultiSamplePlot_"+plotName+"/"+plotName+"_ZJets").c_str());
    TH1D* h8a = (TH1D*) f->Get(("MultiSamplePlot_"+plotName+"/"+plotName+"_Z_ZJets_1jets").c_str());
    TH1D* h8b = (TH1D*) f->Get(("MultiSamplePlot_"+plotName+"/"+plotName+"_Z_ZJets_2jets").c_str());
    TH1D* h8c = (TH1D*) f->Get(("MultiSamplePlot_"+plotName+"/"+plotName+"_Z_ZJets_3jets").c_str());
    TH1D* h8d = (TH1D*) f->Get(("MultiSamplePlot_"+plotName+"/"+plotName+"_Z_ZJets_4jets").c_str());
    TH1D* h8 = (TH1D*) h8a->Clone();
    h8->Add(h8b);
    h8->Add(h8c);
    h8->Add(h8d);

    TH1D* h9 = (TH1D*) f->Get(("MultiSamplePlot_"+plotName+"/"+plotName+"_multijet").c_str());
    
    h1->SetBinContent(h1->GetNbinsX(),h1->GetBinContent(h1->GetNbinsX())+h1->GetBinContent(h1->GetNbinsX()+1));
    h2->SetBinContent(h1->GetNbinsX(),h2->GetBinContent(h1->GetNbinsX())+h2->GetBinContent(h2->GetNbinsX()+1));
    h3->SetBinContent(h1->GetNbinsX(),h3->GetBinContent(h1->GetNbinsX())+h3->GetBinContent(h3->GetNbinsX()+1));
    h4->SetBinContent(h1->GetNbinsX(),h4->GetBinContent(h1->GetNbinsX())+h4->GetBinContent(h4->GetNbinsX()+1));
    h5->SetBinContent(h1->GetNbinsX(),h5->GetBinContent(h1->GetNbinsX())+h5->GetBinContent(h5->GetNbinsX()+1));
    h6->SetBinContent(h1->GetNbinsX(),h6->GetBinContent(h1->GetNbinsX())+h6->GetBinContent(h6->GetNbinsX()+1));
    h7->SetBinContent(h1->GetNbinsX(),h7->GetBinContent(h1->GetNbinsX())+h7->GetBinContent(h7->GetNbinsX()+1));
    h8->SetBinContent(h1->GetNbinsX(),h8->GetBinContent(h1->GetNbinsX())+h8->GetBinContent(h8->GetNbinsX()+1));
    if (h9) h9->SetBinContent(h1->GetNbinsX(),h9->GetBinContent(h1->GetNbinsX())+h9->GetBinContent(h9->GetNbinsX()+1));
        
        h1->SetBinContent(h1->GetNbinsX()+1,0);
        h2->SetBinContent(h2->GetNbinsX()+1,0);
        h3->SetBinContent(h3->GetNbinsX()+1,0);
        h4->SetBinContent(h4->GetNbinsX()+1,0);
        h5->SetBinContent(h5->GetNbinsX()+1,0);
        h6->SetBinContent(h6->GetNbinsX()+1,0);
        h7->SetBinContent(h7->GetNbinsX()+1,0);
        h8->SetBinContent(h8->GetNbinsX()+1,0);
        if (h9) h9->SetBinContent(h9->GetNbinsX()+1,0);
            
            TH1D* added = (TH1D*) h1->Clone();
            added->Sumw2();
            added->Add(h2);
            added->Add(h3);
            added->Add(h4);
            added->Add(h5);
            added->Add(h6);
            added->Add(h7);
            added->Add(h8);
            if (h9) added->Add(h9);
                
                TH1D* hmc = (TH1D*) h1->Clone();
                hmc->Add(hmc,h2,1.,1.);
                hmc->Add(hmc,h3,1.,1.);
                hmc->Add(hmc,h4,1.,1.);
                hmc->Add(hmc,h5,1.,1.);
                hmc->Add(hmc,h6,1.,1.);
                hmc->Add(hmc,h7,1.,1.);
                hmc->Add(hmc,h8,1.,1.);
                if (h9) hmc->Add(hmc,h9,1.,1.);
                    hmc->SetName("hmc");
                    
                    // GET THE DATA TH1
                    
                    TH1D* data = (TH1D*) f->Get(("MultiSamplePlot_"+plotName+"/"+plotName+"_Data").c_str());
                    
                    if (data) data->Scale(data->Integral()/(data->Integral()*btagSF));
                    
                    // GET THE LEGEND
                    
                    TLegend* legend;// = (TLegend*) f->Get(("MultiSamplePlot_"+plotName+"/leg").c_str());
    
    if (data) c->SetBottomMargin(0.3);

    //tmp fix
    
    if (plotName.find("BestJetCombChi2") != string::npos)
      legend = new TLegend(0.7,0.64,0.84,0.89);
    else if (plotName.find("bTagger") != string::npos)
      legend = new TLegend(0.7,0.64,0.94,0.89);
    else
      legend = new TLegend(0.7,0.58,0.94,0.89);
    
    if (data) legend->AddEntry(data,"Data","L E");
    
    h1->SetFillColor(kRed+1);   
    h2->SetFillColor(kRed-7);   
    h3->SetFillColor(kMagenta);   
    h8->SetFillColor(kAzure-2);   
    h7->SetFillColor(kGreen-3);   
    
    legend->AddEntry(h1,"t#bar{t} signal"," F");
    legend->AddEntry(h2,"t#bar{t} other","F");
    legend->AddEntry(h3,"Single-Top","F");

    legend->AddEntry(h7,"W#rightarrowl#nu","F");
    legend->AddEntry(h8,"Z/#gamma*#rightarrowl^{+}l^{-}","F");

    /*datasets_title.push_back("Z/#gamma*#rightarrowl^{+}l^{-}");
     datasets_title.push_back("W#rightarrowl#nu");
     datasets_title.push_back("SingleTop t");
     datasets_title.push_back("SingleTop t");
     datasets_title.push_back("SingleTop tW");
     datasets_title.push_back("SingleTop tW");
     datasets_title.push_back("t#bar{t} other");
     datasets_title.push_back("t#bar{t} signal");
     datasets_title.push_back("Data");*/
    
    
    legend->SetShadowColor(0);
    legend->SetFillColor(0);
    legend->SetLineColor(0);
    legend->SetX1NDC(.75);
    legend->SetX2NDC(.90);  
    legend->SetY1NDC(.60);
    legend->SetY2NDC(.90);
    
    // BUILD THE RATIO PLOT
    
    TH1D* ratio;
    TPad* pad;
    if (data) {
    ratio = (TH1D*) data->Clone();
    ratio->Divide(added);
    
    pad = new TPad("pad", "pad", 0.0, 0.0, 1.0, 1.0);
    pad->SetTopMargin(0.7);
    pad->SetFillColor(0);
    pad->SetFillStyle(0);
    pad->SetGridy(1);
    
    ratio->SetMarkerStyle(20);
    ratio->SetMaximum( 1.5 );
    ratio->SetMinimum(0.5);
    ratio->GetXaxis().SetLabelSize(0.04);
    ratio->GetYaxis().SetLabelSize(0.04);
    ratio->GetXaxis().SetTitle(stack->GetXaxis()->GetTitle());
    ratio->GetYaxis().SetTitle("Ratio");
    //ratio->GetYaxis().SetTitleSize(0.03);
    ratio->SetMarkerSize(0.7);
    ratio->GetYaxis().SetNdivisions(5);
    
    }

    // BUILD THE LUMIERROR PLOT
    
    TH1D* lumiband = (TH1D*) hmc->Clone();
    for (int ilum=0; ilum < lumiband.GetNbinsX()+1; ilum++) {
        
        double l_error = lumiband->GetBinContent(ilum+1)*lumi_error;
        double tt_error = (h1->GetBinContent(ilum+1)+h2->GetBinContent(ilum+1))*ttbar_error;
        double w_error = h7->GetBinContent(ilum+1)*other_error;
        double z_error = h8->GetBinContent(ilum+1)*other_error;
        double st_t_t_error = h3->GetBinContent(ilum+1)*other_error;
        double st_t_tb_error = h4->GetBinContent(ilum+1)*other_error;
        double st_tW_t_error = h5->GetBinContent(ilum+1)*other_error;
        double st_tW_tb_error = h6->GetBinContent(ilum+1)*other_error;
        double btag_error = hmc->GetBinContent(ilum)*ubtagSF;
        // cout << "nTTbar in bin " << ilum+1 << " " <<  << endl;
        
        lumiband->SetBinError(ilum+1,sqrt(pow(l_error,2)+pow(tt_error,2)+pow(w_error,2)+pow(z_error,2)+pow(st_t_t_error,2)+pow(st_t_tb_error,2)+pow(st_tW_t_error,2)+pow(st_tW_tb_error,2)+pow(btag_error,2)));
        //lumiband.GetBinContent(ilum+1)*lumi_error;
    }
    lumiband->SetFillStyle(3004);
    lumiband->SetFillColor(1);
    lumiband->SetMarkerStyle(1);
    
    if (data)
    if (ubtagSF > 0)
        legend->AddEntry( lumiband , "Luminosity+b tagging" , "f");
    else
        legend->AddEntry( lumiband , "Luminosity" , "f");

    
    // DRAW THE PLOT
    
    //stack->SetMaximum( .1*stack->GetMaximum() );
    //stack->SetMinimum( 0.1 );
    //cout << stack->GetMaximum() << endl;
    
    if (plotName.find("bTagger") != string::npos)
      stack->SetMaximum(stack->GetMaximum()*500);
    
    stack->Draw("HIST");
    //stack->GetXaxis()->SetTitle(label);
    stack->GetYaxis()->CenterTitle();
    //stack->GetYaxis()->SetTitle(ylabel);
    if (data) stack->GetXaxis()->SetLabelSize(0);
    else      stack->GetXaxis()->SetLabelSize(0.04);
    stack->GetYaxis().SetLabelSize(0.04);
    if (data) stack->GetXaxis()->SetTitleSize(0);
    if (data) lumiband->Draw("samee2");
    
    if (data) data->SetMarkerStyle(20);
    if (data) data->Draw("SAME:E1");

    legend->Draw();
    
    if (plotName.find("BestJetCombChi2") != string::npos || plotName.find("bTagger") != string::npos) gPad->SetLogy();
        
    if (data) pad->Draw();
    if (data) pad->cd(0);
    if (data) ratio->Draw("e");
        
        // Draw the CMS line
        
        TLatex* latex = new TLatex();
        latex->SetNDC();
        latex->SetTextSize(0.04);
        latex->SetTextAlign(31); // align right
        if (data) latex->DrawLatex(0.38, 0.95, "CMS Preliminary");
        else latex->DrawLatex(0.38, 0.95, "CMS Simulation");
        TLatex* latex2 = new TLatex();
        latex2->SetNDC();
        latex2->SetTextSize(0.04);
        latex2->SetTextAlign(31); // align right
        if (data)latex2->DrawLatex(0.80, 0.95, (lumiPlot + " fb^{-1} at #sqrt{s} = 8 TeV").c_str());
        
        
        // SAVE PLOT
        
        c->Draw();
        
        if (toSave != "NONE") {
            c->SaveAs(toSave.c_str());
            gApplication->Terminate();
        }
    //c->SaveAs(("DATAMCPlot_MLB_BTVStyle.png").c_str());
    //c->SaveAs(("DATAMCPlot_MLB_BTVStyle.pdf").c_str());
    //c->SaveAs(("DATAMCPlot_MLB_BTVStyle.root").c_str());
    
    // SAVE RATIO
    
    /*TFile* ratiofile = new TFile("STACK_RATIO.root","RECREATE");
     
     ratiofile->cd();
     
     ratio->SetName("ratio");
     
     ratio->Write();
     
     ratiofile->Close();
     
     delete ratiofile;*/
    
}
