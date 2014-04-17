{

  TFile* f = new TFile("../TTrees/DATA_6/NTupleAnalyzed.root");

    TH1D* effMC = (TH1D*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Sng_BtagEffAll");
    TH1D* eff = (TH1D*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_BtagEffMeasuredRR");
    
    TH1D* diff = (TH1D*) eff->Clone();
    diff->Add(effMC,-1);
    diff->Divide(effMC);
    
    //TH1D* diff = (TH1D*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_BtagEffMeasuredRRDiff");
    
    diff->SetTitle("Method bias VS discriminant value");
    diff->GetXaxis()->SetTitle("TCHE discriminant");
    diff->GetYaxis()->SetTitle("Bias");
    
    double min=0;
    double max=15;
    
    diff->GetXaxis()->SetRangeUser(min,max);
    
    for (unsigned int i=1; i<diff->GetNbinsX()+1; i++) {
        
        double a = eff->GetBinContent(i);
        double b = effMC->GetBinContent(i);
        double sa = eff->GetBinError(i);
        double sb = effMC->GetBinError(i);
        
        double bias = (a-b)/b;
        
        cout << "-----" <<endl;
        cout << diff->GetBinLowEdge(i) << endl;
        cout << bias << " " << diff->GetBinContent(i) << endl;
        
        double term1 = pow(1/b,2)*pow(sa,2);
        double term2 = pow((-1/b)-((a-b)/pow(b,2)),2)*pow(sb,2);
        double ubias = sqrt(term1+term2);

        cout << ubias << " " << diff->GetBinError(i) << endl << endl;
        cout << "-----" <<endl;

    }
    
    //diff->Sumw2();
    
    diff->Draw();
   
    TLine* zeroLine = new TLine(min,0,max,0);
    zeroLine->SetLineColor(kBlue);
    zeroLine->SetLineStyle(kDashed);
    zeroLine->Draw();
    
        
}
