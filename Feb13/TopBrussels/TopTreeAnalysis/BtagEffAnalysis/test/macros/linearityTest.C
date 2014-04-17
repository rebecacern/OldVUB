{
    //values around TCHEL
    /*double eb_true[5]={0.754,0.796,0.838,0.879,0.929};
    double ueb_true[5]={0.005,0.005,0.005,0.004,0.003};
    double eb_meas[5]={0.771,0.817,0.857,0.885,0.925};
    double ueb_meas[5]={0.026,0.027,0.027,0.025,0.022};
    */
    //values around TCHEM
    /*double eb_true[5]={0.632,0.666,0.702,0.737,0.772};
    double ueb_true[5]={0.006,0.006,0.005,0.005,0.005};
    double eb_meas[5]={0.648,0.681,0.720,0.751,0.794};
    double ueb_meas[5]={0.025,0.026,0.027,0.028,0.029};*/

  double eb_true[3]={0.341,0.686,0.829};
  double ueb_true[3]={0.002,0.002,0.001};
  double eb_meas[3]={0.344,0.687,0.827};
  double ueb_meas[3]={0.004,0.007,0.007};

    //double cuts[5]={2.58,2.1,1.7,1.35,0.915};
    
    int ref=1;
    int nUse=3;
    
    double xVals[5];
    double uxVals[5];
    double yVals[5];
    double uyVals[5];
    
    for (int i=0;i<nUse;i++) {
        
        double SF_in = eb_true[i]/eb_true[ref];
        double SF_out = eb_meas[i]/eb_meas[ref];
        double uSF_in = sqrt((pow(ueb_true[i],2)/pow(eb_true[ref],2))+((pow(eb_true[i],2)*pow(ueb_true[ref],2))/pow(eb_true[ref],4)));
        double uSF_out = sqrt((pow(ueb_meas[i],2)/pow(eb_meas[ref],2))+((pow(eb_meas[i],2)*pow(ueb_meas[ref],2))/pow(eb_meas[ref],4)));;
        
        cout << "-- Point " << i << " -- " << endl;
        cout << "SF In = " << SF_in << " +- " << uSF_in << endl;
        cout << "SF Out = " << SF_out << " +- " << uSF_out << endl;
        
        xVals[i]=SF_in;
        uxVals[i]=uSF_in;
        yVals[i]=SF_out;
        uyVals[i]=uSF_out;
        
    }
    
    double xVals2[5];
    double uxVals2[5];
    double yVals2[5];
    double uyVals2[5];
    
    //for (int i=nUse-1;i>-1;i--) {
    for (int i=0;i<nUse;i++) {
        
        cout << i << endl;
     
        float a = eb_meas[i];
        float b = eb_true[i];
        float sa = ueb_meas[i];
        float sb = ueb_true[i];     
        float bias = (a-b)/b;
                
        float term1 = pow(1/b,2)*pow(sa,2);
        float term2 = pow((-1/b)-((a-b)/pow(b,2)),2)*pow(sb,2);
        float ubias = sqrt(term1+term2);
        
        //xVals2[i]=cuts[i];
        xVals2[i]=xVals[i];
        uxVals2[i]=0;
        yVals2[i]=bias;
        uyVals2[i]=ubias;
        cout << bias << " " << ubias << endl << endl;

        
    }
    
    TCanvas* c = new TCanvas("c","c",800,600);
    
    c->Divide(2,1);
    
    c->cd(1);
    
    TGraphErrors* g = new TGraphErrors(nUse,xVals,yVals,uxVals,uyVals);
    
    g->SetTitle("bSample Method Linearity");
    g->GetXaxis()->SetTitle("SF_{in}");
    g->GetYaxis()->SetTitle("SF_{out}");
    g->GetYaxis()->SetTitleOffset(1.);
    g->Draw("A*");
    
    g->Fit("pol1","Q");
    
    c->cd(2);
    
    TGraphErrors* h = new TGraphErrors(nUse,xVals2,yVals2,uxVals2,uyVals2);
    
    h->SetTitle("bSample Method Bias");
    //h->GetXaxis()->SetTitle("Discriminant cut");
    h->GetXaxis()->SetTitle("SF_{in}");
    h->GetYaxis()->SetTitle("Rel. Bias (#hat{#epsilon_{b}}-#epsilon_{b}^{true})/#epsilon_{b}^{true}");
    h->GetYaxis()->SetTitleOffset(1.);
    h->GetYaxis()->SetLabelSize(0.03);
    h->Draw("A*");
        
    
    
}
