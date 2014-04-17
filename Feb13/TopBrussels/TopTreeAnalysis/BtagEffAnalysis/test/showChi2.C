
{

  string files[3];
  string lables[8];
  int colors[8];

  //files[3]="pseudoexp/allMC_1000_pb.root";
 
  files[0]="TTrees/DATA_6/NTupleAnalyzed.root";

  colors[0] = 1;
  colors[1] = 2;
  colors[2] = 4;
  colors[3] = 5;

  Double_t cut = 90;

  cout << "Displaying cut value " << cut << endl;

  Double_t cutW = (cut*0.2);

  TCanvas* c = new TCanvas("c","c",800,600);

  c->Divide(2,1);

  for (int i=0; i<(int)(sizeof(files)/sizeof(files[0])); i++) {

    if (files[i] == "")
      continue;

    TFile* f = new TFile(files[i].c_str(),"READ");

    c->cd(1);

    gPad->SetGrid();
    gPad->SetLogx();

    TGraph* g1 = (TGraph*) f->Get("Chi2CutCheck/Graph_MtopSymmetry_trough_diff");

    g1->GetYaxis()->SetLabelSize(0.03);
    g1->GetYaxis()->SetTitleOffset(1.);
    g1->GetYaxis()->SetRangeUser(0,300);
 //   g1->GetXaxis()->SetRangeUser(0,10000);

    g1->Draw("A*");

    double yval1 = 0;

    for (int p=0; p<g1->GetN();p++) {

      Double_t x=0;
      Double_t y=0;

      g1->GetPoint(p,x,y);

      if (x == cut)
	yval1=y+2;

    }

    cout << yval1 << endl;

    TLine* line = new TLine(cut,yval1,cut,60);
    line->SetLineColor(kRed);
    line->SetLineWidth(3);
if (cut > 0)     line->Draw();
    line = new TLine(cut,yval1,cut-cutW,yval1+10);
    line->SetLineColor(kRed);
    line->SetLineWidth(3);
if (cut > 0)     line->Draw();
    line = new TLine(cut,yval1,cut+cutW,yval1+10);
    line->SetLineColor(kRed);
    line->SetLineWidth(3);
if (cut > 0)     line->Draw();

    c->cd(2);

    gPad->SetGrid();
    gPad->SetLogx(0);

    TGraph* g2 = (TGraph*) f->Get("Chi2CutCheck/Graph_NeventsRemoved_FromTopMassPeak");

    g2->GetYaxis()->SetRangeUser(0.5,1.2);
    g2->GetXaxis()->SetLimits(1,1000);

    g2->GetYaxis()->SetLabelSize(0.03);
    g2->GetYaxis()->SetTitleOffset(1.);

    g2->Draw("A*");

    double yval = 0;

    for (int p=0; p<g2->GetN();p++) {

      Double_t x=0;
      Double_t y=0;

      g2->GetPoint(p,x,y);

      if (x == cut || (x > cut && yval == 0) )
	yval=y-(y*0.01);
        
        //cout << x << endl;

    }
      
      cout << yval << endl;

    line = new TLine(cut,0.5,cut,yval);
    line->SetLineColor(kRed);
    line->SetLineWidth(3);
if (cut > 0)     line->Draw();
    line = new TLine(cut,yval,cut-cutW,0.9);
    line->SetLineColor(kRed);
    line->SetLineWidth(3);
if (cut > 0)     line->Draw();
    line = new TLine(cut,yval,cut+cutW,0.9);
    line->SetLineColor(kRed);
    line->SetLineWidth(3);
if (cut > 0)     line->Draw();

        gPad->SetLogx();


    c->Draw();

  }

}
